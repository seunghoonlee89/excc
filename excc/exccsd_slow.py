import imp
import sys
import numpy as np
import os

import pyscf
from pyscf.cc import ccsd
from pyscf.cc import rccsd
from pyscf.cc import ccsd_lambda
from pyscf.cc import rintermediates as imd

from pyscf.mcscf import casci
from pyscf       import lib
from pyscf.lib   import logger

from pyscf.ci.cisd  import CISD
from pyscf.ci.cisd  import contract
from pyscf.ci.cisd  import from_fcivec
from pyscf.ci.cisd  import cisdvec_to_amplitudes
from pyscf.ci.cisd  import amplitudes_to_cisdvec
from pyscf.ci.cisd  import tn_addrs_signs
from pyscf.fci.cistring import str2addr 
from pyscf.tools.dump_mat import dump_rec

from ._fci_ci_to_cc import _c1_to_t1, _c2_to_t2, _c3_to_t3, _c4_to_t4
#from .dmrgscf_b2 import DMRGCI

NUM_ZERO = 1e-8

def kernel(the_exccsd, eris=None, t1=None, t2=None, max_cycle=50,
                      tol=1e-8, tolnormt=1e-6, verbose=None):
    log = logger.new_logger(the_exccsd, verbose)
    if eris is None:
        eris = the_exccsd.ao2mo(mo_coeff = the_exccsd.mo_coeff)

    assert t1 is None and t2 is None

    t1, t2, t1_t3c, t2_t3t4c = the_exccsd.get_ext_amps(eris)
    cput1   = cput0 = (logger.process_clock(), logger.perf_counter())
    e_exccsd = the_exccsd.energy(t1, t2, eris)
    e_prev  = e_exccsd
    log.info('Init E_corr(exCCSD) = %.15g', e_exccsd)

    if isinstance(the_exccsd.diis, lib.diis.DIIS):
        adiis = the_exccsd.diis
    elif the_exccsd.diis:
        adiis = lib.diis.DIIS(the_exccsd, the_exccsd.diis_file, incore=the_exccsd.incore_complete)
        adiis.space = the_exccsd.diis_space
    else:
        adiis = None

    alpha = the_exccsd.iterative_damping
    istep = 0
    
    is_converged = False
    is_diverged  = False 
    is_max_cycle = False

    err_ene = 1.0
    err_amp = 1.0

    while not is_converged and not is_diverged and not is_max_cycle:
        t1_new, t2_new = the_exccsd.update_amps(t1, t2, t1_t3c, t2_t3t4c, eris)

        tmp_vec  = the_exccsd.amplitudes_to_vector(t1, t2)
        tmp_vec -= the_exccsd.amplitudes_to_vector(t1_new, t2_new)
        err_amp  = np.linalg.norm(tmp_vec)
        tmp_vec  = None

        if alpha < 1.0:
            t1_new  = (1-alpha) * t1 + alpha * t1_new
            t2_new *= alpha
            t2_new += (1-alpha) * t2

        t1, t2 = t1_new, t2_new
        t1_new = t2_new = None
        
        t1, t2 = the_exccsd.run_diis(t1, t2, istep, err_amp, err_ene, adiis)

        e_prev, e_exccsd  = e_exccsd, the_exccsd.energy(t1, t2, eris)
        err_ene          = (e_exccsd - e_prev)

        log.info('cycle = %d  E_corr(exCCSD) = %.15g  dE = %.9g  norm(t1,t2) = %.6g', istep+1, e_exccsd, err_ene, err_amp)
        cput1 = log.timer('exCCSD iter', *cput1)

        is_converged = abs(err_ene) < tol and err_amp < tolnormt
        is_max_cycle = istep == max_cycle
        is_diverged = abs(e_exccsd - e_prev) > 100
        if adiis is None and the_exccsd.level_shift > 0. and alpha < 1.:  
            is_converged = abs(e_exccsd - e_prev) < tol or e_exccsd - e_prev > 0
        istep += 1

    log.timer('exCCSD', *cput0)
    return is_converged, e_exccsd, t1, t2

def update_amps(cc, t1, t2, t1_t3c, t2_t3t4c, eris):
    assert(isinstance(eris, ccsd._ChemistsERIs))
    nocc, nvir = t1.shape
    fock = eris.fock
    mo_e_o = eris.mo_energy[:nocc]
    mo_e_v = eris.mo_energy[nocc:] + cc.level_shift

    fov = fock[:nocc,nocc:].copy()
    foo = fock[:nocc,:nocc].copy()
    fvv = fock[nocc:,nocc:].copy()

    Foo = imd.cc_Foo(t1,t2,eris)
    Fvv = imd.cc_Fvv(t1,t2,eris)
    Fov = imd.cc_Fov(t1,t2,eris)

    # Move energy terms to the other side
    Foo[np.diag_indices(nocc)] -= mo_e_o
    Fvv[np.diag_indices(nvir)] -= mo_e_v

    # T1 equation
    t1new  = t1_t3c.copy()
    t1new +=-2*np.einsum('kc,ka,ic->ia', fov, t1, t1)
    t1new +=   np.einsum('ac,ic->ia', Fvv, t1)
    t1new +=  -np.einsum('ki,ka->ia', Foo, t1)
    t1new += 2*np.einsum('kc,kica->ia', Fov, t2)
    t1new +=  -np.einsum('kc,ikca->ia', Fov, t2)
    t1new +=   np.einsum('kc,ic,ka->ia', Fov, t1, t1)
    t1new += fov.conj()
    t1new += 2*np.einsum('kcai,kc->ia', eris.ovvo, t1)
    t1new +=  -np.einsum('kiac,kc->ia', eris.oovv, t1)
    eris_ovvv = np.asarray(eris.get_ovvv())
    t1new += 2*lib.einsum('kdac,ikcd->ia', eris_ovvv, t2)
    t1new +=  -lib.einsum('kcad,ikcd->ia', eris_ovvv, t2)
    t1new += 2*lib.einsum('kdac,kd,ic->ia', eris_ovvv, t1, t1)
    t1new +=  -lib.einsum('kcad,kd,ic->ia', eris_ovvv, t1, t1)
    eris_ovoo = np.asarray(eris.ovoo, order='C')
    t1new +=-2*lib.einsum('lcki,klac->ia', eris_ovoo, t2)
    t1new +=   lib.einsum('kcli,klac->ia', eris_ovoo, t2)
    t1new +=-2*lib.einsum('lcki,lc,ka->ia', eris_ovoo, t1, t1)
    t1new +=   lib.einsum('kcli,lc,ka->ia', eris_ovoo, t1, t1)

    # T2 equation
    t2new = t2_t3t4c.copy()
    tmp2  = lib.einsum('kibc,ka->abic', eris.oovv, -t1)
    tmp2 += np.asarray(eris_ovvv).conj().transpose(1,3,0,2)
    tmp = lib.einsum('abic,jc->ijab', tmp2, t1)
    t2new += tmp + tmp.transpose(1,0,3,2)
    tmp2  = lib.einsum('kcai,jc->akij', eris.ovvo, t1)
    tmp2 += eris_ovoo.transpose(1,3,0,2).conj()
    tmp = lib.einsum('akij,kb->ijab', tmp2, t1)
    t2new -= tmp + tmp.transpose(1,0,3,2)
    t2new += np.asarray(eris.ovov).conj().transpose(0,2,1,3)
    if cc.cc2:
        Woooo2 = np.asarray(eris.oooo).transpose(0,2,1,3).copy()
        Woooo2 += lib.einsum('lcki,jc->klij', eris_ovoo, t1)
        Woooo2 += lib.einsum('kclj,ic->klij', eris_ovoo, t1)
        Woooo2 += lib.einsum('kcld,ic,jd->klij', eris.ovov, t1, t1)
        t2new += lib.einsum('klij,ka,lb->ijab', Woooo2, t1, t1)
        Wvvvv = lib.einsum('kcbd,ka->abcd', eris_ovvv, -t1)
        Wvvvv = Wvvvv + Wvvvv.transpose(1,0,3,2)
        Wvvvv += np.asarray(eris.vvvv).transpose(0,2,1,3)
        t2new += lib.einsum('abcd,ic,jd->ijab', Wvvvv, t1, t1)
        Lvv2 = fvv - np.einsum('kc,ka->ac', fov, t1)
        Lvv2 -= np.diag(np.diag(fvv))
        tmp = lib.einsum('ac,ijcb->ijab', Lvv2, t2)
        t2new += (tmp + tmp.transpose(1,0,3,2))
        Loo2 = foo + np.einsum('kc,ic->ki', fov, t1)
        Loo2 -= np.diag(np.diag(foo))
        tmp = lib.einsum('ki,kjab->ijab', Loo2, t2)
        t2new -= (tmp + tmp.transpose(1,0,3,2))
    else:
        Loo = imd.Loo(t1, t2, eris)
        Lvv = imd.Lvv(t1, t2, eris)
        Loo[np.diag_indices(nocc)] -= mo_e_o
        Lvv[np.diag_indices(nvir)] -= mo_e_v

        Woooo = imd.cc_Woooo(t1, t2, eris)
        Wvoov = imd.cc_Wvoov(t1, t2, eris)
        Wvovo = imd.cc_Wvovo(t1, t2, eris)
        Wvvvv = imd.cc_Wvvvv(t1, t2, eris)

        tau = t2 + np.einsum('ia,jb->ijab', t1, t1)
        t2new += lib.einsum('klij,klab->ijab', Woooo, tau)
        t2new += lib.einsum('abcd,ijcd->ijab', Wvvvv, tau)
        tmp = lib.einsum('ac,ijcb->ijab', Lvv, t2)
        t2new += (tmp + tmp.transpose(1,0,3,2))
        tmp = lib.einsum('ki,kjab->ijab', Loo, t2)
        t2new -= (tmp + tmp.transpose(1,0,3,2))
        tmp  = 2*lib.einsum('akic,kjcb->ijab', Wvoov, t2)
        tmp -=   lib.einsum('akci,kjcb->ijab', Wvovo, t2)
        t2new += (tmp + tmp.transpose(1,0,3,2))
        tmp = lib.einsum('akic,kjbc->ijab', Wvoov, t2)
        t2new -= (tmp + tmp.transpose(1,0,3,2))
        tmp = lib.einsum('bkci,kjac->ijab', Wvovo, t2)
        t2new -= (tmp + tmp.transpose(1,0,3,2))

    eia = mo_e_o[:,None] - mo_e_v
    eijab = lib.direct_sum('ia,jb->ijab',eia,eia)
    t1new /= eia
    t2new /= eijab

    return t1new, t2new

class EXCCSD(ccsd.CCSD):
    def __init__(self, mc, frozen=None, mo_coeff=None, mo_occ=None):
        assert isinstance(mc, casci.CASCI)
        assert isinstance(mc.fcisolver, pyscf.fci.direct_spin1.FCISolver)
        ccsd.CCSD.__init__(self, mc._scf, frozen=frozen, mo_coeff=mo_coeff, mo_occ=mo_occ)
        keys = set(('max_cycle', 'conv_tol', 'iterative_damping',
                    'conv_tol_normt', 'diis', 'diis_space', 'diis_file',
                    'diis_start_cycle', 'diis_start_energy_diff', 'direct',
                    'async_io', 'incore_complete', 'cc2', "ncore", "ncas", "nelec_cas", "_casci", "ci", "sparse_tol"))
        self._keys = set(self.__dict__.keys()).union(keys)

        self._casci  = mc
        self.ncore   = mc.ncore
        self.ncas    = mc.ncas
        self.nelec_cas = (mc.nelecas[0], mc.nelecas[1])
        self.ci      = mc.ci

        self.sparse_tol = None

    def dump_flags(self, verbose=None):
        ncas  = self.ncas
        ncore = self.ncore
        nocc  = self.nocc
        nmo   = self.nmo
        nvir  = nmo - nocc

        ncas_vir = ncore + ncas - nocc
        ncas_occ = ncas - ncas_vir

        log = logger.new_logger(self, verbose)
        log.info('')
        log.info('******** %s ********', self.__class__)
        log.info('CC2        = %g', self.cc2)
        log.info('exCCSD nocc     = %s, nvir     = %s, nmo = %s', nocc, nvir, nmo)
        log.info('exCCSD ncas_occ = %s, ncas_vir = %s, ncas = %s', ncas_occ, ncas_vir, ncas)
        if self.frozen is not None:
            log.info('frozen orbitals %s', self.frozen)
        if self.sparse_tol is not None:
            log.info('sparse_tol             = %g', self.sparse_tol)
        log.info('max_cycle              = %d', self.max_cycle)
        log.info('direct                 = %d', self.direct)
        log.info('conv_tol               = %g', self.conv_tol)
        log.info('conv_tol_normt         = %s', self.conv_tol_normt)
        log.info('diis_space             = %d', self.diis_space)
        log.info('diis_start_cycle       = %d', self.diis_start_cycle)
        log.info('diis_start_energy_diff = %g', self.diis_start_energy_diff)
        log.info('max_memory %d MB (current use %d MB)',
                 self.max_memory, lib.current_memory()[0])
        return self

    def update_amps(self, t1, t2, t1_t3c, t2_t3t4c, eris):
        return update_amps(self, t1, t2, t1_t3c, t2_t3t4c, eris)

    def kernel(self, t1=None, t2=None, eris=None):
        return self.exccsd(t1, t2, eris)

    def exccsd(self, t1=None, t2=None, eris=None):
        assert(self.mo_coeff is not None)
        assert(self.mo_occ   is not None)

        if self.verbose >= logger.WARN:
            self.check_sanity()
        self.dump_flags()

        if eris is None:
            eris = self.ao2mo(self.mo_coeff)

        self.e_hf = getattr(eris, 'e_hf', None)
        if self.e_hf is None:
            self.e_hf = self._scf.e_tot

        self.converged, self.e_corr, self.t1, self.t2 = \
                kernel(self, eris, t1, t2, max_cycle=self.max_cycle,
                       tol=self.conv_tol, tolnormt=self.conv_tol_normt,
                       verbose=self.verbose)
        self._finalize()
        return self.e_corr, self.t1, self.t2

    def get_cisd_amplitudes(self, cicas=None):
        if cicas is None:
            cicas = self.ci
        ncas  = self.ncas
        ncore = self.ncore
        nocc  = self.nocc
        nmo   = self.nmo
        nvir  = nmo - nocc
        nelec_a, nelec_b = self._scf.mol.nelec
        nelec_cas_a, nelec_cas_b = self.nelec_cas

        ncas_vir = ncore + ncas - nocc
        ncas_occ = ncas - ncas_vir

        cisd_cas_vec = from_fcivec(cicas, ncas, (nelec_cas_a, nelec_cas_b))
        c0cas, c1cas, c2cas = cisdvec_to_amplitudes(cisd_cas_vec, ncas, ncas_occ)

        c0 = c0cas

        c1 = np.zeros((nocc, nvir))
        c1[ncore:ncore+ncas_occ, 0:ncas_vir] = c1cas

        c2 = np.zeros((nocc, nocc, nvir, nvir))
        c2[ncore:ncore+ncas_occ, ncore:ncore+ncas_occ, 0:ncas_vir, 0:ncas_vir] = c2cas

        return c0, c1, c2

    def get_ext_amps(self, eris):
        """
           exccsd_slow version saves c3, c4, t3, t4 amplitudes in memory
        """
        ncas  = self.ncas
        ncore = self.ncore
        nocc  = self.nocc
        nmo   = self.nmo
        nvir  = nmo - nocc

        ncas_vir = ncore + ncas - nocc
        ncas_occ = ncas - ncas_vir

        t1addrs, t1signs = tn_addrs_signs(ncas, ncas_occ, 1)
        index1 = np.argsort(t1addrs)
        t2addrs, t2signs = tn_addrs_signs(ncas, ncas_occ, 2)
        index2 = np.argsort(t2addrs)
        t3addrs, t3signs = tn_addrs_signs(ncas, ncas_occ, 3)
        index3 = np.argsort(t3addrs)

        fcivec  = self.ci
        c0      = fcivec[0, 0]
        assert np.abs(c0) > NUM_ZERO
        c1a     = fcivec[t1addrs, 0] * t1signs / c0
        c2aa    = fcivec[t2addrs, 0] * t2signs / c0 
        c2ab    = fcivec[t1addrs[:,None], t1addrs] \
                  * np.einsum("i,j->ij", t1signs, t1signs) / c0
        c3aaa   = fcivec[t3addrs, 0] * t3signs / c0
        c3aab   = fcivec[t2addrs[:,None], t1addrs] \
                  * np.einsum("i,j->ij", t2signs, t1signs) / c0
        c4aaab  = fcivec[t3addrs[:,None], t1addrs] \
                  * np.einsum('i,j->ij', t3signs, t1signs) / c0
        c4aabb  = fcivec[t2addrs[:,None], t2addrs] \
                  * np.einsum('i,j->ij', t2signs, t2signs) / c0

        t1a_cas               = _c1_to_t1(c1a[index1].copy(), ncas_occ, ncas_vir)
        t2aa_cas,  t2ab_cas   = _c2_to_t2(t1a_cas, c2aa[index2].copy(), \
                                    c2ab[index1,:][:,index1].copy(), ncas_occ, ncas_vir)
        t3aaa_cas, t3aab_cas  = _c3_to_t3(t1a_cas, t2aa_cas, t2ab_cas, c3aaa[index3].copy(), \
                                    c3aab[index2,:][:,index1].copy(), \
                                    ncas_occ, ncas_vir)
        t4aaab_cas,t4aabb_cas = _c4_to_t4(t1a_cas, t2aa_cas, t2ab_cas, t3aaa_cas, t3aab_cas, \
                                    c4aaab[index3,:][:,index1], c4aabb[index2,:][:,index2], \
                                    ncas_occ, ncas_vir)

        nc  = ncore
        nco = ncore + ncas_occ
        ncov= ncore + ncas_occ + ncas_vir

        eris_ovov_cas = np.asarray(eris.ovvo)[ncore:, :ncas_vir, :ncas_vir, ncore:].transpose(0,1,3,2)
        t1_t3c_cas  = np.einsum('kcld,iklacd->ia', eris_ovov_cas, t3aaa_cas)
        t1_t3c_cas += np.einsum('kcld,iklacd->ia', eris_ovov_cas, t3aab_cas)
        t1_t3c_cas += np.einsum('kcld,ilkadc->ia', eris_ovov_cas, t3aab_cas)
        t1_t3c_cas += np.einsum('kcld,klicda->ia', eris_ovov_cas, t3aab_cas)
        t1_t3c_cas *= 0.5

        eris_ooov_cas = np.asarray(eris.ovoo)[ncore:, :ncas_vir, ncore:, ncore:].transpose(2,3,0,1)
        eris_ovvv = np.asarray(eris.get_ovvv())
        eris_ovvv_cas = eris_ovvv[ncore:, :ncas_vir, :ncas_vir, :ncas_vir]
        fov_cas = np.asarray(eris.fock)[nc:nco, nco:ncov]

        tmp1  = t4aaab_cas + t4aaab_cas.transpose(0,1,3,2,4,5,7,6)
        tmp1 += t4aabb_cas.transpose(0,2,1,3,4,6,5,7)
        tmp1 += t4aabb_cas.transpose(1,2,0,3,5,6,4,7)
        t2_t3t4c_cas = 0.5 * np.einsum('kcld,klijcdab->ijab', eris_ovov_cas, tmp1)

        t3tmp = t3aab_cas + t3aab_cas.transpose(0,2,1,3,5,4)
        tmp2  = np.einsum('kdlc,id->kilc', eris_ovov_cas, t1a_cas)
        tmp2 += eris_ooov_cas.conj()
        t2_t3t4c_cas -= np.einsum('kilc,lkjcab->ijab', tmp2, t3tmp)
        tmp2  = np.einsum('kdlc,jd->kjlc', eris_ovov_cas, t1a_cas)
        tmp2 += eris_ooov_cas.conj()
        t2_t3t4c_cas -= np.einsum('kjlc,likcab->ijab', tmp2, t3tmp)
   
        tmp2  = -np.einsum('kcld,lb->kcbd', eris_ovov_cas, t1a_cas)
        tmp2 += eris_ovvv_cas.conj().copy()
        t2_t3t4c_cas += np.einsum('kcbd,kijcad->ijab', tmp2, t3tmp)
        tmp2  = -np.einsum('kcld,la->kcad', eris_ovov_cas, t1a_cas)
        tmp2 += eris_ovvv_cas.conj().copy()
        t2_t3t4c_cas += np.einsum('kcad,kijcdb->ijab', tmp2, t3tmp)
    
        tmp2  = 2*np.einsum('kcld,ld->kc', eris_ovov_cas, t1a_cas)
        tmp2 +=  -np.einsum('kdlc,ld->kc', eris_ovov_cas, t1a_cas)
        tmp2 += np.asarray(fov_cas).conj().copy()
        t2_t3t4c_cas +=  np.einsum('kc,kijcab->ijab', tmp2, t3tmp)

        t1 = np.zeros((nocc, nvir))
        t2 = np.zeros((nocc, nocc, nvir, nvir))
        t1_t3c   = np.zeros((nocc, nvir))
        t2_t3t4c = np.zeros((nocc, nocc, nvir, nvir))
        t1[nc:, :ncas_vir] = t1a_cas
        t2[nc:, nc:, :ncas_vir, :ncas_vir] = t2ab_cas
        t1_t3c[nc:, :ncas_vir] = t1_t3c_cas
        t2_t3t4c[nc:, nc:, :ncas_vir, :ncas_vir] = t2_t3t4c_cas
        return t1, t2, t1_t3c, t2_t3t4c 

class REXCCSD(EXCCSD):
    pass

FCIEXCCSD = EXCCSD

#class DMRGEXCCSD(FCIEXCCSD):
#    def __init__(self, mc, frozen=None, mo_coeff=None, mo_occ=None):
#        assert isinstance(mc, casci.CASCI)
#        assert isinstance(mc.fcisolver, DMRGCI)
#        ccsd.CCSD.__init__(self, mc._scf, frozen=frozen, mo_coeff=mo_coeff, mo_occ=mo_occ)
#        keys = set(('max_cycle', 'conv_tol', 'iterative_damping',
#                    'conv_tol_normt', 'diis', 'diis_space', 'diis_file',
#                    'diis_start_cycle', 'diis_start_energy_diff', 'direct',
#                    'async_io', 'incore_complete', 'cc2', "ncore", "ncas", "nelec_cas",
#                    "ci", "t1_ext", "t2_ext"))
#        self._keys = set(self.__dict__.keys()).union(keys)
#
#        self._casci  = mc
#        self.ncore   = mc.ncore
#        self.ncas    = mc.ncas
#        self.nelec_cas = (mc.nelecas[0], mc.nelecas[1])
#        assert mc.nelecas[0] == mc.nelecas[1]   # for restricted case
#        self.ci     = None 
#
#        self.t1_ext = None
#        self.t2_ext = None
#
#        self.sparse_tol = None
#
#    def get_ext_amps(self):
#        ncas  = self.ncas
#        ncore = self.ncore
#        nocc  = self.nocc
#        nmo   = self.nmo
#        nvir  = nmo - nocc
#        nelec_cas = self.nelec_cas
#
#        ncas_vir = ncore + ncas - nocc
#        ncas_occ = ncas - ncas_vir
#        ncas_vir2 = ncas_vir*(ncas_vir-1)//2
#        ncas_occ2 = ncas_occ*(ncas_occ-1)//2
#
#        reorder = self._casci.fcisolver.dmrg_args['reorder']
#
#        def get_cisd_vec_cas(vals, dets):
#            c0   = np.zeros(1)
#            c1a  = np.zeros((ncas_occ*ncas_vir))
#            c2ab = np.zeros((ncas_occ*ncas_vir,ncas_occ*ncas_vir))
#            c2aa = np.zeros((ncas_occ2*ncas_vir2))
#
#            def convert_phase_to_Fermi_vacuum(h_a, h_b, val):
#                n_perm = 0 
#                ranka = len(h_a)
#                if ranka > 0:
#                    n_perm += ranka*ncas_occ - sum(h_a) - ranka*(ranka+1)//2
#                rankb = len(h_b)
#                if rankb > 0:
#                    n_perm += rankb*ncas_occ - sum(h_b) - rankb*(rankb+1)//2
#                return val * ( 1 - ( (n_perm & 1) << 1 ) )
#    
#            for val, det in zip(vals, dets):
#                occa =  det & 1
#                occb = (det & 2) >> 1 
#                # hole, ptcl
#                h_a, h_b, p_a, p_b = [], [], [], []
#                if len(reorder) == 0:
#                    for i in range(ncas_occ):
#                        if occa[i] == 0: h_a.append(i)
#                        if occb[i] == 0: h_b.append(i)
#                    for a in range(ncas_occ,ncas,1):
#                        if occa[a] == 1: p_a.append(a-ncas_occ)
#                        if occb[a] == 1: p_b.append(a-ncas_occ)
#                else:
#                    for i in range(ncas):
#                        org_i = reorder[i]
#                        if occa[i] == 0 and org_i < ncas_occ:  h_a.append(org_i)
#                        if occb[i] == 0 and org_i < ncas_occ:  h_b.append(org_i)
#                        if occa[i] == 1 and org_i >= ncas_occ: p_a.append(org_i-ncas_occ)
#                        if occb[i] == 1 and org_i >= ncas_occ: p_b.append(org_i-ncas_occ)
#                    h_a.sort()
#                    h_b.sort()
#                    p_a.sort()
#                    p_b.sort()
#
#                assert len(h_a) == len(p_a) and len(h_b) == len(p_b)
#                if   len(h_a) == 0 and len(h_b) == 0:
#                    c0[0] = val 
#                elif len(h_a) == 1 and len(h_b) == 0:
#                    idx_a = ncas_occ*(p_a[0]+1)-(h_a[0]+1)
#                    c1a[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
#                elif len(h_a) == 1 and len(h_b) == 1:
#                    idx_a = ncas_occ*(p_a[0]+1)-(h_a[0]+1)
#                    idx_b = ncas_occ*(p_b[0]+1)-(h_b[0]+1)
#                    c2ab[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
#                elif len(h_a) == 2 and len(h_b) == 0:
#                    h2 = h_a[1]*(h_a[1]-1)//2 + h_a[0] + 1 
#                    p2 = p_a[1]*(p_a[1]-1)//2 + p_a[0] + 1 
#                    idx_a = ncas_occ2*p2-h2
#                    c2aa[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
#                else:
#                    continue 
#
#            return c0, c1a, c2aa, c2ab 
#
#        scr = self._casci.fcisolver.dmrg_args["scratch"]
#        self._casci.fcisolver.sample_dets(ncas, nelec_cas)
#        vals = np.load("%s/sample-vals.npy"%(scr))
#        dets = np.load("%s/sample-dets.npy"%(scr))
#
#        c0, c1a, c2aa, c2ab = get_cisd_vec_cas(vals, dets)
#
#        c1a  /= c0
#        c2aa /= c0
#        c2ab /= c0
#
#        if self.verbose >= 5:
#            print('c0')
#            print(c0)
#            print('c1a')
#            print(c1a)
#            print('c2ab')
#            print(c2ab)
#
#        t1a        = _c1_to_t1(c1a, ncore, ncas_vir, ncas_occ, nvir, nocc)
#        t2aa, t2ab = _c2_to_t2(t1a, c2aa, c2ab, ncore, ncas_vir, ncas_occ, nvir, nocc)
#
#        t1a_cas  = t1a[ncore:ncore+ncas_occ, 0:ncas_vir]
#        t2aa_cas = t2aa[ncore:ncore+ncas_occ, ncore:ncore+ncas_occ, 0:ncas_vir, 0:ncas_vir]
#        t2ab_cas = t2ab[ncore:ncore+ncas_occ, ncore:ncore+ncas_occ, 0:ncas_vir, 0:ncas_vir]
#
#        return t1a_cas, (t2aa_cas, t2ab_cas)
#
#class DMRGREXCCSD(DMRGEXCCSD):
#    pass
#
#class HCIEXCCSD(FCIEXCCSD):
#    def __init__(self, mc, frozen=None, mo_coeff=None, mo_occ=None):
#        assert isinstance(mc, casci.CASCI)
#        assert isinstance(mc.fcisolver, pyscf.cornell_shci.shci.SHCI)
#
#        ccsd.CCSD.__init__(self, mc._scf, frozen=frozen, mo_coeff=mo_coeff, mo_occ=mo_occ)
#        keys = set(('max_cycle', 'conv_tol', 'iterative_damping',
#                    'conv_tol_normt', 'diis', 'diis_space', 'diis_file',
#                    'diis_start_cycle', 'diis_start_energy_diff', 'direct',
#                    'async_io', 'incore_complete', 'cc2', "ncore", "ncas", "nelec_cas",
#                    "ci", "t1_ext", "t2_ext"))
#        self._keys = set(self.__dict__.keys()).union(keys)
#
#        self.ncore   = mc.ncore
#        self.ncas    = mc.ncas
#        self.nelec_cas = (mc.nelecas[0], mc.nelecas[1])
#        self.ci      = mc.ci
#
#        self.t1_ext = None
#        self.t2_ext = None
#
#        self.output = 'output.dat'
#
#    @staticmethod
#    def read_hci_output(output):
#        vals, orb_as, orb_bs = [], [], []
#        readline = False
#        with open(output) as f:
#            for line in f:
#                if readline:
#                    if line == "----------------------------------------\n":
#                        break 
#                    val, orb_a, orb_b, _ = line.split('|')
#                    vals.append(float(val.split()[-1]))
#                    orb_as.append(np.array([int(i) for i in orb_a.split()]))
#                    orb_bs.append(np.array([int(i) for i in orb_b.split()]))
#                if line == "Excite Lv         Coef      Det (Reordered orb)\n":
#                    readline = True 
#        return vals, orb_as, orb_bs
#        
#    def get_ext_amps(self):
#        ncas  = self.ncas
#        ncore = self.ncore
#        nocc  = self.nocc
#        nmo   = self.nmo
#        nvir  = nmo - nocc
#
#        ncas_vir = ncore + ncas - nocc
#        ncas_occ = ncas - ncas_vir
#        ncas_vir2 = ncas_vir*(ncas_vir-1)//2
#        ncas_occ2 = ncas_occ*(ncas_occ-1)//2
#
#        def get_cisd_vec_cas(vals, orb_as, orb_bs):
#            c0   = np.zeros(1)
#            c1a  = np.zeros((ncas_occ*ncas_vir))
#            c2ab = np.zeros((ncas_occ*ncas_vir,ncas_occ*ncas_vir))
#            c2aa = np.zeros((ncas_occ2*ncas_vir2))
#
#            def convert_phase_to_Fermi_vacuum(h_a, h_b, val):
#                n_perm = 0 
#                ranka = len(h_a)
#                if ranka > 0:
#                    n_perm += ranka*ncas_occ - sum(h_a) - ranka*(ranka+1)//2
#                rankb = len(h_b)
#                if rankb > 0:
#                    n_perm += rankb*ncas_occ - sum(h_b) - rankb*(rankb+1)//2
#                return val * ( 1 - ( (n_perm & 1) << 1 ) )
# 
#            occ_ref = np.array(list(range(ncas_occ)))
#            for i in range(len(vals)):
#                val   = vals[i]
#                orb_a = orb_as[i]
#                orb_b = orb_bs[i]
#
#                h_a = np.setdiff1d(occ_ref, orb_a)
#                h_b = np.setdiff1d(occ_ref, orb_b)
#                p_a = np.setdiff1d(orb_a, occ_ref) - ncas_occ
#                p_b = np.setdiff1d(orb_b, occ_ref) - ncas_occ
#
#                assert len(h_a) == len(p_a) and len(h_b) == len(p_b)
#                if   len(h_a) == 0 and len(h_b) == 0:
#                    c0[0] = val 
#                elif len(h_a) == 1 and len(h_b) == 0:
#                    idx_a = ncas_occ*(p_a[0]+1)-(h_a[0]+1)
#                    c1a[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
#                elif len(h_a) == 1 and len(h_b) == 1:
#                    idx_a = ncas_occ*(p_a[0]+1)-(h_a[0]+1)
#                    idx_b = ncas_occ*(p_b[0]+1)-(h_b[0]+1)
#                    c2ab[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
#                elif len(h_a) == 2 and len(h_b) == 0:
#                    h2 = h_a[1]*(h_a[1]-1)//2 + h_a[0] + 1 
#                    p2 = p_a[1]*(p_a[1]-1)//2 + p_a[0] + 1 
#                    idx_a = ncas_occ2*p2-h2
#                    c2aa[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
#                else:
#                    continue 
#
#            return c0, c1a, c2aa, c2ab 
#
#        vals, orb_as, orb_bs = HCIEXCCSD.read_hci_output(self.output)
#        c0, c1a, c2aa, c2ab = get_cisd_vec_cas(vals, orb_as, orb_bs)
#
#        c1a  /= c0
#        c2aa /= c0
#        c2ab /= c0
#
#        t1a        = _c1_to_t1(c1a, ncore, ncas_vir, ncas_occ, nvir, nocc)
#        t2aa, t2ab = _c2_to_t2(t1a, c2aa, c2ab, ncore, ncas_vir, ncas_occ, nvir, nocc)
#
#        t1a_cas  = t1a[ncore:ncore+ncas_occ, 0:ncas_vir]
#        t2aa_cas = t2aa[ncore:ncore+ncas_occ, ncore:ncore+ncas_occ, 0:ncas_vir, 0:ncas_vir]
#        t2ab_cas = t2ab[ncore:ncore+ncas_occ, ncore:ncore+ncas_occ, 0:ncas_vir, 0:ncas_vir]
#
#        return t1a_cas, (t2aa_cas, t2ab_cas)
#
#class HCIREXCCSD(HCIEXCCSD):
#    pass
#
