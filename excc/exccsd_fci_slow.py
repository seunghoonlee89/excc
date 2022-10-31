import numpy as np

import pyscf
from pyscf.cc import ccsd
from pyscf.cc import rintermediates as imd

from pyscf.mcscf import casci
from pyscf       import lib
from pyscf.lib   import logger

from pyscf.ci.cisd  import tn_addrs_signs
from ._fci_ci_to_cc import _c1_to_t1, _c2_to_t2, _c3_to_t3, _c4_to_t4

from scipy import linalg as la
from scipy import optimize as opt
from scipy.sparse import linalg as spla

NUM_ZERO = 1e-8

def kernel(the_exccsd, eris=None, t1=None, t2=None, max_cycle=50,
                      tol=1e-8, tolnormt=1e-6, verbose=None):
    log = logger.new_logger(the_exccsd, verbose)
    if eris is None:
        eris = the_exccsd.ao2mo(mo_coeff = the_exccsd.mo_coeff)

    assert t1 is None and t2 is None

    t1, t2, t1_t3c, t2_t3t4c = the_exccsd.get_ext_t_amps(eris)
    cput1   = cput0 = (logger.process_clock(), logger.perf_counter())

    #dbg
    #t1, t2 = the_exccsd.get_init_guess(eris)
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

    if the_exccsd.imag_tevol:
        l = the_exccsd.l 
        dl = the_exccsd.dl 
        fac = the_exccsd.fac
        order = the_exccsd.order
        thresh = the_exccsd.thresh
        real = the_exccsd.real

        max_cycle = int(l / dl + 0.1 * dl)
        if not real:
            t1 = np.array(t1, dtype=np.complex128)
            t2 = np.array(t2, dtype=np.complex128)

        normdt = 1.
        t = 0.0
        for istep in range(max_cycle):
            dt1,dt2 = the_exccsd.RK(t1, t2, t1_t3c, t2_t3t4c, eris, dl, fac=fac, order=order)
            tmpvec = the_exccsd.amplitudes_to_vector(dt1, dt2)
            normdt_prev, normdt = normdt, np.linalg.norm(tmpvec)
            t1 += dl * dt1
            t2 += dl * dt2
            tmpvec = the_exccsd.amplitudes_to_vector(t1, t2)
            normt = np.linalg.norm(tmpvec)
            e_prev, e_exccsd = e_exccsd, the_exccsd.energy(t1, t2, eris)
            t += dl
            if real:
                log.info('time={:.2f},Ecorr={:.8f},dE={:.4e},norm(dt1,dt2)={:.4e},normt={:.4f}'.format(t,e_exccsd,e_exccsd-e_prev,normdt,normt))
            else:
                tmpvec = the_exccsd.amplitudes_to_vector(t1.imag, t2.imag)
                normi = np.linalg.norm(tmpvec)
                log.info('time={:.2f},Ecorr={:.8f},dE={:.4e},norm(dt1,dt2)={:.4e},normt={:.4f},normi={:.4e}'.format(t,e_exccsd,e_exccsd-e_prev,normdt,normt,normi))
            tmpvec = None
            if abs(e_exccsd-e_prev) < tol and normdt < tolnormt:
                is_converged = True 
                break
            #if normdt > normdt_prev:
            #    is_converged = False 
            #    break
            #if e_exccsd > e_prev:
            #    is_converged = False 
            #    break
            if abs(e_exccsd) > thresh:
                break
    elif the_exccsd.minres:
        cycle = [0] 
        x0 = the_exccsd.amplitudes_to_vector(t1, t2) 
        def f_res(x):
            t1, t2 = the_exccsd.vector_to_amplitudes(x)
            e_exccsd = the_exccsd.energy(t1, t2, eris)
            t1, t2 = the_exccsd.update_amps(t1, t2, t1_t3c, t2_t3t4c, eris)
            res = the_exccsd.amplitudes_to_vector(t1, t2) 
            norm = max_abs(res)
            log.info("      cycle = %5d , E = %15.8g , norm(res) = %15.5g", cycle[0],
                     e_exccsd, norm)
            cycle[0] += 1
            return res 

        if the_exccsd.precond == 'finv':
            def mop(x):
                return the_exccsd.precond_finv(x, eris)
            M = spla.LinearOperator((x0.shape[-1], x0.shape[-1]), matvec=mop)
        elif the_exccsd.precond == 'diag':
            def mop(x):
                return the_exccsd.precond_diag(x, eris)
            M = spla.LinearOperator((x0.shape[-1], x0.shape[-1]), matvec=mop)
        else:
            M = None
    
        if the_exccsd.method == 'krylov':
            inner_m = the_exccsd.inner_m
            outer_k = the_exccsd.outer_k
            res = opt.root(f_res, x0, method='krylov',
                           options={'fatol': tolnormt, 'tol_norm': safe_max_abs, 
                                    'disp': True, 'maxiter': max_cycle // inner_m,
                                    'line_search': 'wolfe',
                                    'jac_options': {'rdiff': 1e-6, 'inner_maxiter': 100, 
                                                    'inner_inner_m': inner_m, 'inner_tol': tolnormt * 0.5,
                                                    'outer_k': outer_k, 'inner_M': M}
                                   })
        elif the_exccsd.method == 'df-sane':
            res = opt.root(f_res, x0, method='df-sane',
                           options={'fatol': tolnormt, 'disp': True, 'maxfev': max_cycle,
                                    'fnorm': max_abs})
        else:
            raise ValueError

        is_converged = res.success
        t1, t2 = the_exccsd.vector_to_amplitudes(res.x)
        e_exccsd = the_exccsd.energy(t1, t2, eris)
    else: 
        while not is_converged and not is_diverged and not is_max_cycle:
            t1_new, t2_new = the_exccsd.update_amps(t1, t2, t1_t3c, t2_t3t4c, eris)
    
            tmp_vec  = the_exccsd.amplitudes_to_vector(t1, t2)
            tmp_vec -= the_exccsd.amplitudes_to_vector(t1_new, t2_new)
            err_amp_prev, err_amp  = err_amp, np.linalg.norm(tmp_vec)
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
                is_converged = abs(e_exccsd - e_prev) < tol or err_amp - err_amp_prev > 0 #or e_exccsd - e_prev > 0 
            istep += 1

    log.info('max_memory %d MB (current use %d MB)',
             the_exccsd.max_memory, lib.current_memory()[0])

    log.timer('exCCSD', *cput0)
    return is_converged, e_exccsd, t1, t2

def update_amps(cc, t1, t2, t1_t3c, t2_t3t4c, eris, fac=1.0):
    assert(isinstance(eris, ccsd._ChemistsERIs))
    nocc, nvir = t1.shape
    fock = eris.fock
    if not cc.imag_tevol and not cc.minres:
        mo_e_o = eris.mo_energy[:nocc]
        mo_e_v = eris.mo_energy[nocc:] + cc.level_shift

    fov = fock[:nocc,nocc:].copy()
    foo = fock[:nocc,:nocc].copy()
    fvv = fock[nocc:,nocc:].copy()

    Foo = imd.cc_Foo(t1,t2,eris)
    Fvv = imd.cc_Fvv(t1,t2,eris)
    Fov = imd.cc_Fov(t1,t2,eris)

    if not cc.imag_tevol and not cc.minres:
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
        if not cc.imag_tevol and not cc.minres:
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

    if not cc.imag_tevol and not cc.minres:
        eia = mo_e_o[:,None] - mo_e_v
        eijab = lib.direct_sum('ia,jb->ijab',eia,eia)
        t1new /= eia
        t2new /= eijab
    elif cc.imag_tevol:
        t1new *= -fac 
        t2new *= -fac

    return t1new, t2new

def get_ext_t_amps(the_exccsd, eris):
    """
       exccsd_slow version saves c3, c4, t3, t4 amplitudes in memory
    """
    ncas  = the_exccsd.ncas
    ncore = the_exccsd.ncore
    nocc  = the_exccsd.nocc
    nmo   = the_exccsd.nmo
    nvir  = nmo - nocc

    ncas_vir = ncore + ncas - nocc
    ncas_occ = ncas - ncas_vir

    c1a, c2aa, c2ab, c3aaa, c3aab, c4aaab, c4aabb = the_exccsd.get_ext_c_amps()

    t1a_cas               = _c1_to_t1(c1a.copy(), ncas_occ, ncas_vir)
    t2aa_cas,  t2ab_cas   = _c2_to_t2(t1a_cas, c2aa.copy(), c2ab.copy(), \
                                      ncas_occ, ncas_vir)
    t3aaa_cas, t3aab_cas  = _c3_to_t3(t1a_cas, t2aa_cas, t2ab_cas, c3aaa.copy(), \
                                      c3aab.copy(), ncas_occ, ncas_vir)
    t4aaab_cas,t4aabb_cas = _c4_to_t4(t1a_cas, t2aa_cas, t2ab_cas, t3aaa_cas, t3aab_cas, \
                                      c4aaab, c4aabb, ncas_occ, ncas_vir)

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

def max_abs(x):
    """
    Equivalent to np.max(np.abs(x)), but faster.
    """
    if np.iscomplexobj(x):
        return np.abs(x).max()
    else:
        return max(np.max(x), abs(np.min(x)))

def safe_max_abs(x):
    if np.isfinite(x).all():
        return max(np.max(x), abs(np.min(x)))
    else:
        return 1e+12

def precond_finv(cc, vec, eris, tol=1e-8):
    """
    Fock inversion as preconditioner.
    """
    vec = vec.ravel()
    t1, t2 = cc.vector_to_amplitudes(vec)
    nocc, nvir = t1.shape
    mo_e_o = eris.mo_energy[:nocc]
    mo_e_v = eris.mo_energy[nocc:] + cc.level_shift
    
    eia = mo_e_o[:, None] - mo_e_v
    eia[eia > -tol] = -tol
    t1 /= eia
    
    for i in range(nocc):
        t2[i] /= lib.direct_sum('a, jb -> jab', eia[i], eia)
    return cc.amplitudes_to_vector(t1, t2)

def precond_diag(cc, vec, eris):
    """
    Diagonal elements as preconditioner, works not well.
    """
    vec = vec.ravel()
    t1, t2 = cc.vector_to_amplitudes(vec)
    nocc, nvir = t1.shape

    fock = eris.fock
    
    mo_e_o = eris.mo_energy[:nocc]
    mo_e_v = eris.mo_energy[nocc:] + cc.level_shift
    eia = mo_e_o[:, None] - mo_e_v
    eia -= np.einsum('iaai -> ia', eris.ovvo)
    t1 /= eia
    
    eijab = lib.direct_sum('ia,jb->ijab', eia, eia)
    tmp = 0.5 * np.einsum('abab -> ab', eris.vvvv)
    for i in range(nocc):
        eijab[i, i] -= tmp
    tmp = 0.5 * np.einsum('ijij -> ij', eris.oooo)
    for i in range(nvir):
        eijab[:, :, i, i] -= tmp
    t2 /= eijab
    return cc.amplitudes_to_vector(t1, t2)


class EXCCSD(ccsd.CCSD):
    def __init__(self, mc, frozen=None, mo_coeff=None, mo_occ=None):
        assert isinstance(mc, casci.CASCI)
        assert isinstance(mc.fcisolver, pyscf.fci.direct_spin1.FCISolver)
        ccsd.CCSD.__init__(self, mc._scf, frozen=frozen, mo_coeff=mo_coeff, mo_occ=mo_occ)
        keys = set(('max_cycle', 'conv_tol', 'iterative_damping',
                    'conv_tol_normt', 'diis', 'diis_space', 'diis_file',
                    'diis_start_cycle', 'diis_start_energy_diff', 'direct',
                    'async_io', 'incore_complete', 'cc2', "ncore", "ncas",
                    "nelec_cas", "_casci", "ci", "sparse_tol", "ext_orbs"))
        self._keys = set(self.__dict__.keys()).union(keys)
        self._casci  = mc
        self.ncore   = mc.ncore
        self.ncas    = mc.ncas
        self.nelec_cas = (mc.nelecas[0], mc.nelecas[1])
        self.ci      = mc.ci
        self.sparse_tol = None
        self.ext_orbs = [] 

        # parameters for imaginary time evolution
        self.imag_tevol = False
        self.l   = 100. 
        self.dl  = 0.1
        self.fac = 1. 
        self.order = 4 
        self.thresh = 100. 
        self.real = True 
        keys = keys.union(["imag_tevol", "l", "dl", \
                           "fac", "order", "thresh", "real"])

        # parameters for min res
        self.minres = False 
        self.method = "krylov" 
        self.precond = "finv" 
        self.inner_m = 10 
        self.outer_k = 6 
        keys = keys.union(["minres", "method", "precond", "inner_m", "outer_k"])

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
        log.info('imag_tevol             = %r', self.imag_tevol)
        if self.imag_tevol:
            log.info('l                      = %f', self.l)
            log.info('dl                     = %f', self.dl)
            log.info('fac                    = %f', self.fac)
            log.info('order                  = %d', self.order)
            log.info('thresh                 = %f', self.thresh)
            log.info('real                   = %r', self.real)
        log.info('min res solver         = %r', self.minres)
        if self.minres:
            log.info('method                 = %s', self.method)
            log.info('precond                = %s', self.precond)
            log.info('inner_m                = %d', self.inner_m)
            log.info('outer_k                = %d', self.outer_k)
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

    def get_ext_c_amps(self):
        ncas  = self.ncas
        ncas_occ = self.nocc - self.ncore
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

#        print('c0')
#        print(c0)
#        print('c1a/c0')
#        print(c1a[index1])
#        print('c2ab/c0')
#        print(c2ab[index1,:][:,index1])

        return c1a[index1], c2aa[index2], c2ab[index1,:][:,index1], \
               c3aaa[index3], c3aab[index2,:][:,index1], \
               c4aaab[index3,:][:,index1], c4aabb[index2,:][:,index2]

    def compute_derivative(self, t1, t2, t1_t3c, t2_t3t4c, eris, fac=1.0):
        return update_amps(self, t1, t2, t1_t3c, t2_t3t4c, eris, fac=fac)

    def RK(self, t1, t2, t1_t3c, t2_t3t4c, eris, h, fac=1.0, order=1):
        dt11, dt21 = self.compute_derivative(t1, t2, t1_t3c, t2_t3t4c, eris, fac=fac)
        if order == 1:
            dt1, dt2 = dt11, dt21
        else:
            t1_, t2_ = t1 + dt11 * h * 0.5, t2 + dt21 * h * 0.5
            dt12, dt22 = self.compute_derivative(
                 t1_, t2_, t1_t3c, t2_t3t4c, eris, fac=fac)
            t1_, t2_ = t1 + dt12 * h * 0.5, t2 + dt22 * h * 0.5
            dt13, dt23 = self.compute_derivative(
                 t1_, t2_, t1_t3c, t2_t3t4c, eris, fac=fac)
            t1_, t2_ = t1 + dt13 * h, t2 + dt23 * h
            dt14, dt24 = self.compute_derivative(
                 t1_, t2_, t1_t3c, t2_t3t4c, eris, fac=fac)
            dt1 = (dt11 + dt12 * 2.0 + dt13 * 2.0 + dt14) / 6.0
            dt2 = (dt21 + dt22 * 2.0 + dt23 * 2.0 + dt24) / 6.0
        return dt1, dt2

    get_ext_t_amps = get_ext_t_amps
    precond_finv = precond_finv
    precond_diag = precond_diag

class EXRCCSD(EXCCSD):
    pass

FCIEXCCSD = EXCCSD

