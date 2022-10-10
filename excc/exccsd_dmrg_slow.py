import numpy as np

import pyscf
from pyscf.cc import ccsd
from pyscf.cc import rintermediates as imd

from pyscf.mcscf import casci
from pyscf       import lib
from pyscf.lib   import logger

from pyscf.ci.cisd  import tn_addrs_signs
from ._fci_ci_to_cc import _c1_to_t1, _c2_to_t2, _c3_to_t3, _c4_to_t4
from .exccsd_fci_slow import FCIEXCCSD
from pyscf import dmrgscf

NUM_ZERO = 1e-8

class DMRGEXCCSD(FCIEXCCSD):
    def __init__(self, mc, frozen=None, mo_coeff=None, mo_occ=None):
        assert isinstance(mc, casci.CASCI)
        assert isinstance(mc.fcisolver, dmrgscf.DMRGCI)
        ccsd.CCSD.__init__(self, mc._scf, frozen=frozen, mo_coeff=mo_coeff, mo_occ=mo_occ)
        keys = set(('max_cycle', 'conv_tol', 'iterative_damping',
                    'conv_tol_normt', 'diis', 'diis_space', 'diis_file',
                    'diis_start_cycle', 'diis_start_energy_diff', 'direct',
                    'async_io', 'incore_complete', 'cc2', "ncore", "ncas",
                    "nelec_cas", "_casci", "ci", "sparse_tol"))
        self._keys = set(self.__dict__.keys()).union(keys)

        self._casci  = mc
        self.ncore   = mc.ncore
        self.ncas    = mc.ncas
        self.nelec_cas = (mc.nelecas[0], mc.nelecas[1])
        assert mc.nelecas[0] == mc.nelecas[1]   # for restricted case
        self.ci      = None 
        self.sparse_tol = None

    def get_ext_c_amps(self):
        ncas  = self.ncas
        ncore = self.ncore
        nocc  = self.nocc
        nmo   = self.nmo
        nvir  = nmo - nocc

        ncas_vir  = ncore + ncas - nocc
        ncas_occ  = ncas - ncas_vir
        ncas_vir2 = ncas_vir * (ncas_vir-1) // 2
        ncas_occ2 = ncas_occ * (ncas_occ-1) // 2
        ncas_vir3 = ncas_vir * (ncas_vir-1) * (ncas_vir-2) // 6 
        ncas_occ3 = ncas_occ * (ncas_occ-1) * (ncas_occ-2) // 6 

        dS = ncas_occ * ncas_vir
        dD = ncas_occ2 * ncas_vir2
        dT = ncas_occ3 * ncas_vir3

        #reorder = self._casci.fcisolver.dmrg_args['reorder']
        reorder = [] 
        def get_cisdtq_vec_cas(vals, dets):
            c0     = np.zeros(1)
            c1a    = np.zeros((dS))
            c2ab   = np.zeros((dS,dS))
            c2aa   = np.zeros((dD))
            c3aaa  = np.zeros((dT))
            c3aab  = np.zeros((dD,dS))
            c4aaab = np.zeros((dT,dS))
            c4aabb = np.zeros((dD,dD))

            def convert_phase_to_Fermi_vacuum(h_a, h_b, val):
                n_perm = 0 
                ranka = len(h_a)
                if ranka > 0:
                    n_perm += ranka*ncas_occ - sum(h_a) - ranka*(ranka+1)//2
                rankb = len(h_b)
                if rankb > 0:
                    n_perm += rankb*ncas_occ - sum(h_b) - rankb*(rankb+1)//2
                return val * ( 1 - ( (n_perm & 1) << 1 ) )

            S = lambda i, a: ncas_occ * (a + 1) - (i + 1)
            D = lambda i, j, a, b: ncas_occ2 * (b * (b - 1) // 2 + a + 1) \
                                             - (j * (j - 1) // 2 + i + 1)
            T = lambda i, j, k, a, b, c: \
                    ncas_occ3 * (c * (c - 1) * (c - 2) // 6 + b * (b - 1) // 2 + a + 1) \
                              - (k * (k - 1) * (k - 2) // 6 + j * (j - 1) // 2 + i + 1)
    
            for val, det in zip(vals, dets):
                occa =  det & 1
                occb = (det & 2) >> 1 
                # hole, ptcl
                h_a, h_b, p_a, p_b = [], [], [], []
                if len(reorder) == 0:
                    for i in range(ncas_occ):
                        if occa[i] == 0: h_a.append(i)
                        if occb[i] == 0: h_b.append(i)
                    for a in range(ncas_occ,ncas,1):
                        if occa[a] == 1: p_a.append(a-ncas_occ)
                        if occb[a] == 1: p_b.append(a-ncas_occ)
                else:
                    for i in range(ncas):
                        org_i = reorder[i]
                        if occa[i] == 0 and org_i < ncas_occ:  h_a.append(org_i)
                        if occb[i] == 0 and org_i < ncas_occ:  h_b.append(org_i)
                        if occa[i] == 1 and org_i >= ncas_occ: p_a.append(org_i-ncas_occ)
                        if occb[i] == 1 and org_i >= ncas_occ: p_b.append(org_i-ncas_occ)
                    h_a.sort()
                    h_b.sort()
                    p_a.sort()
                    p_b.sort()

                assert len(h_a) == len(p_a) and len(h_b) == len(p_b)
                if   len(h_a) == 0 and len(h_b) == 0:
                    c0[0] = val 
                elif len(h_a) == 1 and len(h_b) == 0:
                    idx_a = S(*h_a, *p_a)
                    c1a[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 2 and len(h_b) == 0:
                    idx_a = D(*h_a, *p_a)
                    c2aa[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 1 and len(h_b) == 1:
                    idx_a = S(*h_a, *p_a)
                    idx_b = S(*h_b, *p_b)
                    c2ab[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 3 and len(h_b) == 0:
                    idx_a = T(*h_a, *p_a)
                    c3aaa[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 2 and len(h_b) == 1:
                    idx_a = D(*h_a, *p_a)
                    idx_b = S(*h_b, *p_b)
                    c3aab[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 3 and len(h_b) == 1:
                    idx_a = T(*h_a, *p_a)
                    idx_b = S(*h_b, *p_b)
                    c4aaab[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 2 and len(h_b) == 2:
                    idx_a = D(*h_a, *p_a)
                    idx_b = D(*h_b, *p_b)
                    c4aabb[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                else:
                    continue 

            assert np.abs(c0) > NUM_ZERO

            #print('c0')
            #print(c0)
            #print('c1a/c0')
            #print(c1a/c0)
            #print('c2ab/c0')
            #print(c2ab/c0)

            return c1a/c0, c2aa/c0, c2ab/c0, c3aaa/c0, c3aab/c0, \
                   c4aaab/c0, c4aabb/c0 

        import os
        scr = self._casci.fcisolver.scratchDirectory
        fvals = "%s/node0/sample-vals.npy"%(scr)
        fdets = "%s/node0/sample-dets.npy"%(scr)
        if os.path.isfile(fvals) and os.path.isfile(fdets):
            vals = np.load(fvals)
            dets = np.load(fdets)
        else:
            assert False

        return get_cisdtq_vec_cas(vals, dets)

class DMRGREXCCSD(DMRGEXCCSD):
    pass


