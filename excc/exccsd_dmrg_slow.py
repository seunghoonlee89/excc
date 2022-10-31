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
        #assert isinstance(mc.fcisolver, dmrgscf.DMRGCI)
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
        assert mc.nelecas[0] == mc.nelecas[1]   # for restricted case
        self.ci      = None 
        self.sparse_tol = None
        self.ext_orbs = [] 

        # parameters for imaginary time evolution
        self.imag_tevol = False
        self.l   = 100. 
        self.dl  = 1. 
        self.fac = 1. 
        self.order = 4 
        self.thresh = 100. 
        self.real = True 
        self._keys = self._keys.union(["imag_tevol", "l", "dl", \
                                       "fac", "order", "thresh", "real"])

        # parameters for min res
        self.minres = False 
        self.method = "krylov" 
        self.precond = "finv" 
        self.inner_m = 10 
        self.outer_k = 6 
        self._keys = self._keys.union(["minres", "method", "precond", "inner_m", "outer_k"])

    def get_ext_c_amps(self):
        if self.sparse_tol is not None:
            tol = self.sparse_tol
        else:
            tol = 0. 

        if len(self.ext_orbs) != 0:
            ext_orbs = list(self.ext_orbs)
            ext_t = lambda p: p in ext_orbs 
            def check_ext(p_a, p_b, h_a, h_b):
                flag_ext = True
                for p in h_a:
                    flag_ext = flag_ext and ext_t(p)  
                for p in h_b:
                    flag_ext = flag_ext and ext_t(p)  
                for p in p_a:
                    flag_ext = flag_ext and ext_t(p)  
                for p in p_b:
                    flag_ext = flag_ext and ext_t(p)  
                return flag_ext
        else:
            ext_orbs = None 

        ncas  = self.ncas
        ncore = self.ncore
        nocc  = self.nocc

        ncas_vir  = ncore + ncas - nocc
        ncas_occ  = ncas - ncas_vir
        ncas_vir2 = ncas_vir * (ncas_vir-1) // 2
        ncas_occ2 = ncas_occ * (ncas_occ-1) // 2
        ncas_vir3 = ncas_vir * (ncas_vir-1) * (ncas_vir-2) // 6 
        ncas_occ3 = ncas_occ * (ncas_occ-1) * (ncas_occ-2) // 6 

        dS = ncas_occ * ncas_vir
        dD = ncas_occ2 * ncas_vir2
        dT = ncas_occ3 * ncas_vir3

        import os
        scr = self._casci.fcisolver.scratchDirectory
        freorder = "%s/node0/orbital_reorder.npy"%(scr)
        if os.path.isfile(freorder):
            reorder = np.load(freorder)
            if np.array_equal(reorder, np.arange(ncas)):
                reorder = [] 
        else:
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

##            #dbg
#            max_val = np.zeros((100))
#            max_val2 = np.zeros((100))

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
                if np.abs(val) < tol:
                    continue 
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
                    if ext_orbs is not None:
                        if not check_ext(p_a, p_b, h_a, h_b):
                            continue
                    idx_a = T(*h_a, *p_a)
                    c3aaa[idx_a] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 2 and len(h_b) == 1:
                    if ext_orbs is not None:
                        if not check_ext(p_a, p_b, h_a, h_b):
                            continue
                    idx_a = D(*h_a, *p_a)
                    idx_b = S(*h_b, *p_b)
                    c3aab[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
                elif len(h_a) == 3 and len(h_b) == 1:
                    if ext_orbs is not None:
                        if not check_ext(p_a, p_b, h_a, h_b):
                            continue
                    idx_a = T(*h_a, *p_a)
                    idx_b = S(*h_b, *p_b)
                    c4aaab[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)

##                    #dbg
#                    if np.abs(val) > max_val[0]:
#                        max_val[0] = np.abs(val)
#                        max_val.sort()

                elif len(h_a) == 2 and len(h_b) == 2:
                    if ext_orbs is not None:
                        if not check_ext(p_a, p_b, h_a, h_b):
                            continue
                    idx_a = D(*h_a, *p_a)
                    idx_b = D(*h_b, *p_b)
                    c4aabb[idx_a, idx_b] = convert_phase_to_Fermi_vacuum(h_a, h_b, val)
##                    #dbg
#                    if np.abs(val) > max_val2[0]:
#                        max_val2[0] = np.abs(val)
#                        max_val2.sort()
                else:
                    continue 

#            print('c4aaab')
#            for v in max_val[::-1]:
#                print(v)
#
#            print('c4aabb')
#            for v in max_val2[::-1]:
#                print(v)
#            #exit() 


            assert np.abs(c0) > NUM_ZERO
            #print('c0')
            #print(c0)
            #print('c1a/c0')
            #print(c1a/c0)
            #print('c2ab/c0')
            #print(c2ab/c0)

            return c1a/c0, c2aa/c0, c2ab/c0, c3aaa/c0, c3aab/c0, \
                   c4aaab/c0, c4aabb/c0 

        fvals = "%s/node0/sample-vals.npy"%(scr)
        fdets = "%s/node0/sample-dets.npy"%(scr)
        if os.path.isfile(fvals) and os.path.isfile(fdets):
            vals = np.load(fvals)
            dets = np.load(fdets)
        else:
            assert False

        return get_cisdtq_vec_cas(vals, dets)

    def exclude(self, r, typ='aaa'):
        assert typ == 'aaa' or typ == 'baa'

        if self.sparse_tol is not None:
            tol = self.sparse_tol
        else:
            tol = 0. 

        ncas  = self.ncas
        ncore = self.ncore
        nocc  = self.nocc
        ncas_vir  = ncore + ncas - nocc
        ncas_occ  = ncas - ncas_vir

        import os
        scr = self._casci.fcisolver.scratchDirectory
        freorder = "%s/node0/orbital_reorder.npy"%(scr)
        if os.path.isfile(freorder):
            reorder = np.load(freorder)
            if np.array_equal(reorder, np.arange(ncas)):
                reorder = [] 
        else:
            reorder = [] 

        fvals = "%s/node0/sample-vals.npy"%(scr)
        fdets = "%s/node0/sample-dets.npy"%(scr)
        if os.path.isfile(fvals) and os.path.isfile(fdets):
            vals = np.load(fvals)
            dets = np.load(fdets)
        else:
            assert False

        ncount = 0
        for val, det in zip(vals, dets):
            if np.abs(val) < tol:
                continue 
            occa =  det & 1
            occb = (det & 2) >> 1 
            # hole, ptcl
            h_a, h_b, p_a, p_b = [], [], [], []
            if len(reorder) == 0:
                for i in range(ncas_occ):
                    if occa[i] == 0: h_a.append(i + ncore)
                    if occb[i] == 0: h_b.append(i + ncore)
                for a in range(ncas_occ, ncas, 1):
                    if occa[a] == 1: p_a.append(a - ncas_occ)
                    if occb[a] == 1: p_b.append(a - ncas_occ)
            else:
                for i in range(ncas):
                    org_i = reorder[i]
                    if occa[i] == 0 and org_i < ncas_occ:  h_a.append(org_i + ncore)
                    if occb[i] == 0 and org_i < ncas_occ:  h_b.append(org_i + ncore)
                    if occa[i] == 1 and org_i >= ncas_occ: p_a.append(org_i - ncas_occ)
                    if occb[i] == 1 and org_i >= ncas_occ: p_b.append(org_i - ncas_occ)
                h_a.sort()
                h_b.sort()
                p_a.sort()
                p_b.sort()

            assert len(h_a) == len(p_a) and len(h_b) == len(p_b)

            if typ == 'aaa' and len(h_a) == 3 and len(h_b) == 0:
                ncount += 1 
                i, j, k = h_a
                a, b, c = p_a 
                r[i, j, k, a, b, c] = 0. 
                r[i, k, j, a, b, c] = 0. 
                r[j, i, k, a, b, c] = 0. 
                r[j, k, i, a, b, c] = 0. 
                r[k, j, i, a, b, c] = 0. 
                r[k, i, j, a, b, c] = 0. 

                r[i, j, k, a, c, b] = 0. 
                r[i, k, j, a, c, b] = 0. 
                r[j, i, k, a, c, b] = 0. 
                r[j, k, i, a, c, b] = 0. 
                r[k, j, i, a, c, b] = 0. 
                r[k, i, j, a, c, b] = 0. 

                r[i, j, k, b, a, c] = 0. 
                r[i, k, j, b, a, c] = 0. 
                r[j, i, k, b, a, c] = 0. 
                r[j, k, i, b, a, c] = 0. 
                r[k, j, i, b, a, c] = 0. 
                r[k, i, j, b, a, c] = 0. 

                r[i, j, k, b, c, a] = 0. 
                r[i, k, j, b, c, a] = 0. 
                r[j, i, k, b, c, a] = 0. 
                r[j, k, i, b, c, a] = 0. 
                r[k, j, i, b, c, a] = 0. 
                r[k, i, j, b, c, a] = 0. 

                r[i, j, k, c, a, b] = 0. 
                r[i, k, j, c, a, b] = 0. 
                r[j, i, k, c, a, b] = 0. 
                r[j, k, i, c, a, b] = 0. 
                r[k, j, i, c, a, b] = 0. 
                r[k, i, j, c, a, b] = 0. 

                r[i, j, k, c, b, a] = 0. 
                r[i, k, j, c, b, a] = 0. 
                r[j, i, k, c, b, a] = 0. 
                r[j, k, i, c, b, a] = 0. 
                r[k, j, i, c, b, a] = 0. 
                r[k, i, j, c, b, a] = 0. 
            elif typ == 'baa' and len(h_a) == 2 and len(h_b) == 1:
                ncount += 1 
                i = h_b
                a = p_b 
                j, k = h_a
                b, c = p_a 

                r[i, j, k, a, b, c] = 0.
                r[i, k, j, a, b, c] = 0.

                r[i, j, k, a, c, b] = 0.
                r[i, k, j, a, c, b] = 0.
            else:
                continue 
        print(typ, '= ', ncount)

    def update_mask(self, r, typ='aaa'):
        assert typ == 'aaa' or typ == 'baa'
        r.fill(False)

        if self.sparse_tol is not None:
            tol = self.sparse_tol
        else:
            tol = 0. 

        ncas  = self.ncas
        ncore = self.ncore
        nocc  = self.nocc
        ncas_vir  = ncore + ncas - nocc
        ncas_occ  = ncas - ncas_vir

        import os
        scr = self._casci.fcisolver.scratchDirectory
        freorder = "%s/node0/orbital_reorder.npy"%(scr)
        if os.path.isfile(freorder):
            reorder = np.load(freorder)
            if np.array_equal(reorder, np.arange(ncas)):
                reorder = [] 
        else:
            reorder = [] 

        fvals = "%s/node0/sample-vals.npy"%(scr)
        fdets = "%s/node0/sample-dets.npy"%(scr)
        if os.path.isfile(fvals) and os.path.isfile(fdets):
            vals = np.load(fvals)
            dets = np.load(fdets)
        else:
            assert False

        ncount = 0
        for val, det in zip(vals, dets):
            if np.abs(val) < tol:
                continue 
            occa =  det & 1
            occb = (det & 2) >> 1 
            # hole, ptcl
            h_a, h_b, p_a, p_b = [], [], [], []
            if len(reorder) == 0:
                for i in range(ncas_occ):
                    if occa[i] == 0: h_a.append(i)
                    if occb[i] == 0: h_b.append(i)
                for a in range(ncas_occ, ncas, 1):
                    if occa[a] == 1: p_a.append(a - ncas_occ)
                    if occb[a] == 1: p_b.append(a - ncas_occ)
            else:
                for i in range(ncas):
                    org_i = reorder[i]
                    if occa[i] == 0 and org_i < ncas_occ:  h_a.append(org_i)
                    if occb[i] == 0 and org_i < ncas_occ:  h_b.append(org_i)
                    if occa[i] == 1 and org_i >= ncas_occ: p_a.append(org_i - ncas_occ)
                    if occb[i] == 1 and org_i >= ncas_occ: p_b.append(org_i - ncas_occ)
                h_a.sort()
                h_b.sort()
                p_a.sort()
                p_b.sort()

            assert len(h_a) == len(p_a) and len(h_b) == len(p_b)

            if typ == 'aaa' and len(h_a) == 3 and len(h_b) == 0:
                ncount += 1 
                i, j, k = h_a
                a, b, c = p_a 
                r[i, j, k, a, b, c] = True 
                r[i, k, j, a, b, c] = True 
                r[j, i, k, a, b, c] = True 
                r[j, k, i, a, b, c] = True 
                r[k, j, i, a, b, c] = True 
                r[k, i, j, a, b, c] = True 

                r[i, j, k, a, c, b] = True 
                r[i, k, j, a, c, b] = True 
                r[j, i, k, a, c, b] = True 
                r[j, k, i, a, c, b] = True 
                r[k, j, i, a, c, b] = True 
                r[k, i, j, a, c, b] = True 

                r[i, j, k, b, a, c] = True 
                r[i, k, j, b, a, c] = True 
                r[j, i, k, b, a, c] = True 
                r[j, k, i, b, a, c] = True 
                r[k, j, i, b, a, c] = True 
                r[k, i, j, b, a, c] = True 

                r[i, j, k, b, c, a] = True 
                r[i, k, j, b, c, a] = True 
                r[j, i, k, b, c, a] = True 
                r[j, k, i, b, c, a] = True 
                r[k, j, i, b, c, a] = True 
                r[k, i, j, b, c, a] = True 

                r[i, j, k, c, a, b] = True 
                r[i, k, j, c, a, b] = True 
                r[j, i, k, c, a, b] = True 
                r[j, k, i, c, a, b] = True 
                r[k, j, i, c, a, b] = True 
                r[k, i, j, c, a, b] = True 

                r[i, j, k, c, b, a] = True 
                r[i, k, j, c, b, a] = True 
                r[j, i, k, c, b, a] = True 
                r[j, k, i, c, b, a] = True 
                r[k, j, i, c, b, a] = True 
                r[k, i, j, c, b, a] = True 
            elif typ == 'baa' and len(h_a) == 2 and len(h_b) == 1:
                ncount += 1 
                i = h_b
                a = p_b 
                j, k = h_a
                b, c = p_a 

                r[i, j, k, a, b, c] = True
                r[i, k, j, a, b, c] = True

                r[i, j, k, a, c, b] = True
                r[i, k, j, a, c, b] = True
            else:
                continue 
        print(typ, '= ', ncount)

class DMRGEXRCCSD(DMRGEXCCSD):
    pass

