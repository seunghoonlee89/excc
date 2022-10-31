import ctypes
import numpy as np
from .lib import load_library

libcc = load_library('libcc')

NUM_ZERO = 1e-9

def _c1_to_t1(c1, ncas_occ, ncas_vir):
    t1 = np.zeros((ncas_occ, ncas_vir), dtype=np.float64)
    libcc.c1_to_t1(
        t1.ctypes.data_as(ctypes.c_void_p),
        c1.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(ncas_occ), ctypes.c_int(ncas_vir)
        )
    return t1

def _c2_to_t2(t1, c2aa, c2ab, ncas_occ, ncas_vir, numzero=NUM_ZERO):
    t2aa = np.zeros((ncas_occ,ncas_occ,ncas_vir,ncas_vir),
                        dtype=np.float64)
    t2ab = np.zeros((ncas_occ,ncas_occ,ncas_vir,ncas_vir),
                        dtype=np.float64)
    libcc.c2_to_t2( 
        t2aa.ctypes.data_as(ctypes.c_void_p),
        t2ab.ctypes.data_as(ctypes.c_void_p),
        c2aa.ctypes.data_as(ctypes.c_void_p),
        c2ab.ctypes.data_as(ctypes.c_void_p),
        t1.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(ncas_occ), ctypes.c_int(ncas_vir),
        ctypes.c_double(numzero)
        )
    return t2aa, t2ab

def _c3_to_t3(t1, t2aa, t2ab, c3aaa, c3aab, ncas_occ, ncas_vir, numzero=NUM_ZERO):
    no, nv = ncas_occ, ncas_vir
    t3aaa = np.zeros((no,no,no,nv,nv,nv), dtype=np.float64)
    t3aab = np.zeros((no,no,no,nv,nv,nv), dtype=np.float64)
    libcc.c3_to_t3(
        t3aaa.ctypes.data_as(ctypes.c_void_p),
        t3aab.ctypes.data_as(ctypes.c_void_p),
        c3aaa.ctypes.data_as(ctypes.c_void_p),
        c3aab.ctypes.data_as(ctypes.c_void_p),
        t1.ctypes.data_as(ctypes.c_void_p),
        t2aa.ctypes.data_as(ctypes.c_void_p),
        t2ab.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(no), ctypes.c_int(nv),
        ctypes.c_double(numzero))
    return t3aaa, t3aab

def _c4_to_t4(t1, t2aa, t2ab, t3aaa, t3aab, c4aaab, c4aabb,
              ncas_occ, ncas_vir, numzero=NUM_ZERO):
    no, nv = ncas_occ, ncas_vir
    t4aaab = np.zeros((no,no,no,no,nv,nv,nv,nv), dtype=np.float64)
    t4aabb = np.zeros((no,no,no,no,nv,nv,nv,nv), dtype=np.float64)
    libcc.c4_to_t4(
        t4aaab.ctypes.data_as(ctypes.c_void_p),
        t4aabb.ctypes.data_as(ctypes.c_void_p),
        c4aaab.ctypes.data_as(ctypes.c_void_p),
        c4aabb.ctypes.data_as(ctypes.c_void_p),
        t1.ctypes.data_as(ctypes.c_void_p),
        t2aa.ctypes.data_as(ctypes.c_void_p),
        t2ab.ctypes.data_as(ctypes.c_void_p),
        t3aaa.ctypes.data_as(ctypes.c_void_p),
        t3aab.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(no),ctypes.c_int(nv),
        ctypes.c_double(numzero))
    return t4aaab, t4aabb

def _c1_to_t1u(c1, ncas_occ, ncas_vir):
    c1a, c1b = c1
    noa, nob = ncas_occ
    nva, nvb = ncas_vir
    t1a = np.zeros((noa, nva), dtype=np.float64)
    t1b = np.zeros((nob, nvb), dtype=np.float64)
    libcc.c1_to_t1(
        t1a.ctypes.data_as(ctypes.c_void_p),
        c1a.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(noa), ctypes.c_int(nva)
        )
    libcc.c1_to_t1(
        t1b.ctypes.data_as(ctypes.c_void_p),
        c1b.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(nob), ctypes.c_int(nvb)
        )
    return t1a, t1b

def _c2_to_t2u(t1, c2, ncas_occ, ncas_vir, numzero=NUM_ZERO):
    t1a, t1b = t1
    c2aa, c2ab, c2bb = c2
    noa, nob = ncas_occ
    nva, nvb = ncas_vir
    t2aa = np.zeros((noa,noa,nva,nva),dtype=np.float64)
    t2ab = np.zeros((noa,nob,nva,nvb),dtype=np.float64)
    t2bb = np.zeros((nob,nob,nvb,nvb),dtype=np.float64)
    libcc.c2_to_t2u( 
        t2aa.ctypes.data_as(ctypes.c_void_p),
        t2ab.ctypes.data_as(ctypes.c_void_p),
        t2bb.ctypes.data_as(ctypes.c_void_p),
        c2aa.ctypes.data_as(ctypes.c_void_p),
        c2ab.ctypes.data_as(ctypes.c_void_p),
        c2bb.ctypes.data_as(ctypes.c_void_p),
        t1a.ctypes.data_as(ctypes.c_void_p),
        t1b.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(noa),ctypes.c_int(nva),
        ctypes.c_int(nob),ctypes.c_int(nvb),
        ctypes.c_double(numzero)
        )
    return t2aa, t2ab, t2bb

def _c3_to_t3u(t1, t2, c3, ncas_occ, ncas_vir, numzero=NUM_ZERO):
    c3aaa, c3aab, c3abb, c3bbb = c3
    t1a, t1b = t1
    t2aa, t2ab, t2bb = t2
    noa, nob = ncas_occ
    nva, nvb = ncas_vir
    t3aaa = np.zeros((noa,noa,noa,nva,nva,nva), dtype=np.float64)
    t3aab = np.zeros((noa,noa,nob,nva,nva,nvb), dtype=np.float64)
    t3abb = np.zeros((noa,nob,nob,nva,nvb,nvb), dtype=np.float64)
    t3bbb = np.zeros((nob,nob,nob,nvb,nvb,nvb), dtype=np.float64)
    libcc.c3_to_t3u(
        t3aaa.ctypes.data_as(ctypes.c_void_p),
        t3aab.ctypes.data_as(ctypes.c_void_p),
        t3abb.ctypes.data_as(ctypes.c_void_p),
        t3bbb.ctypes.data_as(ctypes.c_void_p),
        c3aaa.ctypes.data_as(ctypes.c_void_p),
        c3aab.ctypes.data_as(ctypes.c_void_p),
        c3abb.ctypes.data_as(ctypes.c_void_p),
        c3bbb.ctypes.data_as(ctypes.c_void_p),
        t1a.ctypes.data_as(ctypes.c_void_p),
        t1b.ctypes.data_as(ctypes.c_void_p),
        t2aa.ctypes.data_as(ctypes.c_void_p),
        t2ab.ctypes.data_as(ctypes.c_void_p),
        t2bb.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(noa),ctypes.c_int(nva),
        ctypes.c_int(nob),ctypes.c_int(nvb),
        ctypes.c_double(numzero))
    return t3aaa, t3aab, t3abb, t3bbb

def _c4_to_t4u(t1, t2, t3, c4, ncas_occ, ncas_vir, numzero=NUM_ZERO):
    c4aaaa, c4aaab, c4aabb, c4abbb, c4bbbb = c4
    t1a, t1b = t1
    t2aa, t2ab, t2bb = t2
    t3aaa, t3aab, t3abb, t3bbb = t3
    noa, nob = ncas_occ
    nva, nvb = ncas_vir
    t4aaaa = np.zeros((noa,noa,noa,noa,nva,nva,nva,nva), dtype=np.float64)
    t4aaab = np.zeros((noa,noa,noa,nob,nva,nva,nva,nvb), dtype=np.float64)
    t4aabb = np.zeros((noa,noa,nob,nob,nva,nva,nvb,nvb), dtype=np.float64)
    t4abbb = np.zeros((noa,nob,nob,nob,nva,nvb,nvb,nvb), dtype=np.float64)
    t4bbbb = np.zeros((nob,nob,nob,nob,nvb,nvb,nvb,nvb), dtype=np.float64)
    libcc.c4_to_t4u(
        t4aaaa.ctypes.data_as(ctypes.c_void_p),
        t4aaab.ctypes.data_as(ctypes.c_void_p),
        t4aabb.ctypes.data_as(ctypes.c_void_p),
        t4abbb.ctypes.data_as(ctypes.c_void_p),
        t4bbbb.ctypes.data_as(ctypes.c_void_p),
        c4aaaa.ctypes.data_as(ctypes.c_void_p),
        c4aaab.ctypes.data_as(ctypes.c_void_p),
        c4aabb.ctypes.data_as(ctypes.c_void_p),
        c4abbb.ctypes.data_as(ctypes.c_void_p),
        c4bbbb.ctypes.data_as(ctypes.c_void_p),
        t1a.ctypes.data_as(ctypes.c_void_p),
        t1b.ctypes.data_as(ctypes.c_void_p),
        t2aa.ctypes.data_as(ctypes.c_void_p),
        t2ab.ctypes.data_as(ctypes.c_void_p),
        t2bb.ctypes.data_as(ctypes.c_void_p),
        t3aaa.ctypes.data_as(ctypes.c_void_p),
        t3aab.ctypes.data_as(ctypes.c_void_p),
        t3abb.ctypes.data_as(ctypes.c_void_p),
        t3bbb.ctypes.data_as(ctypes.c_void_p),
        ctypes.c_int(noa),ctypes.c_int(nva),
        ctypes.c_int(nob),ctypes.c_int(nvb),
        ctypes.c_double(numzero))
    return t4aaaa, t4aaab, t4aabb, t4abbb, t4bbbb
                    
