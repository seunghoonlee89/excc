import ctypes
import numpy as np
from .lib import load_library

libcc = load_library('libcc')

NUM_ZERO = 1e-16

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
                    
