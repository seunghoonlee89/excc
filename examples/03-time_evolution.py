#!/usr/bin/env python

'''
A simple example to run externally corrected (ex)-CCSD calculation
using DMRG-solver and imaginary time evolution.

Author:
    Seunghoon Lee
'''

import pyscf
mol = pyscf.M(
    atom  = '''
    Be    0.0000000    0.0000000    1.
    Be    0.0000000    0.0000000   -1.
    ''',
    basis = 'sto-3g',
    spin  = 0
    )
mol.build()

myhf = mol.RHF().run()

from pyscf import dmrgscf, lib
import os
dmrgscf.settings.BLOCKEXE = os.popen("which block2main").read().strip()
dmrgscf.settings.MPIPREFIX = ''

mycas = myhf.CASCI(10, 8)
#mycas.canonicalize = False
mycas.fcisolver = dmrgscf.DMRGCI(mol, maxM=1000, tol=1E-15)
#mycas.fcisolver.twopdm = False
mycas.fcisolver.block_extra_keyword.append('noreorder')
mycas.fcisolver.block_extra_keyword.append('copy_mps ZKET')
mycas.fcisolver.block_extra_keyword.append('trans_mps_to_sz')
mycas.fcisolver.block_extra_keyword.append('sample 1e-15')
mycas.fcisolver.block_extra_keyword.append('sample_reference 4 3333000000')
mycas.fcisolver.block_extra_keyword.append('sample_phase 0 1 2 3 4 5 6 7 8 9')
#mycas.fcisolver.runtimeDir = lib.param.TMPDIR
#mycas.fcisolver.scratchDirectory = lib.param.TMPDIR
#mycas.fcisolver.threads = int(os.environ["OMP_NUM_THREADS"])
#mycas.fcisolver.memory = int(mol.max_memory / 1000) # mem in GB
mycas.kernel()
assert mycas.converged

import excc 
myexcc = excc.DMRGEXCCSD(mycas)
myexcc.verbose = 5
myexcc.conv_tol_normt = 1e-13
myexcc.conv_tol = 1e-15
myexcc.imag_tevol = True
myexcc.l = 1000.
myexcc.dl = 0.1
myexcc.kernel()

