#!/usr/bin/env python

'''
A simple example to run externally corrected (ex)-CCSD calculation
using block2 DMRG-solver.

Author:
    Seunghoon Lee
'''

import pyscf
mol = pyscf.M(
    atom  = '''
    N    0.0000000    0.0000000    0.5600041
    N    0.0000000    0.0000000   -0.5600041
    ''',
    basis = 'cc-pvtz',
    spin  = 0
    )
mol.build()

myhf = mol.RHF().run()

from pyscf import dmrgscf, lib
import os
scratch = './nodex'
os.system('mkdir -p %s/node0' % scratch)
dmrgscf.settings.BLOCKEXE = os.popen("which block2main").read().strip()
dmrgscf.settings.BLOCKSCRATCHDIR = scratch
dmrgscf.settings.MPIPREFIX = ''

ncas = 6
na = 3
nb = 3
nelec_cas = (na, nb)
import random
import numpy as np
reorder = list(np.arange(ncas) + 1)
random.shuffle(reorder)
freorder = "%s/node0/orbital_reorder.txt" % (scratch)

phase_b0 = ''
phase_b1 = ''
for idx in reorder:
    phase_b0 += '%d ' % (idx - 1)
    phase_b1 += '%d ' % idx
f = open(freorder, 'w')
f.write(phase_b1)
f.close()

print("Reordering indices: ", reorder) 

max_rank = 4  # extract 0, S, D, T, Q amplitudes from MPS 
vacuum = '3' * nb + '1' * (na - nb) + '0'  * (ncas - na)
vacuum_str = ""
for i in range(ncas):
    vacuum_str += vacuum[reorder[i] - 1]

mycas = myhf.CASCI(ncas, nelec_cas)
mycas.fcisolver = dmrgscf.DMRGCI(mol, maxM=1000, tol=1E-7)
#mycas.fcisolver.twopdm = False
mycas.fcisolver.block_extra_keyword.append('reorder %s' % freorder)
mycas.fcisolver.block_extra_keyword.append('copy_mps ZKET')
mycas.fcisolver.block_extra_keyword.append('trans_mps_to_sz')
mycas.fcisolver.block_extra_keyword.append('sample 1e-9')
mycas.fcisolver.block_extra_keyword.append('sample_reference %d %s' % (max_rank, vacuum_str))
mycas.fcisolver.block_extra_keyword.append('sample_phase %s' % phase_b0)
#mycas.fcisolver.runtimeDir = lib.param.TMPDIR
#mycas.fcisolver.scratchDirectory = lib.param.TMPDIR
#mycas.fcisolver.threads = int(os.environ["OMP_NUM_THREADS"])
#mycas.fcisolver.memory = int(mol.max_memory / 1000) # mem in GB
mycas.kernel()
assert mycas.converged

import excc 
myexcc = excc.DMRGEXCCSD(mycas).kernel()

