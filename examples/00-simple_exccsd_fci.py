#!/usr/bin/env python

'''
A simple example to run externally corrected (ex)-CCSD calculation
using FCI-solver.

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
    spin  = 0,
    verbose = 5
    )
mol.build()

myhf = mol.RHF().run()

mycas = myhf.CASCI(8, 8)
mycas.kernel()
assert mycas.converged

import excc 
# conventional CC solver
myexcc = excc.EXCCSD(mycas)
myexcc.kernel()
e_excc = myexcc.e_tot

# imaginary time evolution solver
myexcc2 = excc.EXCCSD(mycas)
myexcc2.imag_tevol = True 
myexcc2.dl = 0.01
myexcc2.kernel()
e_excc2 = myexcc2.e_tot

# scipy optimize solver
myexcc3 = excc.EXCCSD(mycas)
myexcc3.minres = True
myexcc3.kernel()
e_excc3 = myexcc3.e_tot

print(e_excc)
print(e_excc2)
print(e_excc3)

