#!/usr/bin/env python

'''
A simple example to run externally corrected (ex)-CCSD calculation
using HCI-solver.

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

from pyscf.shciscf import shci
mycas = myhf.CASCI(6, 6)
mycas.fcisolver = shci.SHCI(mol) 
mycas.fcisolver.extraline = ['printbestdeterminants 400']
mycas.fcisolver.nPTiter = 0 
mycas.fcisolver.DoRDM = False 
mycas.kernel()
assert mycas.converged

import excc 
myexcc = excc.HCIEXCCSD(mycas).kernel()

mycas.fcisolver.cleanup_dice_files()

