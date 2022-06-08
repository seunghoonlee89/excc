#!/usr/bin/env python
#
# Author: Seunghoon Lee <seunghoonlee89@gmail.com>
#

'''
A simple example to run externally corrected (ex)-CCSD calculation.
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

mycas = myhf.CASCI(6, 6)
mycas.kernel()
assert mycas.converged

import excc 
myexcc = excc.EXCCSD(mycas).kernel()

