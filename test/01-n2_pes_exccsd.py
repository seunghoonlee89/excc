#!/usr/bin/env python
#
# Author: Seunghoon Lee <seunghoonlee89@gmail.com>
#

'''
A simple example to run TCCSD and TCCSD(T) calculation.
'''

import pyscf
import excc 

def n2(atom, damp):
    mol = pyscf.M(
        atom = atom,
        basis = 'cc-pvtz',
        spin=0)
    mf = mol.RHF().run()
    
    #################################
    # variational HCI wave function #
    #################################
    no_cas = 6 
    ne_cas = 6 
    from pyscf import mcscf, fci
    mc = mcscf.CASCI(mf, no_cas, ne_cas)
    mc.fcisolver = fci.direct_spin1.FCISolver(mol) 
    mc.fcisolver.conv_tol = 1e-15
    mc.fix_spin_(shift=0.05, ss=0.0)
    mc.fcisolver.nroots = 5 
    mc.fcisolver.spin = 0 
    mc.casci()

    def find_n_spin_gs_state(mc, twos):
        nroots = len(mc.ci)
        norb = mc.ncas
        nelec= mc.nelecas

        e_tot_twos = [] 
        count = 0
        for i in range(nroots):
            twos_i = fci.spin_op.spin_square0(mc.ci[i], norb, nelec)[1]-1
            if twos > twos_i-0.5 and twos < twos_i+0.5:
                return mc.ci[i]
                break 
    mc.ci = find_n_spin_gs_state(mc, 0)
   
    #################
    # tailored CCSD #
    #################
    from pyscf import cc 
    myexcc = excc.EXCCSD(mc)
    myexcc.verbose          = 5 
    myexcc.max_memory       = 10000  # 10 g
    myexcc.max_cycle        = 1000
    myexcc.conv_tol         = 1e-6
    myexcc.diis             = False
    myexcc.level_shift      = 0.3
    myexcc.iterative_damping= damp 
    myexcc.kernel()
    E_EXCCSD_FCI = myexcc.e_tot
    #et = myexcc.ccsd_t()
    #E_EXCCSD_t_FCI = E_TCCSD_FCI+et

    return E_EXCCSD_FCI

if __name__ == '__main__':
    bohr2ang = 0.529177249
    bonds= [1.5 + 0.2*i for i in range(40)] 
    damp = [0.5]*10 + [0.1]*10 + [0.05]*20
    atom = ['''
             N     0.000000   0.000000    0.0 
             N     0.000000   0.000000    %10.6f 
            '''%(b * bohr2ang) for b in bonds]

    Eexcc = []
    for i in range(len(bonds)):
        print(bonds[i])
        Eexcc.append(n2(atom[i], damp[i]))

    for b, e in zip(bonds, Eexcc):
        print(b, e)

