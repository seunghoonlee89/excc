import numpy as np

string_i = np.array(['i', 'j', 'k', 'l'])
string_a = np.array(['a', 'b', 'c', 'd'])

idxs_s = [[0]]
idxs_d = [[0,1], \
          [1,0]]
idxs_t = [[0,1,2], [1,2,0], [2,0,1], \
          [0,2,1], [2,1,0], [1,0,2]]
idxs_q = [[0,1,2,3], [1,2,3,0], [2,3,0,1], [3,0,1,2], \
          [0,2,3,1], [2,3,1,0], [3,1,0,2], [1,0,2,3], \
          [0,3,1,2], [3,1,2,0], [1,2,0,3], [2,0,3,1], \
          [0,1,3,2], [1,3,2,0], [3,2,0,1], [2,0,1,3], \
          [0,3,2,1], [3,2,1,0], [2,1,0,3], [1,0,3,2], \
          [0,2,1,3], [2,1,3,0], [1,3,0,2], [3,0,2,1]]
idxs = [[], idxs_s, idxs_d, idxs_t, idxs_q]

pars_d = [1, -1]
pars_s = [1]
pars_t = [1]*3 + [-1]*3
pars_q = [1]*12 + [-1]*12 
pars = [[], pars_s, pars_d, pars_t, pars_q]

#T: ((((i*nocc+(size_t)(j))*nocc+k)*nvir+a)*nvir+b)*nvir+c
#Q: (((((((int64_t)i*nocc+j)*nocc+k)*nocc+l)*nvir+a)*nvir+b)*nvir+c)*nvir+d


#T: ((((i*nocc+j)*nocc+k)*nvir+a)*nvir+b)*nvir+c
#Q: ((((((i*nocc+j)*nocc+k)*nocc+l)*nvir+a)*nvir+b)*nvir+c)*nvir+d

def idx_string(ia, aa, ib, ab):
    na = len(ia)
    nb = len(ib)
    ntot = na + nb
    l_idx = '(' * (2 * ntot - 1)
    if ntot == 4:
        l_idx += '(int64_t)'
    if na > 0:
        l_idx += ia[0]
        for it in range(1, na, 1):
            l_idx += '*nocca+' + ia[it] + ')'
        for it in range(nb):
            l_idx += '*noccb+' + ib[it] + ')'
    else:
        l_idx += ib[0]
        for it in range(1, nb, 1):
            l_idx += '*noccb+' + ib[it] + ')'
    for it in range(na):
        l_idx += '*nvira+' + aa[it] + ')'
    for it in range(nb):
        l_idx += '*nvirb+' + ab[it] + ')'
    return l_idx 

def code_gen(na, nb):
    assert na >= 0 and na < 5
    assert nb >= 0 and nb < 5
    ntot = na + nb
    assert ntot >= 0 and ntot < 5
    idxa, idxb = idxs[na], idxs[nb]
    para, parb = pars[na], pars[nb]
    ia, ib = string_i[:na], string_i[na:na+nb]
    aa, ab = string_a[:na], string_a[na:na+nb]

    header = 't%d' % ntot + 'a' * na + 'b' * nb + '['       
    if na > 0 and nb > 0:
        for it1, (idx1, par1) in enumerate(zip(idxa, para)):
            for it2, (idx2, par2) in enumerate(zip(idxa, para)):
                for it3, (idx3, par3) in enumerate(zip(idxb, parb)):
                    for it4, (idx4, par4) in enumerate(zip(idxb, parb)):
                        line = header
                        line += idx_string(ia[idx1], aa[idx2], ib[idx3], ab[idx4])
                        if par1 * par2 * par3 * par4 > 0:
                            line += '] =  tmp;'
                        else:
                            line += '] = -tmp;'
                        print(line)
    elif na > 0 and nb == 0:
        for it1, (idx1, par1) in enumerate(zip(idxa, para)):
            for it2, (idx2, par2) in enumerate(zip(idxa, para)):
                line = header
                line += idx_string(ia[idx1], aa[idx2], [], [])
                if par1 * par2 > 0:
                    line += '] =  tmp;'
                else:
                    line += '] = -tmp;'
                print(line)
    elif na == 0 and nb > 0:
        for it3, (idx3, par3) in enumerate(zip(idxb, parb)):
            for it4, (idx4, par4) in enumerate(zip(idxb, parb)):
                line = header
                line += idx_string([], [], ib[idx3], ab[idx4])
                if par3 * par4 > 0:
                    line += '] =  tmp;'
                else:
                    line += '] = -tmp;'
                print(line)


print('aa')
code_gen(2, 0) 
print('bb')
code_gen(0, 2) 
print('ab')
code_gen(1, 1) 

print('aaa')
code_gen(3, 0) 
print('bbb')
code_gen(0, 3) 
print('aab')
code_gen(2, 1) 
print('abb')
code_gen(1, 2) 

print('aaaa')
code_gen(4, 0) 
print('bbbb')
code_gen(0, 4) 
print('aaab')
code_gen(3, 1) 
print('abbb')
code_gen(1, 3) 
print('aabb')
code_gen(2, 2) 
