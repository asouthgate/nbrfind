import ctypes
import numpy as np
from ctypes import cdll
from ctypes import c_int,  c_void_p, c_float, POINTER, c_char_p
from ctypes import sizeof
lib = cdll.LoadLibrary('./libuk2wrap.so')

lib.calculate_dist_resume_freeze.restype  = c_void_p
lib.calculate_dist_resume_freeze.argtypes = [c_char_p, c_char_p, POINTER(c_int), POINTER(c_int), c_int]
lib.calculate_dist_resume_nofreeze.restype  = c_void_p
lib.calculate_dist_resume_nofreeze.argtypes = [c_char_p, c_char_p, POINTER(c_int), POINTER(c_int), c_int]



def uk2_cpp(s1, s2, state_triple=None, state_arr=None, freeze=False):
    if state_triple is None:
        state_triple = np.zeros(3,dtype=np.int32)
        state_triple[1] = -len(s1)
        state_triple[2] = len(s2)
        state_arr = np.zeros((12,len(s1)+len(s2)+1),dtype=np.int32)      
        state_arr[0,:] = -2
        state_arr[1,:] = -1
        state_arr[2,:] = -1   
    prev_rowsize = state_arr.shape[1]
    rowsize = len(s1) + len(s2) + 1
    z = rowsize-prev_rowsize
    if z > 0:
        print("RECEIVED", state_arr[:3])
        state_arr = np.append(state_arr, (np.ones((12, z))*-1).astype(np.int32), axis=1)
        if state_triple[0] == 0:
            state_arr[0,-z:] = -2
        state_arr[3:,-z:] = 0
        state_triple[2] = len(s2)
        print("INPUTTING", state_arr[:3])
    state_triplep = state_triple.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
    state_arrp = state_arr.ctypes.data_as(ctypes.POINTER(ctypes.c_int32))
    if freeze: 
        lib.calculate_dist_resume_freeze(bytes(s1,encoding='utf8'), bytes(s2,encoding='utf8'), state_triplep, state_arrp, rowsize)
        print(state_arr)
        return (state_triple[0], state_arr[5,len(s2)], state_arr[11,len(s2)], state_arr[8,len(s2)], state_triple, state_arr)
    else: 
        lib.calculate_dist_resume_nofreeze(bytes(s1,encoding='utf8'), bytes(s2,encoding='utf8'), state_triplep, state_arrp, rowsize)
        print(state_arr)
        return (state_triple[0], state_arr[5,len(s2)], state_arr[11,len(s2)], state_arr[8,len(s2)])


def get_backtrace_stats(a1,a2):
    # HERE MN IS POSITIONS WITH N
    # M IS MATCH OR MISMATCH POSITIONS WITHOUT N
    snpd = 0
    MN = 0
    M = 0
    for j in range(len(a1)):
        c1 = a1[j]
        c2 = a2[j]
        if "-" not in [c1,c2]:
            # not an indel
            if "N" in [c1,c2]:
                MN += 1
            else:
                M += 1
                if c1 != c2: snpd += 1
    return snpd, MN, M


if __name__ == "__main__":
    import random
    import ukkonen2 as uk2val
    import ukkonen2_withmarrs as uk2valmarr
    import sys
    MINL = 3
    MAXL = 4
    for j in range(100):
        print()
        l1 = random.randint(MINL,MAXL)
        l2 = random.randint(MINL,MAXL)
        s1 = "".join([random.choice("ACGTN") for c in range(l1)])
        s2 = "".join([random.choice("ACGTN") for c in range(l2)])
        print(s1)
        print(s2)
        snpdval, hval, Ltrace = uk2val.ukkonen_lev2(s1, s2, freeze=False, trace=True)
        snpdval2, hval2, NN2, M2TOT = uk2valmarr.ukkonen_lev2(s1, s2, freeze=False)
        M2 = M2TOT-NN2
        a1, a2 = uk2val.backtrace(Ltrace, s1, s2)
        print(a1)
        print(a2)
        d, M_N, M = get_backtrace_stats(a1,a2)
        h, snpd, mn, mtot = uk2_cpp(s1, s2)
        #**** Also check freezing functionality
        randk = random.randint(1, len(s2)-2)
        print(s1,s2,randk, s2[:randk])
        hfreeze, snpdfreeze, mnfreeze, mtotfreeze, st, sa = uk2_cpp(s1, s2[:randk], freeze=True)
        hfreeze, snpdfreeze, mnfreeze, mtotfreeze = uk2_cpp(s1, s2, state_triple=st, state_arr = sa, freeze=False)
        assert hfreeze == h, (hfreeze, h)
        assert snpd == snpdfreeze
        assert mn == mnfreeze, (mn, mnfreeze)
        assert mtot == mtotfreeze
        #****        
        m = mtot-mn
        print(h,snpd,mtot,m,mn,m)
        sys.stderr.write("testing %d\n" % j)
        assert hval == h, (hval, h)
        assert snpdval2 == snpd
        assert snpdval == d, (snpdval, d)
        assert snpdval == snpd, (snpdval, snpd)
#        assert M == M2, (M, M2)
        assert m == M2, (m, M2)
        assert mn == NN2, (mn, NN2) 
#        assert m == M, (m, M)
        assert mn == M_N, (mn, M_N)
    sys.stderr.write("passed basic test..!\n")
