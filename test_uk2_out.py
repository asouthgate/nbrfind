from Bio import SeqIO
from Bio import pairwise2
from libuk import uk2_cpp
from ukkonen2 import ukkonen_lev2
import sys

msarecs = {r.description:str(r.seq) for r in SeqIO.parse(sys.argv[1], format="fasta")}

def get_aln_stats(a1,a2):
    # HERE MN IS POSITIONS WITH N
    # M IS MATCH OR MISMATCH POSITIONS WITHOUT N
    snpd = 0
    NM = 0
    NN = 0
    for j in range(len(a1)):
        c1 = a1[j]
        c2 = a2[j]
        if "-" not in [c1,c2]:
            # not an indel
            if "N" in [c1,c2]:
                NN += 1
            else:
                NM += 1
                if c1 != c2: snpd += 1
    return snpd, NM, NN

with open(sys.argv[2]) as f:
    for l in f:
        first = l.split(",")[0][1:]
        second = l.split(",")[1][1:]
        mseq1 = msarecs[first].upper()
        mseq2 = msarecs[second].upper()
        h, d, NM, NN = [int(i) for i in l.split(",")[2:]]
        if h > -1:
            valres = ukkonen_lev2(mseq1.replace("-",""),mseq2.replace("-",""))
            evedit, evsnpd = valres[1], valres[0]
            snpdval, NMval, NNval = get_aln_stats(mseq1, mseq2)
            print([h, d], [snpdval], [evedit,evsnpd])
            assert evedit == h
            assert snpdval == d
        else:
            valres = ukkonen_lev2(mseq1.replace("-",""),mseq2.replace("-",""))
            evedit, evsnpd = valres[1], valres[0]
            snpdval, NMval, NNval = get_aln_stats(mseq1, mseq2)
            print([h, d], [snpdval], [evedit,evsnpd])
            assert snpdval >= -d

            
