import sys
import subprocess as sp
from Bio import SeqIO
from Bio import pairwise2
import sys

queries = [r for r in SeqIO.parse(sys.argv[1], "fasta")]
recs = [r for r in SeqIO.parse(sys.argv[2], "fasta")]




com = "./build/src/nbrfind %s %s 100" % (sys.argv[1], sys.argv[2])
print(com)
output = sp.check_output(com, shell=True).decode("ascii")
print(output)
resd = {}
for l in output.split("\n"):
    spl = l.split(",")
    print("?",spl)
    if len(spl) > 3:
        resd[spl[0]+spl[1]] = [int(i) for i in spl[2:4]]

failures = []
for q in queries:
    for r in recs:
        key = ">" + q.description + ">" + r.description
        h = -1*pairwise2.align.globalms(str(q.seq), str(r.seq), 0, -1, -2, -2, score_only=True)
        prevres = resd[key]
        print(h, prevres)
        if prevres[0] != h:
            failures.append(key)
            sys.stderr.write("Failed: %s\n" % key)

print("FAILURES")
print(failures)
print(len(failures))

