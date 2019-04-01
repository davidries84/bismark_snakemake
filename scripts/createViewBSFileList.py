import sys
import os

samples = sys.argv[1:]
outlist = []
for sample in samples:
        base = os.path.basename(sample)
        abspath = os.path.abspath(sample)
        samplename = base.split('.')[0]
        outlist.append(abspath+","+samplename)

print ("\n".join(outlist))
