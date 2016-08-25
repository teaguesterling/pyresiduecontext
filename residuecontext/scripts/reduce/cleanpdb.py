#!/usr/bin/env python

import sys
import pdb

OK_RESIDUES = set()
OK_RESIDUES.update(pdb.aminoAcid3Codes)
OK_RESIDUES.update(pdb.keepPolarH.keys())
OK_RESIDUES.remove('HEM')  # Only keep ions and waters in addition

p = pdb.pdbData(sys.argv[1], ignoreWaters=False)
mode = sys.argv[3]

for idx, line in enumerate(p.rawData):
    p.rawData[idx] = line[:79]

if mode == 'final':
#    p.replaceHETATMwithATOM()
    p.removeApolarHydrogen()
    p.write(sys.argv[2]+'.H')
    p = pdb.pdbData(sys.argv[2]+'.H', ignoreWaters=False)
    p.renameHistidines()
    p.renameCysteines()
    

for idx, resName in enumerate(p.resNames):
    if resName not in OK_RESIDUES:
        p.removeLine(idx)

for idx, alt in enumerate(p.altChars):
    if alt not in (' ', 'A'):
        p.removeLine(idx)


p.write(sys.argv[2])



