#!/usr/bin/env python
import sys
import atoms
from acf import *
from xyz import *

if len(sys.argv) < 5:
    print "usage: acf_translate.py input_filename tx ty tz"
    exit()
else:
    filename=sys.argv[1]
    tx=float(sys.argv[2])
    ty=float(sys.argv[3])
    tz=float(sys.argv[4])

atoms_all=acf_read(filename)

for atoms in atoms_all:
    for iat in range(atoms.nat):
        #atoms.rat[iat][0]=str(tx+float(atoms.rat[iat][0]))
        #atoms.rat[iat][1]=str(ty+float(atoms.rat[iat][1]))
        #atoms.rat[iat][2]=str(tz+float(atoms.rat[iat][2]))
        atoms.rat[iat][0]+=tx
        atoms.rat[iat][1]+=ty
        atoms.rat[iat][2]+=tz

filename_out="trans_"+filename
print "\nWriting translated configurations to %s\n" % filename_out
acf_write(atoms_all,filename_out)
