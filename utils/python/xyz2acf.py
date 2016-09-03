#!/usr/bin/env python
import sys
from xyz import *
from acf import *

if len(sys.argv) < 2:
    print "usage: xyz2acf.py input_filename -qtot qtot"
    exit()
else:
    qtot=0.0
    skip=False
    for i in range(1,len(sys.argv)):
        #print i,sys.argv[i]
        if skip:
            skip=False
            continue
        if sys.argv[i]=='-qtot':
            if i+1<len(sys.argv):
                qtot=float(sys.argv[i+1])
                skip=True
                #print qtot
            else:
                sys.exit("ERROR: improper argument list")
        else:
            filename = sys.argv[i]

#print filename,qtot

atoms_all=xyz_read(filename,'bigdft')

for atoms in atoms_all:
    atoms.qtot=qtot

acf_write(atoms_all,"screen")
