#!/usr/bin/env python
import sys
import atoms
from math import *
from acf import *

if len(sys.argv) < 2:
    print "usage: acf_getvol.py input_filename"
    exit()
else:
    filename = sys.argv[1]

#**************************************************************************
def getvol(cellvec):
    vol=(cellvec[0][0]*cellvec[1][1]*cellvec[2][2]-cellvec[0][0]*cellvec[1][2]*cellvec[2][1]- 
        cellvec[0][1]*cellvec[1][0]*cellvec[2][2]+cellvec[0][1]*cellvec[1][2]*cellvec[2][0]+
        cellvec[0][2]*cellvec[1][0]*cellvec[2][1]-cellvec[0][2]*cellvec[1][1]*cellvec[2][0])
#    print vol
    if vol<0.0: 
        print "ERROR: Negative volume!"
        sys.exit("error")
    return vol
#**************************************************************************
atoms_all = []
atoms_all=acf_read(filename)
if len(atoms_all)==1:
#print atoms.cellvec
    vol=getvol(atoms_all[0].cellvec)
    vol_atom=vol/atoms_all[0].nat
    #print vol_atom
    print "%50s%10.5f" % (filename,vol_atom)
else:
    print "\nATTENTION: The are more than one configuration in ACF file."
