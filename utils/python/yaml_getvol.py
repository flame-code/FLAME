#!/usr/bin/env python
import sys
import atoms
from math import *
from io_yaml import *

if len(sys.argv) < 2:
    print("usage: yaml_getvol.py input_filename")
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
        print("ERROR: Negative volume!")
        sys.exit("error")
    return vol
#**************************************************************************
atoms_all = []
atoms_all=read_yaml(filename)
for atoms in atoms_all:
    vol=getvol(atoms.cellvec)
    vol_atom=vol/atoms.nat
    #print "%50s%10.5f" % (filename,vol_atom)
    #print("%5d%10.5f" % (atoms.nat,vol_atom))
    print("%5d%10.5f%15.6E" % (atoms.nat,vol_atom,atoms.epot/atoms.nat))

#if len(atoms_all)==1:
##print atoms.cellvec
#    vol=getvol(atoms_all[0].cellvec)
#    vol_atom=vol/atoms_all[0].nat
#    #print vol_atom
#    print "%50s%10.5f" % (filename,vol_atom)
#else:
#    print "\nATTENTION: The are more than one configuration in YAML file."
