#!/usr/bin/env python
import sys
import atoms
from latvec2dproj import *
from vasp import *
from ascii import *

if len(sys.argv) < 2:
    print "usage: poscar2ascii.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=poscar_read(filename)
dproj,atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)

atoms.cellvec[0][0]=dproj[0]
atoms.cellvec[1][0]=dproj[1]
atoms.cellvec[1][1]=dproj[2]
atoms.cellvec[2][0]=dproj[3]
atoms.cellvec[2][1]=dproj[4]
atoms.cellvec[2][2]=dproj[5]


#atoms_all=[]
#atoms_all.append(Atoms())
#atoms_all[-1]=copy.copy(atoms)
ascii_write(atoms,"screen")
