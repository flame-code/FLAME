#!/usr/bin/env python
import sys
import atoms
from latvec2dproj import *
from aims import *
from acf import *

if len(sys.argv) < 2:
    print "usage: geometry2acf.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=aims_read(filename)
dproj,atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)

atoms.cellvec[0][0]=dproj[0]
atoms.cellvec[1][0]=dproj[1]
atoms.cellvec[1][1]=dproj[2]
atoms.cellvec[2][0]=dproj[3]
atoms.cellvec[2][1]=dproj[4]
atoms.cellvec[2][2]=dproj[5]


#atoms_all=[]
#atoms_all.append(Atoms())
#atoms_all=copy.copy(atoms)
acf_write_b(atoms,"screen")
