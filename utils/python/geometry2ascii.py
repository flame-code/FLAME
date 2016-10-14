#!/usr/bin/env python
import sys
import atoms
from latvec2dproj import *
from aims import *
from ascii import *

if len(sys.argv) < 2:
    print "usage: geometry2ascii.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=aims_read(filename)
atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)

#atoms_all=[]
#atoms_all.append(Atoms())
#atoms_all[-1]=copy.copy(atoms)
ascii_write(atoms,"screen")
