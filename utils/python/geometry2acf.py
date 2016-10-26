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


if not atoms.boundcond=="free":
    atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)

atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
acf_write(atoms_all,"screen")
