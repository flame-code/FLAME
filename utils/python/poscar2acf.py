#!/usr/bin/env python
import sys
import atoms
from latvec2dproj import *
from vasp import *
from acf import *

if len(sys.argv) < 2:
    print "usage: poscar2acf.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=poscar_read(filename)
if not atoms.boundcond=="free":
    atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)

atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
acf_write(atoms_all,"screen")
