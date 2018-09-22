#!/usr/bin/env python
import sys
import atoms
from ascii import *
from vasp import *

if len(sys.argv) < 2:
    print "usage: ascii2poscar.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=ascii_read(filename)

#atoms_all=[]
#atoms_all.append(Atoms())
#atoms_all[-1]=copy.copy(atoms)
poscar_write(atoms,"screen")
