#!/usr/bin/env python
import sys
import atoms
import copy
from cif import *
from ascii import *

if len(sys.argv) < 2:
    print "usage: ascii2cif.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=ascii_read(filename)

atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)

cif_write(atoms_all[0],"screen")
