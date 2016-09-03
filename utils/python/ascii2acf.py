#!/usr/bin/env python
import sys
import atoms
from acf import *
from ascii import *

if len(sys.argv) < 2:
    print "usage: acf2ascii.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=ascii_read(filename)

atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
acf_write(atoms_all,"screen")
