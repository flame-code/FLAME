#!/usr/bin/env python
import sys
import atoms
from acf import *
from vasp import *

if len(sys.argv) < 2:
    print "usage: acf2poscar.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms_all=acf_read(filename)

if len(atoms_all)==1:
    poscar_write(atoms_all[0],"screen")
else:
    print "\nATTENTION: The are more than one configuration in ACF file."
    prefix=raw_input("Please provide a prefix to generate files enumeratedly: ")
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        fnout="%s%5.5d" % (prefix,nconf)
        poscar_write(atoms,fnout)
