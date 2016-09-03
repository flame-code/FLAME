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

atoms_all=acf_read(filename)

if len(atoms_all)==1:
    ascii_write(atoms_all[0],"screen")
else:
    print "\nATTENTION: The are more than one configuration in ACF file."
    prefix=raw_input("Please provide a prefix to generate files enumeratedly: [Default=tt]")
    if prefix=="": prefix="tt"
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        fnout="%s%5.5d.ascii" % (prefix,nconf)
        ascii_write(atoms,fnout)

