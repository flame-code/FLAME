#!/usr/bin/env python
import sys
import atoms
from acf import *
from xyz import *

if len(sys.argv) < 2:
    print "usage: acf2xyz.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms_all=acf_read(filename)
xyz_write(atoms_all,'bigdft')

