#!/usr/bin/env python
import sys
import atoms
from latvec2dproj import *
from vasp import *
from shengbte import *

if len(sys.argv) < 2:
    print "usage: poscar2shengbte.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms=poscar_read(filename)
sheng_write(atoms,"screen")
