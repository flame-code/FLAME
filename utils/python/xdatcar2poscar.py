#!/usr/bin/env python
import sys
import atoms
from vasp import *
from ascii import *

atoms_all=xdatcar_read()

nconf=0
for atoms in atoms_all:
    nconf+=1
#    fnout="POSCAR_%1.1d" % (nconf)
    fnout="POSCAR_%5.5d" % (nconf)
    poscar_write(atoms,fnout)

