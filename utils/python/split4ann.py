#!/usr/bin/env python
import sys
import atoms
from acf import *
from xyz import *
import random

if len(sys.argv) < 3:
    print "usage: split4ann.py input_filename percentage"
    exit()
else:
    filename = sys.argv[1]
    percent=float(sys.argv[2])/100.0

atoms_all=acf_read(filename)

atoms_train=[]
atoms_valid=[]
for atoms in atoms_all:
    if random.random()<percent:
        atoms_valid.append(Atoms())
        atoms_valid[-1]=copy.copy(atoms)
    else:
        atoms_train.append(Atoms())
        atoms_train[-1]=copy.copy(atoms)

nt=len(atoms_train)
nv=len(atoms_valid)
if nt>0:
    fnout="train_%s" % filename
    acf_write(atoms_train,fnout)
if nv>0:
    fnout="valid_%s" % filename
    acf_write(atoms_valid,fnout)
print "number of configurations selected for training and validation: %6d and %6d   from %s" % (nt,nv,filename)
