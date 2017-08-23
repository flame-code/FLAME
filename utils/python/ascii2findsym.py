#!/usr/bin/env python
import sys, copy
import numpy as np
import math as mt
from atoms import *
from ascii import *
from cellutils import *
#************************************************************************#main program
if len(sys.argv) < 3:
     print "usage: ascii2findsym.py input_filename tolerance"
     exit()
else:
    filename = sys.argv[1]
    tol      = sys.argv[2]
tol = float(tol)
atoms=ascii_read(filename)
atoms.rat = backtocell(atoms.nat,atoms.cellvec,atoms.rat)
f = open("findsym.in","w")
f.write("Input for FINDSYM, from filename %s\n" %(filename))
f.write("%15.7E %s\n" %(tol,"Tolerance"))
f.write("1    form of lattice parameters: to be entered as vectors\n")
f.write(" %15.8f %15.8f %15.8f %s\n" % (atoms.cellvec[0][0],atoms.cellvec[0][1],atoms.cellvec[0][2],"Lattice vector 1"))
f.write(" %15.8f %15.8f %15.8f %s\n" % (atoms.cellvec[1][0],atoms.cellvec[1][1],atoms.cellvec[1][2],"Lattice vector 2"))
f.write(" %15.8f %15.8f %15.8f %s\n" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2],"Lattice vector 3"))
f.write("2    form of primitive lattice vectors\n")
f.write("P    unknown centering\n")
f.write("%5d %s\n" %(atoms.nat,"   number of atoms"))
for itype in range(len(atoms.sat)):
    f.write("%3s" %(atoms.sat[itype]))
f.write("   Atom kinds\n")
ratint = rxyz_cart2int(atoms.cellvec,atoms.rat,atoms.nat)
for iat in range(atoms.nat):
    f.write("%26.16E%26.16E%26.16E%2s%5d\n" %(ratint[iat][0],ratint[iat][1],ratint[iat][2], "  Reduced Coordinates", iat+1))
