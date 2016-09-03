#!/usr/bin/env python
import sys
import numpy as np
import math
import copy
import os
cwd=os.getcwd()
path=cwd+"/../../utils/python"
sys.path.insert(1,path)
from atoms import *
from poisson import *
from cube import *
#*****************************************************************************************
if len(sys.argv) < 2:
    print "usage: cube_y2z.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms_old,poisson_old=cube_read(filename)

atoms=copy.deepcopy(atoms_old)
poisson=copy.deepcopy(poisson_old)
poisson.ngpy=poisson_old.ngpz
poisson.ngpz=poisson_old.ngpy
poisson.rho=np.zeros([poisson.ngpz,poisson.ngpy,poisson.ngpx])
for igpz in range(poisson.ngpz):
    for igpy in range(poisson.ngpy):
        for igpx in range(poisson.ngpx):
            poisson.rho[igpz][igpy][igpx]=poisson_old.rho[igpy][igpz][igpx]

poisson.hy=poisson_old.hz
poisson.hz=poisson_old.hy

#for i in range(3):
#    atoms.cellvec[1]=

for iat in range(atoms.nat):
    #print "%14.6E%14.6E" % (atoms_old.rat[iat][1],atoms.rat[iat][1])
    #print "%14.6E%14.6E%14.6E" % (atoms_old.rat[iat][0],atoms_old.rat[iat][1],atoms_old.rat[iat][2])
    atoms.rat[iat][1]=atoms_old.rat[iat][2]
    atoms.rat[iat][2]=atoms_old.rat[iat][1]
    #print "%14.6E%14.6E" % (atoms_old.rat[iat][1],atoms.rat[iat][1])
    #print "%14.6E%14.6E%14.6E" % (atoms_old.rat[iat][0],atoms_old.rat[iat][1],atoms_old.rat[iat][2])

frmt1="%5d%24.15E%24.15E%24.15E\n"
frmt2="%14.6E"
fnout="y2z_%s" % filename
print "writing transformed cube file into %s" % fnout
cube_write(fnout,atoms,poisson.ngpx,poisson.ngpy,poisson.ngpz,poisson.rho,poisson.hx,poisson.hy,poisson.hz,frmt1,frmt2)
#*****************************************************************************************
