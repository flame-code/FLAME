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
if len(sys.argv) < 4:
    print "usage: cube_expand.py input_filename ngpz_bot ngpz_top"
    exit()
else:
    filename = sys.argv[1]
    ngpz_bot = int(sys.argv[2])
    ngpz_top = int(sys.argv[3])

atoms_old,poisson_old=cube_read(filename)

#if ngpz_add%2==1:
#    ngpz_add+=1
#    print "***** notice that ngpz_add is odd number so it is increased by one."
#ngpz_add_half=ngpz_add/2

atoms=copy.deepcopy(atoms_old)
poisson=copy.deepcopy(poisson_old)
poisson.ngpz=poisson_old.ngpz+ngpz_bot+ngpz_top
poisson.rho=np.zeros([poisson.ngpz,poisson.ngpy,poisson.ngpx])
for igpz in range(poisson_old.ngpz):
    for igpy in range(poisson.ngpy):
        for igpx in range(poisson.ngpx):
            poisson.rho[igpz+ngpz_bot][igpy][igpx]=poisson_old.rho[igpz][igpy][igpx]

shift=ngpz_bot*poisson.hz
for iat in range(atoms.nat):
    atoms.rat[iat][2]+=shift

frmt1="%5d%24.15E%24.15E%24.15E\n"
frmt2="%14.6E"
fnout="expanded_%s" % filename
print "writing expanded cube file into %s" % fnout
cube_write(fnout,atoms,poisson.ngpx,poisson.ngpy,poisson.ngpz,poisson.rho,poisson.hx,poisson.hy,poisson.hz,frmt1,frmt2)
#*****************************************************************************************
