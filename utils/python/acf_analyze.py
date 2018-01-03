#!/usr/bin/env python
#import random
import sys
#from cStringIO import StringIO
import math
#path="/home/ghasemi/FLAME/utils/python"
#sys.path.append(path)
from acf import *
from cellutils import *


if len(sys.argv) < 2:
    print "usage: acf_analyze.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms_all=acf_read(filename=filename)

if len(atoms_all)>1:
    print "ERROR: only one configuration in a file is allowed: %6d" % len(atoms_all)
    exit()

for atoms in atoms_all:
    ratred=rxyz_cart2int(atoms.cellvec,atoms.rat,atoms.nat)
    rmin=1.E20
    for iat in range(atoms.nat):
        for i in range(-1,2):
            for j in range(-1,2):
                for k in range(-1,2):
                    for jat in range(atoms.nat):
                        if i==0 and j==0 and k==0 and jat==iat: continue
                        sx=ratred[jat][0]+i
                        sy=ratred[jat][1]+j
                        sz=ratred[jat][2]+k
                        dsx=sx-ratred[iat][0]
                        dsy=sy-ratred[iat][1]
                        dsz=sz-ratred[iat][2]
                        dx=atoms.cellvec[0][0]*dsx+atoms.cellvec[1][0]*dsy+atoms.cellvec[2][0]*dsz
                        dy=atoms.cellvec[0][1]*dsx+atoms.cellvec[1][1]*dsy+atoms.cellvec[2][1]*dsz
                        dz=atoms.cellvec[0][2]*dsx+atoms.cellvec[1][2]*dsy+atoms.cellvec[2][2]*dsz
                        r=(dx**2+dy**2+dz**2)**0.5
                        rmin=min(rmin,r)
    break #assumed one configuration in the file.

print "%50s%10.5f" % (filename,rmin)
