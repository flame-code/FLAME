#!/usr/bin/env python
import random
import sys
from cStringIO import StringIO
from acf import *
import math

if len(sys.argv) < 3:
    print "usage: acf_shake.py input_filename amplitude"
    exit()
else:
    filename = sys.argv[1]
    ampl=float(sys.argv[2])
    filename2 ="screen"
    if len(sys.argv)==4:
        filename2 = sys.argv[3]
    if len(sys.argv)==7:
        filename2 = sys.argv[6]
    if len(sys.argv)>3:
        mvx = float(sys.argv[3])
        mvy = float(sys.argv[4])
        mvz = float(sys.argv[5])
    else:
        mvx = 1
        mvy = 1
        mvz = 1 


atoms_all=acf_read(filename=filename)

#print atoms_all[0].nat

random.seed()
for atoms in atoms_all:
    for iat in range(atoms.nat):
        #ttx=random.random()
        #tty=random.random()
        #ttz=random.random()
        #atoms.positions[iat,0]+=(ttx-0.5)*ampl*2.0
        #atoms.positions[iat,1]+=(tty-0.5)*ampl*2.0
        #atoms.positions[iat,2]+=(ttz-0.5)*ampl*2.0
        theta=math.pi*random.random()
        phi=2.0*math.pi*random.random()
        ttx=mvx*ampl*math.sin(theta)*math.cos(phi)
        tty=mvy*ampl*math.sin(theta)*math.sin(phi)
        ttz=mvz*ampl*math.cos(theta)
        atoms.rat[iat][0]+=ttx
        atoms.rat[iat][1]+=tty
        atoms.rat[iat][2]+=ttz

acf_write(atoms_all,filename2)
