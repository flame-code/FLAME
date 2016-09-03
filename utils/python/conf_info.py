#!/usr/bin/env python
import sys
path="/home/ghasemi/Alborz/utils/python"
sys.path.insert(1,path)
from atoms import *
from acf import *

if len(sys.argv) < 2:
    print "usage: acf2ascii.py input_filename"
    exit()
else:
    filename = sys.argv[1]

atoms_all=acf_read(filename)

stat=True
for atoms in atoms_all:
    xmin=1.E100 ; xmax=-1.E100
    ymin=1.E100 ; ymax=-1.E100
    zmin=1.E100 ; zmax=-1.E100
    for iat in range(atoms.nat):
        if atoms.rat[iat][0]<xmin: xmin=atoms.rat[iat][0]
        if atoms.rat[iat][0]>xmax: xmax=atoms.rat[iat][0]
        if atoms.rat[iat][1]<ymin: ymin=atoms.rat[iat][1]
        if atoms.rat[iat][1]>ymax: ymax=atoms.rat[iat][1]
        if atoms.rat[iat][2]<zmin: zmin=atoms.rat[iat][2]
        if atoms.rat[iat][2]>zmax: zmax=atoms.rat[iat][2]
    dx=xmax-xmin
    dy=ymax-ymin
    dz=zmax-zmin
    amargin=1.E100
    amargin=min(amargin,atoms.cellvec[0][0]-dx)
    amargin=min(amargin,atoms.cellvec[1][1]-dy)
    amargin=min(amargin,atoms.cellvec[2][2]-dz)
    print "%6.1f%6.1f%6.1f%6.1f" % (dx,dy,dz,amargin)
    if amargin<6.0:
        stat=False
        #break
print "STAT: %50s   " % filename,stat

