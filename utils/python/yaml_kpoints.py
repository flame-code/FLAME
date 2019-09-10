#!/usr/bin/env python
import sys
import math
import argparse
from io_yaml import *
from latvec2dproj import *

str1 = "This script reads a file in yaml format and calculates k-points."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('dkpt', action='store' ,type=float, help="Density of k-points")
args=parser.parse_args()

bohr2ang=0.52917720859
ang2bohr=1.0/0.52917720859

atoms_all=read_yaml(args.fn_inp)
iconf=0
for atoms in atoms_all:
    iconf+=1
    #rotation
    atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)
    if atoms.units_length_io=="angstrom":
        ax=atoms.cellvec[0][0]*ang2bohr
        bx=atoms.cellvec[1][0]*ang2bohr
        by=atoms.cellvec[1][1]*ang2bohr
        cx=atoms.cellvec[2][0]*ang2bohr
        cy=atoms.cellvec[2][1]*ang2bohr
        cz=atoms.cellvec[2][2]*ang2bohr
    else:
        ax=atoms.cellvec[0][0] 
        bx=atoms.cellvec[1][0] 
        by=atoms.cellvec[1][1]
        cx=atoms.cellvec[2][0]
        cy=atoms.cellvec[2][1]
        cz=atoms.cellvec[2][2]
    pi=math.acos(-1.0)
    a=ax
    b=math.sqrt(bx**2+by**2)
    c=math.sqrt(cx**2+cy**2+cz**2)
    vol=ax*by*cz
    g1x=2.0*pi*(by*cz)/vol  ;g1y=2.0*pi*(-bx*cz)/vol   ;g1z=2.0*pi*(bx*cy-cx*by)/vol
    g2x=0.0                 ;g2y=2.0*pi*(ax*cz)/vol    ;g2z=2.0*pi*(-ax*cy)/vol
    g3x=0.0                 ;g3y=0.0                   ;g3z=2.0*pi*(ax*by)/vol
    g=[];kpt=[]
    g.append(float(math.sqrt(g1x**2.0+g1y**2.0+g1z**2.0)))
    g.append(float(math.sqrt(g2x**2.0+g2y**2.0+g2z**2.0)))
    g.append(float(math.sqrt(g3x**2.0+g3y**2.0+g3z**2.0)))
    
    for i in range(0,3):
        kpt.append(int(g[i]/(args.dkpt*2.0*pi)))
        if kpt[i]==0:
            kpt[i]=1
        d_test=float(g[i]/(kpt[i]*2.0*pi))
        if d_test>=args.dkpt:
            for j in range(1,25):
                if d_test>=args.dkpt:
                    kpt[i]=kpt[i]+j
                    d_test=float(g[i]/(kpt[i]*2.0*pi))
    print "%s iconf: %5.5d %6d %6d %6d" % (args.fn_inp,iconf,kpt[0],kpt[1],kpt[2]) 
