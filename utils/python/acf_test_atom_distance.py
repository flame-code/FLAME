#!/usr/bin/env python
import sys
import atoms
import math
#import numpy as np
from acf import *
if len(sys.argv) < 2:
    print ""
    print "usage: acf_test_atom_distance.py input_filename dmin output_filename dmax(optional)"
    print " dmin :  minimum distance between 2 atoms."
    print " dmax :  maximum distance between atoms and their nearest neighbour(for free BC). "
    print ""
    exit()
else:
    filename = sys.argv[1]
    dmin=float(sys.argv[2])
    if len(sys.argv) >= 4:
        ofilename=sys.argv[3]
    else :
        ofilename="screen"

if len(sys.argv) == 5:
    dmax=float(sys.argv[4])
else :
    dmax = 5
atoms_all=acf_read(filename)
atoms_all_sel =[]

nconf=-1
nconf2=-1
for atoms in atoms_all:
    nconf+=1
    test = 1
    for iat in range(int(atoms_all[nconf].nat)):
        if (test<1) :
            break
        for jat in range(iat+1,int(atoms_all[nconf].nat)):
            ttx=atoms.rat[iat][0]-atoms.rat[jat][0]
            tty=atoms.rat[iat][1]-atoms.rat[jat][1]
            ttz=atoms.rat[iat][2]-atoms.rat[jat][2]
            distance=math.sqrt(ttx*ttx+tty*tty+ttz*ttz)
            if distance<dmin:
                print "too close atoms: iconf,iat,jat,distance:  " ,nconf+1,iat+1,jat+1,distance
                test=0
                break

    if atoms_all[nconf].boundcond=='free':
        for iat in range(int(atoms_all[nconf].nat)):
            min_iat = 5
            if (test<1) :
                break
            for jat in range(int(atoms_all[nconf].nat)):
                ttx=atoms.rat[iat][0]-atoms.rat[jat][0]
                tty=atoms.rat[iat][1]-atoms.rat[jat][1]
                ttz=atoms.rat[iat][2]-atoms.rat[jat][2]
                distance=math.sqrt(ttx*ttx+tty*tty+ttz*ttz)
                if distance< min_iat :
                    if  iat!=jat:
                         min_iat = distance
            if  min_iat >dmax:
                print "too far atom: iconf,iat,distance:  " ,nconf+1,iat+1, min_iat
                test=0
                break
    if (test>0) :
        nconf2+=1
        atoms_all_sel.append(Atoms())
        atoms_all_sel[-1]=copy.copy(atoms)

acf_write(atoms_all_sel,ofilename)
print "total number of conf = ",nconf+1,"       the number of files deleted = ",nconf-nconf2
