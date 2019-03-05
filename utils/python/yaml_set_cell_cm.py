#!/usr/bin/env python
import argparse
import atoms
import numpy as np
from io_yaml import *

str1 = "This script adds margin to a structure."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file")
parser.add_argument('margin', type=float, help="""Warning: If margin is not set, default
        value (=3.0 Ang) will be used""", nargs='?', const=1, default=3.0)
args=parser.parse_args()
margin=args.margin

atoms_all=read_yaml(args.fn_inp)
nconf=-1
for atoms in atoms_all:
    nconf+=1
    columns =zip (*atoms_all[nconf].rat)
    minimum= []
    maximum= []
    for i, co in enumerate (columns):
        minimum.append(min(co))
        maximum.append(max(co))
    minratx = minimum[0]
    minraty = minimum[1]
    minratz = minimum[2]

    maxratx = maximum[0]
    maxraty = maximum[1]
    maxratz = maximum[2]
#    rat2=np.array(atoms_all[nconf].rat)
#    minratx= rat2.min(axis=0)[0]
#    minraty= rat2.min(axis=0)[1]
#    minratz= rat2.min(axis=0)[2]
#
#    maxratx= rat2.max(axis=0)[0]
#    maxraty= rat2.max(axis=0)[1]
#    maxratz= rat2.max(axis=0)[2]

    if atoms_all[nconf].boundcond=='free':
        atoms_all[nconf].cellvec[0][0]=maxratx-minratx+2*margin
        atoms_all[nconf].cellvec[1][1]=maxraty-minraty+2*margin
        atoms_all[nconf].cellvec[2][2]=maxratz-minratz+2*margin
        for i in range(int(atoms_all[nconf].nat)):
            atoms_all[nconf].rat[i][0]+=margin-minratx
            atoms_all[nconf].rat[i][1]+=margin-minraty
            atoms_all[nconf].rat[i][2]+=margin-minratz
    elif atoms_all[nconf].boundcond=='wire':
        atoms_all[nconf].cellvec[1][1]=maxraty-minraty+2*margin
        atoms_all[nconf].cellvec[2][2]=maxratz-minratz+2*margin
        for i in range(int(atoms_all[nconf].nat)):
            atoms_all[nconf].rat[i][0]=atoms_all[nconf].rat[i][0]%atoms_all[nconf].cellvec[0][0]
            atoms_all[nconf].rat[i][1]+=margin-minraty
            atoms_all[nconf].rat[i][2]+=margin-minratz
    elif atoms_all[nconf].boundcond=='slab':
        atoms_all[nconf].cellvec[2][2]=maxratz-minratz+2*margin
        for i in range(int(atoms_all[nconf].nat)):
            atoms_all[nconf].rat[i][0]=atoms_all[nconf].rat[i][0]%atoms_all[nconf].cellvec[0][0]
            atoms_all[nconf].rat[i][1]=atoms_all[nconf].rat[i][1]%atoms_all[nconf].cellvec[1][1]
            atoms_all[nconf].rat[i][2]+=margin-minratz
    elif atoms_all[nconf].boundcond=='bulk':
        for i in range(int(atoms_all[nconf].nat)):
            atoms_all[nconf].rat[i][0]=atoms_all[nconf].rat[i][0]%atoms_all[nconf].cellvec[0][0]
            atoms_all[nconf].rat[i][1]=atoms_all[nconf].rat[i][1]%atoms_all[nconf].cellvec[1][1]
            atoms_all[nconf].rat[i][2]=atoms_all[nconf].rat[i][2]%atoms_all[nconf].cellvec[2][2]
    else :
        print "ERROR: unknown boundary condition"

write_yaml(atoms_all,args.fn_out)
