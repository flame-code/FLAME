#!/usr/bin/env python
import argparse
import atoms
import copy
from io_yaml import *

str1 = "This script reads a file in the yaml dictionary format and scales it based on given scale factor."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('scale', type=float, help='the scale factor')
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in yaml format")
args=parser.parse_args()

scale = args.scale
atoms_all=read_yaml(args.fn_inp)

atoms_all_out=[]

#print atoms_all[0].rat[0][0]

for iconf in range(len(atoms_all)):
    atoms=Atoms()
    atoms=copy.deepcopy(atoms_all[iconf])
    atoms.cellvec[0][0]=atoms_all[iconf].cellvec[0][0]*scale
    atoms.cellvec[1][0]=atoms_all[iconf].cellvec[1][0]*scale
    atoms.cellvec[1][1]=atoms_all[iconf].cellvec[1][1]*scale
    atoms.cellvec[2][0]=atoms_all[iconf].cellvec[2][0]*scale
    atoms.cellvec[2][1]=atoms_all[iconf].cellvec[2][1]*scale
    atoms.cellvec[2][2]=atoms_all[iconf].cellvec[2][2]*scale
    for iat in range(atoms_all[iconf].nat):
        atoms.rat[iat][0]=atoms_all[iconf].rat[iat][0]*scale
        atoms.rat[iat][1]=atoms_all[iconf].rat[iat][1]*scale
        atoms.rat[iat][2]=atoms_all[iconf].rat[iat][2]*scale
        #print atoms.rat[iat][0]/atoms_all[0].rat[iat][0]
    atoms_all_out.append(Atoms())
    atoms_all_out[-1]=copy.deepcopy(atoms)
    #print scale

for iconf in range(len(atoms_all_out)):
    if atoms_all_out[iconf].boundcond=="free":
        xcmref=0.0
        ycmref=0.0
        zcmref=0.0
        for iat in range(atoms_all[0].nat):
            xcmref+=atoms_all[0].rat[iat][0]
            ycmref+=atoms_all[0].rat[iat][1]
            zcmref+=atoms_all[0].rat[iat][2]
        xcmref=xcmref/atoms_all[0].nat
        ycmref=ycmref/atoms_all[0].nat
        zcmref=zcmref/atoms_all[0].nat
    elif atoms_all_out[iconf].boundcond=="bulk":
        xcmref=0.0 #atoms_all_out[iconf].cellvec[0][0]/2.0
        ycmref=0.0 #atoms_all_out[iconf].cellvec[1][1]/2.0
        zcmref=0.0 #atoms_all_out[iconf].cellvec[2][2]/2.0
    elif atoms_all_out[iconf].boundcond=="slab":
        zcmref=0.0
        for iat in range(atoms_all[0].nat):
            zcmref+=atoms_all[0].rat[iat][2]
        xcmref=atoms_all_out[iconf].cellvec[0][0]/2.0
        ycmref=atoms_all_out[iconf].cellvec[1][1]/2.0
        zcmref=zcmref/atoms_all[0].nat
    else:
        print "ERROR: unknow boundary conditions"
    xcm=0.0
    ycm=0.0
    zcm=0.0
    for iat in range(atoms_all[iconf].nat):
        xcm+=atoms_all_out[iconf].rat[iat][0]
        ycm+=atoms_all_out[iconf].rat[iat][1]
        zcm+=atoms_all_out[iconf].rat[iat][2]
    xcm/=atoms_all_out[0].nat
    ycm/=atoms_all_out[0].nat
    zcm/=atoms_all_out[0].nat
    if atoms_all_out[iconf].boundcond=="bulk":
        xcm=0.0
        ycm=0.0
        zcm=0.0
    for iat in range(atoms_all[iconf].nat):
        atoms_all_out[iconf].rat[iat][0]+=xcmref-xcm
        atoms_all_out[iconf].rat[iat][1]+=ycmref-ycm
        atoms_all_out[iconf].rat[iat][2]+=zcmref-zcm

write_yaml(atoms_all_out,args.fn_out)
