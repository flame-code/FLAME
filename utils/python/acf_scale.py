#!/usr/bin/env python
import sys
import atoms
from acf import *
import copy

if len(sys.argv) < 6:
    print "usage: acf_scale.py filename scale_min scale_max nconf_min nconf_max"
    exit()
else:
    filename=sys.argv[1]
    scale_min=float(sys.argv[2])
    scale_max=float(sys.argv[3])
    nconf_min=int(sys.argv[4])
    nconf_max=int(sys.argv[5])

atoms_all=acf_read(filename)
if len(atoms_all)>1:
    print "WARNING: There are more than one configuration in input file,"
    print "         scaling is applied only to the first one."
atoms_all_out=[]

#print atoms_all[0].rat[0][0]
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

ds=(1.0-scale_min)/nconf_min
for i in range(nconf_min):
    scale=scale_min+float(i)*ds
    atoms=Atoms()
    atoms=copy.deepcopy(atoms_all[0])
    for iat in range(atoms_all[0].nat):
        atoms.rat[iat][0]=atoms_all[0].rat[iat][0]*scale
        atoms.rat[iat][1]=atoms_all[0].rat[iat][1]*scale
        atoms.rat[iat][2]=atoms_all[0].rat[iat][2]*scale
        #print atoms.rat[iat][0]/atoms_all[0].rat[iat][0]
    atoms_all_out.append(Atoms())
    atoms_all_out[-1]=copy.deepcopy(atoms)
    #print scale

atoms_all_out.append(Atoms())
atoms_all_out[-1]=copy.deepcopy(atoms_all[0])

ds=(scale_max-1.0)/nconf_max
for i in range(nconf_max):
    scale=1.0+float(i+1)*ds
    #atoms=Atoms()
    atoms=copy.deepcopy(atoms_all[0])
    for iat in range(atoms_all[0].nat):
        atoms.rat[iat][0]=atoms_all[0].rat[iat][0]*scale
        atoms.rat[iat][1]=atoms_all[0].rat[iat][1]*scale
        atoms.rat[iat][2]=atoms_all[0].rat[iat][2]*scale
        #print scale,iat
    atoms_all_out.append(Atoms())
    atoms_all_out[-1]=copy.deepcopy(atoms)
    #print scale

for iconf in range(len(atoms_all_out)):
    xcm=0.0
    ycm=0.0
    zcm=0.0
    for iat in range(atoms_all[0].nat):
        xcm+=atoms_all_out[iconf].rat[iat][0]
        ycm+=atoms_all_out[iconf].rat[iat][1]
        zcm+=atoms_all_out[iconf].rat[iat][2]
    xcm/=atoms_all_out[0].nat
    ycm/=atoms_all_out[0].nat
    zcm/=atoms_all_out[0].nat
    for iat in range(atoms_all[0].nat):
        atoms_all_out[iconf].rat[iat][0]+=xcmref-xcm
        atoms_all_out[iconf].rat[iat][1]+=ycmref-ycm
        atoms_all_out[iconf].rat[iat][2]+=zcmref-zcm

acf_write(atoms_all_out,"screen")
