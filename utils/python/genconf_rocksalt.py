#!/usr/bin/env python
import sys
import numpy as np
import math
import copy
import os
cwd=os.getcwd()
#path=cwd+"/../../Alborz/utils/python"
path="/home/samare/Alborz/utils/python"
sys.path.insert(1,path)
#print sys.path
from atoms import *
from xyz import *
from ascii import *
from acf import *
#*****************************************************************************************
def cal_rbounds(atoms):
    xmin= 1.E20
    ymin= 1.E20
    zmin= 1.E20
    xmax=-1.E20
    ymax=-1.E20
    zmax=-1.E20
    for iat in range(atoms.nat):
        x=atoms.rat[iat][0]
        if x<xmin:
            xmin=x
        elif x>xmax:
            xmax=x
        y=atoms.rat[iat][1]
        if y<ymin:
            ymin=y
        elif y>ymax:
            ymax=y
        z=atoms.rat[iat][2]
        if z<zmin:
            zmin=z
        elif z>zmax:
            zmax=z
    return xmin,ymin,zmin,xmax,ymax,zmax
#*****************************************************************************************
def put_at_center(atoms,amargin,xmin,ymin,zmin):
    for iat in range(atoms.nat):
        atoms.rat[iat][0]=atoms.rat[iat][0]-xmin+amargin
        atoms.rat[iat][1]=atoms.rat[iat][1]-ymin+amargin
        atoms.rat[iat][2]=atoms.rat[iat][2]-zmin+amargin
#*****************************************************************************************
if len(sys.argv) < 5:
    print "usage: gen_rocksalt.py nat nx ny alat filename -bc boundcond"
    exit()
else:
    nat=int(sys.argv[1])
    nx=int(sys.argv[2])
    ny=int(sys.argv[3])
    alat=float(sys.argv[4])
    filename=sys.argv[5]
    boundcond='unknown'
    skip=False
    for i in range(1,len(sys.argv)):
        #print i,sys.argv[i]
        if skip:
            skip=False
            continue
        if sys.argv[i]=='-bc':
            if i+1<len(sys.argv):
                boundcond=sys.argv[i+1]
                skip=True
                #print boundcond
            else:
                sys.exit("ERROR: improper argument list")

#print boundcond
if nat%2!=0:
    print "ERROR: unacceptable number of atoms, nat must be an even number"
    exit()
if nat<2*nx*ny:
    print "ERROR: improper arguments, nat<2*nx*ny"
    exit()
nz=int(math.ceil(float(nat)/float(2*nx*ny)))
#print nz
#print type(nz)
#exit()

atoms=Atoms()
atoms.nat=0
for iz in range(nz):
    if atoms.nat>=nat: break
    for iy in range(ny):
        if atoms.nat>=nat: break
        for ix in range(nx):
            if atoms.nat>=nat: break
            atoms.nat+=1
            atoms.rat.append([])
            x=ix*alat
            y=iy*alat*0.5
            z=iz*alat*0.5
            atoms.rat[-1].append(x)
            atoms.rat[-1].append(y)
            atoms.rat[-1].append(z)
            ii=(3-(-1)**(2*ix+0+iy+iz))/2
            if ii==1:
                sat='Na'
            else:
                sat='Cl'
            atoms.sat.append(sat)
            atoms.bemoved.append("TTT")
            atoms.nat+=1
            atoms.rat.append([])
            x=ix*alat+alat*0.5
            atoms.rat[-1].append(x)
            atoms.rat[-1].append(y)
            atoms.rat[-1].append(z)
            #print (3-(-1)**(2*ix+1+iy+iz))/2
            ii=(3-(-1)**(2*ix+1+iy+iz))/2
            if ii==1:
                sat='Na'
            else:
                sat='Cl'
            atoms.sat.append(sat)
            atoms.bemoved.append("TTT")

xmin,ymin,zmin,xmax,ymax,zmax=cal_rbounds(atoms)
if boundcond=='free':
    amargin=4.0
    atoms.cellvec[0][0]=xmax-xmin+2.0*amargin
    atoms.cellvec[1][1]=ymax-ymin+2.0*amargin
    atoms.cellvec[2][2]=zmax-zmin+2.0*amargin
    #print atoms.cellvec[0][0],atoms.cellvec[1][1],atoms.cellvec[2][2]
    put_at_center(atoms,amargin,xmin,ymin,zmin)
elif boundcond=='bulk':
    atoms.cellvec[0][0]=nx*alat
    atoms.cellvec[1][1]=ny*alat*0.5
    atoms.cellvec[2][2]=nz*alat*0.5
    amargin=(nx*alat-xmax+xmin)*0.5
    put_at_center(atoms,amargin,xmin,ymin,zmin)
else:
    print 'ERROR: improper boundary conditions'

atoms.boundcond=boundcond
#print atoms.nat,len(atoms.sat),len(atoms.rat)

atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
print "writing to file %20s" % filename
#xyz_write(atoms_all,'bigdft')
#ascii_write(atoms,filename)
acf_write(atoms_all,filename)
#*****************************************************************************************
