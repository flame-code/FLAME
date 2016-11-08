#!/usr/bin/env python
#This program converts POSCAR or geometry.in format to ascii or acf format by applying the rotation matrix.
from math import *
import math as mt
import numpy as np
import sys,copy,argparse
from atoms import *
from vasp import *
from ascii import *
from acf import *
from aims import *
#!*****************************************************************************************
def reportminmax(nat,rat,cellvec):
    xmin=1.E10 ; xmax=-1.E10
    ymin=1.E10 ; ymax=-1.E10
    zmin=1.E10 ; zmax=-1.E10
    for  iat in range(nat):
        xmin=min(xmin,rat[iat][0]) ; xmax=max(xmax,rat[iat][0])
        ymin=min(ymin,rat[iat][1]) ; ymax=max(ymax,rat[iat][1])
        zmin=min(zmin,rat[iat][2]) ; zmax=max(zmax,rat[iat][2])
    
    xdiff=cellvec[0][0]-xmax
    ydiff=cellvec[1][1]-ymax
    zdiff=cellvec[2][2]-zmax    
#!*****************************************************************************************
def latvec2dproj(cellvec,rxyz,nat):
# This subroutine will convert the lattice vector representation of the
# periodic cell (vec1,vec2,vec3) into the projective representation (dxx,dyx,dyy,dzx,dzy,dzz)
# The cell will thus be rotated. The rotational matrix is stored in rotmat as an operator rotmat
# and the atomic position rxyz are transformed into the new coordination sizstem as well
    eps=1E-6
    axes=np.zeros((3,3))
    axes[0][0]=1.0 ; axes[1][1]=1.0 ; axes[2][2]=1.0

    #Calculating dxx
    dproj=np.zeros(6)
    dproj[0]=mt.sqrt(cellvec[0][0]*cellvec[0][0]+cellvec[0][1]*cellvec[0][1]+cellvec[0][2]*cellvec[0][2])
    
    #Calculate the first rotation to align the first axis in x direction---------------------
    rotmat1=np.zeros((3,3))
    for i in range(0,3):
        rotmat1[i][i]=1.0

    tempvec=np.zeros(3)
    tempvec[0]=1.0
    crossp=np.zeros(3)
    axe=np.zeros(3)
    #tempvec is the x-unit vector
    crossp=np.cross(cellvec[:][0],tempvec)

    norotation1=-1
    if (abs(crossp[0])<eps*0.1) and (abs(crossp[1])<eps*0.1) and (abs(crossp[2])<eps*0.1):  #1001
        norotation1=1
    if norotation1<0:
        norm=mt.sqrt(crossp[0]*crossp[0]+crossp[1]*crossp[1]+crossp[2]*crossp[2])
        axe=crossp/norm
        alpha=mt.acos(np.dot(tempvec,cellvec[0][:])/dproj[0])
        rotmat1=rotation(alpha,axe)
        cellvec=matmul(rotmat1,cellvec)
    #1001
    if (cellvec[0][1]>eps) or (cellvec[0][2]>eps):
        print "Error in 1. rotation",cellvec[0][1],cellvec[0][2]
        sys.exit("error")
    #Calculate the second rotation to align the second axis in xy plane------------------------
    rotmat2=np.zeros((3,3))
    for i in range(3):
        rotmat2[i][i]=1.0

    tempvec=cellvec[1][:]
    tempvec[0]=0.0
    axe=np.cross(tempvec,axes[1][:])
    
    norotation2=-1
    if (abs(axe[0])<eps*0.1) and (abs(axe[1])<eps*0.1) and (abs(axe[2])<eps*0.1) and (tempvec[1]>0.0):
        norotation2=1
    if norotation2<0:
        norm=mt.sqrt(axe[0]*axe[0]+axe[1]*axe[1]+axe[2]*axe[2])
        axe=axe/norm
        crossp=np.cross(axe,cellvec[1][:])
    norotation3=-1
    if (abs(crossp[0])<eps*0.1) and (abs(crossp[1])<eps*0.1) and (crossp[2]<0.0):
        norotation3=1
    if norotation3<0:
        alpha=angle_between(axes[1][:],tempvec[:])
        rotmat2=rotation(alpha,axe)
        cellvec=matmul(rotmat2,cellvec)
    #1002 continue
    if cellvec[1][2]>eps:
        print "Error in 2. rotation"
        sys.exit("error") 
    if cellvec[2][2]<0.0:
        print "ERROR: in orientation of the cell"
        sys.exit("error") 
    
    #The total rotational matrix:-------------------------------------------------------------
    rotmat=matmul(rotmat2,rotmat1)
    #Apply rotation on all atoms
    rxyzo = copy.deepcopy(rxyz)
    for iat in range(nat):
        rxyz[iat][0]=rotmat[0][0]*rxyzo[iat][0]+rotmat[1][0]*rxyzo[iat][1]+rotmat[2][0]*rxyzo[iat][2]
        rxyz[iat][1]=rotmat[0][1]*rxyzo[iat][0]+rotmat[1][1]*rxyzo[iat][1]+rotmat[2][1]*rxyzo[iat][2]
        rxyz[iat][2]=rotmat[0][2]*rxyzo[iat][0]+rotmat[1][2]*rxyzo[iat][1]+rotmat[2][2]*rxyzo[iat][2]
    #Calculate all other elements of dproj
    dproj[1]=cellvec[1][0]
    dproj[2]=cellvec[1][1]
    dproj[3]=cellvec[2][0]
    dproj[4]=cellvec[2][1]
    dproj[5]=cellvec[2][2]
    return (dproj,cellvec,rxyz)
#!*****************************************************************************************
def rotation(angle,axe):
#This subroutine will calculate the rotational matrix rotmat for a
#3-dim vector around an axis 'axe' by the angle 'angle'.
    rotator=np.zeros((3,3))
#    !Define Rotation Matrix
    rotator[0][0]=mt.cos(angle)+(axe[0]**2)*(1.0-mt.cos(angle))
    rotator[1][0]=axe[0]*axe[1]*(1.0-mt.cos(angle))-axe[2]*mt.sin(angle)
    rotator[2][0]=axe[0]*axe[2]*(1.0-mt.cos(angle))+axe[1]*mt.sin(angle)
    
    rotator[0][1]=axe[1]*axe[0]*(1.0-mt.cos(angle))+axe[2]*mt.sin(angle)
    rotator[1][1]=mt.cos(angle)+(axe[1]**2)*(1.0-mt.cos(angle))
    rotator[2][1]=axe[1]*axe[2]*(1.0-mt.cos(angle))-axe[0]*mt.sin(angle)
    
    rotator[0][2]=axe[2]*axe[0]*(1.0-mt.cos(angle))-axe[1]*mt.sin(angle)
    rotator[1][2]=axe[2]*axe[1]*(1.0-mt.cos(angle))+axe[0]*mt.sin(angle)
    rotator[2][2]=mt.cos(angle)+(axe[2]**2)*(1.0-mt.cos(angle))
    
    return rotator
#!*****************************************************************************************
def matmul(a,b):
    rows_a = len(a); cols_a = len(a[0])
    rows_b = len(b); cols_b = len(b[0])
    if cols_a != rows_b:
      print "Cannot multiply the two matrices. Incorrect dimensions."
      return
    c=[[0 for row in range(cols_b)] for col in range(rows_a)]
    for i in range(rows_a):
        # iterate through columns of b
        for j in range(cols_b):
            # iterate through rows of b
            for k in range(cols_a):
                c[j][i] += a[k][i]*b[j][k]
    return c
#!*****************************************************************************************
def unit_vector(vector):
    return vector/np.linalg.norm(vector)
#!******************************
def angle_between(v1, v2):
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))
#!*****************************************************************************************
str1="A tool to transform POSCAR or Aims format files to ascii or acf format." 
parser=argparse.ArgumentParser(description=str1)
parser.add_argument('fn_input',action="store",type=str,help="fn_input is the name of the input file")
parser.add_argument('fn_output',action="store",type=str,help="fn_output is the name of the output file")
args=parser.parse_args()

if args.fn_input == 'POSCAR':
    atoms=poscar_read(args.fn_input)
elif 'geometry.in' in args.fn_input:
    atoms=aims_read(args.fn_input)

dproj,atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)
reportminmax(atoms.nat,atoms.rat,atoms.cellvec)

atoms.cellvec[0][0]=dproj[0]
atoms.cellvec[1][0]=dproj[1]
atoms.cellvec[1][1]=dproj[2]
atoms.cellvec[2][0]=dproj[3]
atoms.cellvec[2][1]=dproj[4]
atoms.cellvec[2][2]=dproj[5]

if "ascii" in args.fn_output:
    ascii_write(atoms,args.fn_output)
elif "acf" in args.fn_output:
    atoms_all=[]
    atoms_all.append(Atoms())
    atoms_all[-1]=copy.copy(atoms)
    acf_write(atoms_all,args.fn_output)
else:
    print 'ERROR : please define the type of output file (ascii or acf)'
