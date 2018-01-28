#!/usr/bin/env python
import sys
import atoms
import numpy as np
import math
from acf import *

rotaxis=[]
if len(sys.argv) < 5:
    print ""
    print "usage: acf_rotate.py input_filename axis(x y z) angle(degrees)"
    print ""
    print "example: acf_rotate.py posinp.acf 1 0 0 45 "
    print ""
    exit()
else:
    filename = sys.argv[1]
    rotaxis.append(float(sys.argv[2]))
    rotaxis.append(float(sys.argv[3]))
    rotaxis.append(float(sys.argv[4]))
    rotangle = float(sys.argv[5])
rotangle = rotangle/180*math.pi 
t1=math.sqrt(rotaxis[0]**2 + rotaxis[1]**2 + rotaxis[2]**2)
rotaxis[0]=rotaxis[0]/t1
rotaxis[1]=rotaxis[1]/t1
rotaxis[2]=rotaxis[2]/t1
o=(3,3)
one=np.zeros(o)
pp= np.zeros(o)   
s1= np.zeros(o) 
s2= np.zeros(o) 
s3= np.zeros(o) 
a=  np.zeros(o)

one[0][0]=1.0
one[1][1]=1.0
one[2][2]=1.0
s1[1][2]= 1.0
s1[2][1]=-1.0
s2[0][2]=-1.0
s2[2][0]= 1.0
s3[0][1]= 1.0
s3[1][0]=-1.0

for i in range (3):
    for j in range (3):
        pp[i][j]=rotaxis[i]*rotaxis[j]
        
a=one
for i in range (3):
    for j in range (3):
        a[i][j]=a[i][j]+(1.0-math.cos(rotangle))*(pp[i][j]-one[i][j])
        a[i][j]=a[i][j]-rotaxis[0]*math.sin(rotangle)*s1[i][j]
        a[i][j]=a[i][j]-rotaxis[1]*math.sin(rotangle)*s2[i][j]
        a[i][j]=a[i][j]-rotaxis[2]*math.sin(rotangle)*s3[i][j]

atoms_all=acf_read(filename)
nconf=-1
sumx=0 
sumy=0
sumz=0
for atoms in atoms_all:
    nconf+=1
    for i in range(int(atoms_all[nconf].nat)):
        sumx += atoms_all[nconf].rat[i][0]
        sumy += atoms_all[nconf].rat[i][1]
        sumz += atoms_all[nconf].rat[i][2]
    cm_x = sumx / (atoms_all[nconf].nat)
    cm_y = sumy / (atoms_all[nconf].nat)
    cm_z = sumz / (atoms_all[nconf].nat)
#    print "cm_x =%20.15f    cm_y =%20.15f   cm_z =%20.15f  " %(cm_x ,cm_y,cm_z)
    rat2=np.array(atoms_all[nconf].rat)
    for i in range(int(atoms_all[nconf].nat)):
        rat2[i][0] = rat2[i][0]-cm_x
        rat2[i][1] = rat2[i][1]-cm_y
        rat2[i][2] = rat2[i][2]-cm_z
    for i in range(int(atoms_all[nconf].nat)):
        atoms_all[nconf].rat[i][0]=a[0][0]*rat2[i][0] + a[0][1]*rat2[i][1] + a[0][2]*rat2[i][2] +cm_x
        atoms_all[nconf].rat[i][1]=a[1][0]*rat2[i][0] + a[1][1]*rat2[i][1] + a[1][2]*rat2[i][2] +cm_y
        atoms_all[nconf].rat[i][2]=a[2][0]*rat2[i][0] + a[2][1]*rat2[i][1] + a[2][2]*rat2[i][2] +cm_z
acf_write(atoms_all,"screen")

