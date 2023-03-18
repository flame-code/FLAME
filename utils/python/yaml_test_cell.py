#!/usr/bin/env python
import sys
import atoms
import math
import numpy as np
from acf import *
from io_yaml import *
from cellutils import *

if len(sys.argv) < 3:
    print("")
    print("usage: yaml_test_cell.py input_filename output_filename")
    print("")
    exit()
else:
    filename = sys.argv[1]
    ofilename=sys.argv[2]
    #dmin=float(sys.argv[2])
    #if len(sys.argv) >= 4:
    #    ofilename=sys.argv[3]
    #else :
    #    ofilename="screen"

if len(sys.argv) == 5:
    dmax=float(sys.argv[4])
else :
    dmax = 5
atoms_all=read_yaml(filename)
atoms_all_sel =[]

#for i in range(-10,11):
#    tt1=0.1*i
#    print "%10.0f" % (math.acos(tt1)/math.pi*180.0)

ratiotol=0.1
angmintol=30.0

iconf=-1
nconf_excl=0
for atoms in atoms_all:
    iconf+=1
    bad_conf=False
    #print atoms.cellvec[0][0],atoms.cellvec[0][1],atoms.cellvec[0][2]
    #print atoms.cellvec[1][0],atoms.cellvec[1][1],atoms.cellvec[1][2]
    #print atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2]
    #print "---------"
    a=(atoms.cellvec[0][0]**2+atoms.cellvec[0][1]**2+atoms.cellvec[0][2]**2)**0.5
    b=(atoms.cellvec[1][0]**2+atoms.cellvec[1][1]**2+atoms.cellvec[1][2]**2)**0.5
    c=(atoms.cellvec[2][0]**2+atoms.cellvec[2][1]**2+atoms.cellvec[2][2]**2)**0.5
    ratio=1.E10
    if a/b<ratio: ratio=a/b
    if b/a<ratio: ratio=b/a
    if a/c<ratio: ratio=a/c
    if c/a<ratio: ratio=c/a
    if c/b<ratio: ratio=c/b
    if b/c<ratio: ratio=b/c
    if ratio<ratiotol: bad_conf=True
    #print ratio
    xvec=np.zeros(3) ; xvec[0]=1.0
    yvec=np.zeros(3) ; yvec[1]=1.0
    zvec=np.zeros(3) ; zvec[2]=1.0
    cv=np.zeros((3,3))
    cv[0][0]=atoms.cellvec[0][0] ; cv[0][1]=atoms.cellvec[0][1] ; cv[0][2]=atoms.cellvec[0][2]
    cv[1][0]=atoms.cellvec[1][0] ; cv[1][1]=atoms.cellvec[1][1] ; cv[1][2]=atoms.cellvec[1][2]
    cv[2][0]=atoms.cellvec[2][0] ; cv[2][1]=atoms.cellvec[2][1] ; cv[2][2]=atoms.cellvec[2][2]
    #axe=np.zeros(3)
    #cvp=np.cross(xvec,yvec)
    cvp=np.zeros((3,3))
    cvp[2][:]=np.cross(cv[0][:],cv[1][:])
    cvp[0][:]=np.cross(cv[1][:],cv[2][:])
    cvp[1][:]=np.cross(cv[2][:],cv[0][:])
    tt1=np.dot(cv[0][:],cvp[0][:])
    tt2=(np.dot(cv[0][:],cv[0][:]))**0.5
    tt3=(np.dot(cvp[0][:],cvp[0][:]))**0.5
    ang_a=math.acos(tt1/(tt2*tt3))/math.pi*180.0
    tt1=np.dot(cv[1][:],cvp[1][:])
    tt2=(np.dot(cv[1][:],cv[1][:]))**0.5
    tt3=(np.dot(cvp[1][:],cvp[1][:]))**0.5
    ang_b=math.acos(tt1/(tt2*tt3))/math.pi*180.0
    tt1=np.dot(cv[2][:],cvp[2][:])
    tt2=(np.dot(cv[2][:],cv[2][:]))**0.5
    tt3=(np.dot(cvp[2][:],cvp[2][:]))**0.5
    ang_c=math.acos(tt1/(tt2*tt3))/math.pi*180.0
    if ang_a>90.0: ang_a=180.0-ang_a
    if ang_b>90.0: ang_b=180.0-ang_b
    if ang_c>90.0: ang_c=180.0-ang_c
    ang_a=90.0-ang_a
    ang_b=90.0-ang_b
    ang_c=90.0-ang_c
    #print ang_a,ang_b,ang_c
    if ang_a<angmintol or ang_b<angmintol or ang_c<angmintol:
        print("%8d%5.0f%5.0f%5.0f  F" % (iconf,ang_a,ang_b,ang_c))
    else:
        print("%8d%5.0f%5.0f%5.0f  T" % (iconf,ang_a,ang_b,ang_c))
    if ang_a<angmintol or ang_b<angmintol or ang_c<angmintol: bad_conf=True
    if not bad_conf:
        atoms_all_sel.append(Atoms())
        atoms_all_sel[-1]=copy.deepcopy(atoms)
        nconf_excl+=1

nconf=iconf+1
if nconf!=len(atoms_all): print("ERROR: nconf!=len(atoms_all)")

if len(atoms_all_sel)>0: write_yaml(atoms_all_sel,ofilename)
print("total number of conf = ",nconf,"       the number of configurations deleted = ",nconf-nconf_excl)
