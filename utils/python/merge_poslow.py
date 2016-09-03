#!/usr/bin/python
import os
import sys
import math
import atoms
from acf import *
#*****************************************************************************************
def test_identical(atoms1,atoms2,etol,dist,dtol,i1=-1,i2=-1):
    #print "%d" % len(dist)
    if abs(atoms1.epot-atoms2.epot)<etol:
        same=True
    else:
        same=False
    if same and len(dist)>0:
        #print "sadfasdfdasf"
        if i1==-1 or i2==-1: print "ERROR: i1,i2= %7d%7d" % (i1,i2)
        #print "DIST %19.10E" % dist[i1][i2]
        if dist[i1][i2]>dtol:
            same=False
    return same
#*****************************************************************************************
def merge_poslows(poslow,poslow_all,etol,dist=[],dtol=-1.0):
    #if len(poslow_all)==0:
    #    for atoms in poslow:
    #        poslow_all.append(atoms)
    #    return poslow_all
    poslow_t=[]
    ind_poslow_t=[]
    n1=len(poslow)
    n2=len(poslow_all)
    i1=0
    i2=0
    for i in range(n1+n2):
        print i,i1,i2,n1,n2 #,poslow[i1].epot,poslow_all[i2].epot
        if not i1==n1 and (i2==n2 or poslow[i1].epot<poslow_all[i2].epot):
            #if i==0 or not abs(poslow[i1].epot-poslow_t[-1].epot)<etol:
            if i!=0: same=test_identical(poslow[i1],poslow_t[-1],etol,dist,dtol,i1,ind_poslow_t[-1])
            if i==0 or not same:
                poslow_t.append(poslow[i1])
                ind_poslow_t.append(i1)
                #if len(dist)>0: ind_poslow_t.append(i1)
            i1+=1
        else:
            #if i==0 or not abs(poslow_all[i2].epot-poslow_t[-1].epot)<etol:
            if i!=0: same=test_identical(poslow_all[i2],poslow_t[-1],etol,dist,dtol)
            if i==0 or not same:
                poslow_t.append(poslow_all[i2])
                #next line is just to avoid error message "index out of range" when switching to
                #the other if condition in next iterations.
                ind_poslow_t.append(i2)
            i2+=1
        #print i,i1,i2,n1,n2 #,poslow[i1].epot,poslow_all[i2].epot
    return poslow_t
#*****************************************************************************************
def read_distall():
    dist=[]
    f=open("distall","r")
    nn=-1
    nconf=0
    for line in f.readlines():
        if not int(line.split()[1])==nn:
            nconf+=1
            dist.append([])
            for iconf in range(nconf-1):
                dist[-1].append(dist[iconf][nconf-1])
            dist[-1].append(0.0)
            nn=int(line.split()[1])
        dist[-1].append(float(line.split()[4]))
    f.close()
    nconf+=1
    dist.append([])
    for iconf in range(nconf-1):
        dist[-1].append(dist[iconf][nconf-1])
    dist[-1].append(0.0)
    #for i in range(4):
    #    for j in range(4):
    #        print "%19.10E" % dist[i][j],
    #    print
    return dist
#*****************************************************************************************
#*****************************************************************************************
if len(sys.argv) < 3:
    print "usage: merge_poslow.py output_filename tolerance [-dist dtol]" 
    exit() 
#print "%s" % sys.argv[1]
#print "%s" % sys.argv[2]
#exit()
if len(sys.argv)>=3:
    output_file = sys.argv[1]
    etol=float(sys.argv[2])
    #print "%s" % sys.argv[1]
    #print "%s" % sys.argv[2]
if len(sys.argv)==5:
    if sys.argv[3]=='-dist':
        #filename=sys.argv[2]
        dtol=float(sys.argv[4])
if len(sys.argv)==4 or len(sys.argv)>5:
    print "ERROR: wrong command option."
    exit()

etol*=27.211385
poslow=[]
poslow_all=[]
f=open('list_poslow','r')
nfiles=0
for line in f.readlines():
    nfiles+=1
    if len(sys.argv)==5:
        if nfiles>1:
            print "ERROR: when fingerprints are used, only one file is allowed in list_poslow"
            exit()
        dist=read_distall()
    filename=line.split()[0]
    poslow=acf_read(filename)
    if len(sys.argv)==5:
        poslow_all=merge_poslows(poslow,poslow_all,etol,dist,dtol)
    else:
        poslow_all=merge_poslows(poslow,poslow_all,etol)
f.close()

acf_write(poslow_all,output_file)
#*****************************************************************************************
