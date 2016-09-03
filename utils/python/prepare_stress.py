#!/usr/bin/env python
import numpy as np
from numpy import linalg as LA
import sys

if len(sys.argv) < 2:
    print "usage: prepare_stress.py filename"
    print "WARNING: Notice that the components of reference "
    print "         stress tensor must be in 1st line."
    exit()
else:
    filename = sys.argv[1]

f = open (filename,"r")
n=0
#read line into array 
for line in f.readlines():
    n+=1
    arr = []
    #add a new sublist
    #arr.append([])
    #loop over the elemets, split by whitespace
    for i in line.split():
        #convert to integer and append to the last element of the list
        arr.append(float(i))
    if n==1:
        sref=np.array([[arr[0],0.0,0.0],[arr[1],arr[2],0.0],[arr[3],arr[4],arr[5]]])
        continue
    else:
        sdef=np.array([[arr[0],0.0,0.0],[arr[1],arr[2],0.0],[arr[3],arr[4],arr[5]]])
        s=(sdef-sref)*0.1
        #for i in range(3):
        #    for j in range(3):
        #        if abs(s[i][j])<1.E-3:
        #            s[i][j]=0.0
        #print "%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f" % (s[0][0],s[1][0],s[1][1],s[2][0],s[2][1],s[2][2])
        #print "%9.3f%9.3f%9.3f%9.3f%9.3f%9.3f" % (s[0][0],s[1][0],s[1][1],s[2][0],s[2][1],s[2][2])
        print "%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f" % (s[0][0],s[1][0],s[1][1],s[2][0],s[2][1],s[2][2])

f.closed
