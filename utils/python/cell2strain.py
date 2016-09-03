#!/usr/bin/env python
import numpy as np
from numpy import linalg as LA
import sys

if len(sys.argv) < 2:
    print "usage: cell2strain.py filename"
    print "WARNING: Notice that the components of reference cell"
    print "         vectors must be in 1st line."
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
        #reference cell vector
        href=np.array([[arr[0],0.0,0.0],[arr[1],arr[2],0.0],[arr[3],arr[4],arr[5]]])
        hrefinv=LA.inv(href)
    else:
        #deformed cell vectors
        hdef=np.array([[arr[0],0.0,0.0],[arr[1],arr[2],0.0],[arr[3],arr[4],arr[5]]])
        e=np.dot(hrefinv,hdef-href)
        #e=np.dot(hdef-href,hrefinv)
        e[1][0]/=2.0
        e[2][0]/=2.0
        e[2][1]/=2.0
        print "%14.8f%14.8f%14.8f%14.8f%14.8f%14.8f" % (e[0][0],e[1][0],e[1][1],e[2][0],e[2][1],e[2][2])

f.closed
