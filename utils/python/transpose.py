#!/usr/bin/env python
import sys

if len(sys.argv) < 2:
    print "usage: transpose.py input_filename"
    exit()
else:
    filename = sys.argv[1]

nrow=0
a=[]
f=open(filename,"r")
#read line into a 
for line in f.readlines():
    a.append([])
    ncol=0
    for i in line.split():
        a[nrow].append(i)
        ncol+=1
    nrow+=1
f.closed

at=map(list,zip(*a))
#at=[]
#transpose a to at
#for i in range(len(a)):
#    for j in range(len(a[i])):
#        at

for row in at:
    str_out=""
    for col in row:
        str_out+= "%7s" % col
    print str_out
