#!/usr/bin/env python
import sys
import math

if len(sys.argv) < 3:
    print "usage: find_kpt.py input_filename(in acf format) dkpt"
    exit()
else:
    filename=sys.argv[1]
    dkpt=float(sys.argv[2])

f=open(filename,"r")
iline=0
for line in f.readlines():
    iline+=1
    if iline==9:
        ax=float(line.split()[0])
        bx=float(line.split()[1])
        by=float(line.split()[2])
    elif iline==10:
        cx=float(line.split()[0])
        cy=float(line.split()[1])
        cz=float(line.split()[2])
f.closed

pi=math.acos(-1.0)
a=ax
b=math.sqrt(bx**2+by**2)
c=math.sqrt(cx**2+cy**2+cz**2)
vol=ax*by*cz
g1x=2.0*pi*(by*cz)/vol  ;g1y=2.0*pi*(-bx*cz)/vol   ;g1z=2.0*pi*(bx*cy-cx*by)/vol
g2x=0.0                 ;g2y=2.0*pi*(ax*cz)/vol    ;g2z=2.0*pi*(-ax*cy)/vol
g3x=0.0                 ;g3y=0.0                   ;g3z=2.0*pi*(ax*by)/vol
g=[];kpt=[]
g.append(float(math.sqrt(g1x**2.0+g1y**2.0+g1z**2.0)))
g.append(float(math.sqrt(g2x**2.0+g2y**2.0+g2z**2.0)))
g.append(float(math.sqrt(g3x**2.0+g3y**2.0+g3z**2.0)))

for i in range(0,3):
    kpt.append(int(g[i]/(dkpt*2.0*pi)))
    if kpt[i]==0:
        kpt[i]=1
    d_test=float(g[i]/(kpt[i]*2.0*pi))
    if d_test>=dkpt:
        for j in range(1,25):
            if d_test>=dkpt:
                kpt[i]=kpt[i]+j
                d_test=float(g[i]/(kpt[i]*2.0*pi))
print "%s %6d %6d %6d" % (filename,kpt[0],kpt[1],kpt[2]) 
