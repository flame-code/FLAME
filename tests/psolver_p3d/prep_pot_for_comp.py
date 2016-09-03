#!/usr/bin/env python
import sys
import numpy as np
import math
#*****************************************************************************************
#preparing potential cube file obtained by P3D psolver in Alborz
f1=open("pot_p3d.cube","r")
f2=open("pot_p3d.txt","w")
iline=0
for line in f1.readlines():
    iline+=1
    if iline<8: continue
    for i in line.split():
        if abs(float(i))<0.0000005:
            ii=0.0
        else:
            ii=float(i)
        f2.write("%13.6f" % ii)
    f2.write("\n")
f1.close()
f2.close()
#preparing potential cube file obtained from analatical solution
f1=open("pot_analytic.cube","r")
f2=open("pot_analytic.txt","w")
iline=0
for line in f1.readlines():
    iline+=1
    if iline<8: continue
    for i in line.split():
        if abs(float(i))<0.0000005:
            ii=0.0
        else:
            ii=float(i)
        f2.write("%13.6f" % ii)
    f2.write("\n")
f1.close()
f2.close()
#*****************************************************************************************
