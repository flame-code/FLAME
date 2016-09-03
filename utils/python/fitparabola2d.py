#!/usr/bin/env python
import numpy as numpy
from numpy import linalg as LA
import sys

if len(sys.argv) < 2:
    print "usage: fitparabola2d.py input_filename"
    exit()
else:
    filename = sys.argv[1]

if len(sys.argv) == 3:
    zmaxtol=float(sys.argv[2])
else:
    zmaxtol=1.E100

arr = []
f = open (filename,"r")
#read line into array 
for line in f.readlines():
    #add a new sublist
    arr.append([])
    #loop over the elemets, split by whitespace
    for i in line.split():
        #convert to integer and append to the last element of the list
        arr[-1].append(float(i))

f.closed

npt=len(arr)
print "\nThere are %d data points in file %s.\n" % (npt,filename)

xt = numpy.array([row[:][0] for row in arr])
yt = numpy.array([row[:][1] for row in arr])
zt = numpy.array([row[:][2] for row in arr])

xtt=[]
ytt=[]
ztt=[]
xmin=1.E50 ; xmax=-1.E50
ymin=1.E50 ; ymax=-1.E50
zmin=1.E50 ; zmax=-1.E50
for i in range(npt):
    if xt[i] < xmin: xmin=xt[i]
    if xt[i] > xmax: xmax=xt[i]
    if yt[i] < ymin: ymin=yt[i]
    if yt[i] > ymax: ymax=yt[i]
    if zt[i] < zmin: zmin=zt[i]
    if zt[i] > zmax: zmax=zt[i]
    if not zt[i] > zmaxtol:
        xtt.append(xt[i])
        ytt.append(yt[i])
        ztt.append(zt[i])

print "data bounds:"
print "xmin,xmax %15.5f%15.5f" % (xmin,xmax)
print "ymin,ymax %15.5f%15.5f" % (ymin,ymax)
print "zmin,zmax %15.5f%15.5f" % (zmin,zmax)

x=numpy.array(xtt)
y=numpy.array(ytt)
z=numpy.array(ztt)
np=len(x)
print "\nThere are %d data points after filtering.\n" % np
#print "len(x) %d %s" % (len(x),zmaxtol)
#exit()

var=["numpy.ones(np)","x","y","x**2","x*y","y**2"]
v = numpy.array([eval(i) for i in var])

coefficients, residues, rank, singval = numpy.linalg.lstsq(v.T, z)

print "Fitting process completed."
print "f(x,y)=(%s) + (%s) x + (%s) y + " % (coefficients[0], coefficients[1], coefficients[2])
print "       (%s) x^2 + (%s) xy + (%s) y^2\n" % (coefficients[3], coefficients[4], coefficients[5])

print "residues= ", residues[0], "\n"
print "rank= ", rank, "\n"
#print singval

for i in range(np):
    tt=0.0
    k=0
    for j in var:
        tt+=coefficients[k]*eval(j)[i]
        k+=1
        #print "ERR= %s" % eval(j)[i]
    tt=tt-z[i]
    #print "ERR= %10.4f" % tt

#res_my=0.0
#for i in range(np):
#    tt=coefficients[0]+coefficients[1]*x[i]+coefficients[2]*y[i]+coefficients[3]*x[i]**2+ \
#        coefficients[4]*x[i]*y[i]+coefficients[5]*y[i]**2 - z[i]
#    res_my+=tt**2
#
#print "MY: residue: %f\n" % res_my

print "--------------------------------------------------"
a = [[2*coefficients[3], coefficients[4]], [coefficients[4], 2*coefficients[5]]]
#print a
b = numpy.array(a)
w, v = LA.eig(b)
print "Eigenvalues: ", w
print "Eigenvectors: ", v
print "--------------------------------------------------"
print "To plot in gnuplot use the following: \n"
print "f(x,y)=c0_xy + c1_x * x + c1_y * y + c2_x *x*x + c1_xy *x*y + c2_y * y*y"
print "c0_xy=%20.10f" % float(coefficients[0])
print "c1_x =%20.10f" % float(coefficients[1])
print "c1_y =%20.10f" % float(coefficients[2])
print "c2_x =%20.10f" % float(coefficients[3])
print "c1_xy=%20.10f" % float(coefficients[4])
print "c2_y =%20.10f" % float(coefficients[5])
print "%7s%s%18s" % ("splot '",filename, "' u 1:2:3 , f(x,y)")
print "--------------------------------------------------"
