#!/usr/bin/env python
#
#
from pylab import *
import pylab as pl
from numpy.linalg import *
from numpy import dot,cross,pi
import numpy as np
import itertools
import os,sys

########################

print "============================"
print

print "Reading lattice vectors from geometry.in ..."

latvec = []
pl.rc('font', family='serif')
for line in file("geometry.in"):
    line = line.split("#")[0]
    words = line.split()
    if len(words) == 0:
        continue
    if words[0] == "lattice_vector":
        if len(words) != 4:
            raise Exception("geometry.in: Syntax error in line '"+line+"'")
        latvec += [ array(map(float,words[1:4])) ]

if len(latvec) != 3:
    raise Exception("geometry.in: Must contain exactly 3 lattice vectors")

latvec = asarray(latvec)

print "Lattice vectors:"
for i in range(3):
    print latvec[i,:]
print

#Calculate reciprocal lattice vectors                                                                                                
rlatvec = []
volume = (np.dot(latvec[0,:],np.cross(latvec[1,:],latvec[2,:])))
rlatvec.append(array(2*pi*cross(latvec[1,:],latvec[2,:])/volume))
rlatvec.append(array(2*pi*cross(latvec[2,:],latvec[0,:])/volume))
rlatvec.append(array(2*pi*cross(latvec[0,:],latvec[1,:])/volume))
rlatvec = asarray(rlatvec)

#rlatvec = inv(latvec) Old way to calculate lattice vectors
print "Reciprocal lattice vectors:"
for i in range(3):
    print rlatvec[i,:]
print

########################

print "Reading information from control.in ..."
species = []

max_spin_channel = 1
band_segments = []
band_totlength = 0.0 # total length of all band segments

for line in file("control.in"):
    words = line.split("#")[0].split()
    nline = " ".join(words)

    if nline.startswith("spin collinear"):
        max_spin_channel = 2

    if nline.startswith("output band "):
        if len(words) < 9 or len(words) > 11:
            raise Exception("control.in: Syntax error in line '"+line+"'")
        PLOT_BANDS = True
        start = array(map(float,words[2:5]))
        end = array(map(float,words[5:8]))
        length = norm(dot(rlatvec,end) - dot(rlatvec,start))
        band_totlength += length
        npoint = int(words[8])
        startname = ""
        endname = ""
        if len(words)>9:
            startname = words[9]
        if len(words)>10:
            endname = words[10]
        band_segments += [ (start,end,length,npoint,startname,endname) ]

#######################
prev_end = band_segments[0][0]
distance = band_totlength/30.0 # distance between line segments that do not coincide

iband = 0
xpos = 0.0
labels = [ (0.0,band_segments[0][4]) ]
merged2 = []
for start,end,length,npoint,startname,endname in band_segments:
    iband += 1

    if any(start != prev_end):
        xpos += distance
        labels += [ (xpos,startname) ]

    xvals = xpos+linspace(0,length,npoint)
    xpos = xvals[-1]

    labels += [ (xpos,endname) ]
    
    prev_end = end
    prev_endname = endname
    
    for spin in range(1,max_spin_channel+1):
        fname = "band%i%03i.out"%(spin,iband)
        idx = []
        kvec = []
        band_energies = []
        band_occupations = []
        for line in file(fname):
            words = line.split()
            idx += [ int(words[0]) ]
            kvec += [ map(float,words[1:4]) ]
            band_occupations += [ map(float,words[4::2]) ]
            band_energies += [ map(float,words[5::2]) ]
        assert(npoint) == len(idx)
        merged1 = list(itertools.chain(*band_energies))
        band_energies = asarray(band_energies)


    merged2.append(merged1)    
merged = list(itertools.chain(*merged2))
merged.sort()
valb = max([n for n in merged if n<0]) 
conb = min([n for n in merged if n>0])
egap = conb - valb
print
print "the top of the valence band :", valb
print "the bottom of the conduction band:", conb
print "The band gap of this structure is equal to : %10.3f (eV)" % egap
print "============================"
