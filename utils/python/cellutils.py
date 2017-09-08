#!/usr/bin/env python
import sys, copy
import numpy as np
import math as mt
from atoms import *
#**********************************************************************************#
    #!This function will convert the cartesian coordinates into the internal coordinates 
def rxyz_cart2int(cellvec,ratcart,nat):
    ratint= np.zeros((nat,3))
    cellvecinv = invertmat(cellvec)
    for iat in range(nat):
        ratint[iat][0] = cellvecinv[0][0]*ratcart[iat][0]+cellvecinv[1][0]*ratcart[iat][1]+cellvecinv[2][0]*ratcart[iat][2]
        ratint[iat][1] = cellvecinv[0][1]*ratcart[iat][0]+cellvecinv[1][1]*ratcart[iat][1]+cellvecinv[2][1]*ratcart[iat][2]
        ratint[iat][2] = cellvecinv[0][2]*ratcart[iat][0]+cellvecinv[1][2]*ratcart[iat][1]+cellvecinv[2][2]*ratcart[iat][2]
    return ratint
#********************************************************************#
def latvec2dist_ang(cellvec):
    dist_ang = [0]*6
    dist_ang_tmp = [0]*6
    dist_ang[0] = mt.sqrt(cellvec[0][0]**2+cellvec[0][1]**2+cellvec[0][2]**2)
    dist_ang[1] = mt.sqrt(cellvec[1][0]**2+cellvec[1][1]**2+cellvec[1][2]**2)
    dist_ang[2] = mt.sqrt(cellvec[2][0]**2+cellvec[2][1]**2+cellvec[2][2]**2)
    dist_ang_tmp[3] = np.dot(cellvec[1][:],cellvec[2][:])/dist_ang[1]/dist_ang[2]
    dist_ang_tmp[4] = np.dot(cellvec[2][:],cellvec[0][:])/dist_ang[2]/dist_ang[0]
    dist_ang_tmp[5] = np.dot(cellvec[0][:],cellvec[1][:])/dist_ang[0]/dist_ang[1]
    dist_ang[3]=180.0/mt.pi*mt.acos(max(min(dist_ang_tmp[3],1.0),-1.0))
    dist_ang[4]=180.0/mt.pi*mt.acos(max(min(dist_ang_tmp[4],1.0),-1.0))
    dist_ang[5]=180.0/mt.pi*mt.acos(max(min(dist_ang_tmp[5],1.0),-1.0))

    for i in range(3,6):
        if np.isnan(dist_ang[i]):
            print '%d %.10f %.10f %.10' %(i,dist_ang_tmp[i],dist_ang[i],np.acos(dist_ang[i]))
    return dist_ang
#********************************************************************#
def invertmat(mat):
    matinv = np.zeros((3,3))
    a = copy.deepcopy(mat)
    div = a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) - a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0]) + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0])
    div = 1.0 / div
    matinv[0][0] =  (a[1][1]*a[2][2] - a[1][2]*a[2][1]) *  div
    matinv[0][1] = -(a[0][1]*a[2][2] - a[0][2]*a[2][1]) *  div
    matinv[0][2] =  (a[0][1]*a[1][2] - a[0][2]*a[1][1]) *  div
    matinv[1][0] = -(a[1][0]*a[2][2] - a[1][2]*a[2][0]) *  div
    matinv[1][1] =  (a[0][0]*a[2][2] - a[0][2]*a[2][0]) *  div
    matinv[1][2] = -(a[0][0]*a[1][2] - a[0][2]*a[1][0]) *  div
    matinv[2][0] =  (a[1][0]*a[2][1] - a[1][1]*a[2][0]) *  div
    matinv[2][1] = -(a[0][0]*a[2][1] - a[0][1]*a[2][0]) *  div
    matinv[2][2] =  (a[0][0]*a[1][1] - a[0][1]*a[1][0]) *  div
    #matinv = np.matrix('1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0')
    #matinv = np.linalg.inv(mat)
    return matinv
#********************************************************************#
    #!This function will transform back all atoms into the periodic cell
    #!defined by the 3 lattice vectors in latvec=[v1.v2.v3]
def backtocell(nat,cellvec,rat):
    #!To really be on the safe side, the translation vector can be shortened by  a factor eps in order
    #!to get the atom into the cell.
    dist = [0]*6
    eps=1.0e-15
    count=0.0
    rat2 = copy.copy(rat)
    necessary= True
    while necessary:
        necessary = False
        count = count + 1.0
        # generate 3 normal vectors of the 3 planes
        nvec = nveclatvec(cellvec)
        for iat in range(nat):
           # !3 planes through origin (xy,yz,zx)
            for i in range(3):
                dist[i] = np.dot(rat[iat][:],nvec[i][:])
                if dist[i] < (-abs(dist[i])*eps):
                    k1 = (i+1)%3 + 1
                    if k1==3:
                        k1 = 0
                    rat2[iat][:] = np.array(rat[iat][:]) + np.array(cellvec[k1][:])
                    necessary = True
                k2 = (i+1)%3 + 1
                if k2==3:
                    k2 = 0
                #3 planes on top/side/back (xy,yz,zx)
                dist[i+3] = np.dot(np.array(rat[iat][:])-np.array(cellvec[k2][:]),nvec[i][:])
                if dist[i+3] > (abs(dist[i+3])*eps):
                    k3 = (i+1)%3 + 1
                    if k3 ==3:
                        k3 = 0
                    rat2[iat][:] = np.array(rat[iat][:]) - np.array(cellvec[k3][:])
                    necessary = True
        if count > 1.0e+6:
            print "Too many iterations in back-to-cell"
        else:
            break
    return rat2
#********************************************************************#
    #!Will calculate the normalized normal vector to the 3 planes of the cell
def nveclatvec(cellvec):
    nvec = np.zeros((3,3))
    for i in range(0,3):
        a = cellvec[i][:]
        j = i%3 + 1
        if j==3: 
            j=0
        b = cellvec[j][:]
        crossp = np.cross(a,b)
        norm = mt.sqrt(crossp[0]*crossp[0]+crossp[1]*crossp[1]+crossp[2]*crossp[2])
        nvec[i][:] = crossp[:]/norm
    return nvec
