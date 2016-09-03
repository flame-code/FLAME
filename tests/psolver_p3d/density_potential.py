#!/usr/bin/env python
import sys
import numpy as np
import math
import os
cwd=os.getcwd()
path=cwd+"/../../utils/python"
sys.path.insert(1,path)
from atoms import *
from cube import *
#*****************************************************************************************
def get_potential(cell,ngpx,ngpy,ngpz,hx,hy,hz,a,b,c):
    #pot = [[[] for j in range()] for i in range(3)]
    pot=np.empty([ngpz,ngpy,ngpx])
    pi=4.0*math.atan(1.0)
    z0=cell[2]/2.0
    for igpz in range(ngpz):
        z=(igpz-0)*hz
        tz=math.exp(-(z-z0)**2/c**2)
        for igpy in range(ngpy):
            y=(igpy-0)*hy
            tyz=tz*math.sin(b*math.sin(2.0*pi*y/cell[1]))
            for igpx in range(ngpx):
                x=(igpx-0)*hx
                pot[igpz][igpy][igpx]=tyz*math.sin(a*math.sin(2.0*pi*x/cell[0]))
    return pot
#*****************************************************************************************
def get_density(cell,ngpx,ngpy,ngpz,hx,hy,hz,a,b,c):
    rho=np.empty([ngpz,ngpy,ngpx])
    pi=4.0*math.atan(1.0)
    z0=cell[2]/2.0
    t1=1.0/(2.0*c**4*cell[0]**2*cell[1]**2*pi)
    for igpz in range(ngpz):
        z=(igpz-0)*hz
        for igpy in range(ngpy):
            y=(igpy-0)*hy
            for igpx in range(ngpx):
                x=(igpx-0)*hx
                t2=math.exp(-(z - z0)**2/c**2)
                t3=(2*b*c**4*cell[0]**2*pi**2*math.cos(b*math.sin((2*pi*y)/cell[1]))*
                    math.sin((2*pi*y)/cell[1])*math.sin(a*math.sin((2*pi*x)/cell[0])) +
                    (2*a*c**4*cell[1]**2*pi**2*math.cos(a*math.sin((2*pi*x)/cell[0]))*
                    math.sin((2*pi*x)/cell[0]) + (cell[1]**2*(cell[0]**2*(c**2 - 2*(z - z0)**2) +
                    a**2*c**4*pi**2*(1 + math.cos((4*pi*x)/cell[0]))) +
                    b**2*c**4*cell[0]**2*pi**2*(1 + math.cos((4*pi*y)/cell[1])))*
                    math.sin(a*math.sin((2*pi*x)/cell[0])))*math.sin(b*math.sin((2*pi*y)/cell[1])))
                rho[igpz][igpy][igpx]=t1*t2*t3
    return rho
#*****************************************************************************************

if len(sys.argv) < 7:
    print "usage: density_potential.py cellx celly cellz ngpx ngpy ngpz"
    exit()
else:
    cell=[]
    cell.append(float(sys.argv[1]))
    cell.append(float(sys.argv[2]))
    cell.append(float(sys.argv[3]))
    ngpx=int(sys.argv[4])
    ngpy=int(sys.argv[5])
    ngpz=int(sys.argv[6])

hx=cell[0]/ngpx
hy=cell[1]/ngpy
hz=cell[2]/ngpz
a=1.0
b=1.0
c=1.0
rho=get_density(cell,ngpx,ngpy,ngpz,hx,hy,hz,a,b,c)
pot=get_potential(cell,ngpx,ngpy,ngpz,hx,hy,hz,a,b,c)
atoms=Atoms()
atoms.nat=1
atoms.sat.append('H')
atoms.rat.append([])
atoms.rat[-1].append(float(cell[0]/2.0))
atoms.rat[-1].append(float(cell[1]/2.0))
atoms.rat[-1].append(float(cell[2]/2.0))
atoms.qat.append(float(0.0))
frmt1="%5d%24.15E%24.15E%24.15E\n"
frmt2=" %22.15f "
cube_write('rho.cube',atoms,ngpx,ngpy,ngpz,rho,hx,hy,hz,frmt1,frmt2)
frm1="%5d%13.6f%13.6f%13.6f\n"
#frmt2=" %13.6f "
frmt2=" %12.5E "
cube_write('pot_analytic.cube',atoms,ngpx,ngpy,ngpz,pot,hx,hy,hz,frmt1,frmt2)
#*****************************************************************************************
