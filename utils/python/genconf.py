#!/usr/bin/env python
import sys
import numpy as np
import math
import copy
import os
import random
from set_rcov import setrcov
import pattern
cwd=os.getcwd()
#path=cwd+"/../../Alborz/utils/python"
path="/home/samare/Alborz/utils/python"
sys.path.insert(1,path)
#print sys.path
from atoms import *
from xyz import *
from ascii import *
from acf import *
#*****************************************************************************************
# reading input via arguments:------------------------------------------------------------
#if len(sys.argv) < 5:
#    print "usage: gen_conf.py nat nat_molecule ['typ_1','typ_2',...,'typ_natmolecule'] dboud filename -bc boundcond"
#    exit()
#else:
#    nat=int(sys.argv[1])
#    nat_molecule=int(sys.argv[2])
#    typ_n=[]
#    typ_n.append(int(sys.argv[3]))
#    print typ_n[0]
#    dbond=float(sys.argv[4])
#    filename=sys.argv[5]
#    boundcond='unknown'
#    skip=False
#    for i in range(1,len(sys.argv)):
#        #print i,sys.argv[i]
#        if skip:
#            skip=False
#            continue
#        if sys.argv[i]=='-bc':
#            if i+1<len(sys.argv):
#                boundcond=sys.argv[i+1]
#                skip=True
#                #print boundcond
#            else:
#                sys.exit("ERROR: improper argument list")
#********************************************************************************************i
def Max_val(a,b):
    if a > b:
        return a
    else:
        return b
#*********************************************************************************************
# Reading input via input.conf file
if len(sys.argv) < 2:
    print "usage: genconf.py output_filename.acf"
    exit()
else:
    filename = sys.argv[1]

f=open("input.conf","r")
atoms_all=[]
iline=0
for line in f.readlines():
    iline+=1
    if iline==1:
        atoms=Atoms()
        m=Atoms()
        comment1=line
    elif iline==2:
        atoms.nmolecule=int(line.split()[0])
        m.nat=int(line.split()[1])
        atoms.nat=int(atoms.nmolecule*m.nat)
    elif iline==3:
        comment3=line
    elif iline==4:
        for i in line.split():
            m.sat.append(i)
        print m.sat
        m.sat=[] 
    elif iline==5:
        comment4=line
    elif iline==6:
        atoms.pattern=int(line.split()[0])
    elif iline==9:
        comment5=line
    elif iline>9 and iline<m.nat+10:
        m.rat.append([])
        icol=0
        #loop over the elemets, split by whitespace
        for i in line.split():
            icol+=1
            if icol==1:
                m.sat.append(i)
            elif icol<5:
                m.rat[-1].append(float(i))
    elif iline==atoms.nat+10:
        iline=0
        atoms_all.append(Atoms())
        print atoms_all[-1].nat
        atoms_all[-1]=copy.copy(m)
        #print atoms_all[-1].nat
        print len(atoms_all[-1].rat)
        #print
        #atoms.kill()
        #print len(atoms_all[-1].rat)
        #print
for iat in range(m.nat):
    print m.rat[iat][0]
    #print m.rat[iat][1]
f.closed
print atoms.nmolecule
#return atoms_all
#********************************************************************************************
atoms.rat.append(m.rat)
atoms.sat.append(m.sat)
#atoms.pattern=random.randrange(1,3)
for im in range(atoms.nmolecule):
    x_cm=y_cm=z_cm=0.0
    d=[0.0, 0.0, 0.0]
    for iat in range(m.nat):
        print m.sat[iat], setrcov(m.sat[iat])
        #comparing with "rcov"s the following line is determined 
        rcov_max=setrcov(m.sat[0])
        x_cm+=m.rat[iat][0]/m.nat  
        y_cm+=m.rat[iat][1]/m.nat
        z_cm+=m.rat[iat][2]/m.nat
        d[iat]=math.sqrt((m.rat[iat][0]-x_cm)**2+(m.rat[iat][1]-y_cm)**2+(m.rat[iat][2]-z_cm)**2)
        print d[iat]
        m.rat.append([])
        m.bemoved.append("TTT")
    #print x_cm, y_cm, z_cm
    pattern.gen_pattern(x_cm,y_cm,z_cm,rcov_max,atoms.nmolecule,1) 
    #print x_cm, y_cm, z_cm
    for iat in range(m.nat):
        m.rat[-1-iat].append(x_cm+d[iat])
        m.rat[-1-iat].append(y_cm+d[iat])
        m.rat[-1-iat].append(z_cm+d[iat])
        print m.rat[iat+1][0]
        m.bemoved.append("TTT")
    #print m.rat[2][0]
atoms.rat.append(m.rat)
atoms.bemoved.append(m.bemoved)
atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(m)
print "writing to file %20s" % filename
#xyz_write(atoms_all,'bigdft')
#ascii_write(atoms,filename)
acf_write(atoms_all,filename)
###*****************************************************************************************
