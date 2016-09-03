#!/usr/bin/env python
import sys
import atoms
import numpy as np
from acf import *

if len(sys.argv) < 2:
    print "usage: set_cm.py input_filename cm_x cm_y cm_z "
    exit()
elif len(sys.argv) < 5:
    print "usage: cm is set to the center of box "
    filename = sys.argv[1]
else:
    filename = sys.argv[1]
    cm_x_ref = sys.argv[2]
    cm_y_ref = sys.argv[3]
    cm_z_ref = sys.argv[4]
if len(sys.argv) == 3:
    filename2 = sys.argv[2]
elif len(sys.argv) == 6:
    filename2 = sys.argv[5]
else :
    filename2 ="screen"
atoms_all=acf_read(filename)
nconf=-1
for atoms in atoms_all:
    nconf+=1
    if len(sys.argv) < 5:
        cm_x_ref = atoms_all[nconf].cellvec[0][0]/2.0
        cm_y_ref = atoms_all[nconf].cellvec[1][1]/2.0
        cm_z_ref = atoms_all[nconf].cellvec[2][2]/2.0
        print "cm_x =%20.15f    cm_y =%20.15f   cm_z =%20.15f  " %(cm_x_ref ,cm_y_ref,cm_z_ref)
    sumx = 0
    sumy = 0
    sumz = 0
    for i in range(int(atoms_all[nconf].nat)):
        sumx += atoms_all[nconf].rat[i][0]
        sumy += atoms_all[nconf].rat[i][1]
        sumz += atoms_all[nconf].rat[i][2]
    cm_x = sumx / (atoms_all[nconf].nat)
    cm_y = sumy / (atoms_all[nconf].nat)
    cm_z = sumz / (atoms_all[nconf].nat)
    
    diff_x = cm_x_ref - cm_x
    diff_y = cm_y_ref - cm_y
    diff_z = cm_z_ref - cm_z

    for i in range(int(atoms_all[nconf].nat)):
        atoms_all[nconf].rat[i][0]+=diff_x
        atoms_all[nconf].rat[i][1]+=diff_y
        atoms_all[nconf].rat[i][2]+=diff_z

acf_write(atoms_all,filename2)

