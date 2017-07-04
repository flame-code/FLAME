#!/usr/bin/env python
#This program will freez those atoms which are inside the region which
# is inside the structure and about distance equal with "margin" below the outer surface.
import sys
import atoms
import numpy as np
from acf import *

if len(sys.argv) < 2:
    print "usage: acf_fixing_atoms.py input_filename margin"
    exit()
else:
    filename = sys.argv[1]
    if len(sys.argv)==2:
        margin=float(sys.argv[2])
    else:
        margin=3.0
atoms_all=acf_read(filename)
nconf=-1
#-----Finding the minimum and maximum of coordinates:---------------------------------------
for atoms in atoms_all:
    nconf+=1
    columns =zip (*atoms_all[nconf].rat)
    minimum= []
    maximum= []
    for i, co in enumerate (columns):
        minimum.append(min(co))
        maximum.append(max(co))
    minratx = minimum[0] ; minraty = minimum[1] ; minratz = minimum[2]
    maxratx = maximum[0] ; maxraty = maximum[1] ; maxratz = maximum[2]
#-----Determining the borders of region in which the atoms must be freezed:------------------
    lmin_x = minratx+margin ; lmin_y = minraty+margin ; lmin_z = minratz+margin
    lmax_x = maxratx-margin ; lmax_y = maxraty-margin ; lmax_z = maxratz-margin
#-----Fixing process:------------------------------------------------------------------------
    if atoms_all[nconf].boundcond=='free':
        for i in range(int(atoms_all[nconf].nat)):
            if (atoms_all[nconf].rat[i][0] >lmin_x and atoms_all[nconf].rat[i][0] <lmax_x):
                if (atoms_all[nconf].rat[i][1] >lmin_y and atoms_all[nconf].rat[i][1] <lmax_y):
                    if (atoms_all[nconf].rat[i][2] >lmin_z and atoms_all[nconf].rat[i][2] <lmax_z):
                        atoms_all[nconf].bemoved[i]="FFF"
    elif atoms_all[nconf].boundcond=='wire':
        for i in range(int(atoms_all[nconf].nat)):
            if (atoms_all[nconf].rat[i][1] >lmin_y and atoms_all[nconf].rat[i][1] <lmax_y):
                if (atoms_all[nconf].rat[i][2] >lmin_z and atoms_all[nconf].rat[i][2] <lmax_z):
                    atoms_all[nconf].bemoved[i]="FFF"
    elif atoms_all[nconf].boundcond=='slab':
        for i in range(int(atoms_all[nconf].nat)):
            if (atoms_all[nconf].rat[i][2] >lmin_x and atoms_all[nconf].rat[i][2] <lmax_x):
                atoms_all[nconf].bemoved[i]="FFF"
    else :
        print "ERROR: unknown boundary condition"

acf_write(atoms_all,"screen")
