import sys
import argparse
import atoms
import copy
from atoms import *
#*****************************************************************************************
def aims_write(atoms,filename):
    if atoms.boundcond=='bulk':
        ax=atoms.cellvec[0][0] ; ay=atoms.cellvec[0][1] ; az=atoms.cellvec[0][2]
        bx=atoms.cellvec[1][0] ; by=atoms.cellvec[1][1] ; bz=atoms.cellvec[1][2]
        cx=atoms.cellvec[2][0] ; cy=atoms.cellvec[2][1] ; cz=atoms.cellvec[2][2]
    if filename=="screen":
        #print "%d" % atoms.nat
        #print "%s" % atoms.boundcond
        if atoms.boundcond=='bulk':
            print "lattice_vector  %24.15E%24.15E%24.15E" % (ax,ay,az)
            print "lattice_vector  %24.15E%24.15E%24.15E" % (bx,by,bz)
            print "lattice_vector  %24.15E%24.15E%24.15E" % (cx,cy,cz)
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            print "atom  %24.15E%24.15E%24.15E%5s" % (x,y,z,atoms.sat[i])
    else:
        f= open(filename,"w")
        #f.write("%d\n" % atoms.nat)
        if atoms.boundcond=='bulk':
            f.write("lattice_vector  %24.15E%24.15E%24.15E\n" % (ax,ay,az))
            f.write("lattice_vector  %24.15E%24.15E%24.15E\n" % (bx,by,bz))
            f.write("lattice_vector  %24.15E%24.15E%24.15E\n" % (cx,cy,cz))
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            f.write("atom  %24.15E%24.15E%24.15E%5s\n" % (x,y,z,atoms.sat[i]))
            #f.write(" %12.8f%12.8f%12.8f\n" % (x,y,z))
        f.close()
#*****************************************************************************************
def aims_read(filename):
    f=open(filename,'r')
    atoms=[]
    iline=0
    for line in f.readlines():
        if '#' in line.strip(): continue
        if not line.strip(): continue
        iline+=1
        if iline==1:
            atoms=Atoms()
            atoms.boundcond="bulk"
            atoms.nat=0
            atoms.cellvec[0][0]=float(line.split()[1])
            atoms.cellvec[0][1]=float(line.split()[2])
            atoms.cellvec[0][2]=float(line.split()[3])
        elif iline==2:
            atoms.cellvec[1][0]=float(line.split()[1])
            atoms.cellvec[1][1]=float(line.split()[2])
            atoms.cellvec[1][2]=float(line.split()[3])
        elif iline==3:
            atoms.cellvec[2][0]=float(line.split()[1])
            atoms.cellvec[2][1]=float(line.split()[2])
            atoms.cellvec[2][2]=float(line.split()[3])
        elif (iline>3):
            if 'atom' in str(line.split()[0]):
                atoms.rat.append([])
                icol=0
                #loop over the elemets, split by whitespace
                for i in line.split():
                    icol+=1
                    if icol>1 and icol<5:
                        atoms.rat[-1].append(float(i))
                    elif icol==5:
                        atoms.sat.append(str(i))
                atoms.nat+=1
            else:
                break
    f.closed
    return atoms
#*****************************************************************************************
