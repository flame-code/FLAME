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
    i_lattice_vector=-1
    for line in f.readlines():
        if '#' in line.strip(): continue
        if not line.strip(): continue
        iline+=1
        if iline==1:
            atoms=Atoms()
            atoms.nat=0
        str_lattice_vector=line.strip()[0:14]
        str_atom=line.strip()[0:4]
        if str_lattice_vector=='lattice_vector':
            i_lattice_vector+=1
            atoms.cellvec[i_lattice_vector][0]=float(line.split()[1])
            atoms.cellvec[i_lattice_vector][1]=float(line.split()[2])
            atoms.cellvec[i_lattice_vector][2]=float(line.split()[3])
        elif str_atom=='atom':
            atoms.rat.append([])
            icol=0
            #loop over the elemets, split by whitespace
            for i in line.split():
                icol+=1
                if icol>1 and icol<5:
                    atoms.rat[-1].append(float(i))
                elif icol==5:
                    atoms.sat.append(str(i))
            atoms.bemoved.append("TTT")
            atoms.nat+=1
            #print "%16.8f%16.8f%16.8f" % (atoms.rat[-1][0],atoms.rat[-1][1],atoms.rat[-1][2])
    if i_lattice_vector==-1:
        atoms.boundcond="free"
	atoms.cellvec,atoms.rat=set_cell(atoms)
    else:
        atoms.boundcond="bulk"
    f.closed
    return atoms
#*****************************************************************************************
def set_cell(atoms):
    amargin=3.0
    xmin=sys.float_info.max ; xmax=-sys.float_info.max
    ymin=sys.float_info.max ; ymax=-sys.float_info.max
    zmin=sys.float_info.max ; zmax=-sys.float_info.max
    for iat in range(atoms.nat):
        if atoms.rat[iat][0]>xmax:
            xmax=atoms.rat[iat][0]
        if atoms.rat[iat][0]<xmin:
            xmin=atoms.rat[iat][0]
        if atoms.rat[iat][1]>ymax:
            ymax=atoms.rat[iat][1]
        if atoms.rat[iat][1]<ymin:
            ymin=atoms.rat[iat][1]
        if atoms.rat[iat][2]>zmax:
            zmax=atoms.rat[iat][2]
        if atoms.rat[iat][2]<zmin:
            zmin=atoms.rat[iat][2]
    #---------------------------------------------------------------------------
    cvxx=xmax-xmin+2*amargin
    cvyy=ymax-ymin+2*amargin
    cvzz=zmax-zmin+2*amargin
    cellvec=[[cvxx,0.0,0.0],[0.0,cvyy,0.0],[0.0,0.0,cvzz]]
    rat=[]
    for iat in range(atoms.nat):
        rat.append([])
        rat[-1].append(atoms.rat[iat][0]-xmin+amargin)
        rat[-1].append(atoms.rat[iat][1]-ymin+amargin)
        rat[-1].append(atoms.rat[iat][2]-zmin+amargin)
    return cellvec,rat
#*****************************************************************************************
def get_input_geometry(iline,lines,nat,has_unit_cell):
    atoms=Atoms()
    if has_unit_cell:
        atoms.boundcond="bulk"
        atoms.cellvec[0][0]=float(lines[iline+2].split()[1])
        atoms.cellvec[0][1]=float(lines[iline+2].split()[2])
        atoms.cellvec[0][2]=float(lines[iline+2].split()[3])
        atoms.cellvec[1][0]=float(lines[iline+3].split()[1])
        atoms.cellvec[1][1]=float(lines[iline+3].split()[2])
        atoms.cellvec[1][2]=float(lines[iline+3].split()[3])
        atoms.cellvec[2][0]=float(lines[iline+4].split()[1])
        atoms.cellvec[2][1]=float(lines[iline+4].split()[2])
        atoms.cellvec[2][2]=float(lines[iline+4].split()[3])
        istart=iline+7
    else:
        atoms.boundcond="free"
        istart=iline+4
    #print atoms.cellvec
    atoms.nat=nat
    for iat in range(nat):
        atoms.sat.append(lines[istart+iat].split()[3])
        #print atoms.sat[-1]
        atoms.rat.append([])
        atoms.rat[-1].append(float(lines[istart+iat].split()[4]))
        atoms.rat[-1].append(float(lines[istart+iat].split()[5]))
        atoms.rat[-1].append(float(lines[istart+iat].split()[6]))
        atoms.bemoved.append("TTT")
    return atoms
#*****************************************************************************************
def get_updated_geometry(iline,lines,nat,has_unit_cell):
    atoms=Atoms()
    if has_unit_cell:
        atoms.boundcond="bulk"
        atoms.cellvec[0][0]=float(lines[iline+2].split()[1])
        atoms.cellvec[0][1]=float(lines[iline+2].split()[2])
        atoms.cellvec[0][2]=float(lines[iline+2].split()[3])
        atoms.cellvec[1][0]=float(lines[iline+3].split()[1])
        atoms.cellvec[1][1]=float(lines[iline+3].split()[2])
        atoms.cellvec[1][2]=float(lines[iline+3].split()[3])
        atoms.cellvec[2][0]=float(lines[iline+4].split()[1])
        atoms.cellvec[2][1]=float(lines[iline+4].split()[2])
        atoms.cellvec[2][2]=float(lines[iline+4].split()[3])
        istart=iline+6
    else:
        atoms.boundcond="free"
        istart=iline+2
    #print atoms.cellvec
    atoms.nat=nat
    if lines[istart].split()[0]=='|':
        icol_xyz=3
        icol_sat=3
    else:
        icol_xyz=0
        icol_sat=4
    for iat in range(nat):
        atoms.sat.append(lines[istart+iat].split()[icol_sat])
        #print atoms.sat[-1]
        atoms.rat.append([])
        atoms.rat[-1].append(float(lines[istart+iat].split()[icol_xyz+1]))
        atoms.rat[-1].append(float(lines[istart+iat].split()[icol_xyz+2]))
        atoms.rat[-1].append(float(lines[istart+iat].split()[icol_xyz+3]))
        atoms.bemoved.append("TTT")
    return atoms
#*****************************************************************************************
