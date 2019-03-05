import sys
import copy
from atoms import *
#*****************************************************************************************
def xyz_read(filename,style):
    atoms_all=[]
    iline=0
    f=open(filename,"r")
    #read line into array 
    for line in f.readlines():
        if iline==0:
            atoms=Atoms()
        iline+=1
        #print iline
        if iline==1:
            atoms.nat=int(line.split()[0])
            units=line.split()[1]
            atoms.epot=float(line.split()[2])
            #if units=='atomicd0' or units=='atomic':
            atoms.epot=atoms.epot*27.211385 #in bigdft energy is in Ha disregarding units of position
        if iline==2:
            if style=='mine':
                atoms.boundcond=line.split()[0]
                atoms.cellvec[0][0]=float(line.split()[1])
                atoms.cellvec[1][1]=float(line.split()[2])
                atoms.cellvec[2][2]=float(line.split()[3])
                atoms.cellvec[1][0]=float(line.split()[4])
                atoms.cellvec[2][0]=float(line.split()[5])
                atoms.cellvec[2][1]=float(line.split()[6])
            elif style=='bigdft':
                atoms.boundcond='free'
                atoms.cellvec[0][0]=100.0
                atoms.cellvec[1][1]=100.0
                atoms.cellvec[2][2]=100.0
                atoms.cellvec[1][0]=0.0
                atoms.cellvec[2][0]=0.0
                atoms.cellvec[2][1]=0.0
        if iline>2:
            #add a new sublist
            atoms.rat.append([])
            icol=0
            #loop over the elemets, split by whitespace
            for i in line.split():
                icol+=1
                if icol==1:
                    atoms.sat.append(i)
                elif icol<5:
                    if units=='atomicd0' or units=='atomic':
                        atoms.rat[-1].append(float(i)*0.529177+50.0)
                    else:
                        atoms.rat[-1].append(float(i))
                elif icol==5:
                    atoms.bemoved.append(i)
                if style=='bigdft': atoms.bemoved.append('TTT')
            #print "%d" % atoms.nat
        if iline==int(atoms.nat)+2:
            atoms_all.append(Atoms())
            atoms_all[-1]=copy.copy(atoms)
            iline=0
            break
    f.closed
    return atoms_all
#*****************************************************************************************
def xyz_write(atoms_all,frmt):
    iconf=0
    for atoms in atoms_all:
        iconf+=1
        if frmt=='bigdft':
            first="%d  angstroem  struct%5.5d" % (atoms.nat,iconf)
            second=atoms.boundcond
        else:
            first="%d" % atoms.nat
            second=''
        print first
        print second
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            if frmt=='bigdft':
                print "%5s  %23.14E%23.14E%23.14E" % (atoms.sat[i],x,y,z)
            else:
                print "%5s  %23.14E%23.14E%23.14E%5s" % (atoms.sat[i],x,y,z,atoms.bemoved[i])
#*****************************************************************************************
def xyz_write_b(atoms_all,frmt,filename):
    f= open(filename,"w")
    iconf=0
    for atoms in atoms_all:
        iconf+=1
        if frmt=='bigdft':
            first="%d  angstroem  struct%5.5d" % (atoms.nat,iconf)
            second=atoms.boundcond
        else:
            first="%d" % atoms.nat
            second=''
        f.write("%s \n" % (first))
        f.write("%s \n" % (second))
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            if frmt=='bigdft':
                f.write("%5s  %23.14E%23.14E%23.14E \n" % (atoms.sat[i],x,y,z))
            else:
                f.write("%5s  %23.14E%23.14E%23.14E%5s\n" % (atoms.sat[i],x,y,z,atoms.bemoved[i]))
#*****************************************************************************************
