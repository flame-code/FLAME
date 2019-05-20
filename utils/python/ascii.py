import sys
import atoms
from atoms import *
#*****************************************************************************************
def ascii_write(atoms,filename):
    if filename=="screen":
        print "%d" % atoms.nat
        #print "%s" % atoms.boundcond
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1])
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2])
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            print "%24.15E%24.15E%24.15E%5s" % (x,y,z,atoms.sat[i])
    else:
        f= open(filename,"w")
        f.write("%d\n" % atoms.nat)
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1]))
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2]))
        for i in range(atoms.nat):
            x=atoms.rat[i][0]
            y=atoms.rat[i][1]
            z=atoms.rat[i][2]
            f.write("%24.15E%24.15E%24.15E%5s\n" % (x,y,z,atoms.sat[i]))
            #f.write(" %12.8f%12.8f%12.8f\n" % (x,y,z))
        f.close()
#*****************************************************************************************
def ascii_read(filename):
    f=open(filename,"r")
    atoms=Atoms()
    atoms.boundcond="bulk"
    iline=0
    iline_tot=0
    nconf=0
    atoms.nat=0
    for line in f.readlines():
        iline_tot+=1
        tt=str(line).strip()
        if tt[0]=='#': continue
        #print tt[0]
        iline+=1
        if iline==1:
            if len(line.split()) > 1:
               #atoms.epot=float(line.split()[1])*27.211385
                atoms.epot=float(line.split()[1])
            else:
                atoms.epot=0.0
            #print atoms.epot
            pass
        elif iline==2:
            atoms.cellvec[0][0]=float(line.split()[0])
            atoms.cellvec[1][0]=float(line.split()[1])
            atoms.cellvec[1][1]=float(line.split()[2])
        elif iline==3:
            atoms.cellvec[2][0]=float(line.split()[0])
            atoms.cellvec[2][1]=float(line.split()[1])
            atoms.cellvec[2][2]=float(line.split()[2])
        else:
            atoms.nat+=1
            atoms.bemoved.append("TTT")
            atoms.rat.append([])
            icol=0
            #loop over the elemets, split by whitespace
            for i in line.split():
                #print i
                icol+=1
                if icol<4:
                    atoms.rat[-1].append(float(i))
                elif icol<5:
                    atoms.sat.append(i)
    f.closed
    return atoms
#*****************************************************************************************
