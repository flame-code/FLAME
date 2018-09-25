from atoms import *
import copy
#*****************************************************************************************
def acf_read(filename):
    f = open (filename,"r")
    atoms_all=[]
    iline=0
    nconf=0
    for line in f.readlines():
        iline+=1
        if iline==1:
            atoms=Atoms()
            nconf+=1
            #print nconf
            comment1=line
        elif iline==2:
            comment2=line
        elif iline==3:
            comment3=line
        elif iline==7:
            icol=0
            for i in line.split():
                #if icol==2: print i
                if i=="epot=": atoms.epot=float(line.split()[icol+1])
                if i=="qtot=": atoms.qtot=float(line.split()[icol+1])
                icol+=1
        elif iline==8:
            atoms.nat=int(line.split()[0])
            atoms.boundcond=line.split()[1]
        elif iline==9:
            atoms.cellvec[0][0]=float(line.split()[0])
            atoms.cellvec[1][0]=float(line.split()[1])
            atoms.cellvec[1][1]=float(line.split()[2])
        elif iline==10:
            atoms.cellvec[2][0]=float(line.split()[0])
            atoms.cellvec[2][1]=float(line.split()[1])
            atoms.cellvec[2][2]=float(line.split()[2])
        elif iline>10 and iline<atoms.nat+11:
            atoms.rat.append([])
            icol=0
            #loop over the elemets, split by whitespace
            for i in line.split():
                #print i
                icol+=1
                if icol==1:
                    atoms.sat.append(i)
                elif icol<5:
                    atoms.rat[-1].append(float(i))
                elif icol==5:
                    atoms.bemoved.append(i)
            if iline==atoms.nat+10:
                iline=0
                atoms_all.append(Atoms())
                atoms_all[-1]=copy.copy(atoms)
                #print atoms_all[-1].nat
                #print len(atoms_all[-1].rat)
                #print
                #atoms.kill()
                #print atoms_all[-1].nat
                #print len(atoms_all[-1].rat)
                #print
    f.closed
    return atoms_all
#*****************************************************************************************
def acf_7th_line(atoms,label):
    str_7th_line=""
    str_7th_line+=label
    if str(atoms.epot)!=1.E100:
        str_7th_line+="  epot=%24.15E" % atoms.epot
    if atoms.qtot!=0.0:
        str_7th_line+="  qtot=%8.4f" % atoms.qtot
    return str_7th_line
#*****************************************************************************************
def acf_write(atoms_all,filename,labelpatt='none'):
    nconf=0
    if filename=="screen":
        for atoms in atoms_all:
            nconf+=1
            print "#1st line comment"
            print "#2nd line comment"
            print "#3rd line comment"
            print ""
            print ""
            print "c5=bemoved"
            if labelpatt=='none': label=''
            else: label='label= %s%5.5d' % (labelpatt,nconf)
            print acf_7th_line(atoms,label)
            print "%6d  %s  %d" % (int(atoms.nat),atoms.boundcond,nconf)
            print "%24.15E%24.15E%24.15E" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1])
            print "%24.15E%24.15E%24.15E" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2])
            for i in range(int(atoms.nat)):
                x=atoms.rat[i][0]
                y=atoms.rat[i][1]
                z=atoms.rat[i][2]
                atoms.bemoved.append("TTT")
                print "%5s  %23.14E%23.14E%23.14E%5s" % (atoms.sat[i],x,y,z,atoms.bemoved[i])
    else:
        f= open(filename,"w")
        for atoms in atoms_all:
            nconf+=1
            f.write("#1st line comment\n")
            f.write("#2nd line comment\n")
            f.write("#3rd line comment\n")
            f.write("\n")
            f.write("\n")
            f.write("c5=bemoved\n")
            if labelpatt=='none': label=''
            else: label='label= %s%5.5d' % (labelpatt,nconf)
            f.write(acf_7th_line(atoms,label)+"\n")
            f.write("%6d  %s  %d\n" % (int(atoms.nat),atoms.boundcond,nconf))
            f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[0][0],atoms.cellvec[1][0],atoms.cellvec[1][1]))
            f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2]))
            for i in range(int(atoms.nat)):
                x=atoms.rat[i][0]
                y=atoms.rat[i][1]
                z=atoms.rat[i][2]
                f.write("%5s  %23.14E%23.14E%23.14E%5s\n" % (atoms.sat[i],x,y,z,atoms.bemoved[i]))
        f.close()
#*****************************************************************************************
def read_forces(nat,filename):
    """read the force file """
    fd = open(filename, 'r')
    forces=[]
    while True:
        line = fd.readline()
        if not line:
            break
        if "configuration" in line:
            iat=0
            fat=[]
            continue
        fat.append([float(line.split()[0]), float(line.split()[1]), float(line.split()[2])])
        iat+=1
        if iat==nat:
            forces.append(fat)
            del fat
    return forces
    fd.close()
#*****************************************************************************************
