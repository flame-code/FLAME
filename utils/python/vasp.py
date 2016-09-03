from atoms import *
import copy
#*****************************************************************************************
def poscar_util(atoms):
    atom_types=[]
    atom_nn=[]
    rat=[]
    for iat in range(atoms.nat):
        if not atom_types.__contains__(atoms.sat[iat]):
            atom_types.append(atoms.sat[iat])
            atom_nn.append(1)
            rat.append(atoms.rat[iat])
        else:
            ind=atom_types.index(atoms.sat[iat])
            jat=0
            for jj in range(ind+1):
                jat+=atom_nn[jj]
            rat.insert(jat,atoms.rat[iat])
            atom_nn[ind]+=1
    return atom_types,atom_nn,rat
#*****************************************************************************************
def poscar_write(atoms,filename):
    atom_types,atom_nn,rat=poscar_util(atoms)
    ntypes=len(atom_types)
    if filename=="screen":
        str_line=""
        for itype in range(ntypes):
            str_line+=" %s" % atom_types[itype]
        print str_line
        print "1.0"
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[0][0],atoms.cellvec[0][1],atoms.cellvec[0][2])
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[1][0],atoms.cellvec[1][1],atoms.cellvec[1][2])
        print "%24.15E%24.15E%24.15E" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2])
        print str_line
        str_line=""
        for itype in range(ntypes):
            str_line+=" %d" % atom_nn[itype]
        print str_line
        print atoms.coordinates
        for iat in range(atoms.nat):
            x=float(rat[iat][0])
            y=float(rat[iat][1])
            z=float(rat[iat][2])
            print "%24.15E%24.15E%24.15E" % (x,y,z)
    else:
        str_line=""
        for itype in range(ntypes):
            str_line+=" %s" % atom_types[itype]
        f=open(filename,"w")
        f.write("%s\n" % str_line)
        f.write("1.0\n")
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[0][0],atoms.cellvec[0][1],atoms.cellvec[0][2]))
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[1][0],atoms.cellvec[1][1],atoms.cellvec[1][2]))
        f.write("%24.15E%24.15E%24.15E\n" % (atoms.cellvec[2][0],atoms.cellvec[2][1],atoms.cellvec[2][2]))
        f.write("%s\n" % str_line)
        str_line=""
        for itype in range(ntypes):
            str_line+=" %d" % atom_nn[itype]
        f.write("%s\n" % str_line)
        f.write("%s\n" % atoms.coordinates)
        for iat in range(atoms.nat):
            x=float(rat[iat][0])
            y=float(rat[iat][1])
            z=float(rat[iat][2])
            f.write("%24.15E%24.15E%24.15E\n" % (x,y,z))
        f.close()
#*****************************************************************************************
def xdatcar_read(coordinates):
    f=open("XDATCAR","r")
    atoms_all=[]
    iline=0
    nconf=0
    for line in f.readlines():
        iline+=1
        if iline==1:
            atoms=Atoms()
            comment1=line
        elif iline==2:
            scale=float(line)
        elif iline==3:
            atoms.cellvec[0][0]=float(line.split()[0])
            atoms.cellvec[0][1]=float(line.split()[1])
            atoms.cellvec[0][2]=float(line.split()[2])
        elif iline==4:
            atoms.cellvec[1][0]=float(line.split()[0])
            atoms.cellvec[1][1]=float(line.split()[1])
            atoms.cellvec[1][2]=float(line.split()[2])
        elif iline==5:
            atoms.cellvec[2][0]=float(line.split()[0])
            atoms.cellvec[2][1]=float(line.split()[1])
            atoms.cellvec[2][2]=float(line.split()[2])
        elif iline==6:
            ntypes=0
            atom_types=[]
            for i in line.split():
                atom_types.append(i)
                ntypes+=1
        elif iline==7:
            ncol=0
            atoms.nat=0
            atom_nn=[]
            for i in line.split():
                atom_nn.append(int(i))
                atoms.nat+=int(i)
                ncol+=1
            if ncol!=ntypes:
                print "ERROR: in lines 6 and 7, number of columns differ: %3d%3d",ntypes,ncol
            for itype in range(ntypes):
                for iat in range(atom_nn[itype]): 
                    atoms.sat.append(atom_types[itype])
        elif iline>7:
            if (iline-8)%(atoms.nat+1)>=1:
                atoms.rat.append([])
                xred=float(line.split()[0])
                yred=float(line.split()[1])
                zred=float(line.split()[2])
                if(coordinates=="Cartesian"):
                    atoms.rat[-1].append(xred*atoms.cellvec[0][0]+yred*atoms.cellvec[1][0]+zred*atoms.cellvec[2][0])
                    atoms.rat[-1].append(xred*atoms.cellvec[0][1]+yred*atoms.cellvec[1][1]+zred*atoms.cellvec[2][1])
                    atoms.rat[-1].append(xred*atoms.cellvec[0][2]+yred*atoms.cellvec[1][2]+zred*atoms.cellvec[2][2])
                elif(coordinates=="Reduced"):
                    atoms.rat[-1].append(xred)
                    atoms.rat[-1].append(yred)
                    atoms.rat[-1].append(zred)
                else:
                    print "ERROR: unknown coordinates"
                atoms.coordinates=coordinates
            if (iline-7)%(atoms.nat+1)==0:
                atoms_all.append(Atoms())
                atoms_all[-1]=copy.copy(atoms)
                nconf+=1
                atoms.rat=[]
    f.closed
    return atoms_all
#*****************************************************************************************
def poscar_read(filename):
    f=open(filename,"r")
    atoms=[]
    iline=0
    for line in f.readlines():
        iline+=1
        if iline==1:
            atoms=Atoms()
            atoms.boundcond="bulk"
            #nconf+=1
            comment1=line
        elif iline==2:
            comment2=line
        elif iline==3:
            atoms.cellvec[0][0]=float(line.split()[0])
            atoms.cellvec[0][1]=float(line.split()[1])
            atoms.cellvec[0][2]=float(line.split()[2])
        elif iline==4:
            atoms.cellvec[1][0]=float(line.split()[0])
            atoms.cellvec[1][1]=float(line.split()[1])
            atoms.cellvec[1][2]=float(line.split()[2])
        elif iline==5:
            atoms.cellvec[2][0]=float(line.split()[0])
            atoms.cellvec[2][1]=float(line.split()[1])
            atoms.cellvec[2][2]=float(line.split()[2])
        elif iline==6:
            ntypes=0
            atom_types=[]
            for i in line.split():
                atom_types.append(i)
                ntypes+=1
        elif iline==7:
            ncol=0
            atoms.nat=0
            atom_nn=[]
            for i in line.split():
                atom_nn.append(int(i))
                atoms.nat+=int(i)
                ncol+=1
            if ncol!=ntypes:
                print "ERROR: in lines 6 and 7, number of columns differ: %3d%3d",ntypes,ncol
            for itype in range(ntypes):
                for iat in range(atom_nn[itype]): 
                    atoms.sat.append(atom_types[itype])
        elif iline==8:
            if "Direct" in line.strip(): 
                atoms.coordinates="Direct"
            else:
                atoms.coordinates="Cartesian"
        elif iline>8:
            atoms.rat.append([])
            xred=float(line.split()[0])
            yred=float(line.split()[1])
            zred=float(line.split()[2])
            if(atoms.coordinates=="Direct"):
                atoms.rat[-1].append(xred*atoms.cellvec[0][0]+yred*atoms.cellvec[1][0]+zred*atoms.cellvec[2][0])
                atoms.rat[-1].append(xred*atoms.cellvec[0][1]+yred*atoms.cellvec[1][1]+zred*atoms.cellvec[2][1])
                atoms.rat[-1].append(xred*atoms.cellvec[0][2]+yred*atoms.cellvec[1][2]+zred*atoms.cellvec[2][2])
            elif(atoms.coordinates=="Cartesian"):
                atoms.rat[-1].append(xred)
                atoms.rat[-1].append(yred)
                atoms.rat[-1].append(zred)
            else:
                print "ERROR: unknown coordinates"
            if not line.strip(): break
    f.closed
    return atoms
#*****************************************************************************************
