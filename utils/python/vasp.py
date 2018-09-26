from atoms import *
from cellutils import *
import copy
#*****************************************************************************************
def RepresentsInt(s):
    try: 
        int(s)
        return True
    except ValueError:
        return False
#*****************************************************************************************
def poscar_util(atoms):
    atom_types=[]
    atom_nn=[]
    rat=[]
    posred=[]
    for iat in range(atoms.nat):
        if not atom_types.__contains__(atoms.sat[iat]):
            atom_types.append(atoms.sat[iat])
            atom_nn.append(1)
            rat.append(atoms.rat[iat])
#this generation of the drect coordinates will be replaced! All routines reading files should automatically store the reduced coordinates while reading
            try:
                posred.append(atoms.posred[iat])
            except:
                poscart=[]
                poscart.append(atoms.rat[iat]) 
                n = 1
                posred_tmp=rxyz_cart2int(atoms.cellvec,poscart,n)
                posred.append(posred_tmp)
        else:
            ind=atom_types.index(atoms.sat[iat])
            jat=0
            for jj in range(ind+1):
                jat+=atom_nn[jj]
            rat.insert(jat,atoms.rat[iat])
#this generation of the drect coordinates will be replaced! All routines reading files should automatically store the reduced coordinates while reading
            try:
                posred.insert(jat,atoms.posred[iat])
            except:
                poscart=[]
                poscart.append(atoms.rat[jat]) 
                n = 1
                posred_tmp=rxyz_cart2int(atoms.cellvec,poscart,n)
                posred.insert(jat,posred_tmp)
            atom_nn[ind]+=1
    return atom_types,atom_nn,rat,posred
#*****************************************************************************************
def poscar_write(atoms,filename):
    atom_types,atom_nn,rat,posred=poscar_util(atoms)
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
            if atoms.coordinates=="Cartesian":
                x=float(rat[iat][0])
                y=float(rat[iat][1])
                z=float(rat[iat][2])
            elif atoms.coordinates=="Direct":
                x=float(posred[iat][0])
                y=float(posred[iat][1])
                z=float(posred[iat][2])
            f.write("%24.15E%24.15E%24.15E\n" % (x,y,z))
        f.close()
#*****************************************************************************************
def xdatcar_read():
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
            #Here we check if it is a VASP4 or VASP5 format
            for j,i in enumerate(atom_types):
                vasp4=False
                if RepresentsInt(i):
                    #print "Presumably vasp4 format"
                    vasp4=True
                    atom_types[j]=chr(65+j)
            if vasp4:
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
                iline+=1
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
            if line.strip().split()[0][0]=="D" or line.strip().split()[0][0]=="d": 
                atoms.coordinates="Direct"
            elif line.strip().split()[0][0]=="C" or line.strip().split()[0][0]=="c" or line.strip().split()[0][0]=="K" or line.strip().split()[0][0]=="k":
                atoms.coordinates="Cartesian"
            else:
                print "ERROR: unknown coordinates"
        elif iline>8:
            if (iline-8)%(atoms.nat+1)>=1:
                atoms.rat.append([]) 
                atoms.posred.append([]) 
                xred=float(line.split()[0])
                yred=float(line.split()[1])
                zred=float(line.split()[2])
                if(atoms.coordinates=="Cartesian"):
                    atoms.rat[-1].append(xred)
                    atoms.rat[-1].append(yred)
                    atoms.rat[-1].append(zred)
                    poscart=[]
                    poscart.append([xred,yred,zred]) 
                    n = 1
                    posred=rxyz_cart2int(atoms.cellvec,poscart,n)
                    atoms.posred[-1].append(posred[0][0])
                    atoms.posred[-1].append(posred[0][1])
                    atoms.posred[-1].append(posred[0][2])
                elif(atoms.coordinates=="Direct"):
                    atoms.rat[-1].append(xred*atoms.cellvec[0][0]+yred*atoms.cellvec[1][0]+zred*atoms.cellvec[2][0])
                    atoms.rat[-1].append(xred*atoms.cellvec[0][1]+yred*atoms.cellvec[1][1]+zred*atoms.cellvec[2][1])
                    atoms.rat[-1].append(xred*atoms.cellvec[0][2]+yred*atoms.cellvec[1][2]+zred*atoms.cellvec[2][2])
                    atoms.posred[-1].append(xred)
                    atoms.posred[-1].append(yred)
                    atoms.posred[-1].append(zred)
                else:
                    print "ERROR: unknown coordinates"
            if (iline-8)%(atoms.nat+1)==0:
                atoms_all.append(Atoms())
                atoms_all[-1]=copy.copy(atoms)
                nconf+=1
                atoms.rat=[]
                iline=1
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
            scale_factor=float(line.split()[0])
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
            #scaling cell vectors with scale_factor
            atoms.cellvec[0][0]*=scale_factor
            atoms.cellvec[0][1]*=scale_factor
            atoms.cellvec[0][2]*=scale_factor
            atoms.cellvec[1][0]*=scale_factor
            atoms.cellvec[1][1]*=scale_factor
            atoms.cellvec[1][2]*=scale_factor
            atoms.cellvec[2][0]*=scale_factor
            atoms.cellvec[2][1]*=scale_factor
            atoms.cellvec[2][2]*=scale_factor
        elif iline==6:
            ntypes=0
            atom_types=[]
            for i in line.split():
                atom_types.append(i)
                ntypes+=1
            #Here we check if it is a VASP4 or VASP5 format
            for j,i in enumerate(atom_types):
                vasp4=False
                if RepresentsInt(i):
                    #print "Presumably vasp4 format"
                    vasp4=True
                    atom_types[j]=chr(65+j)
            if vasp4:
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
                iline+=1
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
            if line.strip().split()[0][0]=="D" or line.strip().split()[0][0]=="d": 
                atoms.coordinates="Direct"
            elif line.strip().split()[0][0]=="C" or line.strip().split()[0][0]=="c" or line.strip().split()[0][0]=="K" or line.strip().split()[0][0]=="k":
                atoms.coordinates="Cartesian"
            else:
                print "ERROR: unknown coordinates"
        elif iline>8 and iline<=atoms.nat+8:
            atoms.rat.append([])
            atoms.posred.append([])
            xred=float(line.split()[0])
            yred=float(line.split()[1])
            zred=float(line.split()[2])
            atoms.bemoved.append("TTT")
            if(atoms.coordinates=="Direct"):
                x=xred*atoms.cellvec[0][0]+yred*atoms.cellvec[1][0]+zred*atoms.cellvec[2][0]
                y=xred*atoms.cellvec[0][1]+yred*atoms.cellvec[1][1]+zred*atoms.cellvec[2][1]
                z=xred*atoms.cellvec[0][2]+yred*atoms.cellvec[1][2]+zred*atoms.cellvec[2][2]
                posred= np.zeros((1,3))
                posred[0][:]=([xred,yred,zred])
            elif(atoms.coordinates=="Cartesian"):
                x=xred*scale_factor
                y=yred*scale_factor
                z=zred*scale_factor
                poscart=[]
                poscart.append([x,y,z]) 
                n = 1
                posred=rxyz_cart2int(atoms.cellvec,poscart,n)
            else:
                print "ERROR: unknown coordinates"
            atoms.rat[-1].append(x)
            atoms.rat[-1].append(y)
            atoms.rat[-1].append(z)
            atoms.posred[-1].append(posred[0][0])
            atoms.posred[-1].append(posred[0][1])
            atoms.posred[-1].append(posred[0][2])
            if not line.strip(): break
    f.closed
    return atoms
#*****************************************************************************************
