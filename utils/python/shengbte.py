from atoms import *
import numpy as np
import copy
#*****************************************************************************************
def sheng_util(atoms):
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
#**********************************************************************************#
    #!This function will convert the cartesian coordinates into the internal coordinates 
def rxyz_cart2int(cellvec,ratcart,nat):
    ratint= np.zeros((nat,3))
    cellvecinv = invertmat(cellvec)
    for iat in range(nat):
        ratint[iat][0] = cellvecinv[0][0]*ratcart[iat][0]+cellvecinv[1][0]*ratcart[iat][1]+cellvecinv[2][0]*ratcart[iat][2]
        ratint[iat][1] = cellvecinv[0][1]*ratcart[iat][0]+cellvecinv[1][1]*ratcart[iat][1]+cellvecinv[2][1]*ratcart[iat][2]
        ratint[iat][2] = cellvecinv[0][2]*ratcart[iat][0]+cellvecinv[1][2]*ratcart[iat][1]+cellvecinv[2][2]*ratcart[iat][2]
    return ratint
#********************************************************************#
def invertmat(mat):
    matinv = np.zeros((3,3))
    a = copy.deepcopy(mat)
    div = a[0][0]*(a[1][1]*a[2][2] - a[2][1]*a[1][2]) - a[0][1]*(a[1][0]*a[2][2] - a[1][2]*a[2][0]) + a[0][2]*(a[1][0]*a[2][1] - a[1][1]*a[2][0])
    div = 1.0 / div
    matinv[0][0] =  (a[1][1]*a[2][2] - a[1][2]*a[2][1]) *  div
    matinv[0][1] = -(a[0][1]*a[2][2] - a[0][2]*a[2][1]) *  div
    matinv[0][2] =  (a[0][1]*a[1][2] - a[0][2]*a[1][1]) *  div
    matinv[1][0] = -(a[1][0]*a[2][2] - a[1][2]*a[2][0]) *  div
    matinv[1][1] =  (a[0][0]*a[2][2] - a[0][2]*a[2][0]) *  div
    matinv[1][2] = -(a[0][0]*a[1][2] - a[0][2]*a[1][0]) *  div
    matinv[2][0] =  (a[1][0]*a[2][1] - a[1][1]*a[2][0]) *  div
    matinv[2][1] = -(a[0][0]*a[2][1] - a[0][1]*a[2][0]) *  div
    matinv[2][2] =  (a[0][0]*a[1][1] - a[0][1]*a[1][0]) *  div
    #matinv = np.matrix('1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0')
    #matinv = np.linalg.inv(mat)
    return matinv
#*****************************************************************************************
def sheng_write(atoms,filename):
    atom_types,atom_nn,rat=sheng_util(atoms)
    ntypes=len(atom_types)
    xred=rxyz_cart2int(atoms.cellvec,rat,atoms.nat)
    sc1=4
    sc2=4
    sc3=4
    if filename=="screen":
        str_line=""
        print "&allocations"
        print "    %s %d" % ("nelements= ", len(atom_nn))
        print "    %s %d" % ("natoms= ", atoms.nat)
        print "    %s %d %d %d" % ("ngrid(:)= ", 15, 15, 15)
        print "    %s %d" % ("norientations= ", 0)
        print "&end"

        print "&crystal"
        print "    %s %24.15E" % ("lfactor= ", 1.0)
        print "    %s %24.15E%24.15E%24.15E" % ("lattvec(:,1)= ", 0.1*atoms.cellvec[0][0],0.1*atoms.cellvec[0][1],0.1*atoms.cellvec[0][2])
        print "    %s %24.15E%24.15E%24.15E" % ("lattvec(:,2)= ", 0.1*atoms.cellvec[1][0],0.1*atoms.cellvec[1][1],0.1*atoms.cellvec[1][2])
        print "    %s %24.15E%24.15E%24.15E" % ("lattvec(:,3)= ", 0.1*atoms.cellvec[2][0],0.1*atoms.cellvec[2][1],0.1*atoms.cellvec[2][2])
        for itype in range(ntypes):
            str_line+=" %s," % ("'"+atom_types[itype].replace(' ','')+"'")
        print     "%s %s" % ("    elements=  ", str_line)
        str_line=""
        for itype in range(ntypes):
            for icount in range(atom_nn[itype]):
                str_line+=" %d, " % (itype+1)
        print     "%s %s" % ("    types=  ", str_line)
        for iat in range(atoms.nat):
            istr=str(iat+1)
            x=float(xred[iat][0])
            y=float(xred[iat][1])
            z=float(xred[iat][2])
            print  "    %s %24.15E%24.15E%24.15E" % ("positions(:,"+istr+")= ",x,y,z)
        print "    %s %d %d %d" % ("scell(:)= ",sc1,sc2,sc3)
        print "&end"

        print "&parameters"
        print "    T=300"
        print "    scalebroad=0.2"
        print "&end"

        print "&flags"
        print "    onlyharmonic=.FALSE."
        print "    isotopes=.FALSE.    "
        print "    nonanalytic=.FALSE. "
        print "    nanowires=.FALSE.   "
        print "&end"
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
