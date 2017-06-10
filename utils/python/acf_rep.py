#!/usr/bin/env python
import argparse
import atoms
from acf import *
#****************************************************
def replicate(na,nb,nc,nat,cellvec,rat):
    atoms_t=[]
    atoms_tt=[]
    atoms_sat=[]
    atoms_bemoved=[]
    if atoms.boundcond=="bulk":
        for i in range(int(atoms.nat)):
            for n in range(0,int(na)):
                for m in range(0,int(nb)):
                    for l in range(0,int(nc)):
                        atoms_tt.append(atoms.rat[i][0]+n*atoms.cellvec[0][0]+m*atoms.cellvec[1][0]+l*atoms.cellvec[2][0])
                        atoms_tt.append(atoms.rat[i][1]+m*atoms.cellvec[1][1]+l*atoms.cellvec[2][1])
                        atoms_tt.append(atoms.rat[i][2]+l*atoms.cellvec[2][2])
                        atoms_t.append(atoms_tt)
                        atoms_sat.append(atoms.sat[i]) 
                        atoms_bemoved.append(atoms.bemoved[i])
                        atoms_tt=[]
    else:
        print "ERROR: Input structure is not BULK."
    atoms.cellvec[0][0]=atoms.cellvec[0][0]*na 
    atoms.cellvec[1][0]=atoms.cellvec[1][0]*nb
    atoms.cellvec[1][1]=atoms.cellvec[1][1]*nb
    atoms.cellvec[2][0]=atoms.cellvec[2][0]*nc
    atoms.cellvec[2][1]=atoms.cellvec[2][1]*nc
    atoms.cellvec[2][2]=atoms.cellvec[2][2]*nc
    nn=na*nb*nc*atoms.nat
    atoms.nat = nn    
    atoms.rat=atoms_t
    atoms.sat=atoms_sat
    atoms.bemoved=atoms_bemoved
    return atoms
#****************************************************
str1 = "This script replicate an acf file in the cell vector directions."
str1+="\n The n_i are integers and b_i are reciprocal vectors. The output file is also in acf format."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_input', action='store' ,type=str, help="fn_input is the name of the input file")
parser.add_argument('n1', type=int, help='maxmimum number of k-vectors in direction b1')
parser.add_argument('n2', type=int, help='maxmimum number of k-vectors in direction b2')
parser.add_argument('n3', type=int, help='maxmimum number of k-vectors in direction b3')
args=parser.parse_args()

atoms_all=acf_read(args.fn_input)
atoms=atoms_all[-1]
n1 = args.n1
n2 = args.n2
n3 = args.n3
atoms = replicate(n1,n2,n3,atoms.nat,atoms.cellvec,atoms.rat)
atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
acf_write(atoms_all,"screen")
