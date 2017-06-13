#!/usr/bin/env python
import argparse
import atoms
from acf import *
#****************************************************
def replicate(na,nb,nc,atoms):
    atoms_out=Atoms()
    if atoms.boundcond=="bulk":
        for i in range(int(atoms.nat)):
            for n in range(0,int(na)):
                for m in range(0,int(nb)):
                    for l in range(0,int(nc)):
                        x = atoms.rat[i][0]+n*atoms.cellvec[0][0]+m*atoms.cellvec[1][0]+l*atoms.cellvec[2][0]
                        y = atoms.rat[i][1]+n*atoms.cellvec[0][1]+m*atoms.cellvec[1][1]+l*atoms.cellvec[2][1]
                        z = atoms.rat[i][2]+n*atoms.cellvec[0][2]+m*atoms.cellvec[1][2]+l*atoms.cellvec[2][2]
                        atoms_out.rat.append([])
                        atoms_out.rat[-1].append(x)
                        atoms_out.rat[-1].append(y)
                        atoms_out.rat[-1].append(z)
                        atoms_out.sat.append(atoms.sat[i]) 
                        atoms_out.bemoved.append(atoms.bemoved[i])
    else:
        print "ERROR: Input structure is not BULK."
    atoms_out.cellvec[0][0]=atoms.cellvec[0][0]*na 
    atoms_out.cellvec[0][1]=atoms.cellvec[0][1]*na 
    atoms_out.cellvec[0][2]=atoms.cellvec[0][2]*na 
    atoms_out.cellvec[1][0]=atoms.cellvec[1][0]*nb
    atoms_out.cellvec[1][1]=atoms.cellvec[1][1]*nb
    atoms_out.cellvec[1][2]=atoms.cellvec[1][2]*nb
    atoms_out.cellvec[2][0]=atoms.cellvec[2][0]*nc
    atoms_out.cellvec[2][1]=atoms.cellvec[2][1]*nc
    atoms_out.cellvec[2][2]=atoms.cellvec[2][2]*nc
    atoms_out.nat = na*nb*nc*atoms.nat    
    atoms_out.boundcond="bulk"
    atoms_out.epot=atoms.epot
    return atoms_out
#****************************************************
str1 = "This script replicate a structure in the cell vectors direction."
str1+="\n The n_i are integers and a, b and c are the lattice vectors. The format of the input file must be ACF."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_input', action='store' ,type=str, help="fn_input is the name of the input file")
parser.add_argument('n1', type=int, help='scales a by n1 factor')
parser.add_argument('n2', type=int, help='scales b by n2 factor')
parser.add_argument('n3', type=int, help='scales c by n3 factor')
args=parser.parse_args()

atoms_all=acf_read(args.fn_input)
atoms=atoms_all[-1]
n1 = args.n1
n2 = args.n2
n3 = args.n3
atoms = replicate(n1,n2,n3,atoms)
atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
acf_write(atoms_all,"screen")
