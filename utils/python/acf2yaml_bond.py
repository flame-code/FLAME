#!/usr/bin/env python
import atoms
import argparse
from acf import *
from io_yaml import *

str1 = "This script reads acf then multiplies coordinates with scalling_factor and writes it in a file with yaml format."
#RZX
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in acf format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in yaml format")
parser.add_argument("scale_factor",action='store',type=float,help="Bond length will be multiplied by this number")
args=parser.parse_args()

atoms_all=acf_read(args.fn_inp)
for atoms in atoms_all:
    atoms.units_length_io='angstrom'
    atoms.epot=atoms.epot/27.211385
for iconf in range(len(atoms_all)):
    for iat in range(atoms_all[iconf].nat):
        for i in range(3):
            atoms_all[iconf].rat[iat][i]=atoms_all[iconf].rat[iat][i]*args.scale_factor
    atoms_all[iconf].cellvec[0][0]=atoms_all[iconf].cellvec[0][0]*args.scale_factor
    atoms_all[iconf].cellvec[0][1]=atoms_all[iconf].cellvec[0][1]*args.scale_factor
    atoms_all[iconf].cellvec[0][2]=atoms_all[iconf].cellvec[0][2]*args.scale_factor
    atoms_all[iconf].cellvec[1][0]=atoms_all[iconf].cellvec[1][0]*args.scale_factor
    atoms_all[iconf].cellvec[1][1]=atoms_all[iconf].cellvec[1][1]*args.scale_factor
    atoms_all[iconf].cellvec[1][2]=atoms_all[iconf].cellvec[1][2]*args.scale_factor
    atoms_all[iconf].cellvec[2][0]=atoms_all[iconf].cellvec[2][0]*args.scale_factor
    atoms_all[iconf].cellvec[2][1]=atoms_all[iconf].cellvec[2][1]*args.scale_factor
    atoms_all[iconf].cellvec[2][2]=atoms_all[iconf].cellvec[2][2]*args.scale_factor
write_yaml(atoms_all,args.fn_out)
