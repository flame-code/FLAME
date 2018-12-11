#!/usr/bin/env python
import argparse
import atoms
from acf import *
from io_yaml import *
from io_bin import *

str1 = "This script read yaml dictionary and write it in acf format."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in binary format")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)

bohr2ang=0.52917721
for atoms in atoms_all:
    if atoms.units_length_io=='angstrom':
        atoms.cellvec[0][0]=atoms.cellvec[0][0]/bohr2ang
        atoms.cellvec[0][1]=atoms.cellvec[0][1]/bohr2ang
        atoms.cellvec[0][2]=atoms.cellvec[0][2]/bohr2ang
        atoms.cellvec[1][0]=atoms.cellvec[1][0]/bohr2ang
        atoms.cellvec[1][1]=atoms.cellvec[1][1]/bohr2ang
        atoms.cellvec[1][2]=atoms.cellvec[1][2]/bohr2ang
        atoms.cellvec[2][0]=atoms.cellvec[2][0]/bohr2ang
        atoms.cellvec[2][1]=atoms.cellvec[2][1]/bohr2ang
        atoms.cellvec[2][2]=atoms.cellvec[2][2]/bohr2ang
        for iat in range(atoms.nat):
            atoms.rat[iat][0]=atoms.rat[iat][0]/bohr2ang
            atoms.rat[iat][1]=atoms.rat[iat][1]/bohr2ang
            atoms.rat[iat][2]=atoms.rat[iat][2]/bohr2ang

bin_write(atoms_all,args.fn_out,1)
