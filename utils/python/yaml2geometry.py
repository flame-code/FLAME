#!/usr/bin/env python
import atoms
import argparse
from io_yaml import *
from aims import *

str1="Convert yaml conf(s) to geometry format of FHI-aims."
str1+="\n if there is one configuration, default output file name is geometry.in"
parser=argparse.ArgumentParser(description=str1)
str2="prefix of output files, if there are more than one configuration"
str2+="\n otherwise it is ignored. Default is tt"
parser.add_argument("-prefix",action='store',metavar='prefix',default='tt',help=str2)
parser.add_argument('fn_inp',action="store",type=str,help="Name of the input file in yaml format")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)

bohr2ang=0.52917721
for atoms in atoms_all:
    if atoms.units_length_io=='atomic':
        atoms.cellvec[0][0]=atoms.cellvec[0][0]*bohr2ang
        atoms.cellvec[0][1]=atoms.cellvec[0][1]*bohr2ang
        atoms.cellvec[0][2]=atoms.cellvec[0][2]*bohr2ang
        atoms.cellvec[1][0]=atoms.cellvec[1][0]*bohr2ang
        atoms.cellvec[1][1]=atoms.cellvec[1][1]*bohr2ang
        atoms.cellvec[1][2]=atoms.cellvec[1][2]*bohr2ang
        atoms.cellvec[2][0]=atoms.cellvec[2][0]*bohr2ang
        atoms.cellvec[2][1]=atoms.cellvec[2][1]*bohr2ang
        atoms.cellvec[2][2]=atoms.cellvec[2][2]*bohr2ang
        for iat in range(atoms.nat):
            atoms.rat[iat][0]=atoms.rat[iat][0]*bohr2ang
            atoms.rat[iat][1]=atoms.rat[iat][1]*bohr2ang
            atoms.rat[iat][2]=atoms.rat[iat][2]*bohr2ang

if len(atoms_all)==1:
    aims_write(atoms_all[0],"geometry.in")
else:
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        fnout="%s%5.5d.in" % (args.prefix,nconf)
        aims_write(atoms,fnout)
