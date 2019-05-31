#!/usr/bin/env python
import argparse
import atoms
from latvec2dproj import *
from ascii import *
from io_yaml import *

str1 = "This script reads a file in the yaml dictionary format and writes it in the ascii format."
parser = argparse.ArgumentParser(description=str1)
str2="Name of the output file, if there is one configuration, otherwise it is ignored."
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help=str2)
args=parser.parse_args()

atoms_all = read_yaml(args.fn_inp)

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
    if not atoms_all[0].boundcond=="free":
        atoms_all[0].cellvec,atoms_all[0].rat=latvec2dproj(atoms_all[0].cellvec,atoms_all[0].rat,atoms_all[0].nat)
    ascii_write(atoms_all[0],args.fn_out)
else:
    print "\nATTENTION: The are more than one configuration in YAML file. The given name for the output is ignored!"
    prefix=raw_input("Please provide a prefix to generate files enumeratedly: [Default=tt]")
    if prefix=="": prefix="tt"
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        fnout="%s%5.5d.ascii" % (prefix,nconf)
        if not atoms.boundcond=="free":
            atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)
        ascii_write(atoms,fnout)
