#!/usr/bin/env python
import sys
import argparse
import atoms
from latvec2dproj import *
from ascii import *
from io_yaml import *

str1 = "This script read a file in the yaml dictionary format and write it in the ascii format."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in ascii format")
args=parser.parse_args()

atoms_all = read_yaml(args.fn_inp)

if len(atoms_all)==1:
    if not atoms_all[0].boundcond=="free":
        atoms_all[0].cellvec,atoms_all[0].rat=latvec2dproj(atoms_all[0].cellvec,atoms_all[0].rat,atoms_all[0].nat)
    ascii_write(atoms_all[0],args.fn_out)
else:
    print "\nATTENTION: The are more than one configuration in ACF file."
    prefix=raw_input("Please provide a prefix to generate files enumeratedly: [Default=tt]")
    if prefix=="": prefix="tt"
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        fnout="%s%5.5d.ascii" % (prefix,nconf)
        if not atoms.boundcond=="free":
            atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)
        ascii_write(atoms,fnout)
