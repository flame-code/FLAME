#!/usr/bin/env python
import argparse
import atoms
from acf import *
from io_yaml import *

str1 = "This script read yaml dictionary and write it in acf format."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in acf format")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)
for atoms in atoms_all:
    atoms.epot=atoms.epot*27.211385
acf_write(atoms_all,args.fn_out,"yaml2acf")
