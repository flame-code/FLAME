#!/usr/bin/env python
import atoms
import argparse
from ascii import *
from io_yaml import *

str1 = "This script read an ascii file and write it in the yaml format"
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in acf format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in yaml format")
args=parser.parse_args()

atoms=ascii_read(args.fn_inp)
atoms.units_length_io='angstrom'

atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
write_yaml(atoms_all,args.fn_out)
