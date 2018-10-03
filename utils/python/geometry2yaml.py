#!/usr/bin/env python
import argparse
import atoms 
from aims import *
from io_yaml import *

str1 = "This script read a file in FHI-aims (geometry.in) format and write it in Yaml dictionary format."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in geometry.in format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in the Yaml format")
args=parser.parse_args()

atoms = aims_read(args.fn_inp)
atoms.units_length_io='angstrom'

atoms_all=[]
atoms_all.append(Atoms())
atoms_all[-1]=copy.copy(atoms)
write_yaml(atoms_all,args.fn_out)
