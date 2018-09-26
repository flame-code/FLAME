#!/usr/bin/env python
import sys
import atoms
import argparse
from io_yaml import *
from aims import *

str1 = "This script read a file in yaml format (its last conf.) and write it in FHI-aims (geometry.in) format"
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in geometry.in format")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)
aims_write(atoms_all[-1],args.fn_out)
