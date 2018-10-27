#!/usr/bin/env python

import argparse
from io_yaml import *
import atoms
import random

str1 = "This script shuffle a file in yaml format and write it in another yaml file."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in yaml format")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)
random.shuffle(atoms_all)
write_yaml(atoms_all,args.fn_out)
