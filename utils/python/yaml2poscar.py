#!/usr/bin/env python
import sys
import atoms
import argparse
from io_yaml import *
from vasp import *

str1="Convert yaml conf(s) to POSCAR format."
str1+="\n if there is one configuration, default output file name is POSCAR"
parser=argparse.ArgumentParser(description=str1)
str2="prefix of output files, if there are more than one configuration"
str2+="\n otherwise it is ignored. Default is tt"
parser.add_argument("-prefix",action='store',metavar='prefix',default='tt',help=str2)
parser.add_argument('fn_inp',action="store",type=str,help="Name of the input file in yaml format")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)

if len(atoms_all)==1:
    poscar_write(atoms_all[0],"POSCAR")
else:
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        fnout="%s%5.5d" % (args.prefix,nconf)
        poscar_write(atoms,fnout)
