#!/usr/bin/env python
import atoms
import argparse
from io_yaml import *

str1="Decompose a yaml file"
parser=argparse.ArgumentParser(description=str1)
str2="prefix of output files, if there are more than one configuration"
str2+="\n otherwise it is ignored. Default is tt"
parser.add_argument("-prefix",action='store',metavar='prefix',default='tt',help=str2)
parser.add_argument('fn_inp',action="store",type=str,help="Name of the input file in yaml format")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)

##if len(atoms_all)==1:
##    write_yaml(atoms_all[0],"tt00.yaml")
##else:
nconf=0
for atoms in atoms_all:
    atoms_all_conf=[]
    atoms_all_conf.append(Atoms())
    atoms_all_conf[-1]=copy.copy(atoms)
    nconf+=1
    fnout="%s%5.5d.yaml" % (args.prefix,nconf)
    write_yaml(atoms_all_conf,fnout)
