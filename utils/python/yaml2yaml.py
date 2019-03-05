#!/usr/bin/env python
import argparse
import atoms
from io_yaml import *
str1 = "This script selects structures with specific number of atoms from yaml input and writes them in yaml format."
#RZX
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in acf format")
parser.add_argument("S_nat",action='store',type=int,help="Specified number of atoms")
args=parser.parse_args()
atoms_all=read_yaml(args.fn_inp)
S_atoms_all=[]
i = 0
for iconf in range(len(atoms_all)):
    if atoms_all[iconf].nat==args.S_nat:
        S_atoms_all.append(atoms_all[iconf])
        i=i+1
print "number of structures with %d number of atoms is: %d"%(args.S_nat,i)
write_yaml(S_atoms_all,args.fn_out)

