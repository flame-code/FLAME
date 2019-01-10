#!/usr/bin/env python
import argparse
import atoms
from xyz import *
from io_yaml import *

str1 = "This script reads a file in the yaml dictionary format and writes it in xyz format."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in yaml format")
parser.add_argument('all_conf', type=str, help='True or False, True: write all configurations in the output file, default value: True', 
        nargs='?', const=1, default=True)
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)

if len(atoms_all)==1 or args.all_conf==True:
    xyz_write_b(atoms_all,'bigdft', args.fn_out)
else:
    print "\nATTENTION: The are more than one configuration in YAML file. The given name for the output is ignored!"
    prefix=raw_input("Please provide a prefix to generate files enumeratedly: [Default=tt]")
    if prefix=="": prefix="tt"
    nconf=0
    for atoms in atoms_all:
        atoms_all_conf=[]
        atoms_all_conf.append(Atoms())
        atoms_all_conf[-1]=copy.copy(atoms)
        nconf+=1
        fnout="%s%5.5d.xyz" % (prefix,nconf)
        xyz_write_b(atoms_all_conf,'bigdft', fnout)
