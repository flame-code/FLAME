#!/usr/bin/env python
import argparse
import sys
import atoms
from acf import *
from aims import *
#*****************************************************************************************
str1="Convert ACF files to geometry format of FHI-aims."
str1+="\n if there is one configuration, default output file name is geometry.in"
parser=argparse.ArgumentParser(description=str1)
str2="prefix of output files, if there are more than one configuration"
str2+="\n otherwise it is ignored. Default is tt"
parser.add_argument("-prefix",action='store',metavar='prefix',default='tt',help=str2)
parser.add_argument('fn_input',action="store",type=str,help="fn_input is the name of the input file")
args=parser.parse_args()
#lastconf=not args.last
#if args.prefix: print args.prefix
#print args.prefix

atoms_all=acf_read(args.fn_input)

if len(atoms_all)==1:
    aims_write(atoms_all[0],"geometry.in")
else:
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        fnout="%s%5.5d.in" % (args.prefix,nconf)
        aims_write(atoms,fnout)

