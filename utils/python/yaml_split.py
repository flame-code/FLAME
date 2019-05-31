#!/usr/bin/env python
import atoms
import argparse
from io_yaml import *

str1="To split a yaml file to n-file"
parser=argparse.ArgumentParser(description=str1)
str2="prefix of output files."
str2+="\n Default is tt"
parser.add_argument("-prefix",action='store',metavar='prefix',default='tt',help=str2)
parser.add_argument('fn_inp',action="store",type=str,help="Name of the input file in yaml format")
parser.add_argument('nfiles',action="store",type=int,help="Number of files")
args=parser.parse_args()

atoms_all=read_yaml(args.fn_inp)
nconf=len(atoms_all)
mod=nconf%args.nfiles
bb=nconf/args.nfiles

atoms_all_conf=[]
nf=0
for iconf in range(nconf):
    atoms_all_conf.append(Atoms())
    atoms_all_conf[-1]=copy.copy(atoms_all[iconf])
    if (iconf+1)%bb==0 or iconf==nconf-1:
        nf+=1
        fnout="%s%5.5d.yaml" % (args.prefix,nf)
        write_yaml(atoms_all_conf,fnout)
        atoms_all_conf=[]

if mod!=0:
    print """WARNING: The remainder from the division of nconf=%d by nfiles=%d is %d. 
    The extra conf(s) will be written in the %dth file""" %(nconf, args.nfiles, mod, args.nfiles+1)
