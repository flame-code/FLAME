#!/usr/bin/env python
import argparse
import commands
import string
import re
import sys
from io_yaml import *
from io_bin import *
import copy
#*****************************************************************************************
str1 = "This script selects diverse structures."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in acf format")
parser.add_argument('dtol', action='store' ,type=float, help="tolerance for distances")
args=parser.parse_args()

dist=[]
#f=open("tt1","r")
f=open("distall","r")
nn=-1
nconf=0
for line in f.readlines():
    #print "%25s%6d %25s%6d %20.10E%20.10E" % \
    #        (line.split()[0],int(line.split()[1]),line.split()[2], \
    #int(line.split()[3]),float(line.split()[4]),float(line.split()[5]))
    #con1.append(int(line.split()[1]))
    #con2.append(int(line.split()[3]))
    if not int(line.split()[1])==nn:
        #if nconf!=0: print
        nconf+=1
        dist.append([])
        for iconf in range(nconf-1):
            #dist.insert(iconf,dist[iconf][nconf-1])
            dist[-1].append(dist[iconf][nconf-1])
            #print "%20.10E" % dist[nconf-1][iconf],
        dist[-1].append(0.0)
        #print "%20.10E" % dist[-1][-1],
        nn=int(line.split()[1])
        #continue
    dist[-1].append(float(line.split()[4]))
    #print "%20.10E" % dist[-1][-1],
#../../minhop_data/poslow.     1 ../../minhop_data/poslow.     2     2.5533814195E+00    1.6950290032E-04
f.close()
#print
nconf+=1
dist.append([])
for iconf in range(nconf-1):
    #dist.insert(iconf,dist[iconf][nconf-1])
    dist[-1].append(dist[iconf][nconf-1])
    #print "%20.10E" % dist[nconf-1][iconf],
dist[-1].append(0.0)
#print "%20.10E" % dist[-1][-1],

#for iconf in range(nconf):
#    for jconf in range(nconf):
#        print "%20.10E" % dist[iconf][jconf],
#    print

sel=[]
for iconf in range(nconf):
    if iconf==0:
        sel.append(0) #This is the first configuration
        continue
    new=True
    for jconf in range(len(sel)):
        if dist[iconf][sel[jconf]]<args.dtol:
            new=False
            break
    if new==True:
        sel.append(iconf)
if args.fn_inp.endswith('.bin'):
    atoms_all=bin_read(args.fn_inp)
elif args.fn_inp.endswith('.yaml'):
    atoms_all=read_yaml(args.fn_inp)
atoms_all_out=[]
for iconf in range(len(sel)):
    atoms_all_out.append(Atoms)
    atoms_all_out[-1]=copy.deepcopy(atoms_all[sel[iconf]])

if args.fn_inp.endswith('.bin'):
    bin_write(atoms_all_out,args.fn_out,1)
elif args.fn_inp.endswith('.yaml'):
    write_yaml(atoms_all_out,args.fn_out)
