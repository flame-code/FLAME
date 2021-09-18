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
S_atoms_all_free=[]
S_atoms_all_slab=[]
S_atoms_all_bulk=[]
i = 0
i_free = 0
i_slab = 0
i_bulk = 0
for iconf in range(len(atoms_all)):
    if atoms_all[iconf].nat>args.S_nat:
        i=i+1
        if atoms_all[iconf].boundcond == 'free':
            S_atoms_all_free.append(atoms_all[iconf])
            i_free=i_free+1
        if atoms_all[iconf].boundcond == 'slab':
            S_atoms_all_slab.append(atoms_all[iconf])
            i_slab=i_slab+1
        if atoms_all[iconf].boundcond == 'bulk':
            S_atoms_all_bulk.append(atoms_all[iconf])
            i_bulk=i_bulk+1
print("number of total structures with %d number of atoms is: %d"%(args.S_nat,i))
print("number of free structures with %d number of atoms is: %d"%(args.S_nat,i_free))
print("number of slab structures with %d number of atoms is: %d"%(args.S_nat,i_slab))
print("number of bulk structures with %d number of atoms is: %d"%(args.S_nat,i_bulk))
free_name = args.fn_out.strip('.yaml')+'_free.yaml'
slab_name = args.fn_out.strip('.yaml')+'_slab.yaml'
bulk_name = args.fn_out.strip('.yaml')+'_bulk.yaml'
write_yaml(S_atoms_all_free,free_name)
write_yaml(S_atoms_all_slab,slab_name)
write_yaml(S_atoms_all_bulk,bulk_name)

