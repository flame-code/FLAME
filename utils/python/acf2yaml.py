#!/usr/bin/env python
import atoms
import argparse
from acf import *
from io_yaml import *

str1 = "This script read acf and force files and write them in a file with yaml format"
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in acf format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in yaml format")
parser.add_argument("-force",action='store_false',help="if present, force will be written")
args=parser.parse_args()
force=not args.force

atoms_all=acf_read(args.fn_inp)
for atoms in atoms_all:
    atoms.units_length_io='angstrom'
    atoms.epot=atoms.epot/27.211385

if force:
    fn_force="force_%s" % args.fn_inp
    forces=read_forces(atoms.nat,fn_force)
    if len(forces)!=len(atoms_all):
        print "ERROR: len(forces) force_%s does not match %s, number of confs." \
                % (args.fn_inp,args.fn_inp)
        exit()
    for iconf in range(len(forces)):
        for iat in range(atoms.nat):
            atoms_all[iconf].fat.append(forces[iconf][iat])

write_yaml(atoms_all,args.fn_out)
