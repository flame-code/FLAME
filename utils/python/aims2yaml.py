#!/usr/bin/env python

import argparse
import sys
import copy
from atoms import *
from aims import *
from io_yaml import *
from latvec2dproj import *
#*****************************************************************************************
str1="Extract configurations during MD or geometry optimization from FHI-aims output."
parser=argparse.ArgumentParser(description=str1)
parser.add_argument("-last",action='store_false',help="if present, only last configuration is written")
parser.add_argument('fn_input',action="store",type=str,help="Name of the input file : FHI-aims standard output")
parser.add_argument('fn_output',action="store",type=str,help="Name of the output file in yaml format")
args=parser.parse_args()
lastconf=not args.last
#-----------------------------------------------------------------------------------------
try:
    f=open(args.fn_input,'r')
    lines=f.readlines()
    f.close()
except IOError:
    print("ERROR: Failed to open file %s - check file name and try again." % (args.fn_input))
    sys.exit(0)

atoms_all=[]
conf_complete=False
has_unit_cell=False
iskip=0
nat=0
ev2bohr=27.211385/0.52917721
Ehar=27.211385 #eV

for iline,line in enumerate(lines):
    if iskip>0:
        iskip=-1
        continue
    #-------------------------------------------------------
    str_line=str(line)
    if nat==0:
        if 'Number of atoms                   ' in str_line: nat=int(line.split()[5])
        continue
    #-------------------------------------------------------
    #looks for input geometry
    if 'Input geometry' in str_line:
        if 'Unit cell' in lines[iline+1]: has_unit_cell=True
        #print str_line
        atoms=get_input_geometry(iline,lines,nat,has_unit_cell)
        #rotation to lattice and all atom
        ##if not atoms.boundcond=="free":
        ##    atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)

        #for iat in range(atoms.nat):
        iskip=6+nat
        continue
    #-------------------------------------------------------
    #looks for an updated geometry
    if 'Updated atomic structure' in str_line:
        #print str_line
        atoms=get_updated_geometry(iline,lines,nat,has_unit_cell)
        #rotation to lattice and all atom
        ##if not atoms.boundcond=="free":
        ##    atoms.cellvec,atoms.rat=latvec2dproj(atoms.cellvec,atoms.rat,atoms.nat)

        iskip=5+nat
        continue
    #-------------------------------------------------------
    if 'Total energy uncorrected' in str_line:
        epot_uncorrected=float(line.split()[5])
        continue
    #-------------------------------------------------------
    if 'Total energy corrected' in str_line:
        epot_corrected=float(line.split()[5])
        ediff=abs(1000.0*(epot_corrected-epot_uncorrected)/float(nat))
        if ediff>1.0:
            print("WARNING: difference between energy and free energy in meV/atom: %6.3f" % ediff)
        atoms.epot=epot_uncorrected/Ehar
        continue
    #-------------------------------------------------------
    if 'Charged system requested' in str_line:
        qtot = float(line.split()[5].rstrip('.'))
    elif 'Charge =' in str_line:
        qtot = float(line.split()[2].rstrip(':'))
        continue
    #-------------------------------------------------------
    if 'Total atomic forces' in str_line:
        #print atoms.sat[-1]
        for iat in range(nat):
            atoms.fat.append([])
            atoms.fat[-1].append(float(lines[iline+1+iat].split()[2])/ev2bohr)
            atoms.fat[-1].append(float(lines[iline+1+iat].split()[3])/ev2bohr)
            atoms.fat[-1].append(float(lines[iline+1+iat].split()[4])/ev2bohr)
        conf_complete=True
        iskip=nat
    #-------------------------------------------------------
    if conf_complete:
        atoms_all.append(Atoms())
        atoms.qtot = qtot
        atoms.units_length_io='angstrom'
        atoms_all[-1]=copy.copy(atoms)
        conf_complete=False
        del atoms
    #-------------------------------------------------------
#End of reading from input file.
#-----------------------------------------------------------------------------------------
for atoms in atoms_all:
    if atoms.boundcond=='free':
        atoms.cellvec,atoms.rat=set_cell(atoms)

if lastconf:
    atoms_lastconf=[]
    atoms_lastconf.append(Atoms())
    atoms_lastconf[-1]=copy.copy(atoms_all[-1])
    #writing last configuration into file
    write_yaml(atoms_lastconf,args.fn_output)
else:
    #writing all configurations into file
    write_yaml(atoms_all,args.fn_output)
#*****************************************************************************************
