#!/usr/bin/env python
import argparse
import copy
from atoms import *
from io_yaml import *
#*****************************************************************************************
str1="Extract configurations during MD or geometry optimization from VASP output (vasprun.xml)."
parser=argparse.ArgumentParser(description=str1)
parser.add_argument("-last",action='store_false',help="if present, only last configuration will be written")
parser.add_argument('fn_input',action="store",type=str,help="Name of the input file : vasprun.xml")
parser.add_argument('fn_output',action="store",type=str,help="Name of the output file in yaml format")
args=parser.parse_args()
lastconf=not args.last

f=open(args.fn_input,"r")
atoms_all=[]
ntypes=[]
stypes=[]
sat=[]
iline_basis=-1
iline_atomtypes=-1
iline_positions=-1
iline_forces=-1
conf_complete=False
nat=0
ev2bohr=27.211385/0.52917721
Ehar=27.211385 #eV

for line in f.readlines():
    if iline_atomtypes>-1 and iline_atomtypes<7:
        iline_atomtypes+=1
        continue
    #-------------------------------------------------------
    str_line=str(line)
    if nat==0:
        ind=str_line.find('atoms')
        if ind==3: nat=int(line.split()[1])
        continue
    #looks for keyword atomstypes
    if str_line.find('atomtypes')==15:
        iline_atomtypes=0
        #read_conf=True
        continue
    #looks for keyword basis
    if str_line.find('basis')==18:
        iline_basis=0
        #read_conf=True
        continue
    #looks for keyword positions
    if str_line.find('positions')==17:
        iline_positions=0
        continue
    #looks for keyword forces
    if str_line.find('forces')==16:
        iline_forces=0
        continue
    #-------------------------------------------------------
    if iline_atomtypes>-1:
        if str_line.find('set')==5:
            for itype in range(len(ntypes)):
                for iat in range(ntypes[itype]):
                    sat.append(stypes[itype])
                    #print sat[-1]
            iline_atomtypes=-1
            continue
        ntypes.append(int(str_line[11:15]))
        stypes.append(str_line[22:24])
        #print ntypes[-1],stypes[-1]
        continue
    #-------------------------------------------------------
    if iline_basis>-1:
        if iline_basis==0: atoms=Atoms()
        #print iline_basis
        atoms.cellvec[iline_basis][0]=float(line.split()[1])
        atoms.cellvec[iline_basis][1]=float(line.split()[2])
        atoms.cellvec[iline_basis][2]=float(line.split()[3])
        iline_basis+=1
        if iline_basis==3: iline_basis=-1
        continue
    #-------------------------------------------------------
    if iline_positions>-1:
        if iline_positions==0:
            atoms.nat=nat
            atoms.boundcond="free"
        atoms.sat.append(sat[iline_positions].strip())
        #print atoms.sat[-1]
        atoms.rat.append([])
        xred=float(line.split()[1])
        yred=float(line.split()[2])
        zred=float(line.split()[3])
        x=xred*atoms.cellvec[0][0]+yred*atoms.cellvec[1][0]+zred*atoms.cellvec[2][0]
        y=xred*atoms.cellvec[0][1]+yred*atoms.cellvec[1][1]+zred*atoms.cellvec[2][1]
        z=xred*atoms.cellvec[0][2]+yred*atoms.cellvec[1][2]+zred*atoms.cellvec[2][2]
        atoms.rat[-1].append(x)
        atoms.rat[-1].append(y)
        atoms.rat[-1].append(z)
        atoms.bemoved.append("TTT")
        iline_positions+=1
        if iline_positions==nat: iline_positions=-1
        continue
    #-------------------------------------------------------
    if iline_forces>-1:
        #print atoms.sat[-1]
        atoms.fat.append([])
        atoms.fat[-1].append(float(line.split()[1])/ev2bohr)
        atoms.fat[-1].append(float(line.split()[2])/ev2bohr)
        atoms.fat[-1].append(float(line.split()[3])/ev2bohr)
        iline_forces+=1
        if iline_forces==nat: iline_forces=-1
        continue
    #-------------------------------------------------------
    if str_line.find('e_fr_energy')==12:
        e_fr_energy=float(line.split()[2])
    #-------------------------------------------------------
    if str_line.find('e_wo_entrp')==12:
        atoms.epot=float(line.split()[2])/Ehar
        ediff=abs(1000.0*(e_fr_energy-atoms.epot*Ehar)/float(nat))
        if ediff>1.0:
            print "WARNING: Difference between energy and free energy in (meV/atom): %6.3f" % ediff
        conf_complete=True
    #-------------------------------------------------------
    if conf_complete:
        atoms_all.append(Atoms())
        atoms.units_length_io='angstrom'
        atoms_all[-1]=copy.copy(atoms)
        conf_complete=False
    #-------------------------------------------------------
f.close()

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
