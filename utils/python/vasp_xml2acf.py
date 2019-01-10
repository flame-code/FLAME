#!/usr/bin/env python
import sys
from atoms import *
import copy
from acf import *
#*****************************************************************************************
lastconf=False
if len(sys.argv) < 3:
    print "usage: vasp_xml2acf.py vasprun.xml output_filename"
    exit()
elif len(sys.argv)==4:
    if sys.argv[1]=='-last':
        lastconf=True
    else:
        print "ERROR: wrong command option: %s" % sys.argv[1]
        exit()
    finp = sys.argv[2]
    fout = sys.argv[3]
elif len(sys.argv)==3:
    finp = sys.argv[1]
    fout = sys.argv[2]
else:
    print "ERROR: only two command options are allowed."
    exit()

f=open(finp,"r")
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
        atoms.sat.append(sat[iline_positions])
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
        atoms.fat[-1].append(float(line.split()[1]))
        atoms.fat[-1].append(float(line.split()[2]))
        atoms.fat[-1].append(float(line.split()[3]))
        iline_forces+=1
        if iline_forces==nat: iline_forces=-1
        continue
    #-------------------------------------------------------
    if str_line.find('e_fr_energy')==12:
        e_fr_energy=float(line.split()[2])
    #-------------------------------------------------------
    if str_line.find('e_wo_entrp')==12:
        atoms.epot=float(line.split()[2])
        ediff=abs(1000.0*(e_fr_energy-atoms.epot)/float(nat))
        if ediff>1.0:
            print "WARNING: Difference between energy and free energy in (meV/atom): %6.3f" % ediff
        conf_complete=True
    #-------------------------------------------------------
    if conf_complete:
        atoms_all.append(Atoms())
        atoms_all[-1]=copy.copy(atoms)
        conf_complete=False
    #-------------------------------------------------------


#for atoms in atoms_all:
#    print atoms.nat
#    for iat in range(atoms.nat):
#        print "    <v>%17.8f%17.8f%17.8f </v>" % (atoms.rat[iat][0],atoms.rat[iat][1],atoms.rat[iat][2])
f.close()

if lastconf==True:
    atoms_lastconf=[]
    atoms_lastconf.append(Atoms())
    atoms_lastconf[-1]=copy.copy(atoms_all[-1])
    #writing last configuration into file
    acf_write(atoms_lastconf,fout)
else:
    #writing all configurations into file
    acf_write(atoms_all,fout)
    #writing atomic forces
    filename_force="force_"+fout
    f=open(filename_force,"w")
    tt=27.211385/0.52917721
    nconf=0
    for atoms in atoms_all:
        nconf+=1
        f.write("configuration %5.5d\n" % (nconf))
        for iat in range(atoms.nat):
            fx=atoms.fat[iat][0]/tt
            fy=atoms.fat[iat][1]/tt
            fz=atoms.fat[iat][2]/tt
            #f.write("%17.8f%17.8f%17.8f\n" % (fx,fy,fz))
            f.write("%19.10E%19.10E%19.10E\n" % (fx,fy,fz))
    f.close()
#*****************************************************************************************
