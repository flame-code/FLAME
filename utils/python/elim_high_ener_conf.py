#!/usr/bin/env python
import sys
import atoms
from acf import *

if len(sys.argv) < 5:
    print "usage: elim_high_ener_conf.py input_filename output_filename epot_base ener_ubound"
    exit()
else:
    filename_in = sys.argv[1]
    filename_out = sys.argv[2]
    epot_base = float(sys.argv[3])
    ener_ubound = float(sys.argv[4])

atoms_all=acf_read(filename_in)
ind=filename_in.rfind("/")
if ind==-1:
    force_filename_in="force_"+filename_in[ind+1:]
else:
    force_filename_in=filename_in[:ind+1]+"force_"+filename_in[ind+1:]
f=open(force_filename_in,"r")
iline=0
iconf=-1
for line in f.readlines():
    iline+=1
    #print iline,atoms_all[0].nat
    if (iline-1)%(atoms_all[0].nat+1)==0:
        iat=0
        iconf+=1
    else:
        atoms_all[iconf].fat.append([])
        atoms_all[iconf].fat[-1].append(float(line.split()[0]))
        atoms_all[iconf].fat[-1].append(float(line.split()[1]))
        atoms_all[iconf].fat[-1].append(float(line.split()[2]))
f.close()

atoms_all_sel=[]
for atoms in atoms_all:
    #print "%12.6f%12.6f" % (atoms.epot,epot_base+ener_ubound)
    if not atoms.epot>(epot_base+ener_ubound)*float(atoms.nat):
        atoms_all_sel.append(Atoms())
        atoms_all_sel[-1]=copy.copy(atoms)

if len(atoms_all_sel)>0:
    acf_write(atoms_all_sel,filename_out)
nconf_diff=len(atoms_all)-len(atoms_all_sel)
#print "Number of configurations in %s %6d" % (filename_in,len(atoms_all))
#print "Number of configurations too high in energy %6d" % nconf_diff
#print "Number of configurations written in %s %6d" % (filename_out,len(atoms_all_sel))
print "%40s -> total=%6d  elim=%6d  kept=%6d" % \
(filename_in,len(atoms_all),nconf_diff,len(atoms_all_sel))


#writing atomic forces
filename_force="force_"+filename_out
f=open(filename_force,"w")
nconf=0
for atoms in atoms_all_sel:
    nconf+=1
    f.write("configuration %5.5d\n" % (nconf))
    for iat in range(atoms.nat):
        fx=atoms.fat[iat][0]
        fy=atoms.fat[iat][1]
        fz=atoms.fat[iat][2]
        #f.write("%17.8f%17.8f%17.8f\n" % (fx,fy,fz))
        f.write("%19.10E%19.10E%19.10E\n" % (fx,fy,fz))
f.close()

