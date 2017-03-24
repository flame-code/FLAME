#!/usr/bin/env python
import sys
import atoms
import copy
import random
import argparse
from acf import *
#*****************************************************************************************
str1="Generate configuration by mutation: permutation of atoms"
str2="label of output files, it is useful to include it for future grep commands,"
str2+="\n otherwise it is ignored. Default is unknown"
parser=argparse.ArgumentParser(description=str1)
parser.add_argument("-label",action='store',metavar='label',default='unknown',help=str2)
parser.add_argument('fn_input',action="store",type=str,help="fn_input is the name of the input file")
parser.add_argument('fn_output',action="store",type=str,help="fn_output is the name of the output file")
args=parser.parse_args()
#-----------------------------------------------------------------------------------------
try:
    f=open(args.fn_input,'r')
    #lines=f.readlines()
    f.close()
except IOError:
    print "ERROR: Failed to open file %s - check file name and try again." % (args.fn_input)
    sys.exit(0)

print 'Reading configurations from file %s ... ' % args.fn_input, 
atoms_all=acf_read(args.fn_input)
print 'done.'
#-------------------------------------------------------------------------------
random.seed()
atoms_all_out=[]
for atoms in atoms_all:
    imut=0
    itry=0
    atoms_t=copy.copy(atoms)
    while(imut<10):
        itry+=1
        iat=random.randint(0,atoms_t.nat-1)
        jat=random.randint(0,atoms_t.nat-1)
        #print "%5s %5s %s" % (atoms_t.sat[iat],atoms_t.sat[iat],(atoms_t.sat[iat]=="Sn"))
        correct_types_1=(atoms_t.sat[iat]=="Sn" and atoms_t.sat[jat]=="Zn")
        correct_types_2=(atoms_t.sat[iat]=="Zn" and atoms_t.sat[jat]=="Sn")
        #print "%5s %5s %s" % (atoms_t.sat[iat],atoms_t.sat[iat],(correct_types_1 or correct_types_2))
        if(not (correct_types_1 or correct_types_2)): continue
        #print "%5s %5s" % (atoms_t.sat[iat],atoms_t.sat[jat])
        print "Exchanging atomic positions, itry=%6d     %5s:%5d    %5s:%5d " % (itry,atoms_t.sat[iat],iat,atoms_t.sat[jat],jat)
        #print "Exchanging atomic positions, %10f    %10f   %10f " % (atoms_t.rat[iat][0],atoms_t.rat[iat][1],atoms_t.rat[iat][2])
        #print "Exchanging atomic positions, %10f    %10f   %10f " % (atoms_t.rat[jat][0],atoms_t.rat[jat][1],atoms_t.rat[jat][2])
        #print "-------------------------------------------------"
        imut=imut+1
        x=atoms_t.rat[iat][0]
        y=atoms_t.rat[iat][1]
        z=atoms_t.rat[iat][2]
        atoms_t.rat[iat][0]=atoms_t.rat[jat][0]
        atoms_t.rat[iat][1]=atoms_t.rat[jat][1]
        atoms_t.rat[iat][2]=atoms_t.rat[jat][2]
        atoms_t.rat[jat][0]=x
        atoms_t.rat[jat][1]=y
        atoms_t.rat[jat][2]=z
        #print "Exchanging atomic positions, %10f    %10f   %10f " % (atoms_t.rat[iat][0],atoms_t.rat[iat][1],atoms_t.rat[iat][2])
        #print "Exchanging atomic positions, %10f    %10f   %10f " % (atoms_t.rat[jat][0],atoms_t.rat[jat][1],atoms_t.rat[jat][2])
        #print "-------------------------------------------------"
        atoms_all_out.append(Atoms())
        atoms_all_out[-1]=copy.deepcopy(atoms_t)
    del atoms_t

#-------------------------------------------------------------------------------
print 'Writing configurations to file %s ... ' % args.fn_output, 
acf_write(atoms_all_out,args.fn_output,labelpatt=args.label)
print 'done.'
