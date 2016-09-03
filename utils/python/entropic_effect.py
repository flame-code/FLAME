#!/usr/bin/python
import atoms
import math
import sys
amass=((65.38+15.999)/2.)*1822.88848 #atomic mass
boltzmannconstant=3.16680965e-6
#boltzmannconstant=1
#***********************************************************************
if len(sys.argv) < 4:
    print "usage: entropic_effects.py input_filename number_of_atoms epot"
    exit()
else:
    filename = sys.argv[1]
    nat = int(sys.argv[2])
    epot = float(sys.argv[3])
#***********************************************************************
f = open (filename,"r")
iline=0
landa=[]
for line in f.readlines():
    if iline==nat*3-6: break
    landa.append(float(line.split()[1]))
    iline+=1
    #print landa[iline-1]
#***********************************************************************
freq=[]
for i in range(0,3*nat-6):
    #print landa[i]
    freq.append(float(math.sqrt(landa[i]/amass)))
    #print freq
zpe=0.
for i in range(0,3*nat-6):
    if freq[i]<10.: zpe+=freq[i]*0.5
#print zpe
for itemp in range(1,500,2):
    temp=2*itemp
    fenergy=0.
    tt1=boltzmannconstant*temp
    for i in range(0,3*nat-6):
        if freq[i]<10.: 
            tt2=freq[i]/tt1
        fenergy+=math.log(1.- math.exp(-tt2))
    fenergy=zpe+epot/27.2114+tt1*fenergy
    print temp, fenergy*27.2114
f.closed
#********************************************************************** 
