#!/usr/bin/env python
import commands
import string
import re
import sys
#*****************************************************************************************
if len(sys.argv) < 3:
    print "usage: pickdiffconfs.py \"grep -B 6 -A 63 poslow%5.5d ../../minhop_data/poslow.acf %s\" dtol"
    print "       63 is number of atoms plus 3"
    print "       ../../minhop_data/poslow.acf is filename of configurations"
    exit()
else:
    command = sys.argv[1]
    dtol = float(sys.argv[2])

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
        if dist[iconf][sel[jconf]]<dtol:
            new=False
            break
    if new==True:
        sel.append(iconf)


command=command+"\n"
f=open("r1.sh","w")
for iconf in range(len(sel)):
    #print "grep -B 6 -A 63 poslow%5.5d ../../../minhop_data/poslow.acf" % (sel[iconf]+1)
    if iconf==0:
        tt_com=">all.acf"
    else:
        tt_com=">>all.acf"
    #f.write("grep -B 6 -A 63 poslow%5.5d ../../minhop_data/poslow.acf %s\n" % (sel[iconf]+1,tt_com))
    f.write(command % (sel[iconf]+1,tt_com))
f.close()
#1
#88
#213
#for iconf in range(nconf):
#    if dist[0][iconf]>9.0 and dist[87][iconf]>9.0 and dist[212][iconf]>9.0: print "wrong"
