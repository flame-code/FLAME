#!/usr/bin/env python

import os,shutil ,glob, errno
import numpy as np
#*****************************
#dir="./final"
#if os.path.isdir(dir):
#    shutil.rmtree(dir)
#    os.chdir(dir)    
#os.mkdir(dir)
ofile = open("earr_total.txt","w")
#*****************************naming columns of the earr.dat file 
data=np.ndarray(shape=[0],dtype=[('dir','S8'),('ID',int),('energy',float),('col3',float),('col4',int),('col5',int),('col6',float),('col7',float)])
write=np.ndarray(shape=[0],dtype=[('dir','S8'),('ID',int),('energy',float),('col3',float),('col4',int),('col5',int),('col6',float),('col7',float)])
#data=[[None]*7]
#*****************************
def genfield(directory):
    for line in open(directory+'/earr.dat','r'):
        yield directory+' ' + line
#*****************************read all earr.dat files of each run* and sort by energy value
for file in np.sort(glob.glob("run*/")):
    file=file[:-1]
    print file
    #os.chdir(file)
    
    #with open("earr.dat","r") as ifile:
    tmp=np.genfromtxt(genfield(file),skip_header=2,dtype=[('dir','S8'),('ID',int),('energy',float),('col3',float),('col4',int),('col5',int),('col6',float),('col7',float)])
    if tmp.ndim==0:
        tmp=[tmp]
    data=np.concatenate((data,tmp))
    #print data
    #lines = ifile.readlines()
    #print lines
    #os.chdir("..")
data.sort(order='energy')
#*****************************make directory to save final poslow files 
try:
    os.mkdir('final-earr')
except OSError as exc:
    if exc.errno == errno.EEXIST and os.path.isdir('final-earr'):
        pass
    else:
        raise
#*****************************compare and merge according to both energy and SPG values 
counter=0
for row_i in data:
    for row_j in write:
        if abs(row_i['energy']-row_j['energy'])<1e-4 and row_i['col5']==row_j['col5']:
            break
    else:
        counter+=1
        shutil.copyfile('%s/poslow00%03d.ascii'%(row_i['dir'],row_i['ID']),'final-earr/poslow%05d.ascii'%counter)
        write=np.append(write,row_i)
    #if abs(row['energy']-last_energy)<1e-4 and row['col5']==last_spg:
    #    print str(row)+' del'
    #else:
    #    print str(row)
    #last_energy=row['energy']
    #last_spg=row['col5']
#print data
#*****************************save information about final accepted structures
with open('earr_total.txt', 'a') as f:
    np.savetxt(f, write, fmt='%s %6d %.15f %.15f %5d %5d %.5f %5f')
    
#ifile.close(); ofile.close()
#f[f[:, 1].argsort()]
      #  f.sort()
#f[np.argsort(f[:,1])]
#f.sort(key=lambda s:(s['2'],s['5']))
f.closed
shutil.move('earr_total.txt','final-earr')
print "finish"
    
#ifile=open("earr.dat","r")
#ofile.write(data)
#np.ndarray.tofile(ofile,sep="",format="%s")

#with open('workfile', 'r') as f:
#...     read_data = f.read()
#f.closed
