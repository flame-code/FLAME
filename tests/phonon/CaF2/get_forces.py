import numpy as np
#from munkres import Munkres
import yaml
import copy
import sys
import os.path
sys.path.insert(1,'../../../utils/python')
from io_yaml import *

def force_acf(filename):
	if not os.path.isfile(filename):
		print "%s file does not exist" %(filename)
	fin=open(filename,'r')
	commen = fin.readline().strip()
	commen = fin.readline().strip()
	commen = fin.readline().strip()
	commen = fin.readline().strip()
	commen = fin.readline().strip()
	commen = fin.readline().strip()
	commen = fin.readline().strip()
	nat = int(fin.readline().split()[0])
	print "Nat found",nat 
	commen = fin.readline().strip()
	commen = fin.readline().strip()
	force=[]
	for iat in range(nat):
		fcart = [ float(v) for v in fin.readline().split()[1:4] ]
		force.append(fcart)
#	print force
	return force	

    


#f=file("mesh.yaml","r").read()
f=file("disp.yaml","r").read()
datamap = yaml.load(f,Loader=yaml.CLoader)
fout = open('FORCE_SETS','w')
fout.write("%d \n" % (datamap['natom']))
fout.write("%d \n" % (len(datamap['displacements'])))
lat = np.array(datamap['lattice'])
for arg,dat in zip(sys.argv[1:],datamap['displacements']):
    fout.write(" \n" )
    iat=dat['atom']-1
    try:
        redpos=datamap['atoms'][iat]['position']
    except:
        redpos=datamap['points'][iat]['coordinates']
#    print arg,dat['displacement'],dat['atom'],np.add(np.dot(lat.T,redpos),dat['displacement'])
    print "Extracting forces from", arg
    #force = force_acf(arg)
    atoms_all=read_yaml(arg) #'posout.yaml')
    force=atoms_all[0].fat
    fout.write("%d \n" % (dat['atom']))
   # displ=np.add(np.dot(lat.T,redpos),dat['displacement'])
    displ=dat['displacement']
    fout.write("%21.15f %21.15f %21.15f\n"  %(displ[0],displ[1],displ[2]))
    m = 27.211396132/0.529177249
    for fcart in force:
        fout.write("%21.15f %21.15f %21.15f\n"  %(fcart[0]*m,fcart[1]*m,fcart[2]*m))
    

    
#Explicitly compute the displacement of the atom from disp.yaml
    
#f2 = open('plt2','w')
#
#
#a= datamap['phonon'][0]['band']
#band_dict={}
#for band,iband in zip(a,range(len(a))):
#   f1.write("%21.15f  bandindex: %d qvec: %d \n" % (band['frequency'],iband,0))
#   f2.write("%21.15f  bandindex: %d qvec: %d \n" % (band['frequency'],iband,0))
#   band_dict[iband]=iband
#
#band_dict_old=copy.deepcopy(band_dict)
#
##for k in range(1,6):
#for k in range(1,50):
#   a= datamap['phonon'][k-1]['band']
#   b= datamap['phonon'][k]['band']
#   nband=len(a)
#   smat=np.zeros(shape=(nband,nband))
#   for i,ii in zip(a,range(len(a))):
#      for j,jj in zip(b,range(len(b))):
#          smat[ii,jj]=abs(np.dot(np.array(i['eigenvector']).flatten(),np.array(j['eigenvector']).flatten()))
#   m=Munkres()
#   indexes=m.compute(-smat)
#   total = 0.0
#   iindexes=[0]*len(b)
#   for row, column in indexes:
#       value = smat[row][column]
#       total+= value
#       iindexes[row]=column
#   for band,iband in zip(b,range(len(b))):
#      f1.write("%21.15f  bandindex: %d qvec: %d \n" % (band['frequency'],iband,k))
##update dictionary  
#   for i in range(len(b)):
#      band_dict[i]=band_dict_old[iindexes[i]]
#   for i in range(len(b)):
#      f2.write("%21.15f  bandindex: %d qvec: %d \n" % (b[i]['frequency'],band_dict[i],k))
#      print indexes[iband][1]
#   band_dict_old=copy.deepcopy(band_dict)
