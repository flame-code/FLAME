#!/usr/bin/env python
import sys
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pprint import pprint
#************************************************************
if len(sys.argv) < 5:
    print "usage: kpath_aims.py POSCAR tolerance in_number npoints"
    print "  for phonon band:     in_number = 0"
    print "  for electronic band: in_number = 1"
    print "suggested values for npoints : "
    print "  phonon band :       npoints = 100"
    print "  electronic band :   npoints = 20"
    exit()
else:
    file = sys.argv[1]
    tol = float(sys.argv[2])
    bb = int(sys.argv[3])
    np = int(sys.argv[4])

keyword = ['phonon band','output band']
#structure = read_structure(file)
structure = mg.Structure.from_file("POSCAR")   #str(open("rlx_str040.cif").read(),fmt="cif")  #str(open("POSCAR").read(),fmt="poscar")
#from pymatgen.symmetry.finder import SymmetryFinder
#symmetries = SymmetryFinder(structure, 0.001, 5)
finder = SpacegroupAnalyzer(structure,symprec = tol,angle_tolerance=5)
spg_s = finder.get_spacegroup_symbol()
spg_n = finder.get_spacegroup_number()
pg = finder.get_point_group()
pm = finder.find_primitive()

#print "Spacegroup : ", spg_s
print "space_group : ", spg_n
#print "pointgroup : ", pg
pather = HighSymmKpath(structure,symprec = tol,angle_tolerance=5)
kpoints_path = pather.kpath
kk = kpoints_path["kpoints"]
pp = kpoints_path["path"]
count = 1

#pp[1]=['K','X']
#print pp
#print kk
##print pp[1:2]
##print len(pp)


f = open('KPOINTS_PATH', 'w')
#g = open('KPATH_NAME', 'w')
#h = open('KPATH_COORD', 'w')
#p = open('primitive','w')
iold="LLL"

#print len(pp)  
for i_path in range(0, len(pp)):
  for pp_l in pp[i_path:i_path+1]:
  #  print "Line set ",count," ------------------------------- "
    count=count+1
    count2=0
    size=len(pp_l)
    if i_path>0:
	    i = pp_l[0]
	    if j != i: 
                print >> f, str(keyword[bb]), "  ".join(map(str, kk[j])) , "  ", "  ".join(map(str, kk[i])) , "  ", str(np),j , i
    for i in pp_l[:size-1]:
      if iold != i :
          iold=i
      for j in pp_l[count2+1:count2+2]:
          print >> f, str(keyword[bb]),"  ".join(map(str,kk[i])), "  ","  ".join(map(str, kk[j])) , "  ", str(np),i,j
      count2=count2+1
    iold=j
f.close()

#*****************************************************************************************************************************
#from pymatgen.io.vaspio.vasp_input import Kpoints
#
#struct = structure
#
## First brillouin zone
#ibz = HighSymmKpath(struct)
#print("ibz type     : {0}".format(ibz.name))
#
## suggested path
#print("paths in first brillouin zone :")
#for path in ibz.kpath["path"]:
#	print(path)
#
#kpoints = list()
#labels = list()
#for path in ibz.kpath["path"]:
#	for kpts in path:
#		kpoints.append(ibz.kpath["kpoints"][kpts])
#		labels.append(kpts)
#
## print kpoints file
#Kpoints(comment = "band diagram for monoclinic cell, unique axes a",
#       num_kpts = 100,
#       style = Kpoints.supported_modes.Line_mode,
#       coord_type = "Reciprocal",
#       kpts = kpoints,
#       labels = labels,
#        ).write_file("KPOINTS")
