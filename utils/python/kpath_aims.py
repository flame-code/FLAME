#!/usr/bin/env python
import sys
file = sys.argv[1]
import pymatgen as mg
#from pymatgen.io.smartio import read_structure, write_structure, read_mol, write_mol
#structure = read_structure(file)
structure = mg.Structure.from_file("POSCAR")   #str(open("rlx_str040.cif").read(),fmt="cif")  #str(open("POSCAR").read(),fmt="poscar")
#from pymatgen.symmetry.finder import SymmetryFinder
#symmetries = SymmetryFinder(structure, 0.001, 5)
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
finder = SpacegroupAnalyzer(structure,symprec=0.5,angle_tolerance=5)
spg_s = finder.get_spacegroup_symbol()
spg_n = finder.get_spacegroup_number()
pg = finder.get_point_group()
pm =finder.find_primitive()
#spg=symmetries.get_spacegroup_number()
#print "Spacegroup : ", spg_s
print "Int number : ", spg_n
#print "pointgroup : ", pg
from pymatgen.symmetry.bandstructure import HighSymmKpath
pather = HighSymmKpath(structure,symprec=0.5, angle_tolerance=5)
kpoints_path=pather.kpath
from pprint import pprint
kk=kpoints_path["kpoints"]
pp=kpoints_path["path"]
count=1

#pp[1]=['K','X']
#print pp
#print kk
##print pp[1:2]
##print len(pp)


f = open('KPOINTS_PATH', 'w')
#g = open('KPATH_NAME', 'w')
#h = open('KPATH_COORD', 'w')
p = open('primitive','w')
print >> p,  pm
#print >> f,  "k-points along high symmetry lines"
#print >> f, "20  ! 20 intersections" 
#print >> f,  "Line-mode"
#print >> f, "rec"
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
                print >> f, "phonon band", "  ".join(map(str, kk[j])) , "  ", "  ".join(map(str, kk[i])) , "  ", "100",j , i
                #print >> f, "  ".join(map(str, kk[i])) , "                !", i
                #print >> f, " "
    for i in pp_l[:size-1]:
      if iold != i :
#         print >> g, i
#         print >> h, " ".join(map(str, kk[i]))
          iold=i
      for j in pp_l[count2+1:count2+2]:
          print >> f, "phonon band","  ".join(map(str,kk[i])), "  ","  ".join(map(str, kk[j])) , "  ", "100",i,j
         #print >> f, "  ".join(map(str, kk[j])) , "                !", j
      count2=count2+1
#    print >> g, j
#    print >> h, " ".join(map(str, kk[j]))
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
