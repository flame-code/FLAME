#!/usr/bin/env python
import argparse 
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pprint import pprint
#************************************************************
str1 = "This script should be used for calculations by VASP or PHONOPY and gets all paths in the reciprocal lattice of a structure."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('POSCAR', action='store' ,type=str, help="POSCAR is the name of the input file in VASP5 format")
parser.add_argument('tol', type=float, help='the tolerance for analysing space group. A number less than 0.05 is recommended.')
#parser.add_argument('npoints', type=int, help='the number of sampling points. Usually 50 is proper')
args=parser.parse_args()
#************************************************************
tol = args.tol
#np = args.npoints
#structure = read_structure(file)
structure = mg.Structure.from_file("POSCAR")   #str(open("rlx_str040.cif").read(),fmt="cif")  #str(open("POSCAR").read(),fmt="poscar")
#from pymatgen.symmetry.finder import SymmetryFinder
#symmetries = SymmetryFinder(structure, 0.001, 5)
finder = SpacegroupAnalyzer(structure,symprec=tol,angle_tolerance=5)
spg_s = finder.get_spacegroup_symbol()
spg_n = finder.get_spacegroup_number()
pg = finder.get_point_group()
pm =finder.find_primitive()

#print "Spacegroup : ", spg_s
print "SPG (Int number) : ", spg_n
#print "pointgroup : ", pg
pather = HighSymmKpath(structure,symprec=tol, angle_tolerance=5)
kpoints_path=pather.kpath
kk=kpoints_path["kpoints"]
pp=kpoints_path["path"]
count=1

#pp[1]=['K','X']
#print pp
#print kk
##print pp[1:2]
##print len(pp)


f = open('KPOINTS_PATH', 'w')
g = open('KPATH_NAME', 'w')
h = open('KPATH_COORD', 'w')
print >> f,  "k-points along high symmetry lines"
print >> f, "20  ! 20 intersections" 
print >> f,  "Line-mode"
print >> f, "rec"
iold="LLL"

for i_path in range(0, len(pp)):
  for pp_l in pp[i_path:i_path+1]:
  #  print "Line set ",count," ------------------------------- "
    count=count+1
    count2=0
    size=len(pp_l)
    if i_path>0:
	    i = pp_l[0]
	    if j != i: 
                print >> f, "  ".join(map(str, kk[j])) , "                !", j
                print >> f, "  ".join(map(str, kk[i])) , "                !", i
                print >> f, " "
    for i in pp_l[:size-1]:
      if iold != i :
          print >> g, i
          print >> h, " ".join(map(str, kk[i]))
	  iold=i
      print >> f,  "  ".join(map(str,kk[i])), "                !", i
      for j in pp_l[count2+1:count2+2]:
          print >> f, "  ".join(map(str, kk[j])) , "                !", j
      count2=count2+1
      print >> f, " "
    print >> g, j
    print >> h, " ".join(map(str, kk[j]))
    iold=j
f.close()
g.close()
h.close()
#*****************************************************************************************
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
