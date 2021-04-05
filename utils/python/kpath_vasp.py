#!/usr/bin/env python
import argparse 
import pymatgen as mg
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pprint import pprint
#------------------------------------------------------------
str1 = """ This script can be used for bandstructure calculations by VASP or PHONOPY. 
It gets the k-vectors of highest-symmetry paths in the reciprocal lattice of a structure."""
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('POSCAR', action='store', type=str, 
        help="POSCAR is the name of the input file in VASP5 format.")
parser.add_argument('tol', type=float, 
        help="""the tolerance for analysing space group. 
        A float number less than 0.05 is recommended.""")
args=parser.parse_args()
#------------------------------------------------------------
tol = args.tol
structure = Structure.from_file("POSCAR")
spg_analy = SpacegroupAnalyzer(structure)
primitive_standard_structure=spg_analy.get_primitive_standard_structure(international_monoclinic=False)
primitive_standard_structure.to(fmt="poscar", filename="PPOSCAR")
structure = Structure.from_file("PPOSCAR")
finder = SpacegroupAnalyzer(structure,symprec=tol,angle_tolerance=5)
spg_s = finder.get_space_group_symbol()
spg_n = finder.get_space_group_number()
#pg = finder.get_point_group_symbol()

print ("Spacegroup (symbol) : ", spg_s)
print ("Spacegroup (number) : ", spg_n)
#print ("pointgroup : ", pg)

pather = HighSymmKpath(structure,symprec=tol, angle_tolerance=5)
kpoints_path=pather.kpath
kk=kpoints_path["kpoints"]
pp=kpoints_path["path"]
count=1

f = open('KPOINTS_PATH', 'w')
g = open('KPATH_NAME', 'w')
h = open('KPATH_COORD', 'w')
print ("k-points along high symmetry lines", file=f)
print ("20  ! 20 intersections" , file=f)
print ("Line-mode", file=f)
print ("rec", file=f)
iold="LLL"

for i_path in range(0, len(pp)):
  for pp_l in pp[i_path:i_path+1]:
    count=count+1
    count2=0
    size=len(pp_l)
    if i_path>0:
        i = pp_l[0]
        if j != i: 
            print ("  ".join(map(str, kk[j])) , "                !", j, file=f)
            print ("  ".join(map(str, kk[i])) , "                !", i, file=f) 
            print (" ", file=f)
    for i in pp_l[:size-1]:
      if iold != i :
          print (i, file=g)
          print (" ".join(map(str, kk[i])), file=h)
          iold=i
      print ("  ".join(map(str,kk[i])), "                !", i, file=f)
      for j in pp_l[count2+1:count2+2]:
          print ("  ".join(map(str, kk[j])), "                !", j, file=f)
      count2=count2+1
      print (" ", file=f)
    print (j, file=g)
    print (" ".join(map(str, kk[j])), file=h)
    iold=j
f.close()
g.close()
h.close()
