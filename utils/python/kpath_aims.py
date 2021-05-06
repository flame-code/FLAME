#!/usr/bin/env python
import argparse 
import pymatgen as mg
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.symmetry.bandstructure import HighSymmKpath
from pprint import pprint
#------------------------------------------------------------
str1 = """This script can be used for bandstructure calculations by FHI-aims.
It gets the k-vectors of highest-symmetry paths in the reciprocal lattice of a structure."""
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('POSCAR', action='store', type=str, 
        help="POSCAR is the name of the input file in VASP5 format.")
parser.add_argument('tol', type=float, 
        help="""the tolerance for analysing space group. 
        A float number less than 0.05 is recommended.""")
parser.add_argument('npoints', type=int, help="The number of sampling points (usually 50 is fine).")
parser.add_argument('type', type=int, help="""This number determines the type of output:
0 for phonon and 1 for electronic bandstructure calculations""")
args=parser.parse_args()
#------------------------------------------------------------
tol = args.tol
np = args.npoints
tp = args.type
keyword = ['phonon band','output band']
structure = mg.Structure.from_file("POSCAR")
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
iold="LLL"

for i_path in range(0, len(pp)):
  for pp_l in pp[i_path:i_path+1]:
    count=count+1
    count2=0
    size=len(pp_l)
    if i_path>0:
        i = pp_l[0]
        if j != i: 
            print (str(keyword[tp]), "  ".join(map(str, kk[j])) , "  ", "  ".join(map(str, kk[i])) , "  ", str(np),j , i, file=f)
    for i in pp_l[:size-1]:
      if iold != i :
          iold=i
      for j in pp_l[count2+1:count2+2]:
          print (str(keyword[tp]),"  ".join(map(str,kk[i])), "  ","  ".join(map(str, kk[j])) , "  ", str(np),i,j, file=f)
      count2=count2+1
    iold=j
f.close()
