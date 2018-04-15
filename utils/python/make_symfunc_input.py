#!/usr/bin/env python
import sys
import numpy as np
import math

rotaxis=[]
if len(sys.argv) < 3:
    print ""
    print "usage: make_symfunc_input.py atom method "
    print ""
    print "example: make_symfunc_input.py Na new "
    print ""
    exit()
else:
    atom_name = sys.argv[1]
    method = sys.argv[2]
  #  rotaxis.append(float(sys.argv[2]))
list_new = []    
list_new.append(["g02_001:" ,  0.0010 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g02_002:" ,  0.0100 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g02_003:" ,  0.0200 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g02_004:" ,  0.0350 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g02_005:" ,  0.0600 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g02_006:" ,  0.1000 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g02_007:" ,  0.2000 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g02_008:" ,  0.4000 ,  0.0000 ,  0.0000 ,  0.0000 ])
list_new.append(["g05_001:" ,  0.0001 ,  1.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_002:" ,  0.0001 ,  1.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_003:" ,  0.0001 ,  2.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_004:" ,  0.0001 ,  2.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_005:" ,  0.0030 ,  1.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_006:" ,  0.0030 ,  1.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_007:" ,  0.0030 ,  2.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_008:" ,  0.0030 ,  2.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_009:" ,  0.0080 ,  1.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_010:" ,  0.0080 ,  1.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_011:" ,  0.0080 ,  2.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_012:" ,  0.0080 ,  2.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_013:" ,  0.0150 ,  1.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_014:" ,  0.0150 ,  1.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_015:" ,  0.0150 ,  2.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_016:" ,  0.0150 ,  2.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_017:" ,  0.0150 ,  4.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_018:" ,  0.0150 ,  4.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_019:" ,  0.0150 , 16.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_020:" ,  0.0150 , 16.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_021:" ,  0.0250 ,  1.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_022:" ,  0.0250 ,  1.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_023:" ,  0.0250 ,  2.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_024:" ,  0.0250 ,  2.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_025:" ,  0.0250 ,  4.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_026:" ,  0.0250 ,  4.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_027:" ,  0.0250 , 16.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_028:" ,  0.0250 , 16.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_029:" ,  0.0450 ,  1.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_030:" ,  0.0450 ,  1.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_031:" ,  0.0450 ,  2.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_032:" ,  0.0450 ,  2.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_033:" ,  0.0450 ,  4.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_034:" ,  0.0450 ,  4.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_035:" ,  0.0450 , 16.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_036:" ,  0.0450 , 16.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_037:" ,  0.0800 ,  1.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_038:" ,  0.0800 ,  1.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_039:" ,  0.0800 ,  2.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_040:" ,  0.0800 ,  2.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_041:" ,  0.0800 ,  4.0000 , -1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_042:" ,  0.0800 ,  4.0000 ,  1.0000 ,  0.0000 , 0.0000])
list_new.append(["g05_043:" ,  0.0800 , 16.0000 ,  1.0000 ,  0.0000 , 0.0000])

# -------------------------------------------------------------------------------------

factorH=0.31
factor_Li=1.28
if(atom_name=='H') :
    factor=0.31
elif(atom_name=='He') :
    factor=0.28
elif(atom_name=='Li') :
    factor=1.28
elif(atom_name=='Be') :
    factor=0.96
elif(atom_name=='B') :
    factor=0.84
elif(atom_name=='C') :
    factor=0.76
elif(atom_name=='N') :
    factor=0.71
elif(atom_name=='O') :
    factor=0.66
elif(atom_name=='F') :
    factor=0.57
elif(atom_name=='Ne') :
    factor=0.58
elif(atom_name=='Na') :
    factor=1.66
elif(atom_name=='Mg') :
    factor=1.41
elif(atom_name=='Al') :
    factor=1.21
elif(atom_name=='Si') :
    factor=1.11
elif(atom_name=='P') :
    factor=1.07
elif(atom_name=='S') :
    factor=1.05
elif(atom_name=='Cl') :
    factor=1.02
elif(atom_name=='Ar') :
    factor=1.06
elif(atom_name=='K') :
    factor=2.03
elif(atom_name=='Ca') :
    factor=1.76
elif(atom_name=='Sc') :
    factor=1.70
elif(atom_name=='Ti') :
    factor=1.60
elif(atom_name=='V') :
    factor=1.53
elif(atom_name=='Cr') :
    factor=1.39
elif(atom_name=='Mn') :
    factor=1.39
elif(atom_name=='Fe') :
    factor=1.32
elif(atom_name=='Co') :
    factor=1.26
elif(atom_name=='Ni') :
    factor=1.24
elif(atom_name=='Cu') :
    factor=1.32
elif(atom_name=='Zn') :
    factor=1.22
elif(atom_name=='Ga') :
    factor=1.22
elif(atom_name=='Ge') :
    factor=1.20
elif(atom_name=='As') :
    factor=1.19
elif(atom_name=='Se') :
    factor=1.20
elif(atom_name=='Br') :
    factor=1.20
elif(atom_name=='Kr') :
    factor=1.16
elif(atom_name=='Rb') :
    factor=2.20
elif(atom_name=='Sr') :
    factor=1.95
elif(atom_name=='Y') :
    factor=1.90
elif(atom_name=='Zr') :
    factor=1.75
elif(atom_name=='Nb') :
    factor=1.64
elif(atom_name=='Mo') :
    factor=1.54
elif(atom_name=='Tc') :
    factor=1.47
elif(atom_name=='Ru') :
    factor=1.46
elif(atom_name=='Rh') :
    factor=1.42
elif(atom_name=='Pd') :
    factor=1.39
elif(atom_name=='Ag') :
    factor=1.45
elif(atom_name=='Cd') :
    factor=1.44
elif(atom_name=='In') :
    factor=1.42
elif(atom_name=='Sn') :
    factor=1.39
elif(atom_name=='Sb') :
    factor=1.39
elif(atom_name=='Te') :
    factor=1.38
elif(atom_name=='I') :
    factor=1.39
elif(atom_name=='Xe') :
    factor=1.40
elif(atom_name=='Cs') :
    factor=2.44
elif(atom_name=='Ba') :
    factor=2.15
elif(atom_name=='La') :
    factor=2.07
elif(atom_name=='Ce') :
    factor=2.04
elif(atom_name=='Pr') :
    factor=2.03
elif(atom_name=='Nd') :
    factor=2.01
elif(atom_name=='Pm') :
    factor=1.99
elif(atom_name=='Sm') :
    factor=1.98
elif(atom_name=='Eu') :
    factor=1.98
elif(atom_name=='Gd') :
    factor=1.96
elif(atom_name=='Tb') :
    factor=1.94
elif(atom_name=='Dy') :
    factor=1.92
elif(atom_name=='Ho') :
    factor=1.92
elif(atom_name=='Er') :
    factor=1.89
elif(atom_name=='Tm') :
    factor=1.90
elif(atom_name=='Yb') :
    factor=1.87
elif(atom_name=='Lu') :
    factor=1.87
elif(atom_name=='Hf') :
    factor=1.75
elif(atom_name=='Ta') :
    factor=1.70
elif(atom_name=='W') :
    factor=1.62
elif(atom_name=='Re') :
    factor=1.51
elif(atom_name=='Os') :
    factor=1.44
elif(atom_name=='Ir') :
    factor=1.41
elif(atom_name=='Pt') :
    factor=1.36
elif(atom_name=='Au') :
    factor=1.36
elif(atom_name=='Hg') :
    factor=1.32
elif(atom_name=='Tl') :
    factor=1.45
elif(atom_name=='Pb') :
    factor=1.46
elif(atom_name=='Bi') :
    factor=1.48
elif(atom_name=='Po') :
    factor=1.40
elif(atom_name=='At') :
    factor=1.50
elif(atom_name=='Rn') :
    factor=1.50
elif(atom_name=='Fr') :
    factor=2.60
elif(atom_name=='Ra') :
    factor=2.21
elif(atom_name=='Ac') :
    factor=2.15
elif(atom_name=='Th') :
    factor=2.06
elif(atom_name=='Pa') :
    factor=2.00
elif(atom_name=='U') :
    factor=1.96
elif(atom_name=='Np') :
    factor=1.90
elif(atom_name=='Pu') :
    factor=1.87
elif(atom_name=='Am') :
    factor=1.80
elif(atom_name=='Cm') :
    factor=1.69
else :
    print 'ERROR: no covalent radius stored for atomtype=',atom_name


#***********  ionic *******************         
factor_ion_Li= 90.0/100./0.529
if(atom_name=='Li') :
    factor_ion= 90/100./0.529
    hardness =  2.39/27.21
elif(atom_name=='Na') :
    factor_ion=116/100./0.529
    hardness =  2.30/27.21
elif(atom_name=='K') :
    factor_ion=152/100./0.529
    hardness =  1.92/27.21
elif(atom_name=='F') :
    factor_ion=119/100./0.529
    hardness =  7.01/27.21
elif(atom_name=='S') :
    factor_ion=170/100./0.529
    hardness =  4.14/27.21
elif(atom_name=='Cl') :
    factor_ion=167/100./0.529
    hardness =  4.68/27.21
elif(atom_name=='Br') :
    factor_ion=182/100./0.529
    hardness =  4.22/27.21
elif(atom_name=='Pb') :
    factor_ion=133/100./0.529
    hardness =  5.50/27.21
elif(atom_name=='Te') :
    factor_ion=207/100./0.529
    hardness =  3.52/27.21
elif(atom_name=='Zr') :
    factor_ion=86/100./0.529
    hardness =  3.21/27.21
elif(atom_name=='O') :
    factor_ion=126/100./0.529
    hardness =  6.08/27.21
elif(atom_name=='Al') :
    factor_ion=68/100./0.529
    hardness =  2.77/27.21
elif(atom_name=='Ca') :
    factor_ion=114/100./0.529
    hardness =  4.0/27.21

for i in range((len(list_new))):
    #list_new[i][1]*= (factorH/factor)
    list_new[i][1]*= (factor_ion_Li/factor_ion)**2*2 
fact=factor_ion*0.85
print "ionic radii" , fact

#-----------------------------------------------------------------------------------

filename=atom_name
filename+=".ann.input.yaml"
print filename
f= open(filename,"w")
f.write ( '''main:
    nodes: [5, 5]
    rcut: 15.0
    ampl_chi: 1.00
    prefactor_chi: 2.00
    zion: 2.0
    gausswidth_ion: 0.600
    ener_ref: 0.0
''')
f.write ("    gausswidth: %7.4f\n"% fact)
f.write ("    hardness: %7.4f\n"% hardness)
f.write ('''    chi0: CHI0
    spring_const: 1.00
    qinit: 1.0
''')
f.write ("    method:  %s\n"%(method ))
f.write ( "symfunc:\n")
for i in range((len(list_new))):
    if (i<8) :
        f.write( "    %s %11.6f %9.4f %9.4f %9.4f\n" % (list_new[i][0],list_new[i][1], list_new[i][2],list_new[i][3],list_new[i][4]))
    else:
        f.write( "    %s %11.6f %9.4f %9.4f %9.4f %9.4f\n" % (list_new[i][0],list_new[i][1],list_new[i][2],list_new[i][3],list_new[i][4],list_new[i][5]))
f.close()
