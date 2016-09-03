#!/usr/bin/env python
from ase.io.vasp import read_vasp,write_vasp
import random
import sys
from cStringIO import StringIO
import math

if len(sys.argv) < 3:
    print "usage: shake_poscar.py input_filename amplitude"
    exit()
else:
    filename = sys.argv[1]
    ampl=float(sys.argv[2])

atoms=read_vasp(filename=filename)

nat=len(atoms.numbers)

#print atoms.positions[nat-1]
#print atoms.numbers
#print atoms.cell
#print nat

#ampl=0.1
random.seed(1)
for iat in range(nat):
    #ttx=random.random()
    #tty=random.random()
    #ttz=random.random()
    #atoms.positions[iat,0]+=(ttx-0.5)*ampl*2.0
    #atoms.positions[iat,1]+=(tty-0.5)*ampl*2.0
    #atoms.positions[iat,2]+=(ttz-0.5)*ampl*2.0
    theta=math.pi*random.random()
    phi=2.0*math.pi*random.random()
    ttx=ampl*math.sin(theta)*math.cos(phi)
    tty=ampl*math.sin(theta)*math.sin(phi)
    ttz=ampl*math.cos(theta)
    atoms.positions[iat,0]+=ttx
    atoms.positions[iat,1]+=tty
    atoms.positions[iat,2]+=ttz

#write_vasp('POSCAR2',atoms,symbol_count=None,vasp5=True)
output=StringIO()
write_vasp(output,atoms,symbol_count=None,vasp5=True)
print output.getvalue().rstrip()
output.close()
