#!/usr/bin/env python
import argparse
import random
from cStringIO import StringIO
from io_yaml import *
import math

str1 = "This script shakes a structure."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file")
parser.add_argument('ampl', type=float, help='Amplitude for shaking the structure')
parser.add_argument('mvx', type=float, help='default value: 1', nargs='?', const=1, default=1)
parser.add_argument('mvy', type=float, help='default value: 1', nargs='?', const=1, default=1)
parser.add_argument('mvz', type=float, help='default value: 1', nargs='?', const=1, default=1)
args=parser.parse_args()

mvx=args.mvx
mvy=args.mvy
mvz=args.mvz
ampl=args.ampl

atoms_all=read_yaml(args.fn_inp)
#print atoms_all[0].nat
random.seed()

for atoms in atoms_all:
    for iat in range(atoms.nat):
        #ttx=random.random()
        #tty=random.random()
        #ttz=random.random()
        #atoms.positions[iat,0]+=(ttx-0.5)*ampl*2.0
        #atoms.positions[iat,1]+=(tty-0.5)*ampl*2.0
        #atoms.positions[iat,2]+=(ttz-0.5)*ampl*2.0
        theta=math.pi*random.random()
        phi=2.0*math.pi*random.random()
        ttx=mvx*ampl*math.sin(theta)*math.cos(phi)
        tty=mvy*ampl*math.sin(theta)*math.sin(phi)
        ttz=mvz*ampl*math.cos(theta)
        atoms.rat[iat][0]+=ttx
        atoms.rat[iat][1]+=tty
        atoms.rat[iat][2]+=ttz

write_yaml(atoms_all,args.fn_out)
