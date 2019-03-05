#!/usr/bin/env python
import argparse
from atoms import *
from io_yaml import *
import yaml
import numpy as np
import copy
#RZX
#*****************************************************************************************
def bdft_read(filename):
    atoms_all.append(Atoms())
    stream=open(filename,'r')
    docs=yaml.load(stream)
    atoms_all[0].nat = len(docs["posinp"]["positions"])
    atoms_all[0].rat=[[-1 for i in range(3)] for j in range(atoms_all[0].nat)] 
    atoms_all[0].fat=[[-1 for i in range(3)] for j in range(atoms_all[0].nat)] 
    atoms_all[0].sat=['##' for i in range(atoms_all[0].nat)]
    for i in range(atoms_all[0].nat):
        sat = docs["posinp"]["positions"][i].keys()
        atoms_all[0].sat[i]=sat[0]
        for j in range(3):
            atoms_all[0].rat[i][j]=docs["posinp"]["positions"][i][atoms_all[0].sat[i]][j]
            atoms_all[0].fat[i][j]=docs["Atomic Forces (Ha/Bohr)"][i][atoms_all[0].sat[i]][j]
    if (abs(docs["dft"]["elecfield"][1]) > 0) or (abs(docs["dft"]["elecfield"][2]) > 0) or (abs(docs["dft"]["elecfield"][0]) > 0) :
        atoms_all[0].elecfield=[-1 for i in range(3)]
        atoms_all[0].elecfield_present=True
        for i in range(3):
            atoms_all[0].elecfield[i]=docs["dft"]["elecfield"][i]
    if abs(docs["Electric Dipole Moment (AU)"]["norm(P)"]) > 0:
        atoms_all[0].dpm=[-1 for i in range(3)]
        atoms_all[0].dpm_present=True
        for i in range(3):
            atoms_all[0].dpm[i]=docs["Electric Dipole Moment (AU)"]["P vector"][i]
#    if abs(docs["Quadrupole Moment (AU)"]["trace"]) > 0:
#        atoms_all[0].qpm=[[-1 for i in range(3)] for j in range(3)]
#        atoms_all[0].qpm_present=True
#    for i in range(3):
#        for j in range(3):
#            atoms_all[0].qpm[i][j]=docs["Quadrupole Moment (AU)"]["Q matrix"][i][j]
    
    if (docs["posinp"]["units"]=='angstroem'):
        atoms_all[0].units_length_io='angstrom'
    elif (docs["posinp"]["units"]=='Atomic'):
        atoms_all[0].units_length_io='bohr'
    atoms_all[0].epot = docs["Input Hamiltonian"]["Energies"]["Epot"]
    atoms_all[0].qtot = -1.0*(docs["Input Hamiltonian"]["Total electronic charge"]+docs["Total ionic charge"])
    if abs(atoms_all[0].qtot)< 1.E-6:
        atoms_all[0].qtot=0.E+00
    if (docs["Atomic System Properties"]["Boundary Conditions"]=='Periodic'):
        atoms_all[0].boundcond = 'bulk'
    elif (docs["Atomic System Properties"]["Boundary Conditions"]=='Surface'):
        atoms_all[0].boundcond = 'slab'
    elif (docs["Atomic System Properties"]["Boundary Conditions"]=='Free'):
        atoms_all[0].boundcond = 'free'
    atoms_all[0].cellvec[0][0]= docs["Sizes of the simulation domain"]["Angstroem"][0]
    atoms_all[0].cellvec[1][1]= docs["Sizes of the simulation domain"]["Angstroem"][1]
    atoms_all[0].cellvec[2][2]= docs["Sizes of the simulation domain"]["Angstroem"][2]
    atoms_all[0].bemoved=['TTT' for i in range(atoms_all[0].nat)]
    return atoms_all
#********************************************************************************************

str1 = "This script reads BigDFT log file and writes it in yaml format which is acceptable as posinp.yaml for flam."
parser = argparse.ArgumentParser(description=str1)
parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the BigDFT log in yaml format")
parser.add_argument('fn_out', action='store' ,type=str, help="Name of the output file in yaml format")
args=parser.parse_args()

atoms_all=[]
atoms_all=bdft_read(args.fn_inp)
write_yaml(atoms_all,args.fn_out)
