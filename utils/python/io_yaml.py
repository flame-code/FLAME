from atoms import *
import copy
import yaml

############################################
try:
    from yaml import CLoader as Loader
except ImportError:
    from yaml import Loader
############################################
#*****************************************************************************************
def read_yaml(filename):
    stream=open(filename, "r")
    dict_list = list(yaml.load_all(stream, Loader=Loader))
    atoms_all=[]
    for dl in dict_list:
        atoms=dict2atoms(dl)
        atoms_all.append(Atoms())
        atoms_all[-1]=copy.deepcopy(atoms)
        del atoms
    return atoms_all
#*****************************************************************************************
def write_yaml(atoms_all,filename):
    output_stream=open(filename,"w")
    for atoms in atoms_all:
        dict_atoms=atoms2dict(atoms)
        output_stream.write("---\n")
        yaml.dump(dict_atoms,output_stream,default_flow_style=None)
        del dict_atoms
    output_stream.close()
#*****************************************************************************************
def atoms2dict(atoms):
    dict_atoms = dict()
    dict_atoms={'conf':{}}
    if len(atoms.fat)>0:
        dict_atoms['conf']['force']=[]
    if len(atoms.vat)>0:
        dict_atoms['conf']['velocity']=[]
    if atoms.qtot:
        dict_atoms['conf']['qtot']=atoms.qtot
    if atoms.epot:
        dict_atoms['conf']['epot']=atoms.epot
    dict_atoms['conf']['nat']=atoms.nat
    #local variables
    dict_atoms['conf']['coord']=[]
    for iat in range(atoms.nat):
        dict_atoms['conf']['coord'].append([])
        dict_atoms['conf']['coord'][-1]=([atoms.rat[iat][0],atoms.rat[iat][1], \
                atoms.rat[iat][2],atoms.sat[iat],atoms.bemoved[iat]])
        if len(atoms.fat)>0:
            dict_atoms['conf']['force'].append(atoms.fat[iat])
        if len(atoms.vat)>0:
            dict_atoms['conf']['velocity'].append(atoms.vat[iat])
    dict_atoms['conf']['units_length']=atoms.units_length_io
    dict_atoms['conf']['bc']=atoms.boundcond
    dict_atoms['conf']['cell']=[]
    dict_atoms['conf']['cell'].append([])
    dict_atoms['conf']['cell'][-1].append([])
    dict_atoms['conf']['cell'][0]=atoms.cellvec[0]
    dict_atoms['conf']['cell'].append([])
    dict_atoms['conf']['cell'][-1].append([])
    dict_atoms['conf']['cell'][1]=atoms.cellvec[1]
    dict_atoms['conf']['cell'].append([])
    dict_atoms['conf']['cell'][-1].append([])
    dict_atoms['conf']['cell'][2]=atoms.cellvec[2]
    if atoms.elecfield_present:
        dict_atoms['conf']['elecfield']=[-1 for i in range(3)]
        dict_atoms['conf']['elecfield'][0]=atoms.elecfield[0]
        dict_atoms['conf']['elecfield'][1]=atoms.elecfield[1]
        dict_atoms['conf']['elecfield'][2]=atoms.elecfield[2]
    if atoms.dpm_present:
        dict_atoms['conf']['dpm']=[-1 for i in range(3)]
        dict_atoms['conf']['dpm'][0]=atoms.dpm[0]
        dict_atoms['conf']['dpm'][1]=atoms.dpm[1]
        dict_atoms['conf']['dpm'][2]=atoms.dpm[2]
    if atoms.qpm_present:
        dict_atoms['conf']['qpm']=[]
        dict_atoms['conf']['qpm'].append([])
        dict_atoms['conf']['qpm'][-1].append([])
        dict_atoms['conf']['qpm'][0]=atoms.qpm[0]
        dict_atoms['conf']['qpm'].append([])
        dict_atoms['conf']['qpm'][-1].append([])
        dict_atoms['conf']['qpm'][1]=atoms.qpm[1]
        dict_atoms['conf']['qpm'].append([])
        dict_atoms['conf']['qpm'][-1].append([])
        dict_atoms['conf']['qpm'][2]=atoms.qpm[2]
        
    return dict_atoms
#*****************************************************************************************
def dict2atoms(dict_atoms):
    atoms=Atoms()
    atoms.nat             = dict_atoms['conf']['nat']
    atoms.boundcond       = dict_atoms['conf']['bc']

    if dict_atoms['conf'].has_key('units_length'):
        atoms.units_length_io = dict_atoms['conf']['units_length']
    if dict_atoms['conf'].has_key('epot'):
        atoms.epot_present=True
        atoms.epot = dict_atoms['conf']['epot']
    if dict_atoms['conf'].has_key('qtot'):
        atoms.qtot_present=True
        atoms.qtot = dict_atoms['conf']['qtot']
    if dict_atoms['conf'].has_key('force'):
        atoms.fat_present=True
        for i in range(int(atoms.nat)):
            fr = dict_atoms['conf']['force'][i]
            atoms.fat.append(fr)
    if dict_atoms['conf'].has_key('cell'):
        atoms.cell_present=True
        atoms.cellvec[0] = dict_atoms['conf']['cell'][0]
        atoms.cellvec[1] = dict_atoms['conf']['cell'][1]
        atoms.cellvec[2] = dict_atoms['conf']['cell'][2]
    if dict_atoms['conf'].has_key('velocity'):
        atoms.vat_present=True
        for i in range(int(atoms.nat)):
            vx = dict_atoms['conf']['velocity'][i][0]
            vy = dict_atoms['conf']['velocity'][i][1]
            vz = dict_atoms['conf']['velocity'][i][2]
            atoms.vat.append([vx,vy,vz])
    for i in range(int(atoms.nat)):
        x = dict_atoms['conf']['coord'][i][0]
        y = dict_atoms['conf']['coord'][i][1]
        z = dict_atoms['conf']['coord'][i][2]
        atoms.rat.append([x,y,z])
        sat = dict_atoms['conf']['coord'][i][3] 
        atoms.sat.append(sat)
        if len(dict_atoms['conf']['coord'][i])==5:
            atoms.bemoved_present=True
            abm = dict_atoms['conf']['coord'][i][4]
            atoms.bemoved.append(abm)
        #else:
        #    atoms.bemoved.append(' ')

    return atoms
#*****************************************************************************************
