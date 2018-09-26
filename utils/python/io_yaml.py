from atoms import *
import copy
import yaml
#*****************************************************************************************
def read_yaml(filename):
    stream=open(filename, "r")
    dict_list = list(yaml.load_all(stream))
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
        yaml.dump(dict_atoms,output_stream)
        del dict_atoms
    output_stream.close()
#*****************************************************************************************
def atoms2dict(atoms):
    dict_atoms = dict()
    dict_atoms={'conf':{}}
    if len(atoms.fat)>0:
        dict_atoms['conf']['force']=[]
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
    return dict_atoms
#*****************************************************************************************
def dict2atoms(dict_atoms):
    for key, value in dict_atoms.items(): 
        atoms=Atoms()
        atoms.nat       = dict_atoms[key]['nat']
        atoms.units_length_io = dict_atoms[key]['units_length']
        atoms.boundcond = dict_atoms[key]['bc']
        for sub_key in value:
            if sub_key == 'epot':
                atoms.epot = dict_atoms[key]['epot']
            if sub_key == 'qtot':
                atoms.qtot = dict_atoms[key]['qtot']
            if sub_key == 'force':
                for i in range(int(atoms.nat)):
                    fr = dict_atoms[key]['force'][i]
                    atoms.fat.append(fr)
        atoms.cellvec[0] = dict_atoms[key]['cell'][0]
        atoms.cellvec[1] = dict_atoms[key]['cell'][1]
        atoms.cellvec[2] = dict_atoms[key]['cell'][2]
        for i in range(int(atoms.nat)):
            r   = dict_atoms[key]['coord'][i]
            sat = dict_atoms[key]['coord'][i][3] 
            ab  = dict_atoms[key]['coord'][i][4]
            atoms.rat.append(r)
            atoms.bemoved_present=True
            atoms.sat.append(sat)
            atoms.bemoved.append(ab)

    if dict_atoms['conf'].has_key('cell'): atoms.cell_present=True
    if dict_atoms['conf'].has_key('epot'): atoms.epot_present=True
    if dict_atoms['conf'].has_key('force'): atoms.fat_present=True
    if dict_atoms['conf'].has_key('qtot'): atoms.qtot_present=True
    if dict_atoms['conf'].has_key('bemoved'): atoms.bemoved_present=True #CORRECT_IT
    if dict_atoms['conf'].has_key('velocity'): atoms.vat_present=True
    return atoms
#*****************************************************************************************
