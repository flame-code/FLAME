#!/usr/bin/env python
import subprocess
import math
import sys
import random
import os.path
from time import sleep
import yaml
#path="/home/ghasemi/FLAME/utils/python"
#sys.path.insert(1,path)
from io_yaml import *
#*****************************************************************************************
def read_posinps(posinp_dir,fn_list):
    atoms_all=[]
    for filename in fn_list:
        atoms_all_tmp=read_yaml("%s/%s" % (posinp_dir,filename))
        for atoms in atoms_all_tmp:
            atoms_all.append(Atoms())
            atoms_all[-1]=copy.deepcopy(atoms)
        del atoms_all_tmp
    print "number of configurations read %d " % len(atoms_all)
    #for atoms in atoms_all:
    #    atoms.epot/=27.211385
    return atoms_all
#*****************************************************************************************
def read_flame_in(directory):
    stream = open("%s/flame_in.yaml" % directory,"r")
    dict_flame_in = yaml.load(stream)
    stream.close()
    del stream
    return dict_flame_in
#*****************************************************************************************
def read_ann_param(directory,stypat,fn_ann_param):
    ann_param=[]
    for ap in fn_ann_param:
        filename="%s/%s/%s.ann.param.yaml.%5.5d" % (directory,ap[0],stypat,ap[1])
        stream=open(filename,"r")
        ann_param.append(yaml.load(stream))
        stream.close()
        del stream
    return ann_param
#*****************************************************************************************
def write_ann_param(filename,ann_param):
    output_stream = open(filename,"w")
    yaml.dump(ann_param,output_stream,default_flow_style=False)
    output_stream.close()
    del output_stream
#*****************************************************************************************
def read_flame_log(directory):
    stream = open("%s/flame_log.yaml" % directory,"r")
    dict_flame_log = yaml.load(stream)
    stream.close()
    del stream
    return dict_flame_log
#*****************************************************************************************
def report_result(atoms_all,epot_list,iter_list,success_list):
    dict_res={}
    rmse=0.0
    iter_avg=float(iter_list[0])
    iter_min=float(iter_list[0])
    iter_max=float(iter_list[0])
    de_diff_min=1.E10
    de_diff_max=0.0
    for iatoms in range(1,len(atoms_all)):
        if iatoms==0: print "ERROR: iatoms=0"
        natp=atoms_all[iatoms-1].nat
        nat=atoms_all[iatoms].nat
        epot_per_atom=float(atoms_all[iatoms].epot)/float(nat)
        epot_per_atom_prev=float(atoms_all[iatoms-1].epot)/float(natp)
        de_dft=(epot_per_atom-epot_per_atom_prev)
        de_ann=(epot_list[iatoms]/nat-epot_list[iatoms-1]/natp)
        rmse+=(de_dft-de_ann)**2
        de_diff_min=min(de_diff_min,abs(de_dft-de_ann))
        de_diff_max=max(de_diff_max,abs(de_dft-de_ann))
        #print "%12.4f%12.4f%12.4f%12.4f" % \
        #        (epot_per_atom,epot_per_atom_prev,epot_list[iatoms]/nat, \
        #        epot_list[iatoms-1]/natp),
        #print "%4d%8.4f%8.4f" % (nat,de_dft,de_ann)
        iter_avg+=iter_list[iatoms]
        iter_min=min(iter_min,iter_list[iatoms])
        iter_max=max(iter_max,iter_list[iatoms])
    rmse=(rmse/(len(atoms_all)-1))**0.5
    dict_res['rmse']=rmse*27.211385
    dict_res['de_diff_min']=de_diff_min*27.211385
    dict_res['de_diff_max']=de_diff_max*27.211385
    iter_avg/=len(atoms_all)
    dict_res['iter_avg']=int(iter_avg)
    dict_res['iter_min']=int(iter_min)
    dict_res['iter_max']=int(iter_max)
    return dict_res
#*****************************************************************************************
def report_result_form_ener(atoms_all,form_ener_all,form_ener_all_cent,iter_list,success_list,fn_list):
    dict_res={}
    rmse=0.0
    iter_avg=float(iter_list[0])
    iter_min=float(iter_list[0])
    iter_max=float(iter_list[0])
    de_diff_min=1.E10
    de_diff_max=0.0
    for iconf in range(1,len(atoms_all)):
        if iconf==0: print "ERROR: iconf=0"
        natp=atoms_all[iconf-1].nat
        nat=atoms_all[iconf].nat
        epot_per_atom=float(form_ener_all[iconf])/float(nat)
        epot_per_atom_prev=float(form_ener_all[iconf-1])/float(natp)
        #print "AAA %19.10E%19.10E%19.10E%19.10E" % (form_ener_all_cent[iconf-1],form_ener_all_cent[iconf],form_ener_all[iconf-1],form_ener_all[iconf])
        de_dft=(epot_per_atom-epot_per_atom_prev)
        de_ann=(form_ener_all_cent[iconf]/nat-form_ener_all_cent[iconf-1]/natp)
        #print "BBB %2d%8.3f%5d%5d" % (iconf,de_dft-de_ann,iter_list[iconf-1],iter_list[iconf])
        rmse+=(de_dft-de_ann)**2
        de_diff_min=min(de_diff_min,abs(de_dft-de_ann))
        #de_diff_max=max(de_diff_max,abs(de_dft-de_ann))
        if abs(de_dft-de_ann)>de_diff_max:
            de_diff_max=abs(de_dft-de_ann)
            fn_de_diff_max=[fn_list[iconf-1],fn_list[iconf]]
        #print "%12.4f%12.4f%12.4f%12.4f" % \
        #        (epot_per_atom,epot_per_atom_prev,epot_list[iconf]/nat, \
        #        epot_list[iconf-1]/natp),
        #print "%4d%8.4f%8.4f" % (nat,de_dft,de_ann)
        iter_avg+=iter_list[iconf]
        iter_min=min(iter_min,iter_list[iconf])
        iter_max=max(iter_max,iter_list[iconf])
    rmse=(rmse/(len(atoms_all)-1))**0.5
    dict_res['rmse']=rmse*27.211385
    dict_res['de_diff_min']=de_diff_min*27.211385
    dict_res['de_diff_max']=de_diff_max*27.211385
    iter_avg/=len(atoms_all)
    dict_res['iter_avg']=int(iter_avg)
    dict_res['iter_min']=int(iter_min)
    dict_res['iter_max']=int(iter_max)
    dict_res['fn_de_diff_max_1']=fn_de_diff_max[0]
    dict_res['fn_de_diff_max_2']=fn_de_diff_max[1]
    return dict_res
#*****************************************************************************************
def run_flame(directory,exe):
    com="cd %s ; %s >o1" % (directory,exe)
    #print "    %s" % com
    p=subprocess.Popen(com,shell=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE)
    out, err = p.communicate()
    if out!='': print out
    if err!='': print err
#*****************************************************************************************
def get_nat_per_type(atoms_all,types):
    nat_per_type_all=[]
    for atoms in atoms_all:
        nat_per_type=[0] * len(types)
        for iat in range(atoms.nat):
            for i in range(len(types)):
                if atoms.sat[iat]==types[i]:
                    nat_per_type[i]+=1
                    break
        nat_per_type_all.append(nat_per_type)
    return nat_per_type_all
#*****************************************************************************************
def get_form_ener_all(atoms_all,nat_per_type_all,ener_ref):
    form_ener_all=[]
    for iconf in range(len(atoms_all)):
        fe=0.0
        for i in range(len(nat_per_type_all[iconf])):
            fe+=nat_per_type_all[iconf][i]*ener_ref[i]
        fe=atoms_all[iconf].epot-fe
        form_ener_all.append(fe)
        #print "DFT FE %10.3f" % fe
    return form_ener_all
#*****************************************************************************************
def get_ener_ref_single_cent(stypat,run_dir,exe):
    atoms_all_out=[]
    atoms_all_out.append(Atoms())
    atoms_all_out[-1].nat=1
    atoms_all_out[-1].boundcond="free"
    atoms_all_out[-1].cellvec=[[10.0,0.0,0.0],[0.0,10.0,0.0],[0.0,0.0,10.0]]
    atoms_all_out[-1].rat.append([5.0,5.0,5.0])
    atoms_all_out[-1].sat.append(stypat)
    atoms_all_out[-1].bemoved.append("TTT")
    write_yaml(atoms_all_out,"%s/posinp.yaml" % run_dir)
    del atoms_all_out
    run_flame(run_dir,exe)
    dict_flame_log=read_flame_log(run_dir)
    str1='SQNM optimization iterations'
    str2='SQNM FINISHED'
    sqnm_finished=dict_flame_log[str1][-1][str2]
    ener_ref_single=sqnm_finished['epot']
    return ener_ref_single
#*****************************************************************************************

stream = open("input.yaml","r")
dict_runs = yaml.load(stream)
exe=dict_runs['exe']
params_dir=dict_runs['params_dir']
run_dir=dict_runs['run_dir']
posinp_dir=dict_runs['posinp_dir']
ener_ref=dict_runs['ener_ref']

atoms_all=read_posinps(posinp_dir,dict_runs['confs'])
nat_per_type_all=get_nat_per_type(atoms_all,dict_runs['types'])
form_ener_all=get_form_ener_all(atoms_all,nat_per_type_all,ener_ref)
#print form_ener_all

dict_flame_in=read_flame_in('.')

ann_param_all=[]
for stypat in dict_runs['types']:
    ann_param_all.append(read_ann_param(params_dir,stypat,dict_runs['runs']))

output_stream = open("%s/flame_in.yaml" % run_dir,"w")
yaml.dump(dict_flame_in,output_stream,default_flow_style=False)
output_stream.close()
del output_stream

fn_all=[]
for stypat in dict_runs['types']:
    fn_all.append("%s/%s.ann.param.yaml" % (run_dir,stypat))

output_stream_results = open("results.yaml","w")

for irun in range(len(dict_runs['runs'])):
    for itypat in range(len(dict_runs['types'])):
        write_ann_param(fn_all[itypat],ann_param_all[itypat][irun])
    #-------------------------------------------------------
    ener_ref_cent=[]
    for stypat in dict_runs['types']:
        ener_ref_single=get_ener_ref_single_cent(stypat,run_dir,exe)
        ener_ref_cent.append(ener_ref_single)
    #print ener_ref_cent
    #exit()
    #-------------------------------------------------------
    epot_list=[]
    iter_list=[]
    success_list=[]
    iatoms=-1
    form_ener_all_cent=[]
    for iconf in range(len(atoms_all)):
        iatoms+=1
        atoms_all_out=[]
        atoms_all_out.append(Atoms())
        atoms_all_out[-1]=copy.deepcopy(atoms_all[iconf])
        write_yaml(atoms_all_out,"%s/posinp.yaml" % run_dir)
        del atoms_all_out
        sleep(0.05)
        run_flame(run_dir,exe)
        dict_flame_log=read_flame_log(run_dir)
        str1='SQNM optimization iterations'
        str2='SQNM FINISHED'
        sqnm_finished=dict_flame_log[str1][-1][str2]
        epot_list.append(sqnm_finished['epot']) #*27.211385)
        #print sqnm_finished['epot']
        iter_list.append(sqnm_finished['iter'])
        #print sqnm_finished['iter']
        success_list.append(sqnm_finished['success'])
        print dict_runs['runs'][irun][0],dict_runs['runs'][irun][1], \
                dict_runs['confs'][iatoms]," done."
        #print "%r %d %f" % (success,iterations,epot)
        fe=0.0
        for i in range(len(nat_per_type_all[iconf])):
            fe+=nat_per_type_all[iconf][i]*ener_ref_cent[i]
            #print ener_ref_cent[i]
        fe=epot_list[-1]-fe
        #print "ANN FE %10.3f" % fe
        form_ener_all_cent.append(fe)
    #print form_ener_all_cent
    #print ener_ref
    #print ener_ref_cent
    #dict_res=report_result(atoms_all,epot_list,iter_list,success_list)
    dict_res=report_result_form_ener(atoms_all,form_ener_all,form_ener_all_cent,iter_list,success_list,dict_runs['confs'])
    print "run,epoch,RMSE: %7s %3d %8.4f" % \
            (dict_runs['runs'][irun][0],dict_runs['runs'][irun][1],dict_res['rmse'])
    del epot_list
    del iter_list
    del success_list
    del form_ener_all_cent
    dict_res['run']=dict_runs['runs'][irun][0]
    dict_res['epoch']=dict_runs['runs'][irun][1]
    output_stream_results.write(' - ')
    yaml.dump(dict_res,output_stream_results,default_flow_style=None)

del ann_param_all
output_stream_results.close()
del output_stream_results
#*****************************************************************************************
#example for input.yaml
#types: [Mg,O]
#exe: ~/bin/flame
#params_dir: set5
#run_dir: run_dir
#posinp_dir: posinps
#
#confs:
# - MgO13_00001.yaml
# - MgO13_00002.yaml
# - MgO13_00003.yaml
#
#runs:
# - [run0000,15]
# - [run0001,20]
