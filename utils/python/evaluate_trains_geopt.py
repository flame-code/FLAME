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
from acf import *
#*****************************************************************************************
def read_posinps(posinp_dir,fn_list):
    atoms_all=[]
    for filename in fn_list:
        atoms_all_tmp=acf_read("%s/%s" % (posinp_dir,filename))
        for atoms in atoms_all_tmp:
            atoms_all.append(Atoms())
            atoms_all[-1]=copy.deepcopy(atoms)
        del atoms_all_tmp
    print "number of configurations read %d " % len(atoms_all)
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
    dict_res['rmse']=rmse
    dict_res['de_diff_min']=de_diff_min
    dict_res['de_diff_max']=de_diff_max
    iter_avg/=len(atoms_all)
    dict_res['iter_avg']=int(iter_avg)
    dict_res['iter_min']=int(iter_min)
    dict_res['iter_max']=int(iter_max)

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

stream = open("input.yaml","r")
dict_runs = yaml.load(stream)
exe=dict_runs['exe']
params_dir=dict_runs['params_dir']
run_dir=dict_runs['run_dir']
posinp_dir=dict_runs['posinp_dir']

atoms_all=read_posinps(posinp_dir,dict_runs['confs'])

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

output_stream = open("results.yaml","w")

for irun in range(len(dict_runs['runs'])):
    for itypat in range(len(dict_runs['types'])):
        write_ann_param(fn_all[itypat],ann_param_all[itypat][irun])
    epot_list=[]
    iter_list=[]
    success_list=[]
    iatoms=-1
    for atoms in atoms_all:
        iatoms+=1
        atoms_all_out=[]
        atoms_all_out.append(Atoms())
        atoms_all_out[-1]=copy.deepcopy(atoms)
        acf_write(atoms_all_out,"%s/posinp.acf" % run_dir)
        del atoms_all_out
        sleep(0.05)
        run_flame(run_dir,exe)
        dict_flame_log=read_flame_log(run_dir)
        str1='SQNM optimization iterations'
        str2='SQNM FINISHED'
        sqnm_finished=dict_flame_log[str1][-1][str2]
        epot_list.append(sqnm_finished['epot']*27.211385)
        iter_list.append(sqnm_finished['iter'])
        #print sqnm_finished['iter']
        success_list.append(sqnm_finished['success'])
        print dict_runs['runs'][irun][0],dict_runs['runs'][irun][1], \
                dict_runs['confs'][iatoms]," done."
        #print "%r %d %f" % (success,iterations,epot)
    dict_res=report_result(atoms_all,epot_list,iter_list,success_list)
    print "run,epoch,RMSE: %7s %3d %8.4f" % \
            (dict_runs['runs'][irun][0],dict_runs['runs'][irun][1],dict_res['rmse'])
    del epot_list
    del iter_list
    del success_list
    dict_res['run']=dict_runs['runs'][irun][0]
    dict_res['epoch']=dict_runs['runs'][irun][1]
    output_stream.write(' - ')
    yaml.dump(dict_res,output_stream,default_flow_style=None)

del ann_param_all
output_stream.close()
del output_stream
#*****************************************************************************************
#example for input.yaml
#types: [Mg,O]
#exe: ~/bin/flame
#params_dir: set5
#run_dir: run_dir
#posinp_dir: posinps
#
#confs:
# - MgO13_00001.acf
# - MgO13_00002.acf
# - MgO13_00003.acf
#
#runs:
# - [run0000,15]
# - [run0001,20]
