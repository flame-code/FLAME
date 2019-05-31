#!/usr/bin/env python
import argparse
import subprocess
import math
import sys
import random
import os.path
from time import sleep
import yaml
#path="/home/ghasemi/FLAME/utils/python"
#sys.path.insert(1,path)
#from acf import *
#*****************************************************************************************
def read_ann_param(run,stypat):
    #ann_param=[]
    #for irun in range(nrun):
    filename="%s/%s.ann.input.yaml" % (run,stypat)
    stream=open(filename,"r")
    ann_param=yaml.load(stream)
    stream.close()
    del stream
    return ann_param
#*****************************************************************************************
def read_flame_in(run):
    #ann_param=[]
    #for irun in range(nrun):
    while os.path.exists("%s/flame_in.yaml" % (run)) is True: 
        filename="%s/flame_in.yaml" % (run)
        stream=open(filename,"r")
        flame_in=yaml.load(stream)
        stream.close()
        del stream
        return flame_in
#*****************************************************************************************
def read_train_output(run):
    #ann_param=[]
    #for irun in range(nrun):
    while os.path.exists("%s/train_output.yaml" % (run)) is True: 
        filename="%s/train_output.yaml" % (run)
        stream=open(filename,"r")
        ann_param=yaml.load(stream)
        stream.close()
        del stream
        return ann_param
#*****************************************************************************************
str1="AAA"
parser=argparse.ArgumentParser(description=str1)
#parser.add_argument('fn_input',action="store",type=str,help="fn_input is the name of the input file")
#parser.add_argument('fn_output',action="store",type=str,help="fn_output is the name of the output file")
#parser.add_argument('-epoch',nargs='+',help='BBB',required=True)
str2="""Name of the input file in yaml format. In this file three keys should already be included: 
    (atom)types, epochs,runs(name of the training directories)"""
parser.add_argument('fn_input', action="store", type=str, help=str2)
parser.add_argument("-gausswidth", action='store_false', help="if present, gausswidth will be printed out.")
parser.add_argument("-hardness", action='store_false', help="if present, hardness will be printed out.")
parser.add_argument("-chi0", action='store_false', help="if present, chi0 will be printed out.")
parser.add_argument("-aself", action='store_false', help="if present, aself will be printed out.")
parser.add_argument("-chiavg", action='store_false', help="if present, chiavg will be printed out.")
parser.add_argument("-chivar", action='store_false', help="if present, chivar will be printed out.")
parser.add_argument("-rmse", action='store_false', help="if present, rmse will be printed out.")
parser.add_argument("-frmse", action='store_false', help="if present, frmse will be printed out.")
parser.add_argument("-qavg", action='store_false', help="if present, qavg will be printed out.")
parser.add_argument("-qvar", action='store_false', help="if present, qvar will be printed out.")
parser.add_argument("-qratio", action='store_false', help="if present, qratio will be printed out.")
args=parser.parse_args()
gausswidth=not args.gausswidth
hardness=not args.hardness
chi0=not args.chi0
aself=not args.aself
chiavg=not args.chiavg
chivar=not args.chivar
rmse=not args.rmse
frmse=not args.frmse
qavg=not args.qavg
qvar=not args.qvar
qratio=not args.qratio
#-----------------------------------------------------------------------------------------
#print args.epoch
stream = open(args.fn_input,"r")
dict_input = yaml.load(stream)
#exe=dict_runs['exe']
#params_dir=dict_runs['params_dir']
#run_dir=dict_runs['run_dir']
#posinp_dir=dict_runs['posinp_dir']
runs=dict_input['runs']
epochs=dict_input['epochs']
types=dict_input['types']
#print runs
#print epochs
#print types

if len(runs)==0:
    print "ERROR: runs is empty in %s " % args.fn_input
    exit()

if len(epochs)==0:
    print "ERROR: epochs is empty in %s " % args.fn_input
    exit()

#fd=open('crashed_runs', 'w')
#print chi0
for run in runs:
    #epoch=epochs[0]
    for epoch in epochs:
        train_output=read_train_output(run)
        #print type(train_output)
        #print type(train_output['training iterations'])
        if not train_output: continue
        if train_output['training iterations'] is None: continue
        #print len(train_output['training iterations']),epoch+1
        if len(train_output['training iterations'])<epoch+1: continue
    
        if train_output['training iterations'][epoch]!=None:
            print "%s" % run,
            print "%3d" % epoch,
            ann_params={}
            for typ in types:
                ann_params[typ]=read_ann_param(run,typ)
            flame_in=read_flame_in(run)
            if gausswidth:
                for typ in types:
                    print "%5.2f" % ann_params[typ]['main']['gausswidth'],
            if hardness:
                for typ in types:
                    print "%5.2f" % ann_params[typ]['main']['hardness'],
            if chi0:
                for typ in types:
                    print "%5.2f" % ann_params[typ]['main']['chi0'],
            if chiavg:
                for typ in types:
                    str1="chi_%s" % typ
                    print "%6.3f" % train_output['training iterations'][epoch][str1]['chiavg'],
            if chivar:
                for typ in types:
                    str1="chi_%s" % typ
                    print "%6.3f" % train_output['training iterations'][epoch][str1]['chivar'],
            if rmse:
                if train_output['training iterations'][epoch]['valid']:
                    print "%8.3f" % train_output['training iterations'][epoch]['train']['rmse'],
                    print "%8.3f" % train_output['training iterations'][epoch]['valid']['rmse'],
            if frmse:
                if train_output['training iterations'][epoch]['valid']:
                    print "%6.3f" % train_output['training iterations'][epoch]['train']['frmse'],
                    print "%6.3f" % train_output['training iterations'][epoch]['valid']['frmse'],
            if qavg:
                for typ in types:
                    str1="charge_%s" % typ
                    print "%6.3f" % train_output['training iterations'][epoch][str1]['qavg'],
            if qvar:
                for typ in types:
                    str1="charge_%s" % typ
                    print "%6.3f" % train_output['training iterations'][epoch][str1]['qvar'],
            if qratio:
                for typ in types:
                    str1="charge_%s" % typ
                    tt1=train_output['training iterations'][epoch][str1]['qvar']
                    tt2=train_output['training iterations'][epoch][str1]['qavg']
                    #print type(tt1)
                    #print type(tt2)
                    print "%6.3f" % abs(tt1/tt2),
            if aself:
                print "%6.3f%5.2f" % (flame_in['main']['aself'][1], flame_in['main']['aself'][2]),
        #else:
        #    print 'None',
                #exit()
        #print ann_params.keys()
        #filename="run%4.4d/%s.ann.input.yaml" % (irun,stypat)
        #stream=open(filename,"r")
        #ann_param.append(yaml.load(stream))
        #stream.close()
        #del stream
            print ""
        else:
            #print >> fd, run
            continue
        #for epoch in epochs:
        

