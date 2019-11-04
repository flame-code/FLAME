#!/usr/bin/env python
import sys,copy
import argparse
from argparse import ArgumentParser
import subprocess
import os
import textwrap as _textwrap
import numpy as np

#----------------------------------------------------------------------------------------#
#                                                                                        #
#       This script is prepared to read different data produced during training process  #
#       and print them according to the user requests.                                   #
#       I will work on it to sort read data.                                             #
#       S. Faraji 30/06/2017                                                             #
#                                                                                        #
#----------------------------------------------------------------------------------------#

str1="=============================================================================="
str1+="\n Please choose your desired options to print in the output.                 "
str1+="\n NOTICE 1: If you want to have sorted data according to your desire parameter, you should provide a value."
str1+="\n\n NOTICE 2: The symbol of elements in the compound is required."
parser = argparse.ArgumentParser(description=str1)
#parser.add_argument('type', types=str, nargs='+')
parser.add_argument('-c0' , action="store"     ,type=str,   dest="typesin", required=True , help="The symbol of species in compound")
parser.add_argument('-c1' , action="store"     ,type=int,   dest="epocho" , default = 0,    help="Only consider data of this epoch")
parser.add_argument('-c2' , action="store"     ,type=int,   dest="epocha" , default = 0,    help="Consider data larger than this epoch")
parser.add_argument('-c3' , action="store_true",            dest="case3"  , help="Print RMSE of energy")
parser.add_argument('-c4' , action="store"     ,type=float, dest="case4"  , help="Print RMSE of energy lower than this value.")
parser.add_argument('-c5' , action="store_true",            dest="case5"  , help="Print the value of fat_chi/fat_es")
parser.add_argument('-c6' , action="store"     ,type=float, dest="case6"  , help="Print the value of fat_chi/fat_es lower than this value.")
parser.add_argument('-c7' , action="store_true",            dest="case7"  , help="Print RMSE of force")
parser.add_argument('-c8' , action="store"     ,type=float, dest="case8"  , help="Print the RMSE of force lower than this value.")
parser.add_argument('-c9' , action="store_true",            dest="case9"  , default=False, help="average of chi")
parser.add_argument('-c10', action="store_true",            dest="case10" , default=False, help="change of chi")
parser.add_argument('-c11', action="store_true",            dest="case11" , default=False, help="average of charge")
parser.add_argument('-c12', action="store_true",            dest="case12" , default=False, help="change of charge")
parser.add_argument('-c13', action="store_true",            dest="case13" , default=False, help="print amplitude of chi")
parser.add_argument('-c14', action="store_true",            dest="case14" , default=False, help="print prefactor of chi")
parser.add_argument('-c15', action="store_true",            dest="case15" , default=False, help="print atomic reference energy")
parser.add_argument('-c16', action="store_true",            dest="case16" , default=False, help="print electronic gausswidth")
parser.add_argument('-c17', action="store_true",            dest="case17" , default=False, help="print hardness")
parser.add_argument('-c18', action="store_true",            dest="case18" , default=False, help="print chi0")
parser.add_argument('-c19', action="store_true",            dest="case19" , default=False, help="print ionic gausswidth")
parser.add_argument('-c20', action="store_true",            dest="case20" , default=False, help="print z_ionic")
parser.add_argument('-c21', action="store_true",            dest="case21" , default=False, help="print spring constant")
choice = parser.parse_args()
#print "epocho:    ",choice.epocho
#-----------------------------------------------------------------------------------------
if len(sys.argv) < 3:
    print "usage: getgoodparam.py -h"
    exit()

#~~~~~~~Functions to Read data from different files~~~~~~~~~~~~~~~~~~~~~~
def read_input_ann(filename):
    f = open (filename,"r")
    iline=0
    for strline in f.readlines():
        iline+=1
        #if iline==3:
        #    ampl_chi.append(np.float(strline.split()[1]))
        #    prefac_chi.append(np.float(strline.split()[3]))
        if iline==4:
            zion  = float(strline.split()[1])
            gwion = float(strline.split()[3])
            E0    = float(strline.split()[5])
        if iline==5:
            gw   = float(strline.split()[1])
            hard = float(strline.split()[3])
            chi0 = float(strline.split()[5])
        if iline==6:
            spring_cte = float(strline.split()[1])
            break
    f.close
    return (zion,gwion,E0,gw,hard,chi0,spring_cte)

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def read_err_train(filename,epoch):
    f = open (filename,"r")
    lines = f.readlines()
    f.close
    iline = epoch
    error=[float(i) for i in lines[iline+0].split()[0:5]]
    return error

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def read_chi(filename,epoch):
    f = open (filename,"r")
    lines = f.readlines()
    f.close
    iline = epoch
    datachi = [float(i) for i in lines[iline+0].split()[0:5]]
    return datachi

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def read_charge(filename,epoch):
    f = open (filename,"r")
    lines = f.readlines()
    f.close
    iline = epoch
    datacharge = [float(i) for i in lines[iline+0].split()[0:5]]
    return datacharge

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##def prepare_names():
##    typesin = raw_input('\n>>> Please enter the symbole of species in your compound:  ')
##    types = typesin.split()
##    list_app=[]
##    for itype in range(len(types)):
##        list_app.append("data_ann_e0_"  + str(types[itype]))
##        list_app.append("data_ann_gw_"  + str(types[itype])) 
##        list_app.append("data_ann_hard_"+ str(types[itype]))
##        list_app.append("data_ann_chi0_"+ str(types[itype])) 
##        list_app.append("data_chi_"     + str(types[itype])) 
##        list_app.append("data_charge_"  + str(types[itype])) 
##    return list_app
##
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_epoch_all(filename):
    f = open (filename,"r")
    lines = f.readlines()
    f.close
    epoch_all = -1
    for iline,line in enumerate(lines):
        epoch_all = epoch_all + 1
    return epoch_all

#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
def get_header(choice):
    header = {}
    if choice.case3:  header['case3']  = 'RMSE_E'
    if choice.case4:  header['case4']  = 'RMSE_E'
    if choice.case5:  header['case5']  = 'Fchi/Fes'
    if choice.case6:  header['case6']  = 'Fchi/Fes'
    if choice.case7:  header['case7']  = 'RMSE_F'
    if choice.case8:  header['case8']  = 'RMSE_F'
    if choice.case9:  header['case9']  = 'chi_avg'
    if choice.case10: header['case10'] = 'dchi'
    if choice.case11: header['case11'] = 'q'
    if choice.case12: header['case12'] = 'dq'
    if choice.case12: header['case13'] = 'Amp_chi'
    if choice.case12: header['case14'] = 'prefac_chi'
    if choice.case15: header['case15'] = 'ener_ref'
    if choice.case16: header['case16'] = 'gwe'
    if choice.case17: header['case17'] = 'hard'
    if choice.case18: header['case18'] = 'chi0'
    if choice.case18: header['case19'] = 'gwi'
    if choice.case18: header['case20'] = 'z_ion'
    if choice.case18: header['case21'] = 'k'
    #for i,f in header.iteritems():
    #    print f,
    #print header['case3']
    return header

#%%%%%%%% Main Program %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
#typesin = raw_input('\n>>> Please enter the symbole of species in your compound:  ')
header =  get_header(choice)

fout = open("out_result","w")

if choice.case4 or choice.case6 or choice.case8:
    output = open("sorted_output","w")

types = choice.typesin.split()
#print types

f = open('list_run','r')
nrun = []
for line in f.readlines():
    nrun.append([str(line.split()[0])])
f.close()

rundir        = []
pwd           = os.getcwd()
error_train   = [] 
data_ann_e0   = [] 
data_ann_zion = []
data_ann_gwion= []
data_ann_gw   = []
data_ann_hard = []
data_ann_chi0 = []
data_ann_k    = []
data_chi      = []
data_charge   = []
epoch         = []

path_in= os.getcwd()

iepoch = -1
for irun in range(len(nrun)):
    path = "%s/run%4.4d" % (pwd,irun)
    rundir.append("run%4.4d"%(irun))
    os.chdir(path)
    if irun==0:
        epoch_all = get_epoch_all('err_train')
        if choice.epocho:
            epoch.append(choice.epocho)
        else:
            iepoch = choice.epocha
            while iepoch <= epoch_all:
                epoch.append(iepoch)
                iepoch = iepoch + 1
    error_train1 = []
    if choice.epocho:
        error_train1.append(read_err_train('err_train',iepoch))
    else:
        iepoch = choice.epocha
        while iepoch <= epoch_all:
            error_train1.append(read_err_train('err_train',iepoch))
            iepoch = iepoch + 1
    error_train.append(error_train1)
    
    data_chi1        = []
    data_charge1     = []
    data_ann_e01     = []
    data_ann_gw1     = []
    data_ann_gw11    = []
    data_ann_z1      = []
    data_ann_hard1   = []
    data_ann_chi01   = []
    data_ann_spring1 = []

    for itype in range(len(types)):
        inann="input.ann."+types[itype] 
        inchi="chi."+types[itype]
        incharge="charge."+types[itype]

        ann_zion,ann_gwion,ann_e0,ann_gw,ann_hard,ann_chi0,spring_k=read_input_ann(inann)
        data_ann_z1.append(ann_zion)
        data_ann_gw11.append(ann_gwion)
        data_ann_e01.append(ann_e0)
        data_ann_gw1.append(ann_gw)
        data_ann_hard1.append(ann_hard)
        data_ann_chi01.append(ann_chi0)
        data_ann_spring1.append(spring_k)

        data_chi2    = []
        data_charge2 = []
        if choice.epocho:
            iepoch = choice.epocho
            data_chi2.append(read_chi(inchi,iepoch))
            data_charge2.append(read_charge(incharge,iepoch))
        else:
            iepoch = choice.epocha
            while iepoch <= epoch_all:
                data_chi2.append(read_chi(inchi,iepoch))
                data_charge2.append(read_charge(incharge,iepoch))
                iepoch = iepoch + 1
        data_chi1.append(data_chi2)
        data_charge1.append(data_charge2)
    data_chi.append(data_chi1)
    data_charge.append(data_charge1)
    data_ann_e0.append(data_ann_e01)
    data_ann_zion.append(data_ann_z1)
    data_ann_gw.append(data_ann_gw1)
    data_ann_gwion.append(data_ann_gw11)
    data_ann_hard.append(data_ann_hard1)
    data_ann_chi0.append(data_ann_chi01)
    data_ann_k.append(data_ann_spring1)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

fout.write("%s" %('#run_dir'))
fout.write(" %5s" %('epoch'))

if choice.case15:
    for itype in range(len(types)):
        fout.write(" %13s" %(header['case15']+"%s" %('_') +"%s" %(types[itype])))
if choice.case16:
    for itype in range(len(types)):
        fout.write(" %6s" %(header['case16']+"%s" %('_') +"%s" %(types[itype])))
if choice.case17:
    for itype in range(len(types)):
        fout.write(" %8s" %(header['case17']+"%s" %('_') +"%s" %(types[itype])))
if choice.case18:
    for itype in range(len(types)):
        fout.write(" %7s" %(header['case18']+"%s" %('_') +"%s" %(types[itype])))
if choice.case19:
    for itype in range(len(types)):
        fout.write(" %7s" %(header['case19']+"%s" %('_') +"%s" %(types[itype])))
if choice.case20:
    for itype in range(len(types)):
        fout.write(" %7s" %(header['case20']+"%s" %('_') +"%s" %(types[itype])))
if choice.case21:
    for itype in range(len(types)):
        fout.write(" %7s" %(header['case21']+"%s" %('_') +"%s" %(types[itype])))
#if choice.case13:
#    fout.write(" %10s" %(header['case13']))
#if choice.case14:
#    fout.write(" %10s" %(header['case14']))
if choice.case9:
    for itype in range(len(types)):
        fout.write(" %8s" %(header['case9']+"%s" %('_') +"%s" %(types[itype])))
if choice.case10:
    for itype in range(len(types)):
        fout.write(" %6s" %(header['case10']+"%s" %('_') +"%s" %(types[itype])))
if choice.case11:
    for itype in range(len(types)):
        fout.write(" %7s" %(header['case11']+"%s" %('_') +"%s" %(types[itype])))
if choice.case12:
    for itype in range(len(types)):
        fout.write(" %9s" %(header['case12']+"%s" %('_') +"%s" %(types[itype])))
if choice.case3:
    fout.write(" %9s" %(header['case3']))
if choice.case4:
    fout.write(" %9s" %(header['case4']))
if choice.case5:
    fout.write(" %9s" %(header['case5']))
if choice.case6:
    fout.write(" %9s" %(header['case6']))
if choice.case7:
    fout.write(" %9s" %(header['case7']))
if choice.case8:
    fout.write(" %9s" %(header['case8']))
fout.write("%s\n" %(''))

#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

for irun in range(len(nrun)):
    for i,j in enumerate(epoch):
        fout.write("%s" %(rundir[irun]) + "%5d" %(j))
        if choice.case15:
            for itype in range(len(types)):
                fout.write("%15.7f " %(data_ann_e0[irun][itype]))
        if choice.case16:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_ann_gw[irun][itype]))
        if choice.case17:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_ann_hard[irun][itype]))
        if choice.case18:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_ann_chi0[irun][itype]))
        if choice.case19:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_ann_gwion[irun][itype]))
        if choice.case20:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_ann_zion[irun][itype]))
        if choice.case21:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_ann_k[irun][itype]))
        if choice.case13:
            print "Sorry! The script doesn't work for this option now." 
        if choice.case14:
            print "Sorry! The script doesn't work for this option now." 
        if choice.case9:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_chi[irun][itype][i][1]))
        if choice.case10:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_chi[irun][itype][i][4]))
        if choice.case11:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_charge[irun][itype][i][1]))
        if choice.case12:
            for itype in range(len(types)):
                fout.write("%8.3f " %(data_charge[irun][itype][i][4]))
        if choice.case3:
            fout.write("%8.3f " %(error_train[irun][i][1]))
        if choice.case4:
            fout.write("%8.3f " %(error_train[irun][i][1]))
        if choice.case5:
            fout.write("%8.3f " %(error_train[irun][i][2]))
        if choice.case6:
            fout.write("%8.3f " %(error_train[irun][i][2]))
        if choice.case7:
            fout.write("%8.3f " %(error_train[irun][i][4]))
        if choice.case8:
            fout.write("%8.3f " %(error_train[irun][i][4]))
        fout.write("%s\n" %(''))
fout.close()
os.chdir(path_in)
#%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if choice.case4 or choice.case6 or choice.case8:
    #print "HERE"
    #print os.getcwd()
    ff = open("out_result",'r')
    #for line in ff.readlines():
    #    print line
    lines = copy.deepcopy(ff.readlines())
    ff.close()
    data = []
    linetot = 0
    for iline,line in enumerate(lines):
        data_temp = [str(i) for i in lines[iline+0].split()]
        data.append(data_temp)    
        linetot += 1
    #for i,j in enumerate(data[??]):
    #    output.write("%10s" %(j))
#find = []
#for i in range(len(data[0])):
#    if data[0][i]=='RMSE_E':
#        find.append(i)
#    if data[0][i]=='RMSE_F':
#        find.append(i)
#    if data[0][i]=='Fchi/Fes':
#        find.append(i)
#print find
#print data[0]
###
####log1 = False ; log2 = False ; log3 = False
####if choice.case4 and choice.case6 and choice.case8:
####    log1 = True ; log2 = True ; log3 = True
####    for i,j in enumerate(data[??]):
####        if data[??]< choice.case4:
####            if data[??]< choice.case6:
####                if data[??]<choice.case8:
####                    output.write("%10s" %(j))
####if log1==False and log2 == False:
####    if choice.case4 and choice.case6:
####        log1 = True ; log2 = True
####        for i,j in enumerate(data[??]):
####            if data[??]< choice.case4:
####                if data[??]< choice.case6:
####                    output.write("%10s" %(j))
####if log1==False and log3 == False:
####    if choice.case4 and choice.case8:
####        log1 = True ; log3 = True
####        for i,j in enumerate(data[??]):
####            if data[??]< choice.case4:
####                if data[??]< choice.case8:
####                    output.write("%10s" %(j))
####if log2==False and log3 == False:
####    if choice.case8 and choice.case6:
####        log2 = True ; log3 = True
####        for i,j in enumerate(data[??]):
####            for i in range(linetot):
####                if data[??]< choice.case8:
####                    if data[??]< choice.case6:
####                        output.write("%10s" %(j))
####if log1==False:
####    if choice.case4:
####        log1 = True
####        for i,j in enumerate(data[??]):
####            if data[??]<choice.case4:
####                output.write("%10s" %(j))
####if log2==False:
####    if choice.case6:
####        log2 = True
####        for i,j in enumerate(data[??]):
####            if data[??]<choice.case6:
####                output.write("%10s" %(j))
####if log3==False:
####    if choice.case8:
####        log3 = True
####        for i,j in enumerate(data[??]):
####            if data[??]<choice.case6:
####                output.write("%10s" %(j))
#------Giving different options to user:--------------------------
#class MultilineFormatter(argparse.HelpFormatter):
#    def _fill_text(self, text, width, indent):
#        text = self._whitespace_matcher.sub(' ', text).strip()
#        paragraphs = text.split('|n ')
#        multiline_text = ''
#        for paragraph in paragraphs:
#            formatted_paragraph = _textwrap.fill(paragraph, width, initial_indent=indent, subsequent_indent=indent) + '\n\n'
#            multiline_text = multiline_text + formatted_paragraph
#        return multiline_text
#str2="""options:   
#             |n
#             a = epoch
#             |n
#             b = RMSE of energy
#             |n
#             c = fat_chi/fat_es 
#             |n
#             d = RMSE of force
#             |n
#             e = average of chi
#             |n
#             f = change of chi
#             |n
#             g = average of charge
#             |n
#             h = change of charge
#             |n
#             i = amplitude of chi
#             |n
#             j = prefactor of chi
#             |n
#             k = atomic reference energy
#             |n
#             l = gausswidth
#             |n
#             m = hardness
#             |n
#             n = chi0"""
#parser = argparse.ArgumentParser(description=str2,usage='%(prog)s [a1 a2 b c d e f g h i j k l m n]',formatter_class=MultilineFormatter)
