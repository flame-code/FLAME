#!/usr/bin/env python
import argparse
import commands
import string
import re
#*****************************************************************************************
def is_routine(routine):
    str_routine=str(routine).strip()
    str_routine=str_routine.lower()
    i_sub=str_routine.find('subroutine')
    i_fun=str_routine.find('function')
    check_routine=False
    if (i_sub==0 or i_fun==0): check_routine=True
    i_rrb=str_routine.find(')')
    if i_rrb>-1 and i_rrb+1<len(str_routine):
        str_routine=str_routine[i_rrb+1:]
        str_routine=str_routine.strip()
        i_bind=str_routine.find('bind')
        if i_bind==0: check_routine=False
    return check_routine
#*****************************************************************************************
def check_var_define(line):
    tt=str(line).strip()
    tt=tt.lower()
    tt=tt.replace(" ","")
    i_use=tt.find('use')
    i_mod_interface=tt.find('mod_interface')
    i_comment=tt.find('!')
    i_implicit=tt.find('implicit none')
    i_integer=tt.find('integer')
    i_real=tt.find('real')
    i_complex=tt.find('complex')
    i_external=tt.find('external')
    if tt.find('external_')>-1: i_external=-1
    i_character=tt.find('character')
    i_logical=tt.find('logical')
    i_type=tt.find('type(')
    stat_line='unknown'
    if i_use==0:
        if i_mod_interface>i_use:
            stat_line='mod_interface'
        else:
            stat_line='use'
    if i_comment==0:
        if not stat_line=='unknown':
            print 'ERROR: comment or ',stat_line
        else:
            stat_line='comment'
    if i_implicit==0:
        if not stat_line=='unknown':
            print 'ERROR: implicit or ',stat_line
        else:
            stat_line='implicit'
    if (i_integer==0 or i_real==0 or i_complex==0 or i_character==0 or \
            i_logical==0 or i_type==0 or i_external==0):
        if not stat_line=='unknown':
            print 'ERROR: definition or ',stat_line
        else:
            stat_line='definition'
    if stat_line=='unknown': stat_line='command'
    return stat_line
#*****************************************************************************************
def get_arguments(line):
    tt=str(line).strip()
    tt=tt.lower()
    i_parenthesis=tt.find('(')
    if i_parenthesis>0:
        tt=tt.replace(" ","")
        if '!' in tt:
            tt=tt.split("!")[0]
        i_fun=tt.find('function')
        i_res=tt.find('result')
        if (i_fun>-1 and i_res>-1):
            tt=tt.replace("(",",",1)
            tt=tt.replace(")",",",1)
        tt=tt.replace("(",",",1)
        tt=tt.replace(")","",1)
        ttl=tt.split(",")
        if (i_fun>-1 and i_res>-1):
            ttl.pop(0)
            ttl.remove('result')
        if (i_fun>-1 and not i_res>-1):
            ttl.insert(0,str(ttl[0]).replace("function","",1))
            ttl.pop(1)
    else:
        ttl=[]
    return ttl
#*****************************************************************************************
def get_var_list(line):
    tt=str(line).strip()
    tt=tt.lower()
    i_comment=tt.find('!')
    if(i_comment>-1):
        tt_list=list(tt)
        n=len(tt_list)
        for i in range(i_comment,n):
            tt_list[i]=""
        tt="".join(tt_list)
    tt=tt.replace("use","",1)
    ttl_mod=tt.split(",")[0]
    tt=tt.replace(ttl_mod,"",1)
    tt=tt.replace(",","",1)
    tt=tt.replace("only","",1)
    tt=tt.replace(":","",1)
    tt=tt.replace("except_this_one","",1)
    tt=tt.replace("=>","",1)
    tt=tt.replace(" ","")
    if not tt=="":
        ttl=tt.split(",")
    else:
        ttl=[]
    return ttl
#*****************************************************************************************
def get_subroutine_name(line):
    tt=str(line).strip()
    tt=tt.lower()
    i_parenthesis=tt.find('(')
    if i_parenthesis>0:
        ttl=tt.split()[1].split("(")[0]
    else:
        ttl=tt.split()[1]
    procedure=tt.split()[0]
    return procedure,ttl
#*****************************************************************************************
def get_files():
    cmd='find . -not -path "./futile-suite/*" -iname "*.F90" | LC_COLLATE=en_US.UTF-8 sort'
    str_files=commands.getoutput(cmd)
    files=str_files.splitlines()
    exclude_files=['./src/interface_mod.F90',
                   './wrappers/energyandforces_siesta.f90',
                   './wrappers/energyandforces_openmx.f90',
                   './src/MPIfake.F90',
                   './src/task_netsock.F90',
                   './src/constants_minhocao_mod.F90',
                   './src/minhocao_mod.F90',
                   './src/lammps_mod.F90',
                   './src/fsockets.F90',
                   './src/fingerprint_BCM.F90',
                   './src/fingerprint_gaussmol.F90',
                   './src/recompute_kpt.F90', 
                   './src/optimizer_nlbfgs_minhocao.F90',
                   './src/optimizer_sqnm_minhocao_module.F90',
                   './src/optimizer_subs_minhocao.F90', #keep this file excluded: duplicate routines
                   './src/potential_LAMMPS_interface.F90',
                   './src/POSCAR2ascii.F90',         #Do not touch this !!!"
                   './src/potential_abinit.F90',
                   './src/potential_alborz.F90',
                   './src/ascii2POSCAR.F90',         #Do not touch this !!!"
                   './src/binaries.F90',             #Do not touch this !!!"
                   './src/vasp_recompute_cell.F90',    #Do not touch this !!!"
                   './src/vasp_recompute_kpt.F90',     #Do not touch this !!!"     
                   './src/vasp_recompute_kpt_odd.F90',  #Do not touch this !!!"
                   './src/ternaries.F90',              #Do not touch this !!!"
                   './src/convex_hull.F90',           #Do not touch this !!!"
                   './src/potential_BLJ_minhocao.F90',
                   './src/cell_utils.F90',
                   './src/potential_MLJ.F90',
                   './src/potential_corerepulsion.F90',
                   './src/potential_CP2K.F90',
                   './src/potential_DFTB_minhocao.F90',
                   './src/potential_EDIP.F90',
                   './src/potential_PWSCF.F90',
                   './src/PWSCF_restruct.F90',
                   './src/quaternions.F90',
                   './src/expand_poslows.F90',
                   './src/find_symmetry.F90',
                   './src/fingerprint_oganov.F90',
                   './src/fingerprint_oganov_cont.F90',
                   './src/fingerprint_XYZ2SM.F90',
                   './src/potential_main_minhocao.F90',
#                   './src/io_ascii.F90',
#                   './src/io_vasp_minhocao.F90',
                   './src/potential_IPI.F90',
                   './src/potential_LAMMPS.F90',
                   './src/potential_LenoskyMEAM.F90',
                   './src/lenosky_tb/',
                   './src/potential_LenoskyTB_minhocao.F90',
                   './src/potential_LenoskyTB_LJ_minhocao.F90',
                   './src/potential_LJ_voids.F90',
                   './src/potential_MOPAC.F90',
                   './src/potential_MSOCK.F90',
                   './src/msock_slave_template.F90',
                   './src/parser_core_minhocao.F90',
                   './src/potential_SIESTA_minhocao.F90',
                   './src/potential_TERSOFF.F90',
                   './src/potential_TINKER.F90',
                   './src/potential_VASP_minhocao.F90',
                   './src/train_optimizer.F90',
                   './src/ann_train.F90',
                   './src/atoms_mod.F90',
                   './src/barsaddle.F90',
                   './src/ann_mod.F90',
                   './src/io_acf.F90',
                   './src/io_yaml_conf.F90',
                   './src/io_bin.F90',
                   './src/md_util.F90',
                   './src/parser_core.F90',
                   ]
    
    for file in exclude_files:
        if file in files: files.remove(file)
    
    str_interfaced=commands.getoutput('grep -i subroutine ./src/interface_mod.F90')
    routine_interfaced=str_interfaced.splitlines()
    return files
#*****************************************************************************************
def make_one_line(lines,iline):
    if str(lines[iline]).find("(")==-1:
        return lines[iline],iline
    #ii=str(lines[iline]).find('gauss_grid')
    #if ii==-1: return lines[iline],iline
    oneline_arg=""
    for i in range(100):
        tt=str(lines[iline+i])
        i_comment=tt.find('!')
        if i_comment>-1: tt=tt[:i_comment]
        oneline_arg=oneline_arg+tt
        #print oneline_arg.replace("&","")
        oneline_arg=oneline_arg.replace("&","")
        #print tt
        #print oneline_arg
        if len(tt.rstrip())>0:
            if tt.rstrip()[-1]==')': break
        #if i_comment>0: break
    return oneline_arg,iline+i
#*****************************************************************************************
def routines_missing_use_interface(file,lines,routines,iline_routines):
    excluded_routines=['sort_alborz','sort2_alborz',
        'ylm_mathematica',
        'spg_cell_primitive',
        'spg_cell_refine',
        'get_spg',
        'final_netsock',
        'rotate_stresstensor_other',
        'init_lenosky_tb',
        'conf_latforce',
        'alborz_initialize_timing_categories',
        'erf_over_r_taylor',
        'fp_assign',
        'params_echo',
        'params_check',
        'NLB1',
        'NMCSRCH',
        'NMCSTEP',
        'qe_volume',
        'invmat',
        'cell_force',
        'cryst_to_cart',
        'recips'
        ]
    filename_printed=False
    for iroutine in range(len(routines)):
        tt=str(lines[iline_routines[iroutine]]).split()[1]
        i_roundbracket=tt.find('(')
        if i_roundbracket>-1: tt=tt[0:i_roundbracket]
        if tt in excluded_routines: continue
        if iroutine==len(routines)-1:
            nline_routine=len(lines)-iline_routines[iroutine]-1
        else:
            nline_routine=iline_routines[iroutine+1]-iline_routines[iroutine]-1
        uses_modinterface=False
        for iline in range(nline_routine):
	    tt=str(lines[iline_routines[iroutine]+iline]).strip()
            tt=tt.lower()
            tt=tt.replace(" ","")
            i_modinterface=tt.find('usemod_interface')
            i_comment=tt.find('!')
            if i_comment==0: continue
            if i_modinterface==0:
                uses_modinterface=True
                break
        if not uses_modinterface:
            if not filename_printed:
                print \
                "\nRoutines in file %s do not use mod_interface!" % file
                filename_printed=True
            print routines[iroutine]
#*****************************************************************************************
def make_interface_routines(lines,routines,iline_routines,fout):
    for iroutine in range(len(routines)):
        oneline_arg,iline_arg_end=make_one_line(lines,iline_routines[iroutine])
        nline_arg=iline_arg_end-iline_routines[iroutine]+1
        for iline in range(nline_arg):
            fout.write("%s\n" % lines[iline_routines[iroutine]+iline])
        #fout.write("    use mod_defs\n")
        arg_list=get_arguments(oneline_arg)
        if iroutine==len(routines)-1:
            nline_routine=len(lines)-iline_routines[iroutine]-1
        else:
            nline_routine=iline_routines[iroutine+1]-iline_routines[iroutine]-1
        definitions=[]
        for k in range(nline_routine):
            stat_line=check_var_define(lines[iline_routines[iroutine]+k+1])
            if stat_line=='definition':
                definitions.append(lines[iline_routines[iroutine]+k+1])
        stat_implicit_none=False
        for k in range(nline_routine):
            str_line=str(lines[iline_routines[iroutine]+k+1])
            str_line=str_line.lower()
            i_endsub=str_line.find("end subroutine")
            i_endfun=str_line.find("end function")
            if i_endsub>-1 or i_endfun>-1: break
            stat_line=check_var_define(lines[iline_routines[iroutine]+k+1])
            print_line=False
            if stat_line=='use':
                var_list=get_var_list(lines[iline_routines[iroutine]+k+1])
                for m in range(len(definitions)):
                    tt=str(definitions[m]).strip()
                    tt=tt.lower()
                    for l in range(len(var_list)):
                        ttt=r"\b"+var_list[l]+r"\b"
                        p=re.compile(ttt)
                        if bool(p.search(tt)):
                            print_line=True
                            break
                    if print_line:
                        break
            if stat_line=='implicit':
                print_line=True
            if stat_line=='definition':
                tt=str(lines[iline_routines[iroutine]+k+1]).strip()
                tt=tt.lower()
                arg_list.append("parameter")
                for l in range(len(arg_list)):
                    ttt=r"\b"+arg_list[l]+r"\b"
                    p=re.compile(ttt)
                    if bool(p.search(tt)):
                        print_line=True
                        break
            if print_line:
                fout.write("%s\n" % lines[iline_routines[iroutine]+k+1])
        procedure,subroutine_name=get_subroutine_name(routines[iroutine])
        fout.write("end %s %s\n" % (procedure,subroutine_name))
#*****************************************************************************************
str1="Generates src/interface_mod.F90."
parser=argparse.ArgumentParser(description=str1)
str2="if present, it checks routines for using the module mod_interface."
parser.add_argument("-use",action='store_false',help=str2)
args=parser.parse_args()
use_modinterface=not args.use

fout=open('src/interface_mod.F90',"w")
fout.write("!"+"*"*89+"\n")
fout.write("module mod_interface\n")
fout.write("    implicit none\n")
fout.write("interface\n")
files=get_files()
for file in files:
    allcontent=commands.getoutput('cat '+file)
    lines=allcontent.splitlines()
    routines=[]
    iline_routines=[]
    iline=-1
    for line in lines:
        iline+=1
        str_line=str(line).lower()
        i_sub=str_line.find('subroutine')
        #i_fun=str_line.find('function') commented not to add functions to the interface
        i_fun=-1
        if i_sub>-1 or i_fun>-1:
            check_routine=is_routine(line)
            if check_routine:
                routines.append(line)
                iline_routines.append(iline)
    fout.write("! %s :\n" % file)
    if use_modinterface: routines_missing_use_interface(file,lines,routines,iline_routines)
    make_interface_routines(lines,routines,iline_routines,fout)

fout.write("end interface\n")
fout.write("end module mod_interface\n")
fout.write("!"+"*"*89+"\n")
fout.close()
#*****************************************************************************************
