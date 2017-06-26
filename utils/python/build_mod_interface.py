#!/usr/bin/env python
import commands
import string
import re
#*****************************************************************************************
def changeme(routine):
    #"This changes a passed list into this function"
    #mylist.append([1,2,3,4]);
    #print "Values inside the function: ", mylist
    #if 'minimahopping' in routine: print routine
    tt=str(routine).strip()
    tt=tt.lower()
    #if 'minimahopping' in routine: print routine
    i_sub=tt.find('subroutine')
    i_fun=tt.find('function')
    i_end=tt.find('end')
    i_comment=tt.find('!')
    i_stop=tt.find('stop')
    i_write=tt.find('write')
    #if 'minimahopping' in routine: print i_sub
    #if i_end>i_sub: print i_sub,i_end,i_comment
    check_routine=False
    if (i_sub==0 or i_fun==0): check_routine=True
    #after adding the line above, the next if statements are virtually useless
    #if (i_end>-1 and i_end<i_sub): check_routine=False
    #if (i_comment>-1 and i_comment<i_sub): check_routine=False
    #if (i_stop>-1 and i_stop<i_sub): check_routine=False
    #if (i_write>-1 and i_write<i_sub): check_routine=False
    #if (i_end>-1 and i_end<i_fun): check_routine=False
    #if (i_comment>-1 and i_comment<i_fun): check_routine=False
    #if (i_stop>-1 and i_stop<i_fun): check_routine=False
    #if (i_write>-1 and i_write<i_fun): check_routine=False
    return check_routine
#*****************************************************************************************
def check_var_define(line):
    tt=str(line).strip()
    tt=tt.lower()
    #print tt
    #i_sub=tt.find('subroutine')
    i_use=tt.find('use')
    i_mod_interface=tt.find('mod_interface')
    i_comment=tt.find('!')
    #tt=' '.join(tt.translate(None, string.whitespace[:5]).split())
    i_implicit=tt.find('implicit none')
    i_integer=tt.find('integer')
    i_real=tt.find('real')
    i_external=tt.find('external')
    i_character=tt.find('character')
    i_logical=tt.find('logical')
    i_type=tt.find('type')
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
    if (i_integer==0 or i_real==0 or i_character==0 or i_logical==0 or i_type==0 or i_external==0):
        if not stat_line=='unknown':
            print 'ERROR: definition or ',stat_line
        else:
            stat_line='definition'
    if stat_line=='unknown': stat_line='command'
    #print i_use,i_comment,i_implicit,i_integer,i_real,i_character,i_logical,i_type
    return stat_line
#*****************************************************************************************
def get_arguments(line):
    tt=str(line).strip()
    tt=tt.lower()
    #print len(tt.split())
    i_parenthesis=tt.find('(')
    if i_parenthesis>0:
        #i_comment=tt.find('!')
        tt=tt.replace(" ","")
        if '!' in tt:
            tt=tt.split("!")[0]
        #tt=re.sub(r'\s', '', tt)
        #print tt
        i_fun=tt.find('function')
        i_res=tt.find('result')
        if (i_fun>-1 and i_res>-1):
            tt=tt.replace("(",",",1)
            tt=tt.replace(")",",",1)
        tt=tt.replace("(",",",1)
        tt=tt.replace(")","",1)
        #ttl=tt.split()[1].split(",")
        #print tt.split()[1]
        #tt=tt.join([])
        #tt=tt.split()[1]
        #tt+=r"\n"
        ttl=tt.split(",")
        #ttl.append("")
        ttl.pop(0)
        if (i_fun>-1 and i_res>-1):
            ttl.remove('result')
        #print ttl
        #ttl=tt.split()[1].split("(")[1].split(")")[0].split(",")
    else:
        ttl=[]
    #print len(ttl)
    #print ttl
    return ttl
    #ttp=str(tt.split()[1])
    #print ttp.split("(")
#*****************************************************************************************
def get_var_list(line):
    tt=str(line).strip()
    tt=tt.lower()
    #print len(tt.split())
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
    #print tt
    if not tt=="":
        ttl=tt.split(",")
    else:
        ttl=[]
    #i_parenthesis=tt.find('(')
    #if i_parenthesis>0:
    #    ttl=tt.split()[1].split("(")[1].split(")")[0].split(",")
    #else:
    #    ttl=[]
    #print len(ttl)
    #print ttl,len(ttl)
    return ttl
    #ttp=str(tt.split()[1])
    #print ttp.split("(")
#*****************************************************************************************
def get_subroutine_name(line):
    tt=str(line).strip()
    tt=tt.lower()
    #print len(tt.split())
    i_parenthesis=tt.find('(')
    if i_parenthesis>0:
        ttl=tt.split()[1].split("(")[0]
    else:
        ttl=tt.split()[1]
    procedure=tt.split()[0]
    return procedure,ttl
#*****************************************************************************************
#def check_arguments(n,arg_list,line):
#    print lines[i_routines[j]+k+1]
#    tt=str(lines[i_routines[j]+k+1]).strip()
#    tt=tt.lower()
#    #print tt.find(arg_list[0])
#    ttt=r"\b"+arg_list[0]+r"\b"
#    p=re.compile(ttt)
#    print bool(p.search(tt))
#*****************************************************************************************
#commands.getstatusoutput('cat /bin/junk')
str_files=commands.getoutput('find . -not -path "./futile/*" -not -path "./minhocao/*" -iname "*.F90" |sort')
#str_files=commands.getoutput('find . -iname "cg.F90"')
#str_files=commands.getoutput('find . -iname "minhopp.F90"')
#str_files=commands.getoutput('find . -iname "energyandforces_SIESTA.F90"')
#print str_files
files=str_files.splitlines()
#print files.index('./modules/interface_mod.F90')
#print files[40]
#print "Number of fortran files: %5d" % len(files)
exclude_files=['./modules/interface_mod.F90',
               './wrappers/energyandforces_siesta.f90',
               './wrappers/energyandforces_openmx.f90',
               './src/MPIfake.F90',
               './src/task_netsock.F90',
               './modules/constants_minhocao_mod.F90',
               './modules/minhocao_mod.F90',
               './src/potential_LAMMPS_interface.F90',
               './src/minhocao.F90',
               './src/POSCAR2ascii.F90',
               './src/potential_abinit.F90',
               './src/potential_alborz.F90',
               './src/ascii2POSCAR.F90',
               './src/atoms_minhocao.F90',
               './src/optimizer_bfgs_minhocao.F90',
               './src/binaries.F90',
               './src/potential_BLJ_minhocao.F90',
               './src/cell_utils.F90',
               './src/potential_confinement.F90',
               './src/potential_corerepulsion.F90',
               './src/potential_CP2K.F90',
               './src/potential_DFTB_minhocao.F90',
               './src/potential_EDIP.F90',
               './src/envelope.F90',
               './src/potential_PWSCF.F90',
               './src/PWSCF_restruct.F90',
               './src/quaternions.F90',
               './src/expand_poslows.F90',
               './src/find_symmetry.F90',
               './src/fingerprint_atorb.F90',
               './src/fingerprint_BCM.F90',
               './src/fingerprint_gaussmol.F90',
               './src/fingerprint_GOM.F90',
               './src/fingerprint_MOLGOM.F90',
               './src/fingerprint_oganov.F90',
               './src/fingerprint_oganov_cont.F90',
               './src/fingerprint_XYZ2SM.F90',
               './src/convex_hull.F90',
               './src/potential_main_minhocao.F90',
               './src/io_ascii.F90',
               './src/io_vasp_minhocao.F90',
               './src/potential_IPI.F90',
               './src/potential_LAMMPS.F90',
               './src/potential_LenoskyMEAM.F90',
               './src/lenosky_tb/',
               './src/potential_LenoskyTB_minhocao.F90',
               './src/potential_LenoskyTB_LJ_minhocao.F90',
               './src/optimizer_nlbfgs_minhocao.F90',
               './src/potential_LJ_voids.F90',
               './src/dynamics_md_fixlat.F90',
               './src/potential_MLJ.F90',
               './src/optimizer_sqnm_minhocao_module.F90',
               './src/potential_MOPAC.F90',
               './src/potential_MSOCK.F90',
               './src/msock_slave_template.F90',
               './src/cell_niggli.F90',
               './src/parser_minhocao.F90',
               './src/parser_core_minhocao.F90',
               './src/optimizer_bfgs_qe.F90',
               './src/recompute_kpt.F90',
               './src/optimizer_sd_minhocao.F90',
               './src/potential_SIESTA_minhocao.F90',
               './src/spglib_int.F90',
               './src/spher_harm_mathematica.F90',
               './src/optimizer_sqnm_minhocao.F90',
               './src/optimizer_subs_minhocao.F90',
               './src/ternaries.F90',
               './src/potential_TERSOFF.F90',
               './src/potential_TINKER.F90',
               './src/potential_VASP_minhocao.F90',
               './src/vasp_recompute_cell.F90',
               './src/vasp_recompute_kpt.F90',
               './src/vasp_recompute_kpt_odd.F90']

for f in exclude_files:
    if f in files: files.remove(f)

#print "Number of fortran files: %5d" % len(files)
#print files[40]
str_interfaced=commands.getoutput('grep -i subroutine ./modules/interface_mod.F90')
routine_interfaced=str_interfaced.splitlines()
#all_routines=[]
#print "\n"

f=open('modules/interface_mod.F90',"w")
f.write("!***************************************************************************************************\n")
f.write("module mod_interface\n")
f.write("    implicit none\n")
f.write("interface\n")

n=0
for i in files:
    n+=1
    #test=True
    #if not n==55: continue
    allcontent=commands.getoutput('cat '+i)
    lines=allcontent.splitlines()
    routines=[]
    i_routines=[]
    for line in lines:
        #print line
        if 'subroutine' in line or 'function' in line:
            check_routine=changeme(line)
            if check_routine:
                routines.append(line)
                i_routines.append(lines.index(line))
    #str_routines=commands.getoutput('grep -i subroutine '+i)
    #routines=str_routines.splitlines()
    #string.find(routines[0],'subroutine')
    #tt="".join(routines[0])
    #print
    f.write("! %s :\n" % i)
    #for j in range(2):
    for j in range(len(routines)):
        #check_routine=changeme(routines[j]);
        #if not check_routine: continue
        #check=routines[j] in routine_interfaced
        #print check
        f.write("%s\n" % routines[j])
        arg_list=get_arguments(routines[j])
        if j==len(routines)-1:
            nline_routine=len(lines)-i_routines[j]-1
        else:
            nline_routine=i_routines[j+1]-i_routines[j]-1
        definitions=[]
        for k in range(nline_routine):
            stat_line=check_var_define(lines[i_routines[j]+k+1])
            if stat_line=='definition':
                definitions.append(lines[i_routines[j]+k+1])
                #print lines[i_routines[j]+k+1]
        stat_implicit_none=False
        #for k in range(10):
        for k in range(nline_routine):
            stat_line=check_var_define(lines[i_routines[j]+k+1])
            print_line=False
            if stat_line=='use':
                var_list=get_var_list(lines[i_routines[j]+k+1])
                #print tt.find(var_list[0])
                for m in range(len(definitions)):
                    tt=str(definitions[m]).strip()
                    tt=tt.lower()
                    for l in range(len(var_list)):
                        ttt=r"\b"+var_list[l]+r"\b"
                        #print ttt
                        p=re.compile(ttt)
                        #print bool(p.search(tt))
                        if bool(p.search(tt)):
                            print_line=True
                            break
                    if print_line:
                        break
            if stat_line=='implicit':
                print_line=True
            #if stat_line=='comment':
            #    print_line=True
            if stat_line=='definition':
                tt=str(lines[i_routines[j]+k+1]).strip()
                tt=tt.lower()
                #print tt.find(arg_list[0])
                #print arg_list[0]
                for l in range(len(arg_list)):
                    ttt=r"\b"+arg_list[l]+r"\b"
                    #ttt=arg_list[l]
                    #print ttt
                    p=re.compile(ttt)
                    #print bool(p.search(tt))
                    #if bool(re.search(ttt,tt)):
                    if bool(p.search(tt)):
                        print_line=True
                        break
            if print_line:
                f.write("%s\n" % lines[i_routines[j]+k+1])
        procedure,subroutine_name=get_subroutine_name(routines[j])
        #print procedure
        f.write("end %s %s\n" % (procedure,subroutine_name))
                #arg_list=check_arguments(len(arg_list),arg_list,lines[i_routines[j]+k+1])
                #if not stat_implicit_none:
                #    tt=str(lines[i_routines[j]+k+1]).strip()
                #    tt=tt.lower()
                #    if tt=='implicit none':
                #        stat_implicit_none=True
                #    else:
                #        print 'ERROR: implicit none is missing in '
    #print "++++++++++++++++++++++++++++++++++++++++++++++++++++"
    #print routines[0].find('subroutine')
    #print "\n",i,":\n",str_routines
    #if n==3: break

f.write("end interface\n")
f.write("end module mod_interface\n")
f.write("!***************************************************************************************************\n")
f.close()
#*****************************************************************************************
