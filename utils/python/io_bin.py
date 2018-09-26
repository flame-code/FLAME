from array import array
from atomic_data import *
from atoms import *
#*****************************************************************************************
#written similar to chmod coding
def get_ibemoved(bm):
    bml=list(bm)
    if bml[0]=='T':
        i1=1
    else:
        i1=0
    if bml[1]=='T':
        i2=1
    else:
        i2=0
    if bml[2]=='T':
        i3=1
    else:
        i3=0
    ibm=i1+i2*2+i3*4
    return ibm
#*****************************************************************************************
def bin_write_v1(atoms_all,filename):
    wa=[]
    ver=1.0
    wa.append(ver)
    nconf=len(atoms_all)
    wa.append(float(nconf))
    for iconf in range(nconf):
        if atoms_all[iconf].cell_present:
            wa.append(1.0)
        else:
            wa.append(0.0)
        if atoms_all[iconf].epot_present:
            wa.append(1.0)
        else:
            wa.append(0.0)
        if atoms_all[iconf].fat_present:
            wa.append(1.0)
        else:
            wa.append(0.0)
        if atoms_all[iconf].qtot_present:
            wa.append(1.0)
        else:
            wa.append(0.0)
        if atoms_all[iconf].bemoved_present:
            wa.append(1.0)
        else:
            wa.append(0.0)
        if atoms_all[iconf].vat_present:
            wa.append(1.0)
        else:
            wa.append(0.0)
        if atoms_all[iconf].boundcond=='free':
            ibc=0
        elif atoms_all[iconf].boundcond=='wire':
            ibc=1
        elif atoms_all[iconf].boundcond=='slab':
            ibc=2
        elif atoms_all[iconf].boundcond=='bulk':
            ibc=3
        else:
            print "ERROR: unknown boundary conditions in bin_write_v1, boundcond=" % \
                    atoms_all[iconf].boundcond
            exit()
        wa.append(float(ibc))
        wa.append(float(atoms_all[iconf].nat))
        for iat in range(atoms_all[iconf].nat):
            wa.append(float(atoms_all[iconf].rat[iat][0]))
            wa.append(float(atoms_all[iconf].rat[iat][1]))
            wa.append(float(atoms_all[iconf].rat[iat][2]))
        for iat in range(atoms_all[iconf].nat):
            isat=get_atomic_number(atoms_all[iconf].sat[iat])
            wa.append(float(isat))
        if atoms_all[iconf].cell_present:
            wa.append(float(atoms_all[iconf].cellvec[0][0]))
            wa.append(float(atoms_all[iconf].cellvec[0][1]))
            wa.append(float(atoms_all[iconf].cellvec[0][2]))
            wa.append(float(atoms_all[iconf].cellvec[1][0]))
            wa.append(float(atoms_all[iconf].cellvec[1][1]))
            wa.append(float(atoms_all[iconf].cellvec[1][2]))
            wa.append(float(atoms_all[iconf].cellvec[2][0]))
            wa.append(float(atoms_all[iconf].cellvec[2][1]))
            wa.append(float(atoms_all[iconf].cellvec[2][2]))
        if atoms_all[iconf].epot_present:
            wa.append(float(atoms_all[iconf].epot))
        if atoms_all[iconf].fat_present:
            for iat in range(atoms_all[iconf].nat):
                wa.append(float(atoms_all[iconf].fat[iat][0]))
                wa.append(float(atoms_all[iconf].fat[iat][1]))
                wa.append(float(atoms_all[iconf].fat[iat][2]))
        if atoms_all[iconf].qtot_present:
            wa.append(float(atoms_all[iconf].qtot))
        if atoms_all[iconf].bemoved_present:
            for iat in range(atoms_all[iconf].nat):
                ibm=get_ibemoved((atoms_all[iconf].bemoved[iat]))
                wa.append(float(ibm))
        if atoms_all[iconf].vat_present:
            for iat in range(atoms_all[iconf].nat):
                wa.append(float(atoms_all[iconf].vat[iat][0]))
                wa.append(float(atoms_all[iconf].vat[iat][1]))
                wa.append(float(atoms_all[iconf].vat[iat][2]))
    output_file = open(filename,'wb')
    float_array = array('d', wa)
    float_array.tofile(output_file)
    output_file.close()
    del wa
    del float_array
#*****************************************************************************************
def bin_write(atoms_all,filename,iver):
    if iver==1:
        bin_write_v1(atoms_all,filename)
    else:
        print "ERROR: writing to binary ready only for version=1"
#*****************************************************************************************
def bin_read_v1(filename,wa):
    nconf=int(wa[1])
    bm_list=['FFF','TFF','FTF','TTF','FFT','TFT','FTT','TTT']
    atoms_all=[]
    #print 'nconf ',nconf
    iwa=1
    for iconf in range(nconf):
        atoms_all.append(Atoms())
        iwa+=1
        if int(wa[iwa])==1: atoms_all[-1].cell_present=True
        iwa+=1
        if int(wa[iwa])==1: atoms_all[-1].epot_present=True
        iwa+=1
        if int(wa[iwa])==1: atoms_all[-1].fat_present=True
        iwa+=1
        if int(wa[iwa])==1: atoms_all[-1].qtot_present=True
        iwa+=1
        if int(wa[iwa])==1: atoms_all[-1].bemoved_present=True
        iwa+=1
        if int(wa[iwa])==1: atoms_all[-1].vat_present=True
        iwa+=1
        ibc=int(wa[iwa])
        if ibc==0:
            atoms_all[iconf].boundcond='free'
        elif ibc==1:
            atoms_all[iconf].boundcond='wire'
        elif ibc==2:
            atoms_all[iconf].boundcond='slab'
        elif ibc==3:
            atoms_all[iconf].boundcond='bulk'
        else:
            print "ERROR: invalid coding for boundcond in bin_read_v1, ibc=" % ibc
            exit()
        iwa+=1
        atoms_all[-1].nat=int(wa[iwa])
        for iat in range(atoms_all[-1].nat):
            atoms_all[-1].rat.append([])
            iwa+=1
            atoms_all[-1].rat[-1].append(wa[iwa])
            iwa+=1
            atoms_all[-1].rat[-1].append(wa[iwa])
            iwa+=1
            atoms_all[-1].rat[-1].append(wa[iwa])
        for iat in range(atoms_all[-1].nat):
            iwa+=1
            isat=int(wa[iwa])
            sat=get_atomic_symbol(isat)
            atoms_all[-1].sat.append(sat)
        if atoms_all[-1].cell_present:
            iwa+=1
            atoms_all[-1].cellvec[0][0]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[0][1]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[0][2]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[1][0]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[1][1]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[1][2]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[2][0]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[2][1]=wa[iwa]
            iwa+=1
            atoms_all[-1].cellvec[2][2]=wa[iwa]
        if atoms_all[-1].epot_present:
            iwa+=1
            atoms_all[-1].epot=wa[iwa]
        if atoms_all[-1].fat_present:
            for iat in range(atoms_all[-1].nat):
                atoms_all[-1].fat.append([])
                iwa+=1
                atoms_all[-1].fat[-1].append(wa[iwa])
                iwa+=1
                atoms_all[-1].fat[-1].append(wa[iwa])
                iwa+=1
                atoms_all[-1].fat[-1].append(wa[iwa])
        if atoms_all[-1].qtot_present:
            iwa+=1
            atoms_all[-1].qtot=wa[iwa]
        if atoms_all[-1].bemoved_present:
            for iat in range(atoms_all[-1].nat):
                iwa+=1
                ibm=int(wa[iwa])
                if ibm<0 or ibm>7:
                    print "ERROR: invalid coding for ibm in bin_read_v1 , ibm=%d" % ibm
                    exit()
                atoms_all[-1].bemoved.append(bm_list[ibm])
        if atoms_all[-1].vat_present:
            for iat in range(atoms_all[-1].nat):
                atoms_all[-1].vat.append([])
                iwa+=1
                atoms_all[-1].vat[-1].append(wa[iwa])
                iwa+=1
                atoms_all[-1].vat[-1].append(wa[iwa])
                iwa+=1
                atoms_all[-1].vat[-1].append(wa[iwa])
    #print "----------------------------------------------------------------------"
    #print type(atoms_all[0])
    #print "----------------------------------------------------------------------"
    #print atoms_all[1]
    return atoms_all
#*****************************************************************************************
def bin_read(filename):
    input_file = open(filename,'r')
    float_array = array('d')
    float_array.fromstring(input_file.read())
    wa=list(float_array)
    #print wa
    del float_array
    iver=int(wa[0])
    if iver==1:
        atoms_all=bin_read_v1(filename,wa)
    else:
        print "ERROR: writing to binary ready only for version=1"
    del wa
    return atoms_all
#*****************************************************************************************
