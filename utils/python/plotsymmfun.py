#!/usr/bin/env python
import argparse
import math
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import MultipleLocator, FormatStrFormatter , AutoMinorLocator
from matplotlib import rc

rc('font',**{'family':'serif','serif':['Times']})
rc('text', usetex=True)
#----------------------------------------------------------------------------------------
def symmfun_typ2(eta,r,rs,rc):
    fcij = funccutoff(r,rc)
    v2=math.exp(-1.0*eta*((r-rs)**2))*fcij
    return v2
#---------------------------
def symmfun_typ5r(eta,r,rc):
    fcij = funccutoff(r,rc)
    v5r=math.exp(-1.0*eta*(3*r**2))*fcij**3
    return v5r
#---------------------------
def symmfun_typ5a(eta,zeta,lam,r):
    v5a=2.0**(1.0-zeta)*(lam*(math.cos(r*math.pi/180.0))+1.0)**zeta
    return v5a
#---------------------------
def funccutoff(r,rc):
    fcij=0.5*(math.cos(math.pi*r/rc)+1.0)
    return fcij
#---------------------------
def read_symf(filename):
    import yaml
    fd = open(filename, 'r')
    dict_list = yaml.load(fd)
    sym2=[]
    sym5=[]
    rcut=dict_list['main']['rcut']
    for i in range(len(dict_list['symfunc'])):
        if dict_list['symfunc'].keys()[i].startswith('g02'):
            gg2 = dict_list['symfunc'].keys()[i]
            inp = dict_list['symfunc'][gg2].split()
            sym2.append([float(inp[0]), float(inp[1]), float(rcut)])
        if dict_list['symfunc'].keys()[i].startswith('g05'):
            gg5 = dict_list['symfunc'].keys()[i]
            inp = dict_list['symfunc'][gg5].split()
            sym5.append([float(inp[0]), float(inp[1]), float(inp[2]), float(rcut)])
    return sym2, sym5
    fd.close()
#----------------------------------------------------------------------------------------
if __name__ == "__main__":

    str1 = "This script plots Behler's symmetry functions (only type 2 and 5)."
    parser = argparse.ArgumentParser(description=str1)
    parser.add_argument('fn_inp', action='store' ,type=str, help="Name of the input file in yaml format ?.ann.input.yaml")
    args=parser.parse_args()

    sym2, sym5 = read_symf(args.fn_inp)
    eta2=[]; rs2=[]; rc2=[]
    for i in range(len(sym2)):
        eta2.append(sym2[i][0])
        rs2.append(sym2[i][1])
        rc2.append(sym2[i][2])
    
    eta5=[]; zeta5=[]; lam5=[]; rc5=[]
    for i in range(len(sym5)):
        eta5.append(sym5[i][0])
        zeta5.append(sym5[i][1])
        lam5.append(sym5[i][2])
        rc5.append(sym5[i][3])
#----------------------------------------------------------------
    nval=400
    symf2 = []
    radius=max(rc2)
    dr=radius/(nval-1.0)
    #type 2: radial
    for i in range(len(sym2)):
        r1=0.0
        symf2r=[]
        symf2v=[]
        for j in range(nval):
            r1 = r1 + dr
            arr2 = symmfun_typ2(eta2[i],r1,rs2[i],rc2[i])
            symf2r.append(r1)
            symf2v.append(arr2)
        symf2.append([symf2r,symf2v])

    nval=400
    symf5rr = []
    radius=max(rc5)
    dr=radius/(nval-1.0)
    #type 5: radial
    for i in range(len(sym5)):
        r1=0.0
        symf5r=[]
        symf5v=[]
        for j in range(nval):
            r1 = r1 + dr
            arr2 = symmfun_typ5r(eta5[i],r1,rc5[i])
            symf5r.append(r1)
            symf5v.append(arr2)
        symf5rr.append([symf5r,symf5v])

    nang=360
    symf5aa = []
    r=0.0
    #type 5: angular 
    for i in range(len(sym5)):
        symf5t=[]
        symf5vv=[]
        for j in range(nang):
            r = r + 1.0
            arr5 = symmfun_typ5a(eta5[i],zeta5[i],lam5[i],r)
            symf5t.append(j)
            symf5vv.append(arr5)
        symf5aa.append([symf5t,symf5vv])
#plot type 5 angular part----------------------------------------------------------------
    fig5,ax5 = plt.subplots()
    for i in range(len(symf5aa)):
        plt.plot( symf5aa[i][0],  symf5aa[i][1] ,'-', color='#00BFFF', lw=1.5)
    majorLocator = MultipleLocator(30)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(15)
    plt.xticks(np.linspace(0,360,13, endpoint=True))
    plt.xlim(0.0, 360.0)
    plt.xlabel("Angle (degree)", size=20)
    ax5.tick_params(which='both', width=1.3, labelsize=14)
    #plt.ylabel("Symmetry function type 5",size=20)
    plt.ylabel("$G^5$ (angular)",size=20)
    ax5.xaxis.set_major_locator(majorLocator)
    ax5.xaxis.set_major_formatter(majorFormatter)
    ax5.xaxis.set_minor_locator(minorLocator)
    fig5.savefig("typ5.pdf")
#plot type 2 radial----------------------------------------------------------------
    fig2,ax2 = plt.subplots()
    for i in range(len(symf2)):
        plt.plot(symf2[i][0], symf2[i][1],'-', color='#00BFFF', lw=1.5)
    majorLocator = MultipleLocator(2)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(1)
    plt.xticks(np.linspace(0,13,1, endpoint=True))
    plt.xlim(0.0, 13.0)
    plt.xlabel("r (bohr)", size=20)
    ax2.tick_params(which='both', width=1.3, labelsize=14)
    plt.ylabel("$G^2$ (radial)",size=20)
    ax2.xaxis.set_major_locator(majorLocator)
    ax2.xaxis.set_major_formatter(majorFormatter)
    ax2.xaxis.set_minor_locator(minorLocator)
    fig2.savefig("typ2.pdf")
#plot type 5 radial part----------------------------------------------------------------
    fig52,ax52 = plt.subplots()
    for i in range(len(symf5rr)):
        plt.plot(symf5rr[i][0], symf5rr[i][1],'-', color='#00BFFF', lw=1.5)
    majorLocator = MultipleLocator(2)
    majorFormatter = FormatStrFormatter('%d')
    minorLocator = MultipleLocator(1)
    plt.xticks(np.linspace(0,13,1, endpoint=True))
    plt.xlim(0.0, 13.0)
    plt.xlabel("r (bohr)", size=20)
    ax52.tick_params(which='both', width=1.3, labelsize=14)
    plt.ylabel("$G^5$ (radial)",size=20)
    ax52.xaxis.set_major_locator(majorLocator)
    ax52.xaxis.set_major_formatter(majorFormatter)
    ax52.xaxis.set_minor_locator(minorLocator)
    fig52.savefig("typ5r.pdf")
