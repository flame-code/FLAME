#!/usr/bin/env python
import time
import random
import math
import sys
import os

random.seed()
#----------------------------------------------------------------------------------------------------------
#This programe generate random crystall structure like that described in 10.1103/PhysRevLett.116.075503
#a, b, c and gama are random and alpha=beta=90
#This program requaiers an "input" file like following:
#Ti O            #stypat: symbol of all atom types
#15 15             #nat_kinds: number of atoms of each type per cell: nkinds numbers in this line
#30              #nat_cell: number of atoms per cell
#1.5             # minimum distance between atoms
#----------------------------------------------------------------------------------------------------------
#reading input file
typ1=str(os.popen("cat input |awk 'NR<2 {print$1}'").read())
typ2=str(os.popen("cat input |awk 'NR<2 {print$2}'").read())
ntyp1=int(os.popen("cat input |awk '1<NR && NR<3 {print$1}'").read())
ntyp2=int(os.popen("cat input |awk '1<NR && NR<3 {print$2}'").read())
nat=    int(os.popen("cat input |awk '2<NR && NR<4 {print$1}'").read())
md=   float(os.popen("cat input |awk '3<NR && NR<5 {print$1}'").read())
typ1=typ1.split()[0]
typ2=typ2.split()[0]
#----------------------------------------------------------------------------------------------------------
nan=ntyp2; ncat=ntyp1
pi=4.0*math.atan(1.0)
vc=1+0.2*(nat/30+1)
while True:
    vmin=vc*4*(nat/3)*(2.0)**3
    vmax=vc*4*(nat/3)*(2.3)**3
    a=(random.uniform(vmin,vmax))**(1.0/3)#generating a cellvector 
    b=(random.uniform(vmin,vmax))**(1.0/3)#generating b cellvector 
    gama=pi/(random.uniform(1.5,7))#generating gama angel
    #c=random.uniform(vmin/(a*b*math.sin(gama)),vmax/(a*b*math.sin(gama)))
    #c=vmax/(a*b*math.sin(gama))
    c=vmin/(a*b*math.sin(gama))#generating c cellvector 
    cell_vec=[[0,0,0],[0,0,0],[0,0,0]]
    xan=[0 for i in range(nan)];yan=[0 for i in range(nan)];zan=[0 for i in range(nan)]#empty list for anions
    xcat=[0 for i in range(ncat)];ycat=[0 for i in range(ncat)];zcat=[0 for i in range(ncat)]#empty list for cations
    #----------------------------------------------------------------------------------------------------------
    cell_vec[0][0]=a;                     cell_vec[0][1] =00;                    cell_vec[0][2]=0 #a
    cell_vec[1][0]=b*math.cos(gama);      cell_vec[1][1] =b*math.sin(gama);      cell_vec[1][2]=0 #b
    cell_vec[2][0]=00;                    cell_vec[2][1] =00;                    cell_vec[2][2]=c #c
    
    ax=cell_vec[0][0];      ay=cell_vec[0][1];      az=cell_vec[0][2]
    bx=cell_vec[1][0];      by=cell_vec[1][1];      bz=cell_vec[1][2]
    cx=cell_vec[2][0];      cy=cell_vec[2][1];      cz=cell_vec[2][2]
    #----------------------------------------------------------------------------------------------------------
    #generating receprical lattice vectors g1 and g2 and g3
    v=ax*(by*cz-bz*cy)-ay*(bx*cz-bz*cx)+az*(bx*cy-by*cx)
    g1x=(2.0*pi/v)*(by*cz-cy*bz);  g2x=(2.0*pi/v)*(cy*az-ay*cz);  g3x=(2.0*pi/v)*(ay*bz-by*az);
    g1y=(2.0*pi/v)*(bz*cx-cz*bx);  g2y=(2.0*pi/v)*(cz*ax-az*cx);  g3y=(2.0*pi/v)*(az*bx-bz*ax);
    g1z=(2.0*pi/v)*(bx*cy-cx*by);  g2z=(2.0*pi/v)*(cx*ay-ax*cy);  g3z=(2.0*pi/v)*(ax*by-bx*ay);
    #----------------------------------------------------------------------------------------------------------
    #condition for g vector in order to fix distance between planes in range 1.8<d<2.1
    while True: 
        m1=random.uniform(0,10) 
        m2=random.uniform(0,10) 
        m3=random.uniform(0,10) 
        gx=m1*g1x+m2*g2x+m3*g3x; gy=m1*g1y+m2*g2y+m3*g3y; gz=m1*g1z+m2*g2z+m3*g3z
        ampg=(gx**2+gy**2+gz**2)**0.5
        if 1.8<pi/ampg<2.1: break
    #----------------------------------------------------------------------------------------------------------
    i,ccat,can=0,0,0
    x=[0 for v in range(nat*27)];y=[0 for v in range(nat*27)];z=[0 for v in range(nat*27)];ntry=0
    tv="False"
    while i < 27*nat:
        if ntry>400: vc=vc+0.2; tv="True";  break
        t="False"
        rand1=random.uniform(0,1)
        rand2=random.uniform(0,1)
        rand3=random.uniform(0,1)
        x[i]=rand1*cell_vec[0][0]+rand2*cell_vec[1][0]+rand3*cell_vec[2][0]
        y[i]=rand1*cell_vec[0][1]+rand2*cell_vec[1][1]+rand3*cell_vec[2][1]
        z[i]=rand1*cell_vec[0][2]+rand2*cell_vec[1][2]+rand3*cell_vec[2][2]
        #Image in a direction
        x[i+1]=x[i]+ax; y[i+1]=y[i]+ay; z[i+1]=z[i]+az
        x[i+2]=x[i]-ax; y[i+2]=y[i]-ay; z[i+2]=z[i]-az
        #Image in b direction
        for j in range(i+3,i+6):
            x[j]=x[j-3]+bx;  y[j]=y[j-3]+by;  z[j]=z[j-3]+bz
        for j in range(i+6,i+9):
            x[j]=x[j-6]-bx;  y[j]=y[j-6]-by;  z[j]=z[j-6]-bz
        #Image in c direction
        for j in range(i+9,i+18):
            x[j]=x[j-9]+cx;  y[j]=y[j-9]+cy;  z[j]=z[j-9]+cz
        for j in range(i+18,i+27):
            x[j]=x[j-18]-cx;  y[j]=y[j-18]-cy;  z[j]=z[j-18]-cz
    #-------------------------------------------------------------------------------------------------------
        for j1 in range(i):
            for k in range(i+1,j):
                dx=x[j1]-x[k]
                dy=y[j1]-y[k]
                dz=z[j1]-z[k]
                dd0=(dx**2+dy**2+dz**2)**0.5
                if (dd0<md): t="True"; break    
            if t=="True": break    
        if t=="True": ntry=ntry+1;continue;    
    
    
    
        rg=x[i]*gx+y[i]*gy+z[i]*gz
        pw=math.cos(rg)
        if pw > 0.95 and ccat<=(ncat-1):
            xcat[ccat]=x[i] 
            ycat[ccat]=y[i]
            zcat[ccat]=z[i]
            i=j+1
            ccat=ccat+1
        elif pw<-0.95 and can<=(nan-1):
            xan[can]=x[i] 
            yan[can]=y[i]
            zan[can]=z[i]
            i=j+1
            can=can+1
        else:continue
    if tv=="True": continue
    else: break
#----------------------------------------------------------------------------------------------------------
poscur=open("poscur.ascii","w")
poscur.write("%3d\n"%(nat))
poscur.write("%24.15E%24.15E%24.15E\n"%(cell_vec[0][0],cell_vec[1][0],cell_vec[1][1]))
poscur.write("%24.15E%24.15E%24.15E\n"%(cell_vec[2][0],cell_vec[2][1],cell_vec[2][2]))
for i in range(ncat):
        poscur.write("%24.15E%24.15E%24.15E%5s\n"%(xcat[i],ycat[i],zcat[i],typ1))
for i in range(nan):
        poscur.write("%24.15E%24.15E%24.15E%5s\n"%(xan[i],yan[i],zan[i],typ2))
