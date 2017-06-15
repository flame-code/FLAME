#!/usr/bin/env python
import os
import sys
from numpy import *
import numpy as np
import Gnuplot
import matplotlib.pyplot as plt
from pylab import *
#**********************************************read file: phonopy-FHI-aims-band_structure.dat
if len(sys.argv) < 2:
    print "usage: bandstr.py input_filename"
    exit()
else:
    filename = sys.argv[1]
f=open(filename,'r')
start=[]
S1=[]
S2=[]
end=[]
E1=[]
E2=[]
iline=0
for line in f.readlines():
        path=line.strip()
        if path.startswith("# Start"):
            start = line.split()   
            #print line.split()
            #print start
            s1 = start[6]
            if len(s1)>3:
                s1 = "{/Symbol G}"
            #print s1
            S1.append(s1)
            s2 = start[18]
            #print s2
            S2.append(s2)
        if path.startswith("#   End"):
            end = line.split()
            #print line.split()
            e1 =  end[6]
            if len(e1)>3:
                e1 = "{/Symbol G}"
            #print e1
            E1.append(e1)
            e2 = end[18]
            #print e2
            E2.append(e2)
        if path.startswith("# number of"):
            nbranch = line.split()
            #print nbranch[6]
            nb0 = nbranch[6]
            nb1 = int(nb0)
            nb2 = nb1 + 1 
            nb = str(nb2)
            #print nb
            #print line.split()
        # print line.rstrip("=")
        iline+=1
print S1
print S2
print E1
print E2
print nb
f.closed
#**************************************************plot by gnuplot
g = Gnuplot.Gnuplot()
g.xlabel('Wave number')
g.ylabel('Frequency(cm^{-1})')
g('set encoding utf8')
g('set key inside center bottom vertical noreverse enhanced autotitle nobox')
g('set key font "Times-New-Roman, 12"')
g('set key spacing 0.9')
g('set tics font "Times-New-Roman, 12"')
g('unset xtics')
g('set ytics 0,100')
g('set mytics 2')
g('set xrange [0:'+E2[-1]+'] noreverse nowriteback')
#g('set yrange [*:*] noreverse nowriteback')
g('set yrange [:700] writeback')

#For plotting vertical lines to show zone boundaries:
g('set arrow from '+S2[0]+' ,graph(0,0) to ' +S2[0]+' ,graph (1,1) nohead')
for s in E2:
    g('set arrow from '+s+' ,graph(0,0) to ' +s+' ,graph (1,1) nohead')
#g('set xtics("'+S1[0]+'"'+S2[0]+')')
l = []
for j in xrange(len(E1)):
    s = '"'+E1[j]+'" '+E2[j]+''
    l.append(s)
    #print l
ss = ', '.join(l)
#g('set xtics("{/Symbol G}"  0,'+ss+')')
g('set xtics("'+S1[0]+'"'+S2[0]+','+ss+')')
#**********************************************
#g('set terminal pdfcairo enhanced color solid font "Times-New-Roman,12"')
#g('set output "phband.pdf"')
g('set terminal postscript eps enhanced color font "Times-Roman,14"')
g('set output "phband.eps"')
#g('set terminal svg enhanced font "Times-Roman,12"')
#g('set output "phband.svg"')
g('plot for [col=2:'+nb+'] "phonopy-FHI-aims-band_structure.dat" using 1:col with lines notitle lc rgb "#0072bd" lw 2')
#g('set yrange restore')
g('reset')
g('pause -1')
#10E9E1 #1E90FF
#os.system('epstopdf phband.eps')
#***********************************************plot by matplotlib
