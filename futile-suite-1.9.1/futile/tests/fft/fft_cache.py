#!/usr/bin/env python
# -*- coding: us-ascii -*-

import re
import os
import sys

#Values to write in green
start = "\033[0;32m"
end = "\033[m"

#List of tested values of ncache
list_cache = [ 0, 6, 12, 24, 50, 75, 100]
#list_cache = range(0,2000,10)

#Write fft_cache.gnuplot
fft_cache_gnuplot = '''set term wxt 1
set logscale x
set title "FFT Time versus the FFT size for different cache size"
set xlabel "FFT size array"
set ylabel "Time (s)"
plot "< awk '$1 == 0*1024'    fft_cache.dat" using 2:4 title "no cache"    w lp 1, \
     "< awk '$1 == 6*1024'    fft_cache.dat" using 2:4 title "6kB"    w lp 2, \
     "< awk '$1 == 12*1024'   fft_cache.dat" using 2:4 title "12kB"   w lp 3, \
     "< awk '$1 == 24*1024'   fft_cache.dat" using 2:4 title "24kB"   w lp 4, \
     "< awk '$1 == 50*1024'   fft_cache.dat" using 2:4 title "50kB"   w lp 5, \
     "< awk '$1 == 75*1024'   fft_cache.dat" using 2:4 title "75kB"   w lp 6, \
     "< awk '$1 == 100*1024'  fft_cache.dat" using 2:4 title "100kB"  w lp 7


set term wxt 2
set logscale x
set title "Performance/Best Performance versus the FFT size for different cache size"
set xlabel "FFT size array"
set ylabel "Perf/Best"
plot "fft_columns.dat" using 1:($2/$9) title "no cache"    w lp 1, \
     "fft_columns.dat" using 1:($3/$9) title "6k"    w lp 2, \
     "fft_columns.dat" using 1:($4/$9) title "12k"   w lp 3, \
     "fft_columns.dat" using 1:($5/$9) title "20k"   w lp 4, \
     "fft_columns.dat" using 1:($6/$9) title "50k"   w lp 5, \
     "fft_columns.dat" using 1:($7/$9) title "75k"   w lp 6, \
     "fft_columns.dat" using 1:($8/$9) title "100k"  w lp 7

#plot "fft_columns.dat" using 1:($2/$9) title "0k"    w lp 1, \
#     "fft_columns.dat" using 1:($8/$9) title "100k"  w lp 7

set term wxt 3
unset logscale x
set title "Total time for FFT from 50 to 500 versus the cache size"
set xlabel "Cache size (kB)
set ylabel "Total time (s)"
plot "fft_perf.dat" using ($1/1024):2 notitle w lp 
'''
open("fft_cache.gnuplot","w").write(fft_cache_gnuplot)

#Build new file
dico = dict()
caches = set()
lengthes = set()
for line in open("fft_cache.dat").readlines():
    line = line.split()
    if len(line) > 2:
        cache = int(line[0])
        caches.add(cache)
        length = int(line[1])
        lengthes.add(length)
        if dico.get(length) is None:
            dico[length] = dict()
        dico[length][cache] = float(line[-1])

caches = list(caches)
caches.sort()
lengthes = list(lengthes)
lengthes.sort()

fd = open("fft_columns.dat","w")
fd.write("#n1")
for c in caches:
    fd.write(" %d" % c)
fd.write("\n")

performances = dict()
for c in caches:
    performances[c] = 0.

for l in lengthes:
    fd.write("%d " % l)
    best = 1.e10
    for c in caches:
        perf = dico[l].get(c)
        if perf is not None:
            perf = dico[l][c]
            fd.write("%f20.15 " % perf)
            if perf < best:
                best = perf
            if l >= 50 and l <= 500:
                performances[c] += perf
        else:
            fd.write(" * ")
    fd.write("%f20.15\n" % best)
fd.close()

best = caches[0]
fd = open("fft_perf.dat","w")
fd.write("#ncache performances (50 <= x <+ 500) time\n")
for c in caches:
    perf = performances[c] 
    if performances[c] < performances[best]:
        best = c
    fd.write("%s %s\n" % (c,perf))
fd.close()

print ("Use fft_cache.gnuplot to display the results")
print ("The best ncache between tested values for 50 <= n1 <= 500 is ncache=",str(best))
print (start+"Test succeeded"+end)

