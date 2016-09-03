#!/usr/bin/env python
import sys
import atoms
from acf import *
from ascii import *

if len(sys.argv) < 5:
    print "usage: acf2ascii.py -label labelpatt input_filename output_filename"
    exit()
else:
    if sys.argv[1]=='-label':
        labelpatt=sys.argv[2]
    else:
        print "ERROR: wrong command option: %s" % sys.argv[1]
        exit()
    fninp = sys.argv[3]
    fnout = sys.argv[4]

print 'Reading configurations from file %s ... ' % fninp, 
atoms_all=acf_read(fninp)
print 'done.'

print 'Writing configurations to file %s ... ' % fnout, 
acf_write(atoms_all,fnout,labelpatt=labelpatt)
print 'done.'
