#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------

#> @file
## Define a program to run tests to be called by the Makefiles
## @author
##    Copyright (C) 2016-2016 BigDFT group
##    This file is distributed under the terms of the
##    GNU General Public License, see ~/COPYING file
##    or http://www.gnu.org/copyleft/gpl.txt .
##    For the list of contributors, see ~/AUTHORS
#we should think on reducing the compatibility of the argparse to optparse
import UniParse

parser=UniParse.UniParser(description='Regression checker',method='argparse')

parser.option('-f', '--fldiff', dest='fldiff', default="/dev/null",
              help="script file performing fldiff (default: /dev/null)",
              metavar='FILE')
parser.option('-t', '--tols', dest='tols', default="/dev/null",
              help="yaml file containing the tolerances for each run",
              metavar='FILE')
parser.option('-s', '--srcdir', dest='srcdir', default="/",
              help="yaml file containing the tolerances for each run",
              metavar='FILE')
#
args = parser.args()
#print 'Arguments'
print args.__dict__

import os

instr=os.environ.get('F_REGTEST_INSTRUCTIONS')
if instr is None: quit()

import yaml
d_instr=yaml.load(instr)

def run_test(runs):
    for r in runs:
        print 'executing: ',r
        os.system(r)

def get_time(file):
    if os.path.isfile(file):
        return os.path.getmtime(file)
    else:
        return 0.0


fldiff=args.fldiff
tols=args.tols
base='python '+fldiff+' -t '+tols

for test in d_instr:
    label=test.keys()[0]
    specs=test.values()[0]
    binary=specs.get('binary',label)
    output=specs.get('output',label+'.out.yaml')
    report=label+'.report.yaml'
    ref=specs.get('reference',label+'.ref.yaml')
    if not os.path.isfile(ref):
        ref=os.path.join(args.srcdir,ref)
    dorun=get_time(binary) > get_time(output) or get_time(output) == 0.0
    if dorun: run_test(specs['runs'])
    os.system(base+' --label '+label+' -r '+ref+' -d '+output+' --output '+report)
        

