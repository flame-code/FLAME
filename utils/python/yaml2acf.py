#!/usr/bin/env python
import sys
import atoms
from acf import *
from io_yaml import *

if len(sys.argv) < 3:
    print "usage: yaml2acf.py input_filename output_filename"
    exit()
else:
    fn_inp = sys.argv[1]
    fn_out = sys.argv[2]

atoms_all=read_yaml(fn_inp)
acf_write(atoms_all,fn_out,"yaml2acf")

