#!/usr/bin/env python
import sys
#import atoms
import yaml
from dict_diff import *

if len(sys.argv) < 3:
    print "usage: yaml_diff.py filename1 filename2"
    exit()
else:
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]

#print filename1,filename2
stream=open(filename1,"r")
dict1=yaml.load(stream)
stream=open(filename2,"r")
dict2=yaml.load(stream)
if dict1==dict2:
    print "The two files are identical."
else:
    difference=dict_diff(dict1,dict2)
    print difference
