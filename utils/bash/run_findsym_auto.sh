#!/bin/bash
#export ISODATA=/home/maxamsler/Homefolder/VariableCell/ISOTROPY_2/isobyu/
ascii2findsym.py $1 $2
findsym <findsym.in >findsym.out
grep -A 100000 "# CIF file" findsym.out> findsym.out.cif

