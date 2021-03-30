#!/usr/bin/env python
# -*- coding: us-ascii -*-
#----------------------------------------------------------------------------
# Use to perform comparison between references and output files
# 1 - read the output
# 2 - search all floating point expressions
# 3 - replace it to have a comparable text
# 4 - compare each floating point expressions

# Use diff because difflib has some troubles (TD)
# Date: 13/02/2014
#----------------------------------------------------------------------------

#import difflib
import subprocess
import getopt
import os
import re
import sys
import tempfile

#Check the version of python
version = sys.version_info
if  version < (2,3,0):
    sys.stderr.write("Detected version %d.%d.%d\n" % tuple(version[:3]))
    sys.stderr.write("Minimal required version is python 2.3.0: Use the command diff\n")
    os.system("diff %s %s" % (file1,file2))
    sys.exit(1)

#Match the version number ex. 1.1.9
re_version = re.compile("[(]ver[ ]+[0-9.\-a-z]+[)]",re.IGNORECASE)
#Match a floating number (d or D are not permitted)
re_float = re.compile("([- ]?[0-9]+[.][0-9]*([DEde][-+]?[0-9]+)?)")

#Maximum discrepancy between float results (default)
max_discrepancy = 1.1e-10
#Use for non-significant discrepancy (float <= min_digits digits)
min_digits = 5

def usage():
    print ("fldiff.py [--mode=m] [--discrepancy=d] [--mindigits=e] [--help] file1 file2")
    print ("  --mode=[bigdft,neb,psolver,pseudo] to compare 'bigdft', 'NEB', 'PS_Check' or pseudo output files")
    print ("  --discrepancy=%7.1e                maximal discrepancy between results" % max_discrepancy)
    print ("  --mindigits=%d                        minimal number of digits for relevant float comparison" % min_digits)
    print ("  --help                               display this message")
    sys.exit(1)

def n_digits(figure):
    """Return the number of decimals of the given figure."""
    n = 0
    for ch in figure:
        if ch == "e":
            return n
        elif ch in "0123456789":
            n += 1
    return n

#Check arguments
try:
    optlist, args = getopt.getopt(sys.argv[1:],"mde:h",["mode=","discrepancy=","mindigits=","help"])
except getopt.error:
    sys.stderr.write("Error in arguments\n")
    usage()
#By default, all modes are False
bigdft  = False
neb     = False
psolver = False
pseudo  = False
for opt,arg in optlist:
    if opt == "-m" or opt == "--mode":
        bigdft  = (arg == "bigdft")
        neb     = (arg == "neb")
        psolver = (arg == "psolver")
        pseudo  = (arg == "pseudo")
    elif opt == "-d" or opt == "--discrepancy":
        max_discrepancy=float(arg)
    elif opt == "-e" or opt == "--mindigits":
        min_digits=int(arg)
    elif opt == "-h" or opt == "--help":
        usage()
if len(args) != 2:
    sys.stderr.write("Error in arguments\n")
    usage()


#Arguments
file1 = args[0]
file2 = args[1]

#Check if the output is a tty to print in colour
start_fail = "\033[0;31m"
start_success = "\033[0;32m"
start_pass = "\033[0;33m"
end = "\033[m"
#if sys.stdout.isatty():
#    start_fail = "\033[0;31m"
#    start_success = "\033[0;32m"
#    start_pass = "\033[0;33m"
#    end = "\033[m"
#else:
#    start_fail = ""
#    start_success = ""
#    start_pass = ""
#    end = ""

#Define a junk line
if bigdft:
    #Test if the line should not be compared (bigdft output)
    def line_junk(line):
        "True if the line must not be compared"
        return re_version.search(line) \
            or " |" in line \
            or "CPU time" in line \
            or "SP-TIMINGS" in line \
            or "Load" in line \
            or "memory" in line \
            or "MB" in line \
            or "proc" in line \
            or "Processes" in line \
            or "allocation" in line \
            or "~W" in line \
            or "for the array" in line \
            or "WRITE WAVES" in line \
            or "wavefunction written" in line \
            or "READING WAVES" in line \
            or "average CG stepsize" in line \
            or "GPU data" in line \
            or "Exact exchange calculation" in line \
            or "time for" in line \
            or "orthonormalization:" in line \
            or "Hamiltonian application:" in line \
            or "orthoconstraint:" in line \
            or "preconditioning:" in line \
            or "other:" in line \
            or "statistics:" in line \
            or "total time:" in line \
            or "DIIS" in line \
            or "Communication overlap ratio" in line \
            or "Timestamp" in line \
            or "Logfile already existing" in line \
	        or "Gathering the potential" in line \
            or "<BigDFT>" in line \
            or "logfiles" in line \
            or "alpha" in line \
            or "wavefunctions need NO reformatting" in line \
            or "WARNING:" in line \
            or "/logfiles/" in line \
            or "--" in line \
            or "Walltime since initialization" in line \
            or "Max No. of dictionaries used" in line \
            or "Version Number" in line \
            or "RUN TIME" in line
#	    or "GEOPT" in line
elif neb:
    # Test if the line should not be compared (NEB output)
    def line_junk(line):
        "True if the line must not be compared"
        return re_version.search(line) \
            or "datadir" in line \
            or "workdir" in line \
            or "Reading atomic input" in line \
            or "Start job" in line \
            or "--with" in line
elif psolver:
    #Remove some lines (PS_Check)
    def line_junk(line):
        "True if the line must not be compared"
        return "MEMORY" in line \
            or "CPLX" in line \
            or "CPU time" in line \
            or "memory" in line \
            or "allocation" in line \
            or "Energy diff" in line \
            or "original" in line \
            or "Max diff at" in line \
            or "result" in line \
            or "METHOD COUNT" in line \
            or "for the array" in line
elif pseudo:
    #Remove lines containing the word time
    def line_junk(line):
        "True if the line must not be compared"
        return "time" in line \
            or "date" in line \
            or "changes" in line \
            or "nodes" in line
else:
    def line_junk(line):
        "Always False except for Hostname"
        return "Hostname" in line \
            or "Memory Peak of process" in line \
            or "Walltime since initialization" in line \
            or "RUN TIME" in line

#Check the last line
end_line = "Memory Consumption Report" 
#end_line_old = "MEMORY CONSUMPTION REPORT"

#Read the first file
try:
    original1 = open(file1).read().replace('\r','').splitlines(True)
except IOError:
    sys.stderr.write("The file '%s' does not exist!\n" % file1)
    sys.exit(1)
#Read the second file
try:
    original2 = open(file2).read().replace('\r','').splitlines(True)
except IOError:
    sys.stderr.write("The file '%s' does not exist!\n" % file2)
    sys.exit(1)

maximum = 0.0
ns_discrepancy = False #Non significant discrepancy.
context_discrepancy = ""
context_lines = ""

#First we compare the first two lines in the case of an added prefix 
#(as in platine computer of CCRT)
#We detect a pattern
if bigdft:
    pattern = '                             BBBB         i       ggggg    '
else:
    pattern = ''

try:
    p1 = original1[0].index(pattern)
    p2 = original2[0].index(pattern)
except ValueError:
    #First line not found ??
    p1 = -1
    p2 = -1
except IndexError:
    sys.stdout.write(start_fail+"One file is blank!\n"+end)
    sys.exit(1)

if p1 >= 0 and p2 >= 0 and p1 != p2:
    #we try something
    prefix1 = original1[0].replace(pattern,'')[:-1]
    prefix2 = original2[0].replace(pattern,'')[:-1]
    if prefix1 != '':
        #We remove the prefix
        original1 = map(lambda x: x.replace(prefix1,''), original1)
    if prefix2 != '':
        #We remove the prefix
        original2 = map(lambda x: x.replace(prefix2,''), original2)

#Remove post output
end_left = False
for line in original1:
    end_left = end_line in line
    if end_left:
        break
end_right = False
for line in original2:
    end_right = end_line in line
    if end_right:
        break

#if bigdft:
#    #Do not compare if a file is not properly finished
#    if not end_left:
#        print "WARNING: The file '%s' is not properly finished!" % file1
#    if not end_right:
#        print "WARNING: The file '%s' is not properly finished!" % file2
#    if not (end_left and end_right): 
#        start = start_fail
#        message = "failed    < "
#        print "%sMax discrepancy : %s (%s%s)%s" % (start,"Nan",message,max_discrepancy,end)
#        sys.exit(1)

#Remove line_junk before comparing (the line number is wrong)
memory = 0
if bigdft or psolver:
    time = None
#Open 2 temporary files
t1 = tempfile.NamedTemporaryFile(mode='w')
nl_t1 = 0
for line in original1:
    if not line_junk(line):
        nl_t1 += 1
        t1.write(line)
    else:
        #Only for BigDFT and PSolver
        if bigdft or psolver:
            #Keep sum of elapsed time
            if "CPU time/ELAPSED time" in line:
                if time == None:
                    time = float(line.split()[-2])
                else:
                    time += float(line.split()[-2])
            #Test if memory remaining is 0
            if "Remaining Memory" in line:
                memory = int(line.split()[-1])
t1.flush()

t2 = tempfile.NamedTemporaryFile(mode='w')
nl_t2 = 0
for line in original2:
    if not line_junk(line):
        t2.write(line)
        nl_t2 += 1
t2.flush()

#Generate comparison using the unix diff command
compare = iter(subprocess.getoutput("diff -a -b -d %s %s" %(t1.name,t2.name)).splitlines(True))

t1.close()
t2.close()

print ("max_discrepancy=%8.2e   file=%s   ref=%s" % (max_discrepancy,file1,file2))
try:
    line = next(compare)
    EOF = False
    #Display the maximum discrepancy
except StopIteration:
    #Nothing to compare
    print ('%sNothing to compare (all lines are identical)%s' % (start_success,end))
    EOF = True

context_lines = None
while not EOF:
    #A new context is detected
    context = line
    print_context = False
    left = list()
    line = next(compare)
    while line[0] == "<":
        left.append(line)
        try:
            line = next(compare)
        except StopIteration:
            #We have reached the end of file
            EOF = True
            break
    right = list()
    if line[0] == "-":
        line = next(compare)
    while line[0] == ">":
        right.append(line)
        try:
            line = next(compare)
        except StopIteration:
            #We have reached the end of file
            EOF = True
            break
    #We have a new context and two lists to compare
    n1 = len(left)
    n2 = len(right)
    i1 = -1
    i2 = -1
    while i1 < n1-1 and i2 < n2-1:
        i1 += 1
        i2 += 1
        line1 = left[i1]
        line2 = right[i2]
        #We avoid some lines (try before to be more robust)
        if line_junk(line1) or line_junk(line2):
            continue
        floats1 = list()
        for (one,two) in re_float.findall(line1):
            #'d' is not recognised by python
            one = one.lower().replace("d","e")
            floats1.append((float(one), n_digits(one)))
        floats2 = list()
        for (one,two) in re_float.findall(line2):
            #'d' is not recognised by python
            one = one.lower().replace("d","e")
            floats2.append((float(one), n_digits(one)))
        #Replace all floating point by XXX and ' ' characters
        new1 = re_float.sub('XXX',line1[2:]).replace(' ','')
        new2 = re_float.sub('XXX',line2[2:]).replace(' ','')
        if new1 != new2 and i1 == 0:
            #For the first difference, we display the context
            print (context)
        print_context = True
        if new1 != new2:
            print (line1)
            print (line2)
        n = len(floats1)
        if n == len(floats2):
            diff_discrepancy = False
            for i in range(n):
                tt = abs(floats1[i][0]-floats2[i][0])
                if min(floats1[i][1], floats2[i][1]) > min_digits:
                    #Case of significant discrepancy.
                    if maximum < tt:
                        context_discrepancy = " (line %s)" % context.split("c")[0].split(",")[0]
                        context_lines = "\n"+context_discrepancy[1:]+"\n"+line1+line2
                        maximum = max(maximum,tt)
                    if tt > max_discrepancy:
                        diff_discrepancy = True
                else:
                    #Case of discrepancy on non significant values.
                    if tt > max_discrepancy:
                        diff_discrepancy = True
                        ns_discrepancy = True
            if diff_discrepancy and new1 == new2:
                if not print_context:
                    print (context)
                    print_context = True
                print (line1)
                print (line2)
        elif (maximum < 99):
            print ("%s the number of floating point differs" % context[:-1])
            context_discrepancy = " (line %s)" % context.split("c")[0].split(",")[0]
            context_lines = "\n"+context_discrepancy[1:]+"\n"+line1+line2
            maximum = 99
    #Add lines if necessary
    while i1 < n1-1:
        i1 += 1
        if n1 > 0 and not line_junk(left[i1]):
            if not print_context:
                print (context)
            print_context = True
            print (left[i1])
            #Check if supplementary lines are floating numbers
            #TD (16/11/2011): Do no check!
            #floats = list()
            #for (one,two) in re_float.findall(left[i1]):
            #    #'d' is not recognised by python
            #    one = one.lower().replace("d","e")
            #    floats.append((float(one), n_digits(one)))
            #if len(floats) > 0:
            #    maximum = 99
    while i2 < n2-1:
        i2 += 1
        if n2 > 0 and not line_junk(right[i2]):
            if not print_context:
                print (context)
            print_context = True
            print (right[i2])
            #Check if supplementary lines are floating numbers
            #TD (16/11/2011): Do no check!
            #floats = list()
            #for (one,two) in re_float.findall(right[i2]):
            #    #'d' is not recognised by python
            #    one = one.lower().replace("d","e")
            #    floats.append((float(one), n_digits(one)))
            #if len(floats) > 0:
            #    maximum = 99

if context_lines is not None:
    print (context_lines)
else:
    print ()

#Check if the test is OK.
if maximum > max_discrepancy or memory:
    start = start_fail
    if memory:
        message = "failed-memory remaining-(%sB) " % memory
    else:
        message = "failed    < "
elif ns_discrepancy:
    start = start_pass
    message = "passed    < "
else:
    start = start_success
    message = "succeeded < "

if (bigdft or psolver) and time != None:
    print ("%sMax discrepancy %s: %s (%s%s) -- time %7.2f%s " % \
        (start,context_discrepancy,maximum,message,max_discrepancy,time,end))
elif nl_t1 == 0 and nl_t2 == 0:
    #No lines in the output as waited
    print ("%sNo output as reference%s" % (start,end))
else:
    print ("%sMax discrepancy %s: %s (%s%s)%s" % \
        (start,context_discrepancy,maximum,message,max_discrepancy,end))

sys.exit(0)
