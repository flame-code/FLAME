import sys

def adjust_real(r,n1,n2):
    if r<10.0:
        m1=1
    elif r<100.0:
        m1=2
    elif r<1000.0:
        m1=2
    elif r<10000.0:
        m1=4
    else:
        sys.exit("ERROR: only real values smaller than 10000.0 are accpted.")
    if m1>n1: sys.exit("ERROR: number of requested digits is smaller than digits of the value.")
    frmt="%%%1d.%df" % (m1,n2) #%tr(m1+n2+1)"
    str_r=""
    for i in range(n1-m1):
        str_r+="0"
    str_r+=frmt % r
    return str_r
