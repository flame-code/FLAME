#!/usr/bin/env python
import numpy as numpy
from numpy import linalg as LA
import sys

if len(sys.argv) < 3:
    print "usage: fit_elastic_constant.py input_filename"
    exit()
else:
    filename1 = sys.argv[1]
    filename2 = sys.argv[2]

#var=["x11","x21","x22","x31","x32","x33","x41","x42","x43","x44","x51","x52","x53","x54","x55","x61","x62","x63","x64","x65","x66"]
var=["x11","x21","x31","x41","x51","x61","x22","x32","x42","x52","x62","x33","x43","x53","x63","x44","x54","x64","x55","x65","x66"]
#for j in var:
#    eval(j).array([])

x11=[];x21=[];x31=[];x41=[];x51=[];x61=[];x22=[];x32=[];x42=[];x52=[];x62=[];
x33=[];x43=[];x53=[];x63=[];x44=[];x54=[];x64=[];x55=[];x65=[];x66=[]
b=numpy.empty([6])
f=open(filename1,"r")
n=0
#read line into array 
for line in f.readlines():
    a=[]
    #loop over the elemets, split by whitespace
    for i in line.split():
        #convert to integer and append to the last element of the list
        a.append(float(i))
    b[0]=a[0];b[1]=a[2];b[2]=a[5];b[3]=a[1]*2.0;b[4]=a[3]*2.0;b[5]=a[4]*2.0;
    for j in var: eval(j).append(0.0)
    x11[n]=b[0];x21[n]=b[1];x31[n]=b[2];x41[n]=b[3];x51[n]=b[4];x61[n]=b[5];n+=1
    for j in var: eval(j).append(0.0)
    x21[n]=b[0];x22[n]=b[1];x32[n]=b[2];x42[n]=b[3];x52[n]=b[4];x62[n]=b[5];n+=1
    for j in var: eval(j).append(0.0)
    x31[n]=b[0];x32[n]=b[1];x33[n]=b[2];x43[n]=b[3];x53[n]=b[4];x63[n]=b[5];n+=1
    for j in var: eval(j).append(0.0)
    x41[n]=b[0];x42[n]=b[1];x43[n]=b[2];x44[n]=b[3];x54[n]=b[4];x64[n]=b[5];n+=1
    for j in var: eval(j).append(0.0)
    x51[n]=b[0];x52[n]=b[1];x53[n]=b[2];x54[n]=b[3];x55[n]=b[4];x65[n]=b[5];n+=1
    for j in var: eval(j).append(0.0)
    x61[n]=b[0];x62[n]=b[1];x63[n]=b[2];x64[n]=b[3];x65[n]=b[4];x66[n]=b[5];n+=1

f.closed
if n!=len(x11):
    print "ERROR: "
    exit()

f=open(filename2,"r")
z=[]
#read line into array 
for line in f.readlines():
    a=[]
    #loop over the elemets, split by whitespace
    for i in line.split():
        #convert to integer and append to the last element of the list
        a.append(float(i))
    b[0]=a[0];b[1]=a[2];b[2]=a[5];b[3]=a[1];b[4]=a[3];b[5]=a[4];
    z.append(-b[0])
    z.append(-b[1])
    z.append(-b[2])
    z.append(-b[3])
    z.append(-b[4])
    z.append(-b[5])

f.closed


v = numpy.array([eval(i) for i in var])

coeff, residues, rank, singval = numpy.linalg.lstsq(v.T, z)

#print "Fitting process completed."
#print "f(x,y)=(%s) + (%s) x + (%s) y + " % (coeff[0], coeff[1], coeff[2])
#print "       (%s) x^2 + (%s) xy + (%s) y^2\n" % (coeff[3], coeff[4], coeff[5])

print "residues= ", residues[0], "\n"
print "rank= ", rank, "\n"
#print singval

c=numpy.empty([6,6])

var=["x11","x21","x31","x41","x51","x61","x22","x32","x42","x52","x62","x33","x43","x53","x63","x44","x54","x64","x55","x65","x66"]
c[0][0]=coeff[0];                                                                                          
c[1][0]=coeff[1];c[1][1]=coeff[ 6];                                                                        
c[2][0]=coeff[2];c[2][1]=coeff[ 7];c[2][2]=coeff[11];                                                      
c[3][0]=coeff[3];c[3][1]=coeff[ 8];c[3][2]=coeff[12];c[3][3]=coeff[15];                                    
c[4][0]=coeff[4];c[4][1]=coeff[ 9];c[4][2]=coeff[13];c[4][3]=coeff[16];c[4][4]=coeff[18];                  
c[5][0]=coeff[5];c[5][1]=coeff[10];c[5][2]=coeff[14];c[5][3]=coeff[17];c[5][4]=coeff[19];c[5][5]=coeff[20];

print "%8.2f"                          % (c[0][0]                                        )
print "%8.2f%8.2f"                     % (c[1][0],c[1][1]                                )
print "%8.2f%8.2f%8.2f"                % (c[2][0],c[2][1],c[2][2]                        )
print "%8.2f%8.2f%8.2f%8.2f"           % (c[3][0],c[3][1],c[3][2],c[3][3]                )
print "%8.2f%8.2f%8.2f%8.2f%8.2f"      % (c[4][0],c[4][1],c[4][2],c[4][3],c[4][4]        )
print "%8.2f%8.2f%8.2f%8.2f%8.2f%8.2f" % (c[5][0],c[5][1],c[5][2],c[5][3],c[5][4],c[5][5])

print "\n"

c_lambda=(c[1][0]+c[2][0]+c[2][1])/3.0
c_mu=((c[0][0]+c[1][1]+c[2][2])/3.0-c_lambda)/2.0
c_BulkModulus=((c[0][0]+c[1][1]+c[2][2])/3.0-c_lambda)/3.0+c_lambda
c_mut=(c[3][3]+c[4][4]+c[5][5])/3.0
print "K (Bulk Modulus)=    %12.4f" % c_BulkModulus
print "mu=                  %12.4f" % c_mu
print "mut (Shear Modulus)= %12.4f" % c_mut

for i in range(6):
    for j in range(6):
        if i>=j: continue
        c[i][j]=c[j][i]
        
carr = numpy.array(c)
w, v = LA.eig(carr)
print "Eigenvalues: ", w
#print "Eigenvectors: ", v

#for i in range(n):
#    tt=0.0
#    k=0
#    for j in var:
#        tt+=coeff[k]*eval(j)[i]
#        k+=1
#        #print "ERR= %s" % eval(j)[i]
#    tt=tt-z[i]
#    #print "ERR= %10.4f" % tt

#res_my=0.0
#for i in range(n):
#    tt=coeff[0]+coeff[1]*x[i]+coeff[2]*y[i]+coeff[3]*x[i]**2+ \
#        coeff[4]*x[i]*y[i]+coeff[5]*y[i]**2 - z[i]
#    res_my+=tt**2
#
#print "MY: residue: %f\n" % res_my

#print "--------------------------------------------------"
#a = [[2*coeff[3], coeff[4]], [coeff[4], 2*coeff[5]]]
##print a
#b = numpy.array(a)
#w, v = LA.eig(b)
#print "Eigenvalues: ", w
#print "Eigenvectors: ", v
#print "--------------------------------------------------"
#print "To plot in gnuplot use the following: \n"
#print "f(x,y)=c0_xy + c1_x * x + c1_y * y + c2_x *x*x + c1_xy *x*y + c2_y * y*y"
#print "c0_xy=%20.10f" % float(coeff[0])
#print "c1_x =%20.10f" % float(coeff[1])
#print "c1_y =%20.10f" % float(coeff[2])
#print "c2_x =%20.10f" % float(coeff[3])
#print "c1_xy=%20.10f" % float(coeff[4])
#print "c2_y =%20.10f" % float(coeff[5])
#print "%7s%s%18s" % ("splot '",filename1, "' u 1:2:3 , f(x,y)")
#print "--------------------------------------------------"






#x = numpy.array([row[:][0] for row in arr])
#y = numpy.array([row[:][1] for row in arr])
#z = numpy.array([row[:][2] for row in arr])
#np=len(x)
#print "\nThere are %d data points after filtering.\n" % np
#print "len(x) %d %s" % (len(x),zmaxtol)
#exit()

#non-symmetric
#c11 c12 c13 c14 c15 c16
#c21 c22 c23 c24 c25 c26
#c31 c32 c33 c34 c35 c36
#c41 c42 c43 c44 c45 c46
#c51 c52 c53 c54 c55 c56
#c61 c62 c63 c64 c65 c66

#symmetric
#c11 c21 c31 c41 c51 c61
#c21 c22 c32 c42 c52 c62
#c31 c32 c33 c43 c53 c63
#c41 c42 c43 c44 c54 c64
#c51 c52 c53 c54 c55 c65
#c61 c62 c63 c64 c65 c66

#[["c11","c21","c31","c41","c51","c61"],
# ["c21","c22","c32","c42","c52","c62"],
# ["c31","c32","c33","c43","c53","c63"],
# ["c41","c42","c43","c44","c54","c64"],
# ["c51","c52","c53","c54","c55","c65"],
# ["c61","c62","c63","c64","c65","c66"]]

#"x11","x21","x22","x31","x32","x33","x41","x42","x43","x44","x51","x52","x53","x54","x55","x61","x62","x63","x64","x65","x66"
#"x11","x21","x31","x41","x51","x61","x22","x32","x42","x52","x62","x33","x43","x53","x63","x44","x54","x64","x55","x65","x66"

