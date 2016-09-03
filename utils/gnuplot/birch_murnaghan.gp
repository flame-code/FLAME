#! /usr/bin/gnuplot

#set terminal postscript eps enhanced
#set output "plot1.eps"
#set title "plot1"

set key top left box 3 
set xlabel "lattice constant (atomic units)"
set ylabel "Energy (atomic units)"
#-------------------------------------------------------------------------------
filename=system("echo $filename")
#-------------------------------------------------------------------------------
f(x)=a0+a1*x+a2*x**2
a0=1 ; a1=1 ; a2=1
fit f(x) filename u 1:2 via a0,a1,a2
#print "a2 ",a2
#print "min ",-a1/(2*a2)
#plot filename u 1:2 pt 7 ps 2.0 title "Data" w p , f(x) w l lc 1 lt -1 title "Fit"
#pause -1
#-------------------------------------------------------------------------------
au2GPa=29422.9 #conversion of pressure units from atomic units to GPa

#a --> B_0^1
#c --> V_0
#e --> E_0
#b --> B_0
a=3.5;
c=(-a1/(2*a2))**3
b=2*c*a2
e=a0

f1(x)=(c/x**(3))**(2.0/3.0)
f2(x)=f1(x)-1

E(x) = e+ \
      (9.0*c*b/16.0)* (\
      (f2(x)**3)*a + \
      f2(x)**2 * (6.0-4.0*f1(x)) )

#E(x) = a*x**2 + b*x +c

fit E(x) filename u 1:2 via a,b,c,e

#set xr [10:10.4]

print "dB_0/dp ",a
print "Bulk modulus ",b*au2GPa
print "lattice constant ",c**(1.0/3.0)*0.529177
print "energy minimum ",e

plot filename u 1:2 pt 7 ps 2.0 title "Data" w p , E(x) w l lc 1 lt -1 title "Fit"
#plot E(x) lc 1 lt -1 title "Fit" w l

pause -1

