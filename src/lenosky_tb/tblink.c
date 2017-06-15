#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "defines.h" 

void radelmgenralsp();
void eselfgeneral();
void eselfgeneral2();
void hubbardu();
void setuptb();
void setuptb2();
void clssplint();
double passtbcut();
int isnorth();
void clspairenergy();
int passgstride();
double evalPair();

potl_type pplocal;
int mepslocal,mhtblocal;
double paircut;
int usepairpot;

int passgstride()
{
if (mhtblocal == 1 || mhtblocal == 2 || (mhtblocal==3 && FALSE == TRUE)) return(1);
if ((mhtblocal == 3 ) || mhtblocal == 6) return(2);
if (mhtblocal == 4 || mhtblocal == 8 || mhtblocal ==7 || mhtblocal ==11) return(4);
}

double evalPair(arg1,arg2,arg3,arg4,arg5)
int arg1;
double arg2;
double *arg3;
double *arg4;
int arg5;
{
	clssplint(&(pplocal.phi[arg1]),arg2,arg3,arg4,arg5);

}

double passtbcut(usetb)
int *usetb;
{
/*Returns the necessary cutoff for tight binding h's to an external caller*/ 
double cutoff;
int j;

cutoff= 0.0;

for(j=0;j<mhtblocal;j++)
  cutoff = MAX((pplocal.h[j].x[pplocal.h[j].npt]),cutoff);

*usetb = 0;
if (mhtblocal!=0 && mepslocal!=0) *usetb = 1;

return(cutoff);
}

int isnorth()
{
if (mhtblocal==4 || mhtblocal==7 || mhtblocal == 0) return(0) ;
if (mhtblocal==8 || mhtblocal==11) return(1) ;
}

int isxtal()
{
if (mhtblocal==7 || mhtblocal == 11 ) return(1) ;
if ((mhtblocal < 7 && (!(mhtblocal==3 && FALSE== TRUE ) ) )|| mhtblocal==4 || mhtblocal == 8 || mhtblocal == 0) return(0) ;
/* Note here I've assumed that there won't be any crystal field terms for the hydrogen
   fitting*/
}

void radelmgeneralsp(r,radar,dradar,atomtypei,atomtypej)
double *r,*radar,*dradar;
int atomtypei,atomtypej;
{
/*This routine returns h(r)'s for tight binding link: units eV and Angstrom*/
/* Index i denotes
 0 for s,s,sigma
 1 for s,p,sigma
 2 for p,p,sigma
 3 for p,p,pi */
/*Need to return h(*r), but for now returns a dummy argument*/
int i;
double y,der;
if (atomtypei + atomtypej==0) {
for(i=0;i<4;i++)
  {
  clssplint(&(pplocal.h[i]),*r,&y,&der,0);
  radar[i] = y/(*r)/(*r);
  dradar[i] = der/(*r)/(*r) - 2.0 * y/(*r)/(*r)/(*r) ;
  }
}
else if (atomtypei + atomtypej==1) {
for(i=0;i<2;i++) {
  clssplint(&(pplocal.h[4+i]),*r,&y,&der,0);
  radar[i] = y/(*r)/(*r);
  dradar[i] = der/(*r)/(*r) - 2.0 * y/(*r)/(*r)/(*r) ;
}
radar[2] = dradar[2] = radar[3] = dradar[3] = 0.0;
}
else if (atomtypei + atomtypej==2) {
for(i=0;i<1;i++) {
  clssplint(&(pplocal.h[6+i]),*r,&y,&der,0);
  radar[i] = y/(*r)/(*r);
  dradar[i] = der/(*r)/(*r) - 2.0 * y/(*r)/(*r)/(*r) ;
}
radar[1] = dradar[1] = radar[2] = dradar[2] = radar[3] = dradar[3] = 0.0;
}

}

void eselfgeneral(eself)
double *eself;
{
/*Return self-energy numbers, es, and ep in certain order*/
static int firstcall = 1;
eself[0] = pplocal.eps[0];
//eself[0] = -5.670225;
eself[1] = -eself[0];
eself[2] = -eself[0];
eself[3] = -eself[0];

firstcall = 0;
}

void hubbardu(u)
double *u;
{
*u = 0.0;
}

void eselfgeneral2(eself)
double *eself;
{
*eself = pplocal.eps[1];
}

void setuptb(pp,mepstb,mhtb)
potl_type *pp;
int mepstb,mhtb;
{
pplocal = *pp;
mepslocal = mepstb;
mhtblocal = mhtb;
}

void setuptb2(gnecut,usephi)
double gnecut;
int usephi;
{
paircut = gnecut;
usepairpot = usephi;
}

void clssplint(s,x,y,deriv,extype)
double x,*y,*deriv;
spline_type (*s);
int extype; /*Flag to control method of extrapolation on lower bound 0 = linear, 1= 1/r/r and 1/r */
/*modified version of clssplint: returns df/dx in deriv*/
/*also modified to use endpoint deriv information if out of bounds*/
{
        int klo,khi,k;
        double h,b,a;
if(x>(*s).x[1] && x < (*s).x[(*s).npt]) { /*NR part*/
        klo=1;
        khi=(*s).npt;
        while (khi-klo > 1) {
                k=(khi+klo) >> 1;
                if ((*s).x[k] > x) khi=k;
                else klo=k;
        }
        h=(*s).x[khi]-(*s).x[klo];
        if (h == 0.0) {
              printf("clssplint called with x = %lf\n",x);
              };

        a=((*s).x[khi]-x)/h;
        b=(x-(*s).x[klo])/h;
        *y=a*(*s).y[klo]+b*(*s).y[khi]+((a*a*a-a)*(*s).y2[klo]+(b*b*b-b)*(*s).y2[khi])*(h*h)/6.0;
        *deriv = -(1/h)*((*s).y[klo]-(*s).y[khi] + ((h*h)/6.0)*((3*a*a-1)*(*s).y2[klo] -
                    (3*b*b-1)*(*s).y2[khi])    );
}
else /*Extrapolation past endpoints*/
{
if(x<=(*s).x[1])
{
if(extype == 0) {
*deriv = (*s).yp1;
*y = (*s).y[1] + *deriv * (x - (*s).x[1]);
  }
else
  {
/*Fit 1/x and 1/x/x termination*/
double a,b;
double rc,yc,derc;
rc=(*s).x[1];
yc=(*s).y[1];
derc=(*s).yp1;

b = -rc*rc*rc*derc - rc*rc*yc;
a = rc*yc - b/rc;
 
*deriv = -a/x/x - 2.0 * b/x/x/x;
*y = a/x + b/x/x;
  }
}
else
{
*deriv = (*s).ypn;
*y = (*s).y[(*s).npt] + *deriv * (x - (*s).x[(*s).npt]);
}
}
 
}
