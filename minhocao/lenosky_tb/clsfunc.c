#include <stdlib.h>
#include <math.h>
#include <stdio.h>
#include "defines.h"

#define EBCC 0.0 

#define PHIEXTRAPOLATIONTYPE 0 /*Use 0 for linear, 1 for 1/r and 1/r/r */ 

#define NSYMFLAG TRUE /*Always set to TRUE*/ 

/*Simple routine used in neighbor listing*/
void clsatneighb();

/*Routines to make neighbor lists*/
void add_to_neighbor_list();
void make_neighbor_list();

/*Routine to compute cosine theta given side lengths of a triangle, plus derivatives of cosine
  theta with respect to the side lengths of the triangle*/
void clspccos();

/*The main energy routine.*/
void clsenergy();

/*Calculate all the three body terms.  This is called twice: once by clszrhomaxmin 
 (to set xbeta's) and once directly by clsenergy*/
void triloop();

/*Calculate contribution to stillinger-weber or modified embedded atom potential for any
given triangle of three atoms.  Can be called to set xbetas (in which case forces are not 
computed) or to set forces.*/ 
void calctri();

/*Computes betas and xbetas for a given frame.  Calls triloop with set xbeta flag TRUE.
  Also precomputes values of xf so that they are only computed once in the code.*/
void clszrhomaxmin();

/*Read in a potential from the file "coeff.cls"*/
void prmst38c();
void prmst38c_meam();

/*Read spline parameters from an external file.  Called by setup routine, prmst38c */
void clsfread_spline();

/*This is just a modified version of the Numerical Recipes routine spline, which generates 
  spline coefficients given a set of spline parameters.*/
void clsspline(); 

/*Set nonorthogonal lattice vectors gg.  This routine is used when clsenergy is called by an 
  external fitting routine*/
/*Set parameters.  This routine is used when clsenergy is going to be called by an external fitting
  routine*/
void passgg();

void clssetvar();

/*Initialize ranges so that they will be set to ranges of following data*/
void init_ranges();

/*Pass values of range parameters to main fitting code*/
void pass_ranges();

/*Maximum force*/
double gmaxforce = 0.0;

float ran1();

/*Ranges of various splines in potential*/
double clsgnecut;  /*Range of pair potential*/
double clsrhonecut;  /*Range for rho(r) and xrho(r)*/
double clstricut1;  /*Range for f(r) and xf(r) */

int gnat; /*Number of atoms*/
position ga; /*Lattice constants*/
double ge; /*Global energy value*/
static position rsave[CLSMAXATOM]; /*Save array for positions*/

/*FLAGS and potential*/
static int clsusephi;
static int clsusetri; /*TRUE if Stillinger-weber is being used*/
static int clsusec; /*Use pair potential*/
static int clsmrho;  /*Number of embedded atom potentials*/
static int clsmrho2; /*Number of modified embedded atom potentials*/
static int clsmepstb; /*Number of epsilon values for tight binding*/
static int clsmhtb;  /*Number of h(r) splines for tight binding*/
static potl_type clspp;  /*Record containing all the splines for the potential*/

/*gg is lattice vectors for the cell, and gg_inv is the inverse of the gg matrix*/
/*default is the unit matrix*/
/*Note that frame must be called with coordinates along these lattice vectors, but
  forces are returned in ordinary cartesian coordinates*/
static double gg[3][3] = {{1,0,0},{0,1,0},{0,0,1}};
static double gg_inv[3][3] = {{1,0,0},{0,1,0},{0,0,1}};

/*RANGES for potential:*/
static int ranges_not_set=TRUE; /*This flag means that ranges have not been set*/
static double clsumax[CLSMRHOMAX]; /*Maximum argument for U */
static double clsumin[CLSMRHOMAX]; /*Minimum argument for U */
static double clsxumax[CLSMRHOMAX2]; /*Maximum argument for XU */
static double clsxumin[CLSMRHOMAX2]; /*Minimum argument for XU */
static double clsgmax; /*Maximum argument for stillinger weber g*/
static double clsgmin; /*Minimum argument for stillinger weber g*/
static double clsxgmax[CLSMRHOMAX2]; /*Maximum argument for each g in MEAM potential */
static double clsxgmin[CLSMRHOMAX2]; /*Minimum argument for each g in MEAM potential */
static double clsrmin; /*Minimum neighbor distance*/  

static int ncall = 0; /*Total number of func_LOCALtion calls to clsenergy */
/* ---------------------------------------------------------- */


void clssetvar(p,usetri,usephi,mrho,mhtb,mepstb,gnecut,rhonecut,tricut1,mrho2)
/*Set variables for external calls from fitting code */
potl_type *p;
int usetri,usephi,mrho,mhtb,mepstb,mrho2;
double gnecut,rhonecut,tricut1;
{
clspp=*p;
clsusetri=usetri;
clsusephi=usephi;
clsmrho=mrho;
clsmrho2=mrho2;
clsmhtb =mhtb;
clsmepstb = mepstb;
clsgnecut=gnecut;
clsrhonecut=rhonecut;
clstricut1=tricut1;
setuptb(&clspp,clsmepstb,clsmhtb);
}

void clsatneighb(frame,framep,x,y,at,nb)
/*Routine called for neighbor listing: given atoms x and y
  in given frame, returns lesser of the two in *at,
  and index for greater in neighbor list of *at in *nb */
int x,y;
int *at,*nb;
frame_type *frame;
clsframepp_type *framep;
{
int i,at2,nn;

*at = (x<y) ? x : y ;
at2 = x + y - *at ; /* The other atom besides *at */
nn=(*framep).nb[*at];    
*nb = 0; /*If not in list returns clszero*/
for(i=0;i<nn;i++)
  if ((*framep).nbind[*at][i] == at2)
    *nb = i;
}
 
void clsppccos(pg,x,y,z)
/*Computes cosine theta (pg[0]) given a triangle with side lengths x,y,z.
  Theta is the angle across from the side with length z. 
  Also computes the derivatives of cosine theta with respect to x, y, z, 
  and returns them in pg[1], pg[2], and pg[3] respectively.

  Alternately, if OLDPCCOS is set to false this routine merely returns 
  the triangle side length "z", rescaled in an apropriate fashion.  
  In this case the derivatives with respect to x, y, z are returned 
  as before.
*/
double *pg,x,y,z;
{ 
pg[0] = (x/y + y/x - z*z/x/y)/2.0;
pg[1] = (1/y -y/x/x + z*z /x/x/y)/2.0;
pg[2] = (1/x -x/y/y + z*z /y/y/x)/2.0;
pg[3] = -z/x/y;
}

void clsenergy(frame,framep,err,e1)
frame_type *frame;
clsframepp_type *framep;
int *err;
double *e1;
/*Take *frame and calculate upon it, creating data record (*framep).
Energy and forces are returned inside frame...

Sets *err = 1 if there is some problem, otherwise leaves it alone.
*e1 is the single atom energy of the first atom.

If (*frame).nsym is less than (*frame).n then this routine computes the 
energy using the primitive part of cell only.  This option is broken, and
should not be trusted unless carefully tested in a specific case.
 */
{

#define SWTX(a,b,t) {t=a; a=b; b=t;} /*Swap a and b, using t as a temporary storage variable;
       all three variables can be of any type */

int exflag; /*FLAG which should normally be set to zero to extrapolate splines linearly
      past their endpoints */
int i,j,k,m,mm,jj;
position g;
double val,r,y,afact,der,der2;
potl_type *p; /*Potential data type which contains all the splines needed to evaluate the potential */
static int firstcall = 1;

if (firstcall) init_ranges(); firstcall = 0;
ncall++;
ranges_not_set = FALSE; /*Because this call will set the ranges*/


/*Don't use the broken nsym option if lenosky is set true (recommended) */
#if (NSYMFLAG)
/*MODIFICATION because NSYM is broken for 3 body terms*/
(*frame).nsym = (*frame).n ;
#endif

/*Set p (a local variable) equal to the global variable clspp that contains the potential*/
p = &clspp;

/* Set forces, energy, and lattforce (energy derivatives with respect to the nonorthogonal
   lattice vectors) to zero to start the routine) */
for(i=0;i<(*frame).n;i++)
  { (*frame).f[i].x = 0.0; (*frame).f[i].y = 0.0; (*frame).f[i].z = 0.0;}
(*frame).en = 0.0; 
*e1 = 0;
for(i=0;i<3;i++) 
  for(j=0;j<3;j++)
    (*frame).lattforce[i][j] = 0.0;

/* If this test is true, then something is wrong, so return without computing anything */ 
if (clsmhtb==0 && clsmepstb!=0)
  return; 

/*MAKE NEIGHBOR LISTS FOR FRAME*/
make_neighbor_list(frame,framep,err);

/* START OF MAIN CODE TO EVALUATE THE POTENTIAL */

/*EVALUATE THE ARGUMENTS of the embedding potentials "u" and "xu" */
/*These arguments are returned in (*framep).beta, (*framep).xbeta */
/*This routine also precomputes the values of xf so that they are only evaluated once. */
clszrhomaxmin(frame,framep,p);

/*Evaluate the value of the embedded atom potential by finding U(beta).
  Additionally set calculate dU/d(beta) and set beta equal to this value.
  The chain rule will then be used to calculate the derivatives of the 
  embedded atom potential. */
for(j=0;j<(*frame).nsym;j++) /*NOTE presence of nsym*/
  for(mm=0;mm<clsmrho;mm++)
  { /*Final step is to get deriv of u at beta*/
  double t; 
  t = (*framep).beta[mm][j];
  if (t < clsumin[mm]) clsumin[mm] = t;
  if (t > clsumax[mm]) clsumax[mm] = t;
  clssplint(&(*p).u[mm],t,&y,&der,0);
  (*framep).beta[mm][j] = der;
  (*frame).en += y;
if (j==0)
  *e1 += y;
  }


/*Now do the same thing for the modified embedded atom potential. 
  Here xbeta takes the place of beta.  The energy is evaluated,
  and dU/d(xbeta) is calculated.  xbeta is set equal to this value
  so that the chain rule can later be used to evaluate derivatives of the
  potential */
for(j=0;j<(*frame).nsym;j++) /*NOTE presence of nsym*/
  for(mm=0;mm<clsmrho2;mm++)
  { /*Final step is to get deriv of u at beta*/
  double t;
  t = (*framep).xbeta[mm][j];
  if (t < clsxumin[mm]) clsxumin[mm] = t;
  if (t > clsxumax[mm]) clsxumax[mm] = t;
  clssplint(&(*p).xu[mm],t,&y,&der,0);
  (*framep).xbeta[mm][j] = der;
  (*frame).en += y;
if (j==0)
  *e1 += y;
  }

/*Evaluate the pair potential and its derivatives. */
 
for(j=0;j<(*frame).n;j++)
  for(k=0;k<(*framep).nb[j];k++) {
    m=(*framep).nbind[j][k];
    g=(*framep).grad[j][k];
    r=(*framep).r[j][k];

if (clsusephi) {
    clssplint(&((*p).phi[0]),r,&y,&der,PHIEXTRAPOLATIONTYPE);
    (*frame).f[j].x += -der * (g.x);
    (*frame).f[j].y += -der * (g.y);
    (*frame).f[j].z += -der * (g.z);
    (*frame).f[m].x += der * (g.x);
    (*frame).f[m].y += der * (g.y);
    (*frame).f[m].z += der * (g.z);
/*FLAGY*/
/*FLAGY.5*/
if (TRUE) {
position ptemp,ptemp2;
int iloc,jloc;
ptemp.x = der*g.x*r; ptemp.y = der*g.y*r; ptemp.z = der*g.z*r;

REVERT(ptemp,ptemp2);
(*frame).lattforce[0][0] += g.x * ptemp2.x;
(*frame).lattforce[0][1] += g.x * ptemp2.y;
(*frame).lattforce[0][2] += g.x * ptemp2.z;
(*frame).lattforce[1][0] += g.y * ptemp2.x;
(*frame).lattforce[1][1] += g.y * ptemp2.y;
(*frame).lattforce[1][2] += g.y * ptemp2.z;
(*frame).lattforce[2][0] += g.z * ptemp2.x;
(*frame).lattforce[2][1] += g.z * ptemp2.y;
(*frame).lattforce[2][2] += g.z * ptemp2.z;
}


/*Evaluate the forces and "lattforce" (derivatives of the energy with respect to 
  the nonorthogonal lattice vectors) due to the embedded atom potential*/

if (j<(*frame).nsym)
    (*frame).en += y/2;
if (m<(*frame).nsym)
    (*frame).en += y/2;
if (j==0)
  *e1 += y/2;
    } /*if (clsusephi)*/

    if (k<(*framep).nbrho[j])
      for(mm=0;mm<clsmrho;mm++) {
    if (mm==0 && clsusephi == 0) 
      exflag = 1; 
    else 
      exflag = 0; 
    clssplint(&(*p).rho[mm],r,&y,&der2,exflag);
    afact = ((*framep).beta[mm][j] + (*framep).beta[mm][m]) ;
/*FLAGY*/
    (*frame).f[j].x += -der2 * (g.x) * afact;
    (*frame).f[j].y += -der2 * (g.y) * afact;
    (*frame).f[j].z += -der2 * (g.z) * afact;
    (*frame).f[m].x += der2 * (g.x) * afact;
    (*frame).f[m].y += der2 * (g.y) * afact;
    (*frame).f[m].z += der2 * (g.z) * afact;
/*FLAGY.6*/
if (TRUE) {
position ptemp,ptemp2;
int iloc,jloc;
ptemp.x = der2*g.x*r*afact; ptemp.y = der2*g.y*r*afact; ptemp.z = der2*g.z*r*afact;

REVERT(ptemp,ptemp2);
(*frame).lattforce[0][0] += g.x * ptemp2.x;
(*frame).lattforce[0][1] += g.x * ptemp2.y;
(*frame).lattforce[0][2] += g.x * ptemp2.z;
(*frame).lattforce[1][0] += g.y * ptemp2.x;
(*frame).lattforce[1][1] += g.y * ptemp2.y;
(*frame).lattforce[1][2] += g.y * ptemp2.z;
(*frame).lattforce[2][0] += g.z * ptemp2.x;
(*frame).lattforce[2][1] += g.z * ptemp2.y;
(*frame).lattforce[2][2] += g.z * ptemp2.z;
}
      };

/*Evaluate the forces and lattforces (derivatives of the energy with respect to the nonorthogonal
  lattice vectors) resulting from the modified embedded atom potential, only for the func_LOCALtion
  rho(r).  Forces and lattforces due to f(r), g(r), et cetera are computed seperately. */
    if (k<(*framep).nbrho[j])
      for(mm=0;mm<clsmrho2;mm++) {
    exflag = 0;
    clssplint(&(*p).xrho[mm],r,&y,&der2,exflag);
    afact = ((*framep).xbeta[mm][j] + (*framep).xbeta[mm][m]) ;
/*FLAGY*/
    (*frame).f[j].x += -der2 * (g.x) * afact;
    (*frame).f[j].y += -der2 * (g.y) * afact;
    (*frame).f[j].z += -der2 * (g.z) * afact;
    (*frame).f[m].x += der2 * (g.x) * afact;
    (*frame).f[m].y += der2 * (g.y) * afact;
    (*frame).f[m].z += der2 * (g.z) * afact;
/*FLAGY.1*/
if (TRUE) { 
position ptemp,ptemp2;
int iloc,jloc;
ptemp.x = der2*g.x*r*afact; ptemp.y = der2*g.y*r*afact; ptemp.z = der2*g.z*r*afact;

REVERT(ptemp,ptemp2);
(*frame).lattforce[0][0] += g.x * ptemp2.x;
(*frame).lattforce[0][1] += g.x * ptemp2.y;
(*frame).lattforce[0][2] += g.x * ptemp2.z;
(*frame).lattforce[1][0] += g.y * ptemp2.x;
(*frame).lattforce[1][1] += g.y * ptemp2.y;
(*frame).lattforce[1][2] += g.y * ptemp2.z;
(*frame).lattforce[2][0] += g.z * ptemp2.x;
(*frame).lattforce[2][1] += g.z * ptemp2.y;
(*frame).lattforce[2][2] += g.z * ptemp2.z;
}
      };

} /*End of loop over pairs*/

/*If Stillinger-Weber is being used, then the values for f(r) and the derivatives f'(r) need
to be tabulated in a special record "fval" in (*framep).  Here we do this*/  
if (clsusetri)
{ for(j=0;j<(*frame).n;j++)
  for(k=0;k<(*framep).nbtri[j];k++)
      {
      clssplint(&((*p).f),(*framep).r[j][k],&val,&der,0);
      (*framep).fval[j][k][0] = val;
      (*framep).fval[j][k][1] = der;
      }
 
}

/*Compute the energy for Stillinger-Weber; the forces for Stillinger-Weber; and 
  finally the forces due to the f and g splines in the modified embedded atom potential*/
triloop(frame,framep,err,e1,0);

/*Return the total energy */
(*frame).en *= (   (double)(*frame).n / ((double)(*frame).nsym)   ) ;

}

void triloop(frame,framep,err,e1,sxb)
frame_type *frame;
clsframepp_type *framep;
int *err;
double *e1;
int sxb;
/*
This routine performs calculations for all the three body terms in the potential.
It locates each independent triangle and then calls calctri on it.

Because neighbor listing method is fairly complicated there are three separate
types of triangles.  For the final two types, calctri is only called if a given triangle is
not calculated in one of the other cases.
 
Sets *err = 1 if there is some problem, otherwise leaves it alone.
*e1 is the single atom energy of the first atom.

Sxb is the flag to set xbeta's (add tri terms to them).  If it is 1,
then (*frame).en is not changed, and forces are not added.

 */
{
#define SWTX(a,b,t) {t=a; a=b; b=t;}
 
int exflag,i,j,k,m,mm;
position g;
double r,y,afact,der,der2;
potl_type *p;
p = &clspp;

/* ---------------------------- */
if(clsusetri || clsmrho2 != 0)
{
double r1,r2,r3,val,der;
position g1,g2,g3;
int ii,jj,kk,iix,jjx,kkx,tempi;


/*LOOP over triangles: TRIANGLES OF FIRST TYPE*/
for(i=0;i<(*frame).n-2;i++)
  for(jj=0;jj<(*framep).nbtri[i];jj++) {
    j=(*framep).nbind[i][jj];
    for(kk=0;kk<(*framep).nbtri[j];kk++)
      {
      k=(*framep).nbind[j][kk];
      r1=(*framep).r[i][jj];
      r3=(*framep).r[j][kk];
      g1=(*framep).grad[i][jj];
      g3=(*framep).grad[j][kk];
      CALCGRADR(r1,g1,r3,g3,r2,g2);
//      CELLDISTNO(((*frame).p[i]),((*frame).p[k]),((*frame).a),r2);
//      CELLGRADNO(((*frame).p[i]),((*frame).p[k]),((*frame).a),g2);
      calctri(i,j,k,r1,r2,r3,g1,g2,g3,frame,framep,p,e1,sxb);
      }
    }

/*TRIANGLES OF SECOND TYPE*/
for(i=0;i<(*frame).n-2;i++)
  for(jj=0;jj<(*framep).nbtri[i];jj++)
    for(kk=jj+1;kk<(*framep).nbtri[i];kk++) 
      {
      j=(*framep).nbind[i][jj];
      k=(*framep).nbind[i][kk];
{
jjx =jj;
kkx =kk;
if (j>k) {SWTX(j,k,tempi); SWTX(jjx,kkx,tempi);}
}
      r1 = (*framep).r[i][jjx];
      r2 = (*framep).r[i][kkx];
      g1 = (*framep).grad[i][jjx];
      g2 = (*framep).grad[i][kkx];
      CALCGRADR(-r1,g1,r2,g2,r3,g3);

//      CELLDISTNO(((*frame).p[j]),((*frame).p[k]),((*frame).a),r3);
      if(r3 > clstricut1) {/*Don't count a triangle twice */
//CELLGRADNO(((*frame).p[j]),((*frame).p[k]),((*frame).a),g3);
calctri(i,j,k,r1,r2,r3,g1,g2,g3,frame,framep,p,e1,sxb); 
}
      }
/*TRIANGLES OF THIRD TYPE*/
for(k=2;k<(*frame).n;k++)
  for(ii=0;ii<(*framep).nbtrixx[k];ii++)
    for(jj=ii+1;jj<(*framep).nbtrixx[k];jj++)
      {
      i=(*framep).nbindxx[k][ii];
      j=(*framep).nbindxx[k][jj];
{
iix =ii;
jjx =jj;
if (i>j) {SWTX(i,j,tempi); SWTX(iix,jjx,tempi);}
}


//      CELLDISTNO(((*frame).p[i]),((*frame).p[j]),((*frame).a),r1);
//      if(r1>clstricut1) { /*Don't count a triangle twice*/
//      CELLGRADNO(((*frame).p[i]),((*frame).p[j]),((*frame).a),g1);
      CELLDISTNO(((*frame).p[i]),((*frame).p[k]),((*frame).a),r2);
      CELLGRADNO(((*frame).p[i]),((*frame).p[k]),((*frame).a),g2); 
      CELLDISTNO(((*frame).p[j]),((*frame).p[k]),((*frame).a),r3);
      CELLGRADNO(((*frame).p[j]),((*frame).p[k]),((*frame).a),g3);
      CALCGRADR(r2,g2,-r3,g3,r1,g1);

      if(r1>clstricut1) { /*Don't count a triangle twice*/
      calctri(i,j,k,r1,r2,r3,g1,g2,g3,frame,framep,p,e1,sxb);


      }
      }
}
}

void calctri(i,j,k,r1,r2,r3,g1,g2,g3,frame,framep,p,e1,sxb)
/*This routine performs calculations on a triangle given the list in framep.
  i,j,k are the atoms for vertices -- must currently be in ascending order.
  r1,r2,r3 are side lengths; g1,g2,g3 are gradients used in force calculations.
  Returns force and energy from term, by adding them to contents of frame..

  Because terms have the generic form f(rij)f(rik) g(cos(theta jik) )
  there are three terms that can be nonzero for a given triangle:
     centered on i
     centered on j
     centered on k
  
  This routine computes each of these terms seperately and then adds them up
  to get the net forces.

  Sets the array xbeta if the flag sxb == 1 only; in this cases forces are
  not computed. 
*/
clsframepp_type *framep;
frame_type *frame;
int i,j,k;
double r1,r2,r3;
position g1,g2,g3;
potl_type *p;
double *e1;
int sxb;
{
int m1,m2;
int tu0,tu1,tu2;
double *p1,*p2,*p3,*p1a,*p2a,*p3a;
double c1[4],c2[4],c3[4];
double der1,der2,der3,v;
double y,der;
double va,vb,vc,der1a,der1b,der1c,der2a,der2b,der2c;
double derextraa,derextrab,derextrac;
double der3a,der3b,der3c;
double dxbeta_dz;
int mm, sw, farg,swon,garg;
int test;

if (TRUE) {
test = FALSE && ((i==0) || (j==0) || (k==0));
if (test) printf("FLAGX, calctri, i, j, k = %d %d %d\n",i,j,k);
if (test) printf("FLAGX, calctri, r1,r2,r3 = %lf %lf %lf\n",r1,r2,r3);
if (test) printf("FLAGX, calctri, g1 = %lf %lf %lf\n",g1.x,g1.y,g1.z);
if (test) printf("FLAGX, calctri, g2 = %lf %lf %lf\n",g2.x,g2.y,g2.z);
if (test) printf("FLAGX, calctri, g3 = %lf %lf %lf\n",g3.x,g3.y,g3.z);

}

clsatneighb(frame,framep,i,j,&m1,&m2);
p1 = (*framep).fval[m1][m2];
clsatneighb(frame,framep,i,k,&m1,&m2);
p2 = (*framep).fval[m1][m2];
clsatneighb(frame,framep,j,k,&m1,&m2);
p3 = (*framep).fval[m1][m2];
tu0 = (r1<=clstricut1 && r2<=clstricut1);
tu1 = (r1<=clstricut1 && r3<=clstricut1);
tu2 = (r2<=clstricut1 && r3<=clstricut1);
if(tu0)
  clsppccos(c1,r1,r2,r3);
if(tu1)
  clsppccos(c2,r1,r3,r2);
if(tu2)
  clsppccos(c3,r2,r3,r1);

/*LOOP over terms:*/
if (clsusetri==1) sw =1 ; else sw = 0;
if (sxb==1) sw = 0; /*Don't do sw for clszrhomaxmin*/

for(mm=0;mm<sw+clsmrho2;mm++) {  
farg = 2*(mm + (1 - sw)) ; /*Argument for looking up precomputed f*/
swon = (sw==1 && mm==0) ; /*Flag to tell if sw is on & index is on it*/
garg = mm - sw ; /*Index for argument of g*/
p1a = p1 + farg;
p2a = p2 + farg;
p3a = p3 + farg;

v=der1=der2=der3=0.0;
va=vb=vc=der1a=der1b=der1c = 0;
derextraa=derextrab=derextrac=0;
der2a=der2b=der2c =0 ; 
der3a = der3b = der3c =0;

/*Compute three body term with angle CENTERED on i*/
if (tu0) {
  if (swon)
/*NOTE that GAUSSIANG is a flag to use a gaussian instead of the func_LOCALtion g(r).  Normally
  it is set to FALSE*/
    { 
    double t;
    t = c1[0];
    clssplint(&((*p).g),c1[0],&y,&der,0);
    if (t > clsgmax) clsgmax = t;
    if (t < clsgmin) clsgmin = t;
    }
  else
    { double t;
      t = c1[0];
    clssplint(&((*p).xg[garg]),c1[0],&y,&der,0);
    if (t > clsxgmax[garg]) clsxgmax[garg] = t;
    if (t < clsxgmin[garg]) clsxgmin[garg] = t;
    }
  CLSACPOTL(va,/*r1,r2,r3,*/ der1a,der2a,der3a,p1a,p2a,c1,y,der);
  }

if (i==0 && sxb==0 && swon)
  *e1 += va;

if (i<(*frame).nsym && sxb==0 && swon)
  (*frame).en += va;

/*Compute three body term with angle CENTERED on j*/
if (tu1) { 
  if (swon)
    { 
    double t;
    t = c2[0];
    clssplint(&((*p).g),c2[0],&y,&der,0);
    if (t > clsgmax) clsgmax = t;
    if (t < clsgmin) clsgmin = t;
    }
  else
    { double t;
      t = c2[0];
    clssplint(&((*p).xg[garg]),c2[0],&y,&der,0);
    if (t > clsxgmax[garg]) clsxgmax[garg] = t;
    if (t < clsxgmin[garg]) clsxgmin[garg] = t;
    }
  CLSACPOTL(vb,/*r1,r3,r2,*/ der1b,der3b,der2b,p1a,p3a,c2,y,der);
  }
if (j<(*frame).nsym && sxb==0 && swon)
  (*frame).en += vb;

/*Compute three body term with angle CENTERED ON k*/
if (tu2) {
  if (swon)
    { 
    double t;
    t = c3[0];
    clssplint(&((*p).g),c3[0],&y,&der,0);
    if (t > clsgmax) clsgmax = t;
    if (t < clsgmin) clsgmin = t;
    }
  else
    { double t;
      t = c3[0];
    clssplint(&((*p).xg[garg]),c3[0],&y,&der,0);
    if (t > clsxgmax[garg]) clsxgmax[garg] = t;
    if (t < clsxgmin[garg]) clsxgmin[garg] = t;
    }
  CLSACPOTL(vc,/*r2,r3,r1,*/ der2c,der3c,der1c,p2a,p3a,c3,y,der);

  }
if (k<(*frame).nsym && sxb == 0 && swon) 
  (*frame).en += vc;

v=va+vb+vc;
if (swon)
  {/* OLD CODE for SW (before NSYM): if (sxb==0) (*frame).en += v; */ }
else
  if (sxb==1)
  { /*Add to xbeta values*/
  (*framep).xbeta[garg][i] += va;
  (*framep).xbeta[garg][j] += vb;
  (*framep).xbeta[garg][k] += vc; 
  }

if (swon) {
  der1 = der1a + der1b + der1c;
  der2 = der2a + der2b + der2c;
  der3 = der3a + der3b + der3c;
  } 
else {
  double b1,b2,b3;
/*Let b1,b2,b3 = xbeta for atom i,j,k and meam potl # == garg:*/ 
  b1 = (*framep).xbeta[garg][i];
  b2 = (*framep).xbeta[garg][j];
  b3 = (*framep).xbeta[garg][k];
  der1 = der1a*b1 + der1b*b2 + der1c*b3;
  der2 = der2a*b1 + der2b*b2 + der2c*b3;
  der3 = der3a*b1 + der3b*b2 + der3c*b3;
  if (test) printf("der1, der2, der3 = %lf %lf %lf \n",der1,der2,der3);
  if (test) printf("(i,j),(i,k),(j,k)\n");
 }

if (sxb==0) {
/*FLAGY*/
/*Compute force due to partial derivative of r1 = distance (m,j)*/
    if(test) printf("i,j grad = %lf %lf %lf\n",g1.x,g1.y,g1.z);
    (*frame).f[i].x += -der1 * (g1.x);
    (*frame).f[i].y += -der1 * (g1.y);
    (*frame).f[i].z += -der1 * (g1.z);
    (*frame).f[j].x += der1 * (g1.x);
    (*frame).f[j].y += der1 * (g1.y);
    (*frame).f[j].z += der1 * (g1.z);
/*FLAGY.2*/
if (TRUE) {
position ptemp,ptemp2;
int iloc,jloc;
ptemp.x = der1*g1.x*r1; ptemp.y = der1*g1.y*r1; ptemp.z = der1*g1.z*r1;

REVERT(ptemp,ptemp2);
(*frame).lattforce[0][0] += g1.x * ptemp2.x;
(*frame).lattforce[0][1] += g1.x * ptemp2.y;
(*frame).lattforce[0][2] += g1.x * ptemp2.z;
(*frame).lattforce[1][0] += g1.y * ptemp2.x;
(*frame).lattforce[1][1] += g1.y * ptemp2.y;
(*frame).lattforce[1][2] += g1.y * ptemp2.z;
(*frame).lattforce[2][0] += g1.z * ptemp2.x;
(*frame).lattforce[2][1] += g1.z * ptemp2.y;
(*frame).lattforce[2][2] += g1.z * ptemp2.z;
}
/*Compute force due to partial derivative of r2 = distance (m,k)*/
    if(test) printf("i,k grad = %lf %lf %lf\n",g2.x,g2.y,g2.z);
    (*frame).f[i].x += -der2 * (g2.x);
    (*frame).f[i].y += -der2 * (g2.y);
    (*frame).f[i].z += -der2 * (g2.z);
    (*frame).f[k].x += der2 * (g2.x);
    (*frame).f[k].y += der2 * (g2.y);
    (*frame).f[k].z += der2 * (g2.z);
/*FLAGY.3*/
if (TRUE) {
   position ptemp,ptemp2;
   int iloc,jloc;
   ptemp.x = der2*g2.x*r2; ptemp.y = der2*g2.y*r2; ptemp.z = der2*g2.z*r2;
   
  REVERT(ptemp,ptemp2);
  (*frame).lattforce[0][0] += g2.x * ptemp2.x;
  (*frame).lattforce[0][1] += g2.x * ptemp2.y;
  (*frame).lattforce[0][2] += g2.x * ptemp2.z;
  (*frame).lattforce[1][0] += g2.y * ptemp2.x;
  (*frame).lattforce[1][1] += g2.y * ptemp2.y;
  (*frame).lattforce[1][2] += g2.y * ptemp2.z;
  (*frame).lattforce[2][0] += g2.z * ptemp2.x;
  (*frame).lattforce[2][1] += g2.z * ptemp2.y;
  (*frame).lattforce[2][2] += g2.z * ptemp2.z;
  }
/*Compute force due to partial derivative of r3 = distance (j,k)*/
    if (test) printf("j,k grad =%lf %lf %lf\n",g3.x,g3.y,g3.z);
    (*frame).f[j].x += -der3 * (g3.x);
    (*frame).f[j].y += -der3 * (g3.y);
    (*frame).f[j].z += -der3 * (g3.z);
    (*frame).f[k].x += der3 * (g3.x);
    (*frame).f[k].y += der3 * (g3.y);
    (*frame).f[k].z += der3 * (g3.z);
/*FLAGY.4*/
if (TRUE) {
   position ptemp,ptemp2;
   int iloc,jloc;
   ptemp.x = der3*g3.x*r3; ptemp.y = der3*g3.y*r3; ptemp.z = der3*g3.z*r3;
   
  REVERT(ptemp,ptemp2);
  (*frame).lattforce[0][0] += g3.x * ptemp2.x;
  (*frame).lattforce[0][1] += g3.x * ptemp2.y;
  (*frame).lattforce[0][2] += g3.x * ptemp2.z;
  (*frame).lattforce[1][0] += g3.y * ptemp2.x;
  (*frame).lattforce[1][1] += g3.y * ptemp2.y;
  (*frame).lattforce[1][2] += g3.y * ptemp2.z;
  (*frame).lattforce[2][0] += g3.z * ptemp2.x;
  (*frame).lattforce[2][1] += g3.z * ptemp2.y;
  (*frame).lattforce[2][2] += g3.z * ptemp2.z;
  }
}
  } /*End of loop over mm*/
}

void clszrhomaxmin(frame,framep,p)
/*Computes betas (for embedded atom) and xbetas (for modified embedded atom method)
  for the given frame.  These values are just (Sum on j rho(rij)) and  
(Sum on j xrho(rij) + Sum on j,k f(rij) f(rik) g (cos(theta(jik))) ) evaluated on each
site i.

  This routine also computes the values "fval" which are precomputed values of xf and "zval" which are precomputed values of xz and xz'
*/
potl_type *p;
frame_type *frame;
clsframepp_type *framep;
{
int jj,j,k;
double y,der;
int m,mm;
double r;
int exflag;
int err;
double e1;
double hjk,derhjk;

/*Evaluate betas for embedded atom method*/
for(mm=0;mm<clsmrho;mm++) {
for(j=0;j<(*frame).n;j++)
  (*framep).beta[mm][j] = 0.0;
/*Precompute betas for pair potential*/
for(j=0;j<(*frame).n;j++)
  for(k=0;k<(*framep).nbrho[j];k++) {
  m = (*framep).nbind[j][k];
  r = (*framep).r[j][k];
  /*now add rho(r) to beta*/
  if (mm==0 && clsusephi == 0)
    exflag = 1;
  else
    exflag = 0;
  clssplint(&((*p).rho[mm]),r,&y,&der,exflag);
  (*framep).beta[mm][j] += y;
  (*framep).beta[mm][m] += y;
  }  
}    

/*Precompute values of xf */
for (jj=0;jj<clsmrho2;jj++)
 for(j=0;j<(*frame).n;j++)
  for(k=0;k<(*framep).nbtri[j];k++)
      {
      double val,der;
      clssplint(&((*p).xf[jj]),(*framep).r[j][k],&val,&der,0);
      (*framep).fval[j][k][2+2*jj] = val;
      (*framep).fval[j][k][3+2*jj] = der;
      }

/*Compute xbetas for modified embedded atom method*/
/*Part 1 of computing xbetas: compute contribution from xrho*/ 
/*Exflag controls extrapolation method, here exflag = 0 is used for linear extrapolation*/
for(mm=0;mm<clsmrho2;mm++) {
for(j=0;j<(*frame).n;j++)
  (*framep).xbeta[mm][j] = 0.0;
for(j=0;j<(*frame).n;j++)
  for(k=0;k<(*framep).nbrho[j];k++) {
  m = (*framep).nbind[j][k];
  r = (*framep).r[j][k];
  /*now add rho(r) to xbeta*/
  exflag = 0;
  clssplint(&((*p).xrho[mm]),r,&y,&der,exflag);
  (*framep).xbeta[mm][j] += y;
  (*framep).xbeta[mm][m] += y;
  }  
}    

/*Part 2 of computing xbetas: compute contribution from three body terms.  Note set
xbeta FLAG is set to 1 in this call: */
if (clsmrho2!=0)
  triloop(frame,framep,&err,&e1,1); /*Add tri terms to xbetas*/
}

void passgg(ggx)
/*Pass the nonorthogonal lattice vectors gg, when clsenergy is to be called by an external
fitting code*/
double ggx[3][3];
{
int i,j;
for(i=0;i<3;i++)
  for(j=0;j<3;j++)
    gg[i][j] = ggx[i][j];
/*Set gg_inv = inverse of (gg) */
{
double det;
                det = gg[0][0]*gg[1][1]*gg[2][ 2] +
                    gg[1][ 0]*gg[2][1]*gg[0][ 2] +
                    gg[2][ 0]*gg[0][1]*gg[1][ 2] -
                    gg[0][ 0]*gg[2][1]*gg[1][ 2] -
                    gg[1][ 0]*gg[0][1]*gg[2][ 2] -
                    gg[2][ 0]*gg[1][1]*gg[0][ 2];
                gg_inv[0][0] = (1.0/det) * (gg[1][1]*gg[2][2] - gg[2][1]*gg[1][2]);
                gg_inv[0][1] = (1.0/det) * (gg[2][1]*gg[0][2] - gg[0][1]*gg[2][2]);
                gg_inv[0][2] =-(1.0/det) * (gg[1][1]*gg[0][2] - gg[0][1]*gg[1][2]);

                gg_inv[1][0] =-(1.0/det) * (gg[1][0]*gg[2][2] - gg[2][0]*gg[1][2]);
                gg_inv[1][1] =-(1.0/det) * (gg[2][0]*gg[0][2] - gg[0][0]*gg[2][2]);
                gg_inv[1][2] = (1.0/det) * (gg[1][0]*gg[0][2] - gg[0][0]*gg[1][2]);

                gg_inv[2][0] = (1.0/det) * (gg[1][0]*gg[2][1] - gg[2][0]*gg[1][1]);
                gg_inv[2][1] = (1.0/det) * (gg[2][0]*gg[0][1] - gg[0][0]*gg[2][1]);
                gg_inv[2][2] =-(1.0/det) * (gg[1][0]*gg[0][1] - gg[0][0]*gg[1][1]);
}
}

void lmndouble_frame(l,mm,n,frx,fror)
int l,mm,n;
frame_type *frx, *fror;
{
/*Note incoming frame only has p,a,and n set, everything else can be ignored*/
/*Everything is malloc'ed externally*/
/*frx is 8 frames of fror, BUT itype is set to 1 for every atom */
int norig,i,j,k,ind,m;
position disp,pin,pres;

norig = (*fror).n;

(*frx).n = l*mm*n * norig;
((*frx).a).x = l * ((*fror).a).x;
((*frx).a).y = mm * ((*fror).a).y;
((*frx).a).z = n * ((*fror).a).z;

ind = 0;
for(i=0;i<=l-1;i++)
  for(j=0;j<=mm-1;j++)
    for(k=0;k<=n-1;k++)
      for(m=0;m<norig;m++) {
  disp.x = ((*fror).a).x * i;
  disp.y = ((*fror).a).y * j;
  disp.z = ((*fror).a).z * k;

  pin = (*fror).p[m];
  VECT_ADD(disp,pin,pres);
  (*frx).p[ind] = pres;

  ind++;
  }

(*frx).nsym = (*fror).n; /*NSYM is number of atoms in primitive cell*/
//if (lowersym)
//  (*frx).nsym = lowersymnum;
}


void lenoskymeam_(natom,alat,rxyz,ff,erxyz,count_md,glattforce,gg1,n_double,pot_type)
int *natom;
double *alat;
double *rxyz;
double *ff;
double *erxyz;
double *count_md;
double gg1[3][3];
double glattforce[3][3];
int n_double;
int pot_type;
{
static clsframepp_type framepp;
static frame_type frame;
static frame_type doubled_frame;
double e1;
int err;
static int firstcall = 1;
int i,j;

if (firstcall == 1)
  {
    prmst38c_meam(pot_type);
    doubled_frame.p = (position *) malloc (CLSMAXATOM*sizeof(position));
    doubled_frame.f = (position *) malloc (CLSMAXATOM*sizeof(position));
    firstcall = 0;
  }

passgg(gg1);
*count_md =  *count_md + 1.0;

frame.n = *natom;
frame.nsym = *natom; /*No symmetries assumed*/
frame.a.x = alat[0];
frame.a.y = alat[1];
frame.a.z = alat[2];
frame.p = (position *) rxyz;
frame.f = (position *) ff;

lmndouble_frame(n_double,n_double,n_double,&doubled_frame,&frame);

/*Perform the calculations:*/
clsenergy(&doubled_frame,&framepp,&err,&e1);
/*Forces are set already by reference*/
*erxyz = doubled_frame.en / (n_double * n_double * n_double);
for(i=0;i<3;i++)
	for(j=0;j<3;j++)
		glattforce[i][j] = doubled_frame.lattforce[i][j] / (n_double * n_double * n_double);
for(i=0;i<*natom;i++)
	frame.f[i] = doubled_frame.f[i];
}

void prmst38c()
{
int mm;
FILE *f;

f = fopen("coeff.cls","r");

clsgnecut=5.24;
clsrhonecut=clstricut1 = 4.0;

fscanf(f,"%d %d %d %d %d %d",&clsusephi,&clsusetri,&clsmrho,&clsmhtb,&clsmepstb,&clsmrho2);
if (clsusephi)
  {
  clsfread_spline(f,&(clspp.phi[0]));
  clsfread_spline(f,&(clspp.phi[1]));
  clsfread_spline(f,&(clspp.phi[2]));
  }
for(mm=0;mm<clsmrho;mm++) {
  clsfread_spline(f,&(clspp.rho[mm]));
  clsfread_spline(f,&(clspp.u[mm]));
  }
if (clsusetri) {
  clsfread_spline(f,&(clspp.f));
  clsfread_spline(f,&(clspp.g));
  };
for(mm=0;mm<clsmhtb;mm++)
  clsfread_spline(f,&(clspp.h[mm]));
for(mm=0;mm<clsmepstb;mm++)
  fscanf(f,"%lf",&(clspp.eps[mm]));
for(mm=0;mm<clsmrho2;mm++) {
#if (MEAMRHO)
  clsfread_spline(f,&(clspp.xrho[mm]));
#endif
  clsfread_spline(f,&(clspp.xu[mm]));
  clsfread_spline(f,&(clspp.xf[mm]));
  clsfread_spline(f,&(clspp.xg[mm]));
  }

fclose(f);

setuptb(&clspp,clsmepstb,clsmhtb);
setuptb2(clsgnecut,clsusephi);
}


void prmst38c_meam(pot_type)
int pot_type;
{
int mm;
FILE *f;

printf("prmst38c_meam called, reading in meam parameters\n");
if (pot_type == 1)
  f = fopen("coeff.meam.cls","r");
else
	f = fopen("coeff.meam.cls.Ti","r");

if (pot_type == 1) {
clsgnecut=4.5;
clsrhonecut=clstricut1 = 3.5;
}
else
{
	clsgnecut = 5.5;
	clsrhonecut = clstricut1 = 4.41;
}
fscanf(f,"%d %d %d %d %d %d",&clsusephi,&clsusetri,&clsmrho,&clsmhtb,&clsmepstb,&clsmrho2);
if (clsusephi)
  clsfread_spline(f,&(clspp.phi[0]));
for(mm=0;mm<clsmrho;mm++) {
  clsfread_spline(f,&(clspp.rho[mm]));
  clsfread_spline(f,&(clspp.u[mm]));
  }
if (clsusetri) {
  clsfread_spline(f,&(clspp.f));
  clsfread_spline(f,&(clspp.g));
  };
for(mm=0;mm<clsmhtb;mm++)
  clsfread_spline(f,&(clspp.h[mm]));
for(mm=0;mm<clsmepstb;mm++)
  fscanf(f,"%lf",&(clspp.eps[mm]));
for(mm=0;mm<clsmrho2;mm++) {
#if (MEAMRHO)
  clsfread_spline(f,&(clspp.xrho[mm]));
#endif
  clsfread_spline(f,&(clspp.xu[mm]));
  clsfread_spline(f,&(clspp.xf[mm]));
  clsfread_spline(f,&(clspp.xg[mm]));
  }

if(clsusephi)
  clsgnecut = clspp.phi[0].x[clspp.phi[0].npt] ;
if(clsmrho!=0)
  clsrhonecut = clspp.rho[0].x[clspp.rho[0].npt] ;
if(clsusetri)
  clstricut1 = clspp.f.x[clspp.f.npt] ;
else if(clsmrho2!=0)
  clstricut1 = clspp.xf[0].x[clspp.xf[0].npt] ;

fclose(f);

setuptb(&clspp,clsmepstb,clsmhtb);
setuptb2(clsgnecut,clsusephi);
}

void clsfread_spline(f,s)
/*Read spline parameters from an external file*/
FILE *f;
spline_type *s;
{
/*Read s from f, using format common with fwrite_spline(f,s) */
int i;
fscanf(f,"%d",&((*s).npt));
fscanf(f,"%lf %lf",&((*s).yp1),&((*s).ypn));
fscanf(f,"%d %d %d %d",&((*s).iset1),&((*s).isetn),&((*s).ider1),&((*s).idern));
for(i=1;i<=(*s).npt;i++)
  fscanf(f,"%lf %lf %lf",&((*s).x[i]),&((*s).y[i]),&((*s).y2[i]));

clsspline(s);
}

void clsspline(s)
/*This is just a modified version of the Numerical Recipes routine spline, which generates 
  spline coefficients given a set of spline parameters.*/
spline_type (*s);
{
int i,k;
double p,qn,sig,un,u[NSPMAX+1];

if (((*s).yp1) > 0.99e30)
(*s).y2[1]=u[1]=0.0;
else {
(*s).y2[1] = -0.5;
u[1]=(3.0/((*s).x[2]-(*s).x[1]))*(((*s).y[2]-(*s).y[1])/((*s).x[2]-(*s).x[1])-(*s).yp1);
}
for (i=2;i<=(*s).npt-1;i++) {
sig=((*s).x[i]-(*s).x[i-1])/((*s).x[i+1]-(*s).x[i-1]);
p=sig*((*s).y2[i-1])+2.0;
(*s).y2[i]=(sig-1.0)/p;
u[i]=((*s).y[i+1]-(*s).y[i])/((*s).x[i+1]-(*s).x[i]) - ((*s).y[i]-(*s).y[i-1])/((*s).x[i]-(*s).x[i-1]);
u[i]=(6.0*u[i]/((*s).x[i+1]-(*s).x[i-1])-sig*u[i-1])/p;
}
if ((*s).ypn > 0.99e30)
qn=un=0.0;
else {
qn=0.5;
un=(3.0/((*s).x[(*s).npt]-(*s).x[(*s).npt-1]))*((*s).ypn-((*s).y[(*s).npt]-(*s).y[(*s).npt-1])/((*s).x[(*s).npt]-(*s).x[(*s).npt-1]));
}
(*s).y2[(*s).npt]=(un-qn*u[(*s).npt-1])/(qn*((*s).y2[(*s).npt-1])+1.0);
for (k=(*s).npt-1;k>=1;k--)
(*s).y2[k]= ((*s).y2[k]*((*s).y2[k+1]))+u[k];
}


void init_ranges()
/*Initialize ranges so that they will be set to the ranges present in the data set which
follows*/
{
int i;
for(i=0;i<CLSMRHOMAX;i++)
  {
  clsumax[i] = -1e3;
  clsumin[i] = 1e3;
  };
for(i=0;i<CLSMRHOMAX2;i++)
  {
  clsxumax[i] = -1e3;
  clsxumin[i] = 1e3;
  clsxgmax[i] = -1e3;
  clsxgmin[i] = 1e3;
  };
clsgmax = -1e3;
clsgmin = 1e3;
clsrmin = 1e3;
}

void pass_ranges(umax, umin, xumax, xumin, gmax, gmin, xgmax, xgmin, rmin, use_dummy_ranges) 
double *umax;
double *umin;
double *xumax;
double *xumin;
double *gmax;
double *gmin;
double *xgmax;
double *xgmin;
double *rmin;
int use_dummy_ranges;

/*Pass ranges to main fitting code, once they have been accumulated*/
{

int i;

printf("PASS RANGES CALLED WWW\n");

if (ranges_not_set == TRUE) {init_ranges(); ranges_not_set = FALSE; };
 
for(i=0;i<CLSMRHOMAX;i++)
  {
  umax[i] = clsumax[i];
  umin[i] = clsumin[i]; 
  };
for(i=0;i<CLSMRHOMAX2;i++)
  {
  xumax[i] = clsxumax[i];
  xumin[i] = clsxumin[i];
  xgmax[i] = clsxgmax[i];
  xgmin[i] = clsxgmin[i];  
  };
*gmax = clsgmax;
*gmin = clsgmin;
*rmin = clsrmin; /*See line immediately below */
}

void add_to_neighbor_list(frame,framep,r,j,k,err) 
frame_type *frame;
clsframepp_type *framep;
double r;
int j,k;
int *err;
/*On call k is greater than j; this routine adds atom k to the neighbor list of 
  atom i.  r is the distance between j and k, and it must be less than clsgnecut.
  This routine also creates reverse neighbor list. */
{
position g;
CELLGRADNO(((*frame).p[j]),((*frame).p[k]),((*frame).a),g);
(*framep).nbind[j][(*framep).nb[j]] = k;
(*framep).r[j][(*framep).nb[j]] = r;
(*framep).grad[j][(*framep).nb[j]] = g;
(*framep).nb[j]++;
if(r <= clstricut1) {
  /*Add atom j to reversed 'xx' list of atom k*/
  (*framep).nbindxx[k][(*framep).nbtrixx[k]] = j;
  (*framep).nbtrixx[k]++;
  }
if ((*framep).nb[j] >= CLSMAXNEIGHB || (*framep).nbtrixx[k] >= CLSMAXNEIGHB)
  {
  printf("WARNING clsenergy: neighbor limit exceeded (%d %d)\n",(*framep).nb[j], (*framep).nbtrixx[k]);
  *err=1;
  }
if (r==0.0)
  {
  printf("WARNING clsenergy: r= 0.0 error\n");
  *err = 1;
  }
} 

void make_neighbor_list(frame,framep,err)
frame_type *frame;
clsframepp_type *framep;
int *err;
{
int i,j,k,m;
position q1,q2,q3;
int nxcell,nycell,nzcell,nx,ny,nz,mx,my,mz;
/* *************************************************************************************************** */
/* NEIGHBOR LISTING CODE */
/* 
  (*framep) is the record which contains all the neighbor lists.

  There are two neighbor lists, one for atoms with greater index, and one with atoms with lesser
  index:
  (*framep).nb[j] is the number of neighbors of atom j; only atoms with index k>j are stored
      in this list.  The cutoff for this list is clsgnecut.
  (*framep).nbtrixx[j] is the number of atoms of atom j; only atoms with index k<j are stored
      in this list.  The cutoff for this list is clstricut1 which is the cutoff for the three body
      terms in the potential.  

  In the first neighbor list: 
  (*framep).nbind[j][m] stores the index of the mth neighbor of atom j.
  (*framep).r[j][m] stores the distance to the mth neighbor of atom j. 
  (*framep).grad[j][m] stores the direction vector to the mth neighbor of atom j.
  (*framep).nbrho[j] is the number of neighbors of atom j inside the cutoff "clsrhonecut".
  (*framep).nbtri[j] is the number of neighbors of atom j inside the cutoff "clstricut1".

  In the second neighbor list:
  (*framep).nbindxx[j][m] stores the index of the mth neighbor of atom j.
  
  Note that the main neighbor list is sorted: for an atom j, the neighbors are sorted in increasing
  order of distance from j.
*/

/*Set number of neigbors to 0*/
for(j=0;j<(*frame).n;j++)
  {
  (*framep).nb[j] = 0;
  (*framep).nbtrixx[j] = 0; 
  }

/*ACTUAL CREATION OF NEIGHBOR LISTS*/

/*Determine amount of subdivisions for cell (nxcell, nycell, nzcell):*/
{
double det,l1,l2,l3,t1,t2,t3;
position v1,v2,v3;
        det = gg[0][0]*gg[1][1]*gg[2][ 2] +
            gg[1][ 0]*gg[2][1]*gg[0][ 2] +
            gg[2][ 0]*gg[0][1]*gg[1][ 2] -
            gg[0][ 0]*gg[2][1]*gg[1][ 2] -
            gg[1][ 0]*gg[0][1]*gg[2][ 2] -
            gg[2][ 0]*gg[1][1]*gg[0][ 2];
v1.x = gg[0][0]*(*frame).a.x; v1.y = gg[1][0]*(*frame).a.x; v1.z = gg[2][0]*(*frame).a.x;
v2.x = gg[0][1]*(*frame).a.y; v2.y = gg[1][1]*(*frame).a.y; v2.z = gg[2][1]*(*frame).a.y;
v3.x = gg[0][2]*(*frame).a.z; v3.y = gg[1][2]*(*frame).a.z; v3.z = gg[2][2]*(*frame).a.z;

CROSS(v1,v2,q3);
t3 = 1/(LENGTH(q3));
l3 = t3*det*(*frame).a.x*(*frame).a.y*(*frame).a.z;
q3.x *= t3; q3.y *= t3; q3.z *= t3; 

CROSS(v1,v3,q2);
t2 = 1/(LENGTH(q2));
l2 = t2*det*(*frame).a.x*(*frame).a.y*(*frame).a.z;
q2.x *= t2; q2.y *= t2; q2.z *= t2;

CROSS(v2,v3,q1);
t1 = 1/(LENGTH(q1));
l1 = t1*det*(*frame).a.x*(*frame).a.y*(*frame).a.z;
q1.x *= t1; q1.y *= t1; q1.z *= t1;

nxcell = (int) floor(l1 / clsgnecut);
if (nxcell > MAXCELLS) nxcell = MAXCELLS;
if (nxcell < 3) nxcell = 3;

nycell = (int) floor(l2 / clsgnecut);
if (nycell > MAXCELLS) nycell = MAXCELLS;
if (nycell < 3) nycell = 3;

nzcell = (int) floor(l3 / clsgnecut);
if (nzcell > MAXCELLS) nzcell = MAXCELLS;
if (nzcell < 3) nzcell = 3;

}

//printf("%d %d %d cell vects\n",nxcell,nycell,nzcell);

/*Now zero out linked list*/
for(i=0;i<nxcell;i++)
  for(j=0;j<nycell;j++)
    for(k=0;k<nzcell;k++)
      (*framep).cellptr[i][j][k] = -1;  

for(i=0;i<(*frame).n;i++) 
  {
  int cellx, celly, cellz;
/*For each atom i, determine which cell it is in: cellx, celly, cellz*/
  cellx = floor( nxcell * (*frame).p[i].x / (*frame).a.x ); 
  celly = floor( nycell * (*frame).p[i].y / (*frame).a.y );
  cellz = floor( nzcell * (*frame).p[i].z / (*frame).a.z );
  while (cellx < 0) cellx += nxcell;
  while (celly < 0) celly += nycell;
  while (cellz < 0) cellz += nzcell;
  cellx = cellx%nxcell; celly = celly%nycell; cellz=cellz%nzcell;
/*Now add the atom to the list*/
  (*framep).linkptr[i] = (*framep).cellptr[cellx][celly][cellz];
  (*framep).cellptr[cellx][celly][cellz] = i;
  }

for(nx=0;nx<nxcell;nx++)
  for(ny=0;ny<nycell;ny++)
    for(nz=0;nz<nzcell;nz++)
      for(mx=-1;mx<2;mx++)
        for(my=-1;my<2;my++)
          for(mz=-1;mz<2;mz++)
            {
            int ux,uy,uz,atom1,atom2;
            ux = (nx+mx+nxcell)%nxcell;
            uy = (ny+my+nycell)%nycell;
            uz = (nz+mz+nzcell)%nzcell; 
            atom2 = (*framep).cellptr[ux][uy][uz];
            while (atom2!=-1) { 
            atom1 = (*framep).cellptr[nx][ny][nz];
            while (atom1!=-1 ) {
                 double r;
                 position g;
                 int temp;
                 if (atom1 < atom2) 
                   {
                    CELLDISTNO(((*frame).p[atom1]),((*frame).p[atom2]),((*frame).a),r);
                    if (r < clsrmin) clsrmin = r;
                    if (r <= clsgnecut)
                      {
                      add_to_neighbor_list(frame,framep,r,atom1,atom2,err);
                      }
                   }
                 atom1 = (*framep).linkptr[atom1];
                 } /*End while */
                 atom2 = (*framep).linkptr[atom2];
                 } /*End while */
            } /*End for loops */ 
/*END OF ACTUAL CREATION OF NEIGHBOR LISTS*/


/*Do sorting of neighbor lists in order of increasing distance*/
{
        int globalMaxNeighb = 0;

for(j=0;j<(*frame).n;j++) {
        if ((*framep).nb[j] > globalMaxNeighb) globalMaxNeighb = (*framep).nb[j];
      {
      int x,tempi;
      double tempd;
      position temp;

      for(k=0;k<(*framep).nb[j];k++)
for(x=k+1;x<(*framep).nb[j];x++)
  if ((*framep).r[j][k] > (*framep).r[j][x]) {
    SWTX((*framep).r[j][k],(*framep).r[j][x],tempd);
    SWTX((*framep).nbind[j][k],(*framep).nbind[j][x],tempi);
    SWTX((*framep).grad[j][k],(*framep).grad[j][x],temp);
  };
      };

      (*framep).nbrho[j] = 0;
      (*framep).nbtri[j] = 0;
      for(k=0;k<(*framep).nb[j];k++) {
if ((*framep).r[j][k] <= clsrhonecut) (*framep).nbrho[j]++;
if ((*framep).r[j][k] <= clstricut1) (*framep).nbtri[j]++;
}
  };
//printf("GLOBALMAXNEIGHB = %d\n",globalMaxNeighb);
}

/* *************************************************************************************************** */
/* END OF NEIGHBOR LISTING CODE */
}

