#include <math.h>
#include <stdio.h>
#include <stdlib.h> //has qsort()
#include <complex.h>
#include "defines.h"
#include "nrutil.h"

void prmst38c();
void eselfgeneral();
void radelmgeneralsp(double *,double *,double *,int,int);
double evalPair(int,double,double *,double *,int);
void eselfgeneral2();
void dsyevr();
void zheevd();
void lenoskymeam_();

/*Notes on code revision
atomtype 0 = Si
atomtype 1 = hydrogen
silicons come first in list 0...nsilicon-1
hydrogens are nsilicon...nat

hydrogens have 1 orbital 1 electron
silicons have 4 orbitals 4 electrons
*/

#define ZEROTFORCE FALSE /*If TRUE Turns off finite temperature corrections to forces*/

#define SFLAG 1 /*Controls some of the printing*/

#define FERMIBROADENING 300 //1000 // 100 /*Fermi-Dirac BROADENING in Kelvin */

/*Number of eigenvalues to compute = MIN(dimension of matrix * EIGENFRACT + EXTRAEIGEN, dimension of matrix)*/
/*This saves time in the eigenvalue routine*/
#define EIGENFRACT 1.0
#define EXTRAEIGEN  10
#define NC(nc) ((int) ( floor(nc*EIGENFRACT +EXTRAEIGEN) > nc ? nc : floor(nc*EIGENFRACT + EXTRAEIGEN) ))

#define MMATOM 100 /*Maximum number of atoms within cell, also set CLSMAXATOM and CLSMAXNEIGHB in defines.h*/
#define MAXKPT 6000 //1728 /*12x12x12 grid with offset*/
#define KPTEIGENSAVE 4 /*2x2x2 grid with offset -- save this many sets of eigenvectors to avoid recomputation*/

//Disabled in Max's call so he can set them
#define INTERACTIVE_CELLSEARCH 2 /*Controls loop for finding nearest image -- may need to increase in highly nonorthogonal or small cells
                                   also the correct value of this parameter differs for MEAM and TB */
#define INTERACTIVE_TB_OR_MEAM 0 //0 TB 1 Si-MEAM 2 Ti-MEAM

#define INTERACTIVE_NX 2
#define INTERACTIVE_NY 2
#define INTERACTIVE_NZ 2
#define INTERACTIVE_OFFSET_X 1
#define INTERACTIVE_OFFSET_Y 1
#define INTERACTIVE_OFFSET_Z 1
#define INTERACTIVE_DO_KPT 1

static int CELLSEARCH=INTERACTIVE_CELLSEARCH;
static int nx = INTERACTIVE_NX;
static int ny = INTERACTIVE_NY;
static int nz = INTERACTIVE_NZ;
static int offset_x = INTERACTIVE_OFFSET_X;
static int offset_y = INTERACTIVE_OFFSET_Y;
static int offset_z = INTERACTIVE_OFFSET_Z;
static int do_kpt = INTERACTIVE_DO_KPT;

#define SCALE 0.03 // Note 0.03 is good for 16 atom cell and seems OK for 2 atom cell as well

void lenoskytb_();
void main_tbnorth();
void yfdocclocal();
void forcediagonalizeg();
void forcediagonalizeg_complex();
void slatercoupling();
void gammacoupling();
void gammamat();
void gammaenergy();
void gammaenergy_Kpt();
double totalenergy();
double func();
double dfunc();
double newclscall();
void setgg();
double MHM_interface_function();
void generate_Kpts();
int compare_kpts();

typedef struct {
	double eigen,gfocc;
	int weights,which_kpt,index_within_kpt,idummy;
} eigenRecord ;

/*Notes:
totalenergy calls gammaenergy (TB energy)  and pairenergy (pair pot energy)
gammaenergy calls gammamat (for setting up tb matrix), 
                  forcediagonalizeg, (to diagonalize matrix) 
                  yfdocclocal (fermi-dirac occupancy), 
                  then gammamat (for derivs of matrix)
gammamat calls gammacoupling
gammacoupling calls slatercoupling

func and dfunc are interface routines for numerical recipes optimizer

lenoskytb_() is interface routine by which TB routines may be called from external code

main_tbnorth() is a seperate driver routine that can be called as main()
 
newclscall() is an interface routine for the fitting code to call tight-binding
*/

extern double passtbcut();
extern int isnorth();
extern int isxtal(); 

#define Zprintf if (printflag==1 && SFLAG) printf

static double eigenenergy[4*MMATOM*MAXKPT]; /*All the eigenvalues at the gamma point */
static int gweights[4*MMATOM*MAXKPT]; /*weighting*/
static int which_kpt[4*MMATOM*MAXKPT];
static int index_within_kpt[4*MMATOM*MAXKPT];

static double gtbcut_tbnorth; /*This is the tight--binding cutoff, not currently used by the code.*/
static int gusetb_tbnorth; /*This flag tells whether tight binding is to be used*/
static int gusenorth; /*This flag tells whether nonorth. tb is present rather than orthogonal; 
			this version of code does not support nonorthogonal*/
static int gusextal; /*flag to use xtal field terms, or not*/ 
static int printflag = 1; /* This flag controls whether printed output is on or off, in Zprintf*/

extern int passgstride();

static int gdopair = 1;
static int gstride = 4;

static int gncalls = 0;
static int docalc  = 1;

static double fermitemp = FERMIBROADENING; /*Fermi temperature in eV, is set by code also*/
static double gfocc[4*MMATOM*MAXKPT]; /*Array for returning occupancies*/
static double gtoten; /*Contains total energy as a global variable for bands*/
static double gpairen; /*passing variable for pairenergy returns to main()*/
static int freeat[1000]; /*Free atom labels for hydrogen*/

static double efermi;

static int gnat,gnsilicon;
static position glatt[3];
static double gg[3][3],gg_inv[3][3];
static double glattforce[3][3];


/*PUT THIS BACK IN TO COMPILE TBNORTH*/
//void main(void)
//{
//main_tbnorth();
//}

void main_tbnorth(void)
{
/*Read in a standard file from clsman, generated by ``out filename.''*/
/*from stdin.  Printing can be controlled by the value of printflag, see above...*/
int jjj,i,nat,idummy;
double ax,ay,az,fdummy,en;
position p[4*MMATOM],force[4*MMATOM];
position latt[3];

printf("GAMMA POINT only tight binding code\n");

Zprintf("Reading spline potential coeff.cls\n");
prmst38c(); /* Reads potential and does setup for tblink.c */
gtbcut_tbnorth = passtbcut(&gusetb_tbnorth);
gusenorth = isnorth(); 
gusextal = isxtal(); 
gstride = passgstride();
gusenorth = 0;
gusextal = 0;
gstride = 4;
gusetb_tbnorth = 1;
gtbcut_tbnorth = 5.24;
Zprintf("gusetb_tbnorth == %d\n",gusetb_tbnorth);
Zprintf("gtbcut_tbnorth == %lf\n", gtbcut_tbnorth);
Zprintf("gusenorth == %d\n",gusenorth);
Zprintf("gusextal == %d\n",gusextal);
Zprintf("gstride == %d\n",gstride);

if (docalc) {
scanf("%d %d",&nat,&gnsilicon);
Zprintf("##nat = %d\n",nat);

scanf("%lf %lf %lf",&latt[0].x,&latt[0].y,&latt[0].z);
scanf("%lf %lf %lf",&latt[1].x,&latt[1].y,&latt[1].z);
scanf("%lf %lf %lf",&latt[2].x,&latt[2].y,&latt[2].z);
Zprintf("##lattice vectors = \n##%lf %lf %lf\n##%lf %lf %lf\n##%lf %lf %lf\n",latt[0].x,latt[0].y,latt[0].z,latt[1].x,latt[1].y,latt[1].z,latt[2].x,latt[2].y,latt[2].z);

for(i=0;i<nat;i++)
  {
    scanf("%lf %lf %lf",&(p[i ].x),&(p[i ].y),&(p[i ].z));
    Zprintf("##%lf %lf %lf\n",p[i].x,p[i].y,p[i].z);
  }

printf("Calling totalenergy\n");
for(i=0;i<1;i++)
en = totalenergy(p,latt,nat,force,INTERACTIVE_TB_OR_MEAM);

if (1 == 0) {
	int count = 0;
	int jj,kk;
	double errfract, errfractmax = 0.0;
	double errabsolute,errabsolutemax = 0.0;
	double detot = 0.0;
	int seed=10000;
	position dx[MMATOM],dlatt[3];
	srand(seed);
#define RVAL (   (double)rand() / ((double)(RAND_MAX)+(double)(1)) )
#define RANDSTEPSIZE 0.000
#define CELLSTEPSIZE 0.001
#define MAXRANDITER 10
#define RANDOMDIRECTIONEACHTIME 0
for(jj=0;jj<MAXRANDITER;jj++) {
	int ii;
	position p1[MMATOM];
	position f0[MMATOM];
	position f1[MMATOM];
	double en1, en2;
	double de1, de2;

	//if (RANDOMDIRECTIONEACHTIME || count == 0)
	for(ii=0;ii<nat;ii++) {
		if (RANDOMDIRECTIONEACHTIME || jj==0 ) {
		dx[ii].x = (2.0 - 2*RVAL)*RANDSTEPSIZE;
		dx[ii].y = (2.0 - 2*RVAL)*RANDSTEPSIZE;
		dx[ii].z = (2.0 - 2*RVAL)*RANDSTEPSIZE;
		for(kk=0;kk<3;kk++)
		{
			dlatt[kk].x = (2.0 - 2*RVAL)*CELLSTEPSIZE;
			dlatt[kk].y = (2.0 - 2*RVAL)*CELLSTEPSIZE;
			dlatt[kk].z = (2.0 - 2*RVAL)*CELLSTEPSIZE;
		}
		}
	  VECT_ADD(p[ii],dx[ii],p1[ii]);
	}

	de2 = 0.0;
	en1 = totalenergy(p,latt,nat,f0,INTERACTIVE_TB_OR_MEAM);

	printf("LATTFORCE\n");
	printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][0],glattforce[1][0],glattforce[2][0]);
	printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][1],glattforce[1][1],glattforce[2][1]);
	printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][2],glattforce[1][2],glattforce[2][2]);

	de2 = de2 + glattforce[0][0] * dlatt[0].x;
	de2 = de2 + glattforce[0][1] * dlatt[1].x;
	de2 = de2 + glattforce[0][2] * dlatt[2].x;
	de2 = de2 + glattforce[1][0] * dlatt[0].y;
	de2 = de2 + glattforce[1][1] * dlatt[1].y;
	de2 = de2 + glattforce[1][2] * dlatt[2].y;
	de2 = de2 + glattforce[2][0] * dlatt[0].z;
	de2 = de2 + glattforce[2][1] * dlatt[1].z;
	de2 = de2 + glattforce[2][2] * dlatt[2].z;

	//Before applying the change I need to REVERT coordinates
	for(kk=0;kk<nat;kk++) {
		position temp;
		REVERT(p[kk],temp);
		p[kk] = temp;
		REVERT(p1[kk],temp);
	    p1[kk] = temp;
	}

    VECT_ADD(latt[0],dlatt[0],latt[0]);
    VECT_ADD(latt[1],dlatt[1],latt[1]);
    VECT_ADD(latt[2],dlatt[2],latt[2]);

    glatt[0] = latt[0];
    glatt[1] = latt[1];
    glatt[2] = latt[2];
    setgg();
    //After applying the change I need to CONVERT coordinates
    	for(kk=0;kk<nat;kk++) {
    		position temp;
    		CONVERT(p[kk],temp);
    		p[kk] = temp;
    		CONVERT(p1[kk],temp);
    		p1[kk] = temp;
    	}

//	de2 = de2 + 0.5*(DOT(dlatt[0],glattforce[0])) + DOT(dlatt[1],glattforce[1]) + DOT(dlatt[2],glattforce[2])));
	en2 = totalenergy(p1,latt,nat,f1,INTERACTIVE_TB_OR_MEAM);

	printf("LATTFORCE\n");
	printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][0],glattforce[1][0],glattforce[2][0]);
	printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][1],glattforce[1][1],glattforce[2][1]);
	printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][2],glattforce[1][2],glattforce[2][2]);

	de2 = de2 + glattforce[0][0] * dlatt[0].x;
		de2 = de2 + glattforce[0][1] * dlatt[1].x;
		de2 = de2 + glattforce[0][2] * dlatt[2].x;
		de2 = de2 + glattforce[1][0] * dlatt[0].y;
		de2 = de2 + glattforce[1][1] * dlatt[1].y;
		de2 = de2 + glattforce[1][2] * dlatt[2].y;
		de2 = de2 + glattforce[2][0] * dlatt[0].z;
		de2 = de2 + glattforce[2][1] * dlatt[1].z;
		de2 = de2 + glattforce[2][2] * dlatt[2].z;

		de2 = -de2 / 2.0;
//	de2 = de2 + 0.5*(DOT(dlatt[0],glattforce[0])) + DOT(dlatt[1],glattforce[1]) + DOT(dlatt[2],glattforce[2])));

	printf("en1 = %le\n",en1);

	de1 = en1-en2;

	detot = detot + fabs(de1);
	//Now compute estimated change in energy
	for(ii=0;ii<nat;ii++)
		de2 = de2 + 0.5*DOT(dx[ii],f0[ii]) + 0.5*DOT(dx[ii],f1[ii]);

	errfract = fabs(de1 / de2 - 1.0);
	errabsolute = fabs(de1 - de2);

//		printf("%3.18le %3.18le\n",de1,de2);
		printf("%3.18le\n",errfract);
//		printf("%3.18le\n",errabsolute);

	if (errfract > errfractmax)
		errfractmax = errfract;
	if (errabsolute > errabsolutemax)
		errabsolutemax = errabsolute;

		for(ii=0;ii<nat;ii++)
			p[ii] = p1[ii];

		count++;
}
printf("Greatest error fraction = %le\n",errfractmax);
printf("Greatest absolute error = %le\n",errabsolutemax);
printf("Typical change in energy = %le\n",detot / (double)MAXRANDITER);
}


#if (SFLAG)
printf("Returned energy eV = %3.18lf\n",en);

printf("\n");
printf("Forces for atoms  eV/Angstrom (fx, fy, fz):\n");
for(i=0;i<nat;i++)
  printf("%d (%3.18lf %3.18lf %3.18lf)\n",i,force[i].x, force[i].y, force[i].z); 
printf("LATTFORCE\n");
printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][0],glattforce[1][0],glattforce[2][0]);
printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][1],glattforce[1][1],glattforce[2][1]);
printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][2],glattforce[1][2],glattforce[2][2]);
#endif

#if (1==0)
printf("OPTIMIZATION:\n");
glatt[0] = latt[0];
glatt[1] = latt[1];
glatt[2] = latt[2];
setgg();
gnat = nat;
{
double pq[3*MMATOM+10],xi[3*MMATOM+10];
int i,j,k;
double fret;
int iter;
int index;

index = 1;
for(i=0;i<nat;i++)
  {
	position temp;
	REVERT(p[i],temp);
  pq[3*i+1] = temp.x;
  pq[3*i+2] = temp.y;
  pq[3*i+3] = temp.z;
  index = index + 3;
  }

pq[index] = SCALE*glatt[0].x; index++;
pq[index] = SCALE*glatt[0].y; index++;
pq[index] = SCALE*glatt[0].z; index++;
pq[index] = SCALE*glatt[1].x; index++;
pq[index] = SCALE*glatt[1].y; index++;
pq[index] = SCALE*glatt[1].z; index++;
pq[index] = SCALE*glatt[2].x; index++;
pq[index] = SCALE*glatt[2].y; index++;
pq[index] = SCALE*glatt[2].z; index++;

//NOTE DISABLED WITH ZERO FUNCTION CALLS --> change to 500 for relax
frprmn(pq,3*nat+9,0.0,&iter,&fret,5000,0.0);

index = 1;
for(i=0;i<nat;i++)
  {
	position temp;
  temp.x = pq[3*i+1];
  temp.y = pq[3*i+2];
  temp.z = pq[3*i+3];
  CONVERT(temp,p[i]);
  index = index +3;
  }

latt[0].x = pq[index]/SCALE; index++;
latt[0].y = pq[index]/SCALE; index++;
latt[0].z = pq[index]/SCALE; index++;
latt[1].x = pq[index]/SCALE; index++;
latt[1].y = pq[index]/SCALE; index++;
latt[1].z = pq[index]/SCALE; index++;
latt[2].x = pq[index]/SCALE; index++;
latt[2].y = pq[index]/SCALE; index++;
latt[2].z = pq[index]/SCALE; index++;
}

en = totalenergy(p,latt,nat,force,INTERACTIVE_TB_OR_MEAM);

#if (SFLAG)
printf("Returned energy eV = %3.18lf\n",en);
printf("Relative to isolated charge neutral atoms\n");
printf("\n");
printf("Forces for atoms  eV/Angstrom (fx, fy, fz):\n");
for(i=0;i<nat;i++)
  printf("%d (%3.18lf %3.18lf %3.18lf)\n",i,force[i].x, force[i].y, force[i].z); 
printf("LATTFORCE\n");
printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][0],glattforce[1][0],glattforce[2][0]);
printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][1],glattforce[1][1],glattforce[2][1]);
printf("%3.18lf %3.18lf %3.18lf\n",glattforce[0][2],glattforce[1][2],glattforce[2][2]);
printf("LATT\n");
printf("%lf %lf %lf\n",glatt[0].x,glatt[0].y,glatt[0].z);
printf("%lf %lf %lf\n",glatt[1].x,glatt[1].y,glatt[1].z);
printf("%lf %lf %lf\n",glatt[2].x,glatt[2].y,glatt[2].z);
if (0==1) {
FILE *C;

C = fopen("coord.out","w");

printf("Coords for atoms Angstrom (x,y,z): \n");
for(i=0;i<nat;i++)
  {
  printf("%d (%3.18lf %3.18lf %3.18lf)\n",i,p[i].x,p[i].y,p[i].z);
  fprintf(C,"%3.18lf %3.18lf %3.18lf\n",p[i].x,p[i].y,p[i].z);
  }

fclose(C);
}

//For each atom compute distance from (0,0,0) and norm of force
{
double fnorm[10000], rdist[10000];
FILE *fout;

fout = fopen("F.out","w");

for(i=0;i<nat;i++)
  {
    double dx, dy, dz;
    dx = p[i].x; 
    dy = p[i].y; 
    dz = p[i].z; 
    if (dx > ax/2) dx = dx - ax;
    if (dy > ay/2) dy = dy - ay;
    if (dz > az/2) dz = dz - az;
    rdist[i] = sqrt(dx*dx + dy*dy + dz*dz);
    fnorm[i] = sqrt(force[i].x*force[i].x + force[i].y*force[i].y + force[i].z*force[i].z);
  }

for(i=0;i<nat;i++)
  fprintf(fout, "%3.18lf %3.18lf\n",rdist[i], fnorm[i]);
  //fprintf(fout, "%3.18lf %3.18lf\n",rdist[i], fnorm[i]*fnorm[i]*rdist[i]*rdist[i]); 

fclose(fout);

}

#endif
printf("Total number of calls = %d\n",gncalls);
#endif
}
}

void setgg()
{
	gg[0][0] = glatt[0].x; gg[1][0] = glatt[0].y; gg[2][0] = glatt[0].z;
	gg[0][1] = glatt[1].x; gg[1][1] = glatt[1].y; gg[2][1] = glatt[1].z;
	gg[0][2] = glatt[2].x; gg[1][2] = glatt[2].y; gg[2][2] = glatt[2].z;

	{ //FORMULA CHECKED TL 7.9.10 -- also transposing GG is a bad idea, it is transposed correctly.
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

double newpairenergy(coord,latt,nat,force)
position *coord,*latt,*force;
int nat;
/*Given atoms with reciprocal lattice latt, number of atoms nat, at coordinates coord, return energy from
  pair potential, with forces in array force[] */
{
    double pairen;
    int j,m;
    position p1, p2;
    int atomtypej,atomtypem;
    position diff;
    double dist;
    double y,der;
//Alireza start
    double t1,zarib=5.0,alpha=0.5*1.0/(0.7*0.7);
        //double t1,zarib=-2.4,alpha=0.5*1.0/(0.7*0.7);
        //double t1,zarib=-1.8,alpha=0.5*1.0/(0.7*0.7);
//Alireza stop 


    for(j=0;j<3;j++)
    	for(m=0;m<3;m++)
    		glattforce[j][m] = 0.0;

    for(j=0;j<nat;j++)
      force[j].x = force[j].y = force[j].z = 0.0;

    pairen = 0.0;
    for(j=0;j<nat;j++)
      for(m=j; m<nat; m++)
      {
    	p1 = coord[j];
    	p2 = coord[m];

    /*Here want section to compute i,j,k for nearest image, rather than
      running through choices*/
    {
          position d,c,cc,p;
    	  int ii,jj,kk;
    	  position diff_save;
    	  double dist_save;
    	  dist_save = 1e30;
    	  VECT_SUBTRACT(p1,p2,d);
    	  REVERT(d,c);
    	  c.x = floor(c.x);
    	  c.y = floor(c.y);
    	  c.z = floor(c.z);
    	  CONVERT(c,cc);
    	  VECT_SUBTRACT(d,cc,p);
    	  for (ii=-CELLSEARCH;ii<=CELLSEARCH;ii++)
    		  for (jj=-CELLSEARCH;jj<=CELLSEARCH;jj++)
    			  for (kk=-CELLSEARCH;kk<=CELLSEARCH;kk++)
    			  {
    				  position chi,chi2;
    				  chi.x = ii; chi.y = jj; chi.z = kk;
    				  CONVERT(chi,chi2);
    				  VECT_ADD(p,chi2,diff);;
    				  dist = LENGTH(diff);
    				  if (dist < gtbcut_tbnorth && dist > 0.001) {
    					  atomtypej =0;
    					         if (j>=gnsilicon) atomtypej=1;
    					         atomtypem=0;
    					         if (m>=gnsilicon) atomtypem=1;
    					         evalPair(atomtypej+atomtypem, dist, &y, &der, 1);
						 if (m==j)
 							{
							y = y/2.0;
							der = der/2.0;
         						}
                             //Alireza begins
                             if(atomtypej+atomtypem==2) {
                                   t1=zarib*exp(-alpha*dist*dist);
                                   y+=t1;
                                   der+=-2.0*alpha*dist*t1;
                             }
                             //Alireza Ends

                                                 force[j].x += -der * (diff.x) / dist;
    					         force[j].y += -der * (diff.y) / dist;
    					         force[j].z += -der * (diff.z) / dist;
    					         force[m].x += der * (diff.x) / dist;
    					         force[m].y += der * (diff.y) / dist;
    					         force[m].z += der * (diff.z) / dist;
    					         {
    					        	 position ptemp,ptemp2;
    					        	 ptemp.x = der*diff.x;
    					        	 ptemp.y = der*diff.y;
    					        	 ptemp.z = der*diff.z;
    					        	 REVERT(ptemp,ptemp2);
    					         glattforce[0][0] += (diff.x) / dist * ptemp2.x;
    					         glattforce[0][1] += (diff.x) / dist * ptemp2.y;
    					         glattforce[0][2] += (diff.x) / dist * ptemp2.z;
    					         glattforce[1][0] += (diff.y) / dist * ptemp2.x;
    					         glattforce[1][1] += (diff.y) / dist * ptemp2.y;
    					         glattforce[1][2] += (diff.y) / dist * ptemp2.z;
    					         glattforce[2][0] += (diff.z) / dist * ptemp2.x;
    					         glattforce[2][1] += (diff.z) / dist * ptemp2.y;
    					         glattforce[2][2] += (diff.z) / dist * ptemp2.z;
    					         }
    					     pairen += y;
    				  }  } } }
	return(pairen);
}

void forcediagonalizeg(m,eigen,nc,eigenv)
/*Does the eigenvalue problem, with arguments having the same format as the
routine diagonalize.  m is H and is real, nc by nc
matrix.  Eigen is the list of eigenvalues, and eigenv is the list of eigenvectors.*/
double **m;
double *eigen,**eigenv;
int nc;
{
int i,j,k;

int itype,n,lda,ldb,lwork,ierr=0;
static double a[16*MMATOM*MMATOM];
static double work[16*MMATOM*MMATOM + 200*MMATOM];
static double z[16*MMATOM*MMATOM];
static int iwork[16*MMATOM*MMATOM + 200*MMATOM];
static int isuppz[8*MMATOM];
char range;
char jobz,uplo;
static int errcount = 0;
double vl,vu;
int il,iu,ldz;
double abstol;
int m_param;
int liwork;

k=0;
for(i=1;i<=nc;i++)
  for(j=1;j<=nc;j++)
    {
if (j>i)
     a[k] = m[i][j];
else a[k] = m[j][i];
     k++;
    }

itype = 1; /*Problem type*/
jobz = 'V' ; /*Set to 'N' if eigenvectors are not needed.*/
range = 'I' ; /*Find all eigenvectors*/
uplo = 'U' ;
n = lda = ldb =ldz=nc ;
vl = vu = 0.0;
il = 1;
iu = NC(nc);
abstol = 0;
m_param = 0;
/*ierr does not need to be set on entry*/     
/*work is the workspace array, already set*/
/*eigen is also an output*/
lwork = 16*MMATOM*MMATOM+200*MMATOM - 1;
liwork = 16*MMATOM*MMATOM + 200*MMATOM - 1; 

#ifdef HAVE_MKL
dsyevr(&jobz,&range,&uplo,&n,a,&lda,&vl,&vu,&il,&iu,&abstol,&m_param,eigen,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&ierr);
#else
dsyevr_(&jobz,&range,&uplo,&n,a,&lda,&vl,&vu,&il,&iu,&abstol,&m_param,eigen,z,&ldz,isuppz,work,&lwork,iwork,&liwork,&ierr);
#endif

if (ierr != 0  && errcount < 250)
  {   
   printf("TBNORTH WARNING: ierr , errcount == %d   %d from dsygv diagonalize\n",ierr,errcount);
   printf("Will not print this message when errcount exceeds 250\n");
   errcount++ ;
  }

/*Now z array is returned with eigenvectors, normalization and transposition to be determined*/
k=0;
for(i=1;i<=NC(nc);i++)
  for(j=1;j<=nc;j++)
    {
       eigenv[i][j] = z[k];
       k++;
    };

}

void forcediagonalizeg_complex(m,eigen,nc,eigenv)
/*Does the eigenvalue problem, with arguments having the same format as the
routine forcediagonalizeg.  m is H and is Hermitian, nc by nc
matrix.  Eigen is the list of eigenvalues, and eigenv is the list of eigenvectors.*/
double complex **m;
double *eigen;
double complex **eigenv;
int nc;
{
int i,j,k;

int n,lda,lwork, info,lrwork,liwork;
static double complex a[16*MMATOM*MMATOM];
static double complex work[16*MMATOM*MMATOM + 200*MMATOM];
static double rwork[16*MMATOM*MMATOM + 200*MMATOM];
static int iwork[16*MMATOM*MMATOM + 200*MMATOM];
char jobz,uplo;
static int errcount = 0;

k=0;
for(i=1;i<=nc;i++)
  for(j=1;j<=nc;j++)
    {
if (j>i)
     a[k] = m[i][j];
else a[k] = m[j][i];
     k++;
    }

jobz = 'V' ; /*Set to 'N' if eigenvectors are not needed.*/
uplo = 'U' ;
n = lda  = nc;
lwork = 16*MMATOM*MMATOM+200*MMATOM - 1;
lrwork = 16*MMATOM*MMATOM+200*MMATOM - 1;
liwork = 16*MMATOM*MMATOM+200*MMATOM - 1;

//Choose one of these two calls
//zheev(&jobz,&uplo,&n,a,&lda,eigen,work,&lwork,rwork,&info);

//NOTE zheevd is faster:in tests 2 routines agree to 14 digits
#ifdef HAVE_MKL
zheevd(&jobz,&uplo,&n,a,&lda,eigen,work,&lwork,rwork,&lrwork,iwork,&liwork,&info);
#else
zheevd_(&jobz,&uplo,&n,a,&lda,eigen,work,&lwork,rwork,&lrwork,iwork,&liwork,&info);
#endif

if (info != 0 && errcount < 250)
  {
   printf("TBNORTH WARNING: info == %d from dsygv diagonalize\n",info);
   printf("Will not print this message when errcount exceeds 250\n");
   errcount++ ;
  }

/*Now a array is returned with eigenvectors*/
k=0;
for(i=1;i<=NC(nc);i++)
  for(j=1;j<=nc;j++)
    {
       eigenv[i][j] = a[k];
       k++;
    };

}

double newclscall(p,g,a,natom,nsilicon)
position *p, *g, a;
int natom,nsilicon;
{
/*This routine imitates the routine clscall that was used to link to clsman.  Call to this routine is hardwired
  into the program fastfit.c, now compiled as fastfit*/
double energy;
position latt[3];
int i;
latt[0].x = a.x ;
latt[1].y = a.y ;
latt[2].z = a.z ;
latt[0].y = latt[0].z = latt[1].x = latt[1].z = latt[2].x = latt[2].y = 0.0;

gnsilicon = nsilicon;

/*
printf("Newclscall natom == %d\n",natom);
for(i=0;i<natom;i++)
  printf("Newclscall %d (%lf %lf %lf)\n",i,p[i].x,p[i].y,p[i].z);
*/

gtbcut_tbnorth = passtbcut(&gusetb_tbnorth);
gusenorth = isnorth();
gusextal = isxtal();
gstride = passgstride();

gusenorth = 0;
gusextal = 0;
gstride = 4;
gusetb_tbnorth = 1;
gtbcut_tbnorth = 5.24;


fermitemp = FERMIBROADENING; /*Fixed temperature of 1000 kelvin for now*/

gdopair = 0; /*Don't try to evaluate the pair potential*/
energy  = totalenergy(p,latt,natom,g,INTERACTIVE_TB_OR_MEAM);
gdopair = 1; /*Reset flag to its default value*/

/*
printf("Newclscall returning energy = %lf\n",energy);
for(i=0;i<natom;i++)
  printf("Newclscall, f[%d].x, f[%d].y, f[%d].z = %lf %lf %lf\n",i,i,i,g[i].x,g[i].y,g[i].z);
*/

return(energy);

}

/*Interface for Goedecker global minima finder code
 *
 * NOTES ON HOW TO CALL
 * latt0 = first lattice vector
 * latt1 = second lattice vector
 * latt2 = third lattice vector
 * lattforce0 = derivatives with respect to first lattice vector
 * lattforce1 = derivatives with respect to second lattice vector
 * lattforce2 = derivatives with respect to third lattice vector
 *
 *NOTE THAT the "alatt" argument has been removed.
 *
 * Positions (rxyz) and forces (ff) are in cartesian coordinates (Angstroms) and NOT direct coordinates
 *
 * NOTE that parameter CELLSEARCH (see line 31 above) must be set large enough so that all nearby images are summed,
 * however if it is set too large the performance of the code will suffer.
 *
 * parameter tb_or_meam needs to be set to choose potential type
 * 0 = silicon/hydrogen tight binding
 * 1 = silicon meam (all atoms treated as silicon and nsilicon ignored)
 * 2 = titanium meam (all atoms treated as titanium and nsilicon ignored)
 * NOTE that for the MEAM options cellsearch still has to be set, and the correct value is different for each tb_or_meam option
 *
 * K-point TB instructions:
 * do_kptIn = 1 for k-point tight binding, 0 for gamma point
 * CELLSEARCH needs to be set appropriately regardless of k-points
 * nxIn, nyIn, nzIn are the k-point grids along the three reciprocal lattice vectors
 * offset_xIn, offset_yIn, offset_zIn are set to 1 to offset from the origin along each direction, and 0 to not offset
 * The code uses inversion symmetry to reduce the number of k-points (e.g. a 2x2x2 mesh with offsets = 1 will have 4 kpoints)
 * MAXKPT is the maximum number of actual k-points (after taking the inversion symmetry into account)
 * KPTEIGENSAVE is the number of k-points to save the eigenvalues to avoid recomputation.  It affects only computational time, and not accuracy, but
 * if it is set too high there could be memory issues.  If this happens you can either lower MMATOM (the maximum number of atoms) or KPTEIGENSAVE.
 * These routines use the C library complex.h so there may be some issues on older C compilers that don't support this extension.
 * */
void lenoskytb_(nat,rxyz,ff,erxyz,count_md,nsilicon,latt0,latt1,latt2,lattforce0,lattforce1,lattforce2,tb_or_meam,cellsearch,nxIn,nyIn,nzIn,offset_xIn,offset_yIn,offset_zIn,do_kptIn)
int *nat;
double *rxyz;
double *ff;
double *erxyz;
double *count_md; 
int *nsilicon; /*NOTE silicon atoms must be first*/
double *latt0,*latt1,*latt2,*lattforce0,*lattforce1,*lattforce2;
int *cellsearch;
int *tb_or_meam;
int *nxIn;
int *nyIn;
int *nzIn;
int *offset_xIn;
int *offset_yIn;
int *offset_zIn;
int *do_kptIn;
{

CELLSEARCH=*cellsearch;
static int firstcall = 1;
int i,j,k;

if(*nat>100) {
    printf("Due to improper large static arrays nat is limited to 100\n");
    printf("You are running with nat=%d\n",nat);
    printf("Stopping in C program in function lenoskytb_\n");
    exit(0);
}

if (firstcall == 1)
  {
  docalc = 0;
  main_tbnorth();
  firstcall = 0;
  };
*count_md  = *count_md + 1.0;
gnat = *nat;
gnsilicon = *nsilicon;

nx = *nxIn;
ny = *nyIn;
nz = *nzIn;
offset_x = *offset_xIn;
offset_y = *offset_yIn;
offset_z = *offset_zIn;
do_kpt = *do_kptIn;

  {
  /*PRINT SOME WARNINGS*/
  if (*nat > MMATOM && (*tb_or_meam)==0) {
                   printf("WARNING Number of atoms = %d is greater than MMATOM = %d\n",*nat,MMATOM);
                   printf("THIS MAY CAUSE FAILURE OR UNPREDICTABLE BEHAVIOR\n");
                   printf("INCREASE MMATOM in tbnorth.c and RECOMPILE\n");
                   }
  if (*nat > CLSMAXATOM*CELLSEARCH*CELLSEARCH*CELLSEARCH && (*tb_or_meam)>0) {
                   printf("WARNING Number of atoms = %d is greater than CLSMAXATOM*(CELLSEARCH^3) = %d\n",*nat,CLSMAXATOM*CELLSEARCH*CELLSEARCH*CELLSEARCH);
                   printf("THIS MAY CAUSE FAILURE OR UNPREDICTABLE BEHAVIOR\n");
                   printf("INCREASE CLSMAXATOM in defines.h and RECOMPILE\n");
                   }
  if (*nat > CLSMAXNEIGHB) {
                   printf("WARNING Number of atoms = %d is greater than CLSMAXNEIGHB = %d\n",*nat,CLSMAXNEIGHB);
                   printf("THIS MAY CAUSE FAILURE OR UNPREDICTABLE BEHAVIOR\n");
                   printf("INCREASE CLSMAXNEIGHB in defines.h and RECOMPILE\n");
                   }
  if (*nsilicon > *nat ) printf("WARNING nsilicon = %d is greater than number of atoms =%d\n",*nsilicon,*nat);
  if (*nat == 0) printf("WARNING lenoskytb called with zero atoms\n");
  if (*nat < 0) printf("WARNING lenoskytb called with negative number of atoms\n");
  }
*erxyz = MHM_interface_function(rxyz-1,ff-1,latt0-1,latt1-1,latt2-1,lattforce0-1,lattforce1-1,lattforce2-1,*tb_or_meam);
}

/*******************************************************************************/
/*func() and dfunc() are interface for frprmn, numerical recipes minimizer code*/
double func(p)
double *p;
{
int i,j,k,index;
position coord[MMATOM], force[MMATOM];


index = 3*gnat+1;
glatt[0].x = p[index]/SCALE; index++;
glatt[0].y = p[index]/SCALE; index++;
glatt[0].z = p[index]/SCALE; index++;
glatt[1].x = p[index]/SCALE; index++;
glatt[1].y = p[index]/SCALE; index++;
glatt[1].z = p[index]/SCALE; index++;
glatt[2].x = p[index]/SCALE; index++;
glatt[2].y = p[index]/SCALE; index++;
glatt[2].z = p[index]/SCALE; index++;
setgg();

index = 1;
for(i=0;i<gnat;i++)
  {
	position temp;
  temp.x = p[3*i+1];
  temp.y = p[3*i+2];
  temp.z = p[3*i+3];
  CONVERT(temp,coord[i]);
  index = index+3;
  }



return(totalenergy(coord,glatt,gnat,force,INTERACTIVE_TB_OR_MEAM));
}

double dfunc(p,xi)
double *p,*xi;
{
	int i,j,k,index;
	position coord[MMATOM], force[MMATOM];
    double energy;

    index = 3*gnat+1;
    	glatt[0].x = p[index]/SCALE; index++;
    	glatt[0].y = p[index]/SCALE; index++;
    	glatt[0].z = p[index]/SCALE; index++;
    	glatt[1].x = p[index]/SCALE; index++;
    	glatt[1].y = p[index]/SCALE; index++;
    	glatt[1].z = p[index]/SCALE; index++;
    	glatt[2].x = p[index]/SCALE; index++;
    	glatt[2].y = p[index]/SCALE; index++;
    	glatt[2].z = p[index]/SCALE; index++;
    	setgg();

	index = 1;
	for(i=0;i<gnat;i++)
	  {
		position temp;
	  temp.x = p[3*i+1];
	  temp.y = p[3*i+2];
	  temp.z = p[3*i+3];
	  CONVERT(temp,coord[i]);
	  index = index+3;
	  }

	energy = totalenergy(coord,glatt,gnat,force,INTERACTIVE_TB_OR_MEAM);

	index = 1;
	for(i=0;i<gnat;i++)
	  {
		position temp,temp2;
	  temp.x = force[i].x;
	  temp.y = force[i].y;
	  temp.z = force[i].z;

	  CONVERT(temp,temp2);
	  xi[3*i+1] = temp2.x;
	  xi[3*i+2] = temp2.y;
	  xi[3*i+3] = temp2.z;

	  index = index+3;
	  }
	xi[index] = -glattforce[0][0]/SCALE; index++;
	xi[index] = -glattforce[1][0]/SCALE; index++;
	xi[index] = -glattforce[2][0]/SCALE; index++;
	xi[index] = -glattforce[0][1]/SCALE; index++;
	xi[index] = -glattforce[1][1]/SCALE; index++;
	xi[index] = -glattforce[2][1]/SCALE; index++;
	xi[index] = -glattforce[0][2]/SCALE; index++;
	xi[index] = -glattforce[1][2]/SCALE; index++;
	xi[index] = -glattforce[2][2]/SCALE; index++;

	return(energy);

}

double MHM_interface_function(p,xi,latt0,latt1,latt2,lattforce0,lattforce1,lattforce2,tb_or_meam)
double *p,*xi;
double *latt0,*latt1,*latt2,*lattforce0,*lattforce1,*lattforce2;
int tb_or_meam;
{
	int i,j,k,index;
	position coord[MMATOM], force[MMATOM];
    double energy;

    glatt[0].x = latt0[1];
    glatt[0].y = latt0[2];
    glatt[0].z = latt0[3];
    glatt[1].x = latt1[1];
    glatt[1].y = latt1[2];
    glatt[1].z = latt1[3];
    glatt[2].x = latt2[1];
    glatt[2].y = latt2[2];
    glatt[2].z = latt2[3];

    	setgg();

	index = 1;
	for(i=0;i<gnat;i++)
	  {
		position temp;
	  temp.x = p[3*i+1];
	  temp.y = p[3*i+2];
	  temp.z = p[3*i+3];
	  coord[i] = temp;
	  index = index+3;
	  }

	energy = totalenergy(coord,glatt,gnat,force,tb_or_meam);

	index = 1;
	for(i=0;i<gnat;i++)
	  {
	  xi[3*i+1] = force[i].x;
	  xi[3*i+2] = force[i].y;
	  xi[3*i+3] = force[i].z;
	  index = index+3;
	  }

	lattforce0[1] = glattforce[0][0];
	lattforce0[2] = glattforce[1][0];
	lattforce0[3] = glattforce[2][0];
	lattforce1[1] = glattforce[0][1];
	lattforce1[2] = glattforce[1][1];
	lattforce1[3] = glattforce[2][1];
	lattforce2[1] = glattforce[0][2];
	lattforce2[2] = glattforce[1][2];
	lattforce2[3] = glattforce[2][2];

	return(energy);

}
/****************************************************************************/


double totalenergy(coord,latt,nat,force,tb_or_meam)
position *coord,*latt,*force;
int nat;
int tb_or_meam;
//int nx,ny,nz,offset_x,offset_y,offset_z,do_kpt;
/*Calls gammaenergy, at the gamma point only, and calls pair energy*/
{
int i,ii;
double toten,pairen;

gncalls++;

toten = pairen = 0.0;

glatt[0] = latt[0]; glatt[1] = latt[1]; glatt[2] = latt[2];
setgg();

//Force coord to be in primitive cell only
for(ii=0;ii<nat;ii++)
{
	position primitive;
	REVERT(coord[ii],primitive);
	primitive.x = primitive.x - floor(primitive.x);
	primitive.y = primitive.y - floor(primitive.y);
	primitive.z = primitive.z - floor(primitive.z);
	CONVERT(primitive,coord[ii]);
}

if (tb_or_meam == 0) {

if (gdopair)
  pairen = newpairenergy(coord,latt,nat,force);
else
  {
  pairen = 0.0;
  for(i=0;i<nat;i++)
    {force[i].x = force[i].y = force[i].z = 0.0; }
  };

gpairen = pairen;

if (gusetb_tbnorth!=0) {

	if (do_kpt==0)
gammaenergy(coord,latt,nat,force);
	else
gammaenergy_Kpt(coord,latt,nat,force,nx,ny,nz,offset_x,offset_y,offset_z);

toten += gtoten ; /*returned by fdocclocal routine*/

{double eself[4],es; eselfgeneral(eself); 
pairen -= gnsilicon*(2*eself[0] + 2*eself[1]);
eselfgeneral2(&es);
pairen -= (nat-gnsilicon)*es;
 };

  }; /*end if gusetb_tbnorth!=0 */

//printf("Totalenergy = %3.18lf\n",toten+pairen);
return(toten+pairen); }
else
{
	double a[3];
	double erxyz;
	int count_md = 0;
	int i;
	//passgg(gg);
	a[0] = a[1] = a[2] = 1.0;

	//REVERT coordinates
	for(i=0;i<nat;i++)
	{
		position temp;
		REVERT(coord[i],temp);
		coord[i]=temp;
	}
	lenoskymeam_(&nat,a,((double *)coord),((double *)force),&erxyz,&count_md,glattforce,gg,CELLSEARCH,tb_or_meam);
	//CONVERT coordinates
		for(i=0;i<nat;i++)
		{
			position temp;
			CONVERT(coord[i],temp);
			coord[i]=temp;
		}
//printf("CALLING LENOSKY MEAM\n, erxyz = %lf\n",erxyz);

	//Now need to set glattforce somehow

  return (erxyz);
}
}

//int compare_eigenvalues(const eigenRecord *a, const eigenRecord *b)
int compare_eigenvalues(const void *a_void, const void *b_void)
{
    const eigenRecord *a=(const eigenRecord *)a_void;
    const eigenRecord *b=(const eigenRecord *)b_void;
	int retVal;
	double v;

	v = (*a).eigen - (*b).eigen;
	if (v==0) retVal = 0;
	else if (v>0) retVal = 1;
	else if (v<0) retVal = -1;
	return(retVal);
}

//int compare_kpts(const eigenRecord *a, const eigenRecord *b)
int compare_kpts(const void *a_void, const void *b_void)
{
    const eigenRecord *a=(const eigenRecord *)a_void;
    const eigenRecord *b=(const eigenRecord *)b_void;
	int v;
	v = (*a).which_kpt - (*b).which_kpt;
	if (v == 0)
		v = (*a).index_within_kpt - (*b).index_within_kpt;
	return (v);
}

void sort_eigen(n,sort_type,eigen,weights,which_kpt,index_within_kpt,gfocc)
//On entry n is number of eigenvalues, sort_type is 0 to sort by eigenvalues and 1 to sort by which_kpt and index_with_within_kpt
//All arrays start at 0, order is ascending
int n,sort_type;
double *eigen,*gfocc;
int *weights,*which_kpt,*index_within_kpt;
{
	static eigenRecord myArray[4*MMATOM*MAXKPT];
	int i;
	for(i=0;i<n;i++)
	{
		myArray[i].eigen = eigen[i];
		//printf("i, eigen[i], myArray[i].eigen = %d %lf %lf\n",i,eigen[i],myArray[i].eigen);
		myArray[i].weights = weights[i];
		myArray[i].which_kpt = which_kpt[i];
		myArray[i].index_within_kpt = index_within_kpt[i];
		myArray[i].gfocc = gfocc[i];
	}
	if (sort_type==0)
	  qsort(myArray, n, sizeof(eigenRecord), compare_eigenvalues);
	else
	  qsort(myArray, n, sizeof(eigenRecord), compare_kpts);
	for(i=0;i<n;i++)
		{
			eigen[i] = myArray[i].eigen;
			//printf("i, eigen[i], myArray[i].eigen = %d %lf %lf\n",i,eigen[i],myArray[i].eigen);
			weights[i] = myArray[i].weights;
			which_kpt[i] = myArray[i].which_kpt;
			index_within_kpt[i] = myArray[i].index_within_kpt;
			gfocc[i] = myArray[i].gfocc;
		}
}



void generate_Kpts(Kpts,Weights,numK,SumWeight,nx,ny,nz,offset_x,offset_y,offset_z)
position *Kpts;
int *Weights; //both of these start at 0 index
int *numK;
int *SumWeight; //Total weight
int nx,ny,nz; //Number of kpts in x,y,z
int offset_x,offset_y,offset_z; //These are 0 or 1 flags
{
int i,j,k,ii;
position K;
int found;
//loop over points
//if point has k <--> -k symmetry add it with weight = 1
//otherwise add it with weight = 2 only if reflected point is not already in the list
*numK = 0;
for(i=0;i<nx;i++)
	for(j=0;j<ny;j++)
		for(k=0;k<nz;k++)
		{
			K.x = (i+0.5*offset_x - 0.5 *(double) nx) / (double) nx;
			K.y = (j+0.5*offset_y - 0.5 *(double) ny) / (double) ny;
			K.z = (k+0.5*offset_z - 0.5 *(double) nz) / (double) nz;
			found = 0;
			for(ii=0;ii<*numK;ii++)
			{
				position R,diff;
				R = Kpts[ii];
				R.x = -R.x;
				R.y = -R.y;
				R.z = -R.z;
				VECT_SUBTRACT(K,R,diff);
				diff.x = diff.x - floor(diff.x);
				diff.y = diff.y - floor(diff.y);
				diff.z = diff.z - floor(diff.z);
				if (fabs(diff.x) + fabs(diff.y) + fabs(diff.z) < 1e-9)
				{
					found = 1;
					Weights[ii] = Weights[ii] + 1;
					break;
				}
			}
			if (found == 0) {
				Kpts[*numK] = K;
				Weights[*numK] = 1;
				*numK = *numK + 1;
			}
		}
for(ii=0;ii<*numK;ii++)
{
	position K,R;
	K = Kpts[ii];
	K.x = K.x * 2.0 * M_PI;
	K.y = K.y * 2.0 * M_PI;
	K.z = K.z * 2.0 * M_PI;
	//REVERT(K,R); //Convert to reciprocal lattice vector

	//Maybe REVERT is not quite right ... we need to multiply by transpose of gg_inv

	CONVERT_TO_RECIPROCAL_LATTICE(K,R);

	if (1==0) {
		//Test section
		double dot1,dot2,dot3;
		printf("Testing Kpoint %d\n",ii);
		dot1 = DOT(glatt[0],R);
		dot2 = DOT(glatt[1],R);
		dot3 = DOT(glatt[2],R);
		printf("Testing Kpoint %d: %lf %lf %lf\n",ii,dot1,dot2,dot3);
	}


	Kpts[ii] = R;
}
//Return sum of weights
	*SumWeight = 0;
	for(ii=0;ii<*numK;ii++)
		*SumWeight = *SumWeight+ Weights[ii];
}

void gammaenergy(coord,latt,nat,force)
position *coord,*latt;
int nat;
position *force;
/*Returns the band structure contribution to the energy/atom at the gamma point
  for nat atoms with orthogonal positions in cell given by coord,
  for a coordinate system with lattice vectors latt.  Nonorthogonal case is assumed.
  ADDS to the forces in array force[], which are already present from pair term.
*/
{
static double **drem,**rem,**remLatt, **eigenv,**rho,**rho2;
static double *ggocc;

static int firstcall = 1;
int i,j,k,nmatrix;

static double complex **remC, **remLattC,**eigensave;
position K;
int do_K;
int SumWeight = 1;

if (firstcall == 1) {
drem = dmatrix(1,gstride*MMATOM,1,MMATOM);
rem = dmatrix(1,gstride*MMATOM,1,gstride*MMATOM); 
remLatt = dmatrix(1,gstride*MMATOM,1,gstride*MMATOM*9);
remC = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
remLattC = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM*9);
eigenv = dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
rho = dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
rho2 = dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
ggocc = dvector(0,gstride*MMATOM-1);
firstcall = 0;
  }

nmatrix = gstride * gnsilicon + (nat-gnsilicon);

/*Find Hamiltonian and diagonalize it*/
K.x = K.y = K.z = 0.0;
do_K = 0;
gammamat(coord,latt,rem ,remLatt, nat,0,remC,remLattC,K,do_K); /* flag = 0 for H matrix, 1 for S, if present*/

forcediagonalizeg(rem,eigenenergy,nmatrix,eigenv);

for(i=1;i<=nmatrix;i++)
  for(j=1;j<=nmatrix;j++)
    rho[i][j] = 0.0;

if (1==0) {
//OUTPUT EIGENENERGIES

FILE *E;
E = fopen("eigen.out","w");

for(i=1;i<=NC(nmatrix);i++)
  fprintf(E,"%3.18lf\n",eigenenergy[i-1] - efermi);


fclose(E);

}


/*Fermi diract distribution: returns total energy*/
for(i=0;i<NC(nmatrix);i++)
	gweights[i] = 1;
yfdocclocal(&gtoten,fermitemp,eigenenergy,NC(nmatrix),nmatrix,gfocc,gweights,SumWeight);

/*Create density matrix rho*/
/*Need to include extra terms.  We have gfocc[i-1]...gfocc[NC(gstride*nat)-1]
also fermitemp = temp in Kelvin so "beta" = 1 / kT = 11604 / fermitemp
ggocc maximum is 2.0 */

if (ZEROTFORCE==FALSE) {
double term1,  term2,beta;
double totaln;
term1 = 0;
term2 = 0;
beta = 11604.0 / fermitemp;
totaln =0;
for(i=1;i<=NC(nmatrix);i++)
  term1 += gfocc[i-1] * -1.0*(1.0 - gfocc[i-1] / 2.0);
for(i=1;i<=NC(nmatrix);i++)
  term2 += eigenenergy[i-1] * gfocc[i-1] * -1.0*(1.0 - gfocc[i-1] / 2.0);
for(i=1;i<=NC(nmatrix);i++)
  ggocc[i-1] = gfocc[i-1] + beta * gfocc[i-1] * -1.0*(1.0 - gfocc[i-1] / 2.0) * (eigenenergy[i-1] - 1.0*term2 / term1); 
for(i=1;i<=NC(nmatrix);i++)
  totaln += gfocc[i-1];
/*TEST: printf("Total number of electrons= %lf\n",totaln);
*/
}

for(i=1;i<=NC(nmatrix);i++) /*loop was over occupied only*/
  {
  /*For ith eigenvector eigenv[i][1..gstride*nat], compute rho and rho2 */
  /*Use factors (gfocc[i-1] / 2.0) for eigenv[i][]  */
  /*Rho is quadratic in eigenfunctions, so contains this factor (probability) once*/
  double sterm;
  sterm = 0.0;

/*
  for(j=1;j<=nmatrix;j++)
    sterm += 0.5 * eigenv[i][j] * eigenv[i][j] ;
  sterm = 1.0 / sterm ; 

printf("Sterm = %lf\n",sterm);
*/

sterm = 2.0; /*Well behaved eigensolver*/

  for(j=1;j<=nmatrix;j++)
    { 
    double reg1, reg2;
#if (ZEROTFORCE)
    reg1 = eigenv[i][j]*sterm*(gfocc[i-1] / 2.0) ;
#else
    reg1 = eigenv[i][j]*sterm*(ggocc[i-1] / 2.0) ;
#endif
    for(k=j;k<=nmatrix;k++) 
      rho[j][k] += reg1 * eigenv[i][k];
    }
  }

/*
for(i=1;i<=NC(nmatrix);i++)
  for(k=1;k<=nmatrix;k++)
    printf("eigenv[%d][%d] = %lf\n",i,k,eigenv[i][k]);
printf("nmatrix,NC(matrix) = %d %d\n",nmatrix,NC(nmatrix));
*/
/*Take trace rho * H for Hellman Feynman theorem*/
/*Note that other triangle of rho, and rho2 is missing and we do sum with
  factors of two to compensate.*/
for(i=1;i<=3;i++) {
	K.x = K.y = K.z = 0.0;
	do_K = 0;
gammamat(coord,latt,rem,remLatt,nat,i,remC,remLattC,K,do_K);

for(j=1;j<=nmatrix;j++)
  for(k=j;k<=nmatrix;k++) {
int a1,a2;
double temp;
/*MARK*/
if ((k-1)/4>=gnsilicon) a2 = gnsilicon + (k-1) - gnsilicon*4; else a2 = (k-1)/4;
if ((j-1)/4>=gnsilicon) a1 = gnsilicon + (j-1) - gnsilicon*4; else a1 = (j-1)/4; 

/*OLD CODE
a1 = (j-1)/gstride;
a2 = (k-1)/gstride;
*/
{
	int ii,jj;

if (j!=k || 1) {
	       int prefactor;
	       if (j!=k)
		 prefactor = 2;
	       else
 		 prefactor = 1;
	       temp = rho[j][k]*rem[j][k]*2.0 ;
	       for (ii=0;ii<3;ii++)
	    	   for (jj=0;jj<3;jj++)
	    		   glattforce[ii][jj] += rho[j][k]*remLatt[j][(k-1)*9+1+3*ii+jj] * prefactor;
}
//if (j==k) temp = 0.0;
if (i==1) {force[a1].x -= temp; force[a2].x += temp;};
if (i==2) {force[a1].y -= temp; force[a2].y += temp;};
if (i==3) {force[a1].z -= temp; force[a2].z += temp;};
}

}
} /*i loop*/

}

void gammamat(coord,latt,rem,remLatt,nat,flag2,remC,remLattC,K,do_K)
position *coord,*latt;
double **rem,**remLatt;
int nat,flag2;
double complex **remC;
double complex **remLattC;
position K;
int do_K;
/*Given ''nat'' atoms at orthogonal positions ''coord''
within the unit cell, this routine
sets up the entire matrix to be diagonalized at the gamma point.
Lattice translation vectors are also given in orthogonal coordinates.
Entire coupling matrix is returned in rem, must be dimensioned
properly on call.  Note hgen and dhgen are arrays used to store values of splines and their derivatives for each atom pair.*/
{
int ii,jj,i,j;
int nmatrix,indexi,indexj,atomtypei,atomtypej,nsizei,nsizej;
static double **rex,**rexLatt;
static double complex **rexC,**rexLattC;
//static double rex[4][4];
//static double rexLatt[4][40];
static int firstcall = 1;
double eself[4];
double tvect[4]; /*crystal field terms on ith atom*/

if (firstcall==1)
    {
	rex = dmatrix(0,gstride-1,0,gstride-1);
	rexLatt = dmatrix(0,gstride-1,0,gstride*9-1);
	rexC = Complex_dmatrix(0,gstride-1,0,gstride-1);
	rexLattC = Complex_dmatrix(0,gstride-1,0,gstride*9-1);
    }

nmatrix = gstride * gnsilicon + (nat - gnsilicon);

for(i=1;i<=nmatrix;i++)
  for(j=1;j<=nmatrix;j++)
  {
	  int ii,jj;
    rem[i][j] = 0.0;
    if (do_K==0) {
    for(ii=0;ii<3;ii++)
    	for(jj=0;jj<3;jj++)
    		remLatt[i][(j-1)*9+1+3*ii+jj] = 0.0;
    } else
    {
    	remC[i][j] = 0.0;
        for(ii=0;ii<3;ii++)
        	for(jj=0;jj<3;jj++)
        		remLattC[i][(j-1)*9+1+3*ii+jj] = 0.0;
        }
  }

for(i=0;i<nat;i++)
  for(j=i;j<nat;j++) {
atomtypei = 0;
if (i>=gnsilicon) atomtypei=1;

atomtypej = 0;
if (j>=gnsilicon) atomtypej=1;

nsizei = 4-3*atomtypei;
nsizej = 4-3*atomtypej;

if (atomtypei==0) indexi=i*gstride; else indexi=gnsilicon*gstride+(i-gnsilicon);
if (atomtypej==0) indexj=j*gstride; else indexj=gnsilicon*gstride+(j-gnsilicon);

gammacoupling(i,j,coord[i],coord[j],latt,rex,rexLatt,flag2,tvect,atomtypei,atomtypej,rexC,rexLattC,K,do_K);

if (do_K==0) {
  for(ii=0;ii<nsizei;ii++)
    for(jj=0;jj<nsizej;jj++)
    {
    	int i,j;
      rem[1+indexi+ii][1+indexj+jj] += rex[ii][jj];
      for(i=0;i<3;i++)
    	  for(j=0;j<3;j++)
    		  remLatt[1+indexi+ii][(1+indexj+jj)*9-8+3*i+j] += rexLatt[ii][jj*9 + 3*i+j];
    }
  } else
  {
    for(ii=0;ii<nsizei;ii++)
      for(jj=0;jj<nsizej;jj++)
      {
      	int i,j;
        remC[1+indexi+ii][1+indexj+jj] += rexC[ii][jj];
        for(i=0;i<3;i++)
      	  for(j=0;j<3;j++)
      		  remLattC[1+indexi+ii][(1+indexj+jj)*9-8+3*i+j] += rexLattC[ii][jj*9 + 3*i+j];
      }

    }
  } //i,j loops

    eselfgeneral(eself);

    if (do_K == 0) {

    	if (flag2==0) {
    for(i=0;i<gnsilicon;i++)
      for(j=0;j<gstride;j++)
        rem[1+i*gstride+j][1+i*gstride+j] += eself[j];
    {
    double es;
    eselfgeneral2(&es);
    for(i=gnsilicon;i<nat;i++)
      rem[1+gnsilicon*gstride+(i-gnsilicon)][1+gnsilicon*gstride+(i-gnsilicon)] += es;
    }
    } // flag2

    } else {

    	if (flag2==0) {
  for(i=0;i<gnsilicon;i++)
    for(j=0;j<gstride;j++)
      remC[1+i*gstride+j][1+i*gstride+j] += eself[j];
  {
  double es;
  eselfgeneral2(&es);
  for(i=gnsilicon;i<nat;i++)
    remC[1+gnsilicon*gstride+(i-gnsilicon)][1+gnsilicon*gstride+(i-gnsilicon)] += es;
  }
    } //flag2

    }

firstcall = 0;
}

/*Find couplings between atomi and atomj*/
void gammacoupling(atomi,atomj,p1,p2,latt,rem,remLatt,flag2,tvect,atomtypei,atomtypej,remC,remLattC,K,do_K)
position p1,p2,latt[3];
double **rem,**remLatt;
int flag2;
double *tvect;
int atomi,atomj,atomtypei,atomtypej;
double complex **remC;
double complex **remLattC;
position K;
int do_K;
/* p1 and p2 represent positions of atoms within the cell.
To speed this up, I have assumed that latt (coordinate vectors) form a diagonal matrix
*/
{
int i,j,k,ii,jj;
static int firstcall = 1;
position diff;
double dist;
static double **mc,**mc2,**mc3;
static double *vect;
static double *hgen,*dhgen; 
double complex phaseFactor;

if (firstcall==1) {mc = dmatrix(0,gstride-1,0,gstride-1);
                   mc2 = dmatrix(0,gstride-1,0,gstride-1);
                   mc3 = dmatrix(0,gstride-1,0,gstride-1);
                   vect = dvector(0,gstride-1);
                   hgen = dvector(0,4);dhgen=dvector(0,4);};

for(i=0;i<=gstride-1;i++)
  for(j=0;j<=gstride-1;j++)
    {
	  int ii,jj;
    rem[i][j] = 0.0; 
    remC[i][j] = 0.0;
    tvect[i] = 0.0;
    for(ii=0;ii<3;ii++)
    	for(jj=0;jj<3;jj++) {
          remLatt[i][(j)*9+3*ii+jj] = 0.0;
          remLattC[i][(j)*9+3*ii+jj] = 0.0;
    	}
    }

{
	  position d,c,cc,p;
	  int ii,jj,kk;
	  position diff_save;
	  double dist_save;
	  dist_save = 1e30;
	      	  VECT_SUBTRACT(p1,p2,d);
	      	  REVERT(d,c);
	      	  c.x = floor(c.x);
	      	  c.y = floor(c.y);
	      	  c.z = floor(c.z);
	      	  CONVERT(c,cc);
	      	  VECT_SUBTRACT(d,cc,p);
//          printf("CELLSEARCH %d \n",CELLSEARCH);
	  for (ii=-CELLSEARCH;ii<=CELLSEARCH;ii++)
		  for (jj=-CELLSEARCH;jj<=CELLSEARCH;jj++)
			  for (kk=-CELLSEARCH;kk<=CELLSEARCH;kk++)
			  {
				  position chi,chi2;
				  chi.x = ii; chi.y = jj; chi.z = kk;
				  CONVERT(chi,chi2);
				  VECT_ADD(p,chi2,diff);;
				  dist = LENGTH(diff);
				  if (dist <= gtbcut_tbnorth && dist>=0.0001)
				  {
					  int ii,jj;
					  diff.x /= dist;
					  diff.y /= dist;
					  diff.z /= dist;
					       radelmgeneralsp(&dist,hgen,dhgen,atomtypei,atomtypej);
					       if (do_K)
					    	   phaseFactor = cexp(I*dist*DOT(K,diff));
					       if (flag2 == 1)
					       {
					    	   slatercoupling(diff,dist,mc2,vect,2,hgen,dhgen);
					    	   slatercoupling(diff,dist,mc3,vect,3,hgen,dhgen);
					    	   slatercoupling(diff,dist,mc,vect,flag2,hgen,dhgen);
					       }
					       else
					    	   slatercoupling(diff,dist,mc,vect,flag2,hgen,dhgen);
					       for(ii=0;ii<=gstride-1;ii++)
					         for(jj=0;jj<=gstride-1;jj++)
					         {
					        	 if (do_K)
					        		 remC[ii][jj]+=mc[ii][jj] * phaseFactor;
					        	 else
					                 rem[ii][jj] += mc[ii][jj];

					           if (flag2==1) {
					           position ptemp,ptemp2;
					           ptemp.x = ptemp.y = ptemp.z = 0.0;
					               ptemp.x = diff.x*dist;
					        	   ptemp.y = diff.y*dist;
					        	   ptemp.z = diff.z*dist;
					           REVERT(ptemp,ptemp2);
					           if (do_K) {
					        	   remLattC[ii][jj*9+3*0+0] += mc[ii][jj]* ptemp2.x*phaseFactor;
					        	   remLattC[ii][jj*9+3*0+1] += mc[ii][jj]* ptemp2.y *phaseFactor;
					        	   remLattC[ii][jj*9+3*0+2] += mc[ii][jj]* ptemp2.z *phaseFactor;
					        	   remLattC[ii][jj*9+3*1+0] += mc2[ii][jj]* ptemp2.x *phaseFactor;
					        	   remLattC[ii][jj*9+3*1+1] += mc2[ii][jj]* ptemp2.y *phaseFactor;
					        	   remLattC[ii][jj*9+3*1+2] += mc2[ii][jj]* ptemp2.z *phaseFactor;
					        	   remLattC[ii][jj*9+3*2+0] += mc3[ii][jj]* ptemp2.x *phaseFactor;
					        	   remLattC[ii][jj*9+3*2+1] += mc3[ii][jj]* ptemp2.y *phaseFactor;
					        	   remLattC[ii][jj*9+3*2+2] += mc3[ii][jj]* ptemp2.z *phaseFactor;
					           } else {
					           remLatt[ii][jj*9+3*0+0] += mc[ii][jj]* ptemp2.x;
					           remLatt[ii][jj*9+3*0+1] += mc[ii][jj]* ptemp2.y;
					           remLatt[ii][jj*9+3*0+2] += mc[ii][jj]* ptemp2.z;
					           remLatt[ii][jj*9+3*1+0] += mc2[ii][jj]* ptemp2.x;
					           remLatt[ii][jj*9+3*1+1] += mc2[ii][jj]* ptemp2.y;
					           remLatt[ii][jj*9+3*1+2] += mc2[ii][jj]* ptemp2.z;
					           remLatt[ii][jj*9+3*2+0] += mc3[ii][jj]* ptemp2.x;
					           remLatt[ii][jj*9+3*2+1] += mc3[ii][jj]* ptemp2.y;
					           remLatt[ii][jj*9+3*2+2] += mc3[ii][jj]* ptemp2.z;
					           }
					           }
					         }
					       for(ii=0;ii<=gstride-1;ii++)
					         tvect[ii] += vect[ii] ;
				  }	 } }
firstcall = 0;
}

void slatercoupling(u,r,mat,vect,flag2,hgen,dhgen)
position u;
double r;
double **mat; /*Runs from 0..3,0..3*/
double *vect; /*Again from 0 to 3*/
int flag2;
double *hgen,*dhgen;
/*Generates the matrix of couplings for the program, if u is unit vector between atoms
r is distance, returns coupling matrix mat; 
eselfgeneral must be treated as a seperate case. Flag is set to zero for H, 1 for S
      
flag2 is 0 for normal case, 1 for x derivative, 2 for y, 3 for z, where
u = (x,y,z) / sqrt(x^2+y^2+z^2)...

returns on--site energies in vect, or their derivatives, if xtal field terms are set.
Note that right now the onsite energies are added elsewhere and vect is returned as zero.
 */
{
double ess,esx,esy,esz,exs,eys,ezs,exx,eyy,ezz,exy,eyz,exz,eyx,ezy,ezx;
double v0,v1,v2,v3;
double isss,ipps,ippp,disss,dipps,dippp;
double hsss,hsps,hpps,hppp;
double dhsss,dhsps,dhpps,dhppp;
int j,i;
static int firstcall = 1;

v0 = v1 = v2= v3 = 0.0;

hsss=hgen[0];
hsps=hgen[1];
hpps=hgen[2];
hppp=hgen[3];

if (flag2!=0) {
dhsss = dhgen[0];
dhsps = dhgen[1];
dhpps = dhgen[2];
dhppp = dhgen[3];
  }

if (flag2==0) {
ess = hsss;
esx = hsps*u.x;
esy = hsps*u.y;
esz = hsps*u.z;
exx = hpps*u.x*u.x + (1-u.x*u.x)*hppp;
eyy = hpps*u.y*u.y + (1-u.y*u.y)*hppp;
ezz = hpps*u.z*u.z + (1-u.z*u.z)*hppp;
exy = (u.x*u.y) * (hpps - hppp);
eyz = (u.y*u.z) * (hpps - hppp);
exz = (u.x*u.z) * (hpps - hppp);
}

/*A1*/
  
if (flag2!=0) {
/*Use dr/dx = x / r = u.x*/
/*Assume fv (type position) is unit vector along x,y, or z depending on flag2*/
double drdx;
position fv,g;
if(flag2 == 1) {fv.x=1.0; fv.y=fv.z=0.0;};
if(flag2 == 2) {fv.y=1.0; fv.x=fv.z=0.0;};
if(flag2 == 3) {fv.z=1.0; fv.x=fv.y=0.0;};
drdx = u.x * fv.x + u.y * fv.y + u.z * fv.z ;
g.x = fv.x / r - drdx * u.x / r ; /*derivative of u.x in direction of fv*/
g.y = fv.y / r - drdx * u.y / r ;
g.z = fv.z / r - drdx * u.z / r ;

/*Calculate derivatives with respect to r dependence of matrix elements*/
ess = drdx*dhsss;
esx = drdx*dhsps*u.x;
esy = drdx*dhsps*u.y;
esz = drdx*dhsps*u.z;
exx = drdx*(dhpps*u.x*u.x + (1-u.x*u.x)*dhppp);
eyy = drdx*(dhpps*u.y*u.y + (1-u.y*u.y)*dhppp);
ezz = drdx*(dhpps*u.z*u.z + (1-u.z*u.z)*dhppp);
exy = drdx* (u.x*u.y) * (dhpps - dhppp);
eyz = drdx* (u.y*u.z) * (dhpps - dhppp);
exz = drdx* (u.x*u.z) * (dhpps - dhppp);

/*A2*/

/*Now calculate derivatives with respect to x,y,z appearing in u = (x,y,z) / Sqrt(...) */
/*Here g is the derivative of u with respect to either x,y,z (for flag2 = 1, 2, or 3) */
esx += hsps*g.x;
esy += hsps*g.y;
esz += hsps*g.z;
exx += 2.0*hpps*g.x*u.x + (-2.0*g.x*u.x)*hppp;
eyy += 2.0*hpps*g.y*u.y + (-2.0*g.y*u.y)*hppp;
ezz += 2.0*hpps*g.z*u.z + (-2.0*g.z*u.z)*hppp;
exy += (g.x*u.y + u.x*g.y) * (hpps - hppp);
eyz += (g.y*u.z + u.y*g.z) * (hpps - hppp);
exz += (g.x*u.z + u.x*g.z) * (hpps - hppp);
}

eyx = exy;
ezy = eyz;
ezx = exz;
exs = -esx;
eys = -esy;
ezs = -esz;

mat[0][0] = ess; mat[0][1] = exs; mat[0][2] = eys; mat[0][3] = ezs;
mat[1][0] = esx; mat[1][1] = exx; mat[1][2] = eyx; mat[1][3] = ezx; 
mat[2][0] = esy; mat[2][1] = exy; mat[2][2] = eyy; mat[2][3] = ezy; 
mat[3][0] = esz; mat[3][1] = exz; mat[3][2] = eyz; mat[3][3] = ezz; 

vect[0] = v0; vect[1] = v1; vect[2] = v2; vect[3] = v3;

firstcall = 0;

}

/*Returns Fermi Dirac distribution*/
void yfdocclocal(ebond,temp,eval,natob,nel,focc,weights,SumWeight)
double *ebond;
double temp;
double *eval;
int natob;
int nel;
double *focc;
int *weights;
int SumWeight;
/*
Natob is number of eigenvalues, *eval is array of eigenvalues, 
  *focc is array of returned eigenvalue occupancies,
  nel is number of electrons (not occupied states, but electrons)
  *ebond is returned energy
  temp is temperature in degrees Kelvin
Eigenvalues must be sorted on entry.  Number of electrons nel
must be even.
*/

#define FDMAXIT 20000 /*Maximum number of NR iterations*/
#define FDEPSOCC 1e-9 /*Allowed error in number of electrons*/
{
double rsid,sumocc,sumder,delta,beta,fdexp,fder;      
int it,i,nocc;

for(i=0;i<natob;i++)
  focc[i] = 0.0;

//for(i=0;i<natob;i++)
//  printf("weights[i] = %d\n", weights[i]);
//printf("SumWeight = %d\n",SumWeight);

//nocc = nel/2;
//for(i=0;i<nocc;i++)
//  focc[i] = 2.0;
//if (nocc*2!=nel)
//  focc[nocc]=1.0;
//efermi = eval[nocc-1];
//if (nocc*2!=nel)
//  efermi=eval[nocc];

{//Figure out efermi to start
	double total_occ;
	nocc = SumWeight * nel/2;
	total_occ = 0.0;
	for(i=0;i<natob;i++)
	{
		total_occ += weights[i];
		if (total_occ > nocc) {
			efermi = eval[i];
			nocc = i;
			break;
		}
	}
	for(i=0;i<nocc;i++)
		focc[i] = 2.0;
}


beta = 11604.0 / temp; /*Conversion of temp into eV*/

it = 0; 
do {
  it++;
  if (it > FDMAXIT) printf("WARNING:Iteration count exceeded in yfdocc\n");
  sumocc = sumder =0.0;
  for(i=0;i<natob;i++) {
    delta = eval[i] - efermi;
    if (delta > 0.0) {
      fdexp = exp (-delta*beta);
      focc[i] = weights[i] / (double) SumWeight * 2.0*fdexp/(1.0+fdexp);
      sumocc+=focc[i];
      fder = weights[i] / (double) SumWeight * 2.0 * beta * fdexp / (1.0+fdexp) / (1.0 + fdexp);
      sumder += fder;
      }
    else
      {
      fdexp = exp (delta*beta);
      focc[i] = weights[i] / (double) SumWeight * 2.0 / (1.0+fdexp);
      sumocc += focc[i];
      fder = weights[i] / (double) SumWeight * 2.0 * beta*fdexp / (1.0+fdexp) / (1.0 + fdexp);
      sumder += fder;
      }
    } /*for i*/
  rsid = sumocc - (double) nel;
  delta = (-rsid/sumder);
  if (fabs(rsid) > FDEPSOCC)
    efermi -= rsid/sumder;
  } while (fabs(rsid) > FDEPSOCC && it <= FDMAXIT);

*ebond = 0.0;
for(i=0;i<natob;i++)
  *ebond += focc[i] * eval[i] ;
for(i=0;i<natob;i++)
	focc[i] = focc[i] * (SumWeight/ (double) weights[i]);
}

void gammaenergy_Kpt(coord,latt,nat,force,nx,ny,nz,offset_x,offset_y,offset_z)
position *coord,*latt;
int nat;
position *force;
int nx,ny,nz;
int offset_x,offset_y,offset_z;
/*Returns the band structure contribution to the energy/atom at the gamma point
  for nat atoms with orthogonal positions in cell given by coord,
  for a coordinate system with lattice vectors latt.  Nonorthogonal case is assumed.
  ADDS to the forces in array force[], which are already present from pair term.
*/
{
static double **drem,**rem,**remLatt;
static double complex **eigenv,**rho,**rho2;
static double *ggocc;

static int firstcall = 1;
int i,j,k,nmatrix;

static double complex **remC, **remLattC, **eigensave;
static double total_eigen;

static position Kpts[MAXKPT];
static int Weights[MAXKPT];
int NumKPoints;

position K;
int do_K;
int SumWeight = 1;
int ii,jj,kk;

if (firstcall == 1) {
drem = dmatrix(1,gstride*MMATOM,1,MMATOM);
rem = dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
remLatt = dmatrix(1,gstride*MMATOM,1,gstride*MMATOM*9);
remC = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
remLattC = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM*9);
eigensave = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM*KPTEIGENSAVE);
eigenv = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
rho = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
rho2 = Complex_dmatrix(1,gstride*MMATOM,1,gstride*MMATOM);
ggocc = dvector(0,gstride*MMATOM*MAXKPT);
firstcall = 0;
  }

generate_Kpts(Kpts,Weights,&NumKPoints,&SumWeight,nx,ny,nz,offset_x,offset_y,offset_z);

nmatrix = gstride * gnsilicon + (nat-gnsilicon);

for(ii=0;ii<NumKPoints;ii++) {
/*Find Hamiltonian and diagonalize it*/
K = Kpts[ii];
do_K = 1;
gammamat(coord,latt,rem ,remLatt, nat,0,remC,remLattC,K,do_K); /* flag = 0 for H matrix, 1 for S, if present*/


forcediagonalizeg_complex(remC,eigenenergy+ii*NC(nmatrix),nmatrix,eigenv);

if (ii < KPTEIGENSAVE)
  for(i=1;i<=nmatrix;i++)
	for(jj=1;jj<=nmatrix;jj++)
      eigensave[i][jj+ii*nmatrix] = eigenv[i][jj];

for(jj=0;jj<NC(nmatrix);jj++)
{
	gweights[jj+ii*NC(nmatrix)] = Weights[ii];
	gfocc[jj+ii*NC(nmatrix)] = 0.0;
	which_kpt[jj+ii*NC(nmatrix)] = ii;
	index_within_kpt[jj+ii*NC(nmatrix)] = jj;
}
{
	double eigentot;
	eigentot = 0.0;
	for(jj=0;jj<NC(nmatrix)/2;jj++)
		eigentot += 2.0*eigenenergy[jj+ii*NC(nmatrix)];
}
} //ii

sort_eigen(NumKPoints*NC(nmatrix),0,eigenenergy,gweights,which_kpt,index_within_kpt,gfocc);

/*Fermi diract distribution: returns total energy*/
yfdocclocal(&gtoten,fermitemp,eigenenergy,NumKPoints*NC(nmatrix),nmatrix,gfocc,gweights,SumWeight);

sort_eigen(NumKPoints*NC(nmatrix),1,eigenenergy,gweights,which_kpt,index_within_kpt,gfocc);

/*Create density matrix rho*/
/*Need to include extra terms.  We have gfocc[i-1]...gfocc[NC(gstride*nat)-1]
also fermitemp = temp in Kelvin so "beta" = 1 / kT = 11604 / fermitemp
ggocc maximum is 2.0 */

if (ZEROTFORCE==FALSE) {
double term1,  term2,beta;
double totaln;
term1 = 0;
term2 = 0;
beta = 11604.0 / fermitemp;
totaln =0;
for(i=1;i<=NumKPoints*NC(nmatrix);i++)
  term1 += gfocc[i-1] * -1.0*(1.0 - gfocc[i-1] / 2.0);
for(i=1;i<=NumKPoints*NC(nmatrix);i++)
  term2 += eigenenergy[i-1] * gfocc[i-1] * -1.0*(1.0 - gfocc[i-1] / 2.0);
for(i=1;i<=NumKPoints*NC(nmatrix);i++)
  ggocc[i-1] = gfocc[i-1] + beta * gfocc[i-1] * -1.0*(1.0 - gfocc[i-1] / 2.0) * (eigenenergy[i-1] - 1.0*term2 / term1);
for(i=1;i<=NumKPoints*NC(nmatrix);i++)
  totaln += gfocc[i-1];
//printf("Total number of electrons= %lf\n",totaln);
//for(i=0;i<NumKPoints*NC(nmatrix);i++)
//{
//	printf("%d, %3.18lf %lf %lf\n",i,eigenenergy[i],ggocc[i],gfocc[i]);
//}
}

total_eigen = 0.0;
for(kk=0;kk<NumKPoints;kk++)
{

	if (kk < KPTEIGENSAVE)
	  for(i=1;i<=nmatrix;i++)
		for(jj=1;jj<=nmatrix;jj++)
	      eigenv[i][jj] = eigensave[i][jj+kk*nmatrix];
	else {
	//get eigenV
	K = Kpts[kk];
	do_K = 1;
	gammamat(coord,latt,rem ,remLatt, nat,0,remC,remLattC,K,do_K); /* flag = 0 for H matrix, 1 for S, if present*/
	forcediagonalizeg_complex(remC,eigenenergy+kk*NC(nmatrix),nmatrix,eigenv);
	};

for(i=1;i<=nmatrix;i++)
  for(j=1;j<=nmatrix;j++)
    rho[i][j] = 0.0;

for(i=1;i<=NC(nmatrix);i++) /*loop was over occupied only*/
  {
  /*For ith eigenvector eigenv[i][1..gstride*nat], compute rho and rho2 */
  /*Use factors (gfocc[i-1] / 2.0) for eigenv[i][]  */
  /*Rho is quadratic in eigenfunctions, so contains this factor (probability) once*/
  double sterm;
  sterm = 0.0;

sterm = 2.0; /*Well behaved eigensolver*/

  for(j=1;j<=nmatrix;j++)
    {
    double complex reg1, reg2;
#if (ZEROTFORCE)
    reg1 = eigenv[i][j]*sterm*(gfocc[kk*NC(nmatrix)+i-1] / 2.0) ;
#else
    reg1 = eigenv[i][j]*sterm*(ggocc[kk*NC(nmatrix)+i-1] / 2.0) ;
#endif
#define SIGN1 (-1)
#define SIGN2 (1)
#define SIGN3 (-1)
    for(k=j;k<=nmatrix;k++) {
    	double complex temp, temp2;
    temp = creal(eigenv[i][k]) + SIGN1* I*cimag(eigenv[i][k]); //conjugate
    temp2 = creal(reg1) + SIGN2* I*cimag(reg1);
      rho[j][k] += temp2 * temp;
    } //k
    }//j
  }//i

/*Take trace rho * H for Hellman Feynman theorem*/
/*Note that other triangle of rho, and rho2 is missing and we do sum with
  factors of two to compensate.*/
for(i=1;i<=3;i++) {
	K = Kpts[kk];
	do_K = 1;
gammamat(coord,latt,rem,remLatt,nat,i,remC,remLattC,K,do_K);

for(j=1;j<=nmatrix;j++)
  for(k=j;k<=nmatrix;k++) {
int a1,a2;
double temp;
/*MARK*/
if ((k-1)/4>=gnsilicon) a2 = gnsilicon + (k-1) - gnsilicon*4; else a2 = (k-1)/4;
if ((j-1)/4>=gnsilicon) a1 = gnsilicon + (j-1) - gnsilicon*4; else a1 = (j-1)/4;

{
	int ii,jj;

if (j!=k || 1) {
	       int prefactor;
	       if (j!=k)
		 prefactor = 2;
	       else
 		 prefactor = 1;
	       {
	    	   double complex cc;
	    	   cc = creal(remC[j][k]) + SIGN3*I*cimag(remC[j][k]);
	           temp = creal(rho[j][k]*cc*2.0)*Weights[kk]/(double)SumWeight ;
	       }
	       for (ii=0;ii<3;ii++)
	    	   for (jj=0;jj<3;jj++)
	    	   {
	    		   double complex cc;
	    	       cc = creal(rho[j][k]) - I*cimag(rho[j][k]);
	    		   glattforce[ii][jj] += Weights[kk]/(double)SumWeight*creal(cc*remLattC[j][(k-1)*9+1+3*ii+jj] * prefactor);
               }
} //if
//if (j==k) temp = 0.0;
if (i==1) {force[a1].x -= temp; force[a2].x += temp;};
if (i==2) {force[a1].y -= temp; force[a2].y += temp;};
if (i==3) {force[a1].z -= temp; force[a2].z += temp;};
} //

} /*j,k loop*/
} /*i loop*/

} // end of kk loop over Kpoints

}


