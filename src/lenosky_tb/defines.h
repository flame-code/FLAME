/*----------------------------------------------------------------------------
defines.h: file for vector type defintions,etc. 
*/

#define RMIN 1.75059816 

#define TRUE 1
#define FALSE 0

#define MAXCELLS 50

#define GAUSSIANG FALSE /*USE Gaussian g(r) for MEAM fitting*/
#define MEAMFSQ FALSE /*USE Fsq in meam if MEAMRHO is false */
#define MEAMRHO TRUE /*HAVE a rho(r) in MEAM formalism */
#define OLDPCCOS TRUE 
#define FIXTB FALSE
#define MHMAX 11 /*Number of H splines, maximum*/
#define MEPSMAX 10 /*Number of epsilon parameters, maximum*/
#define MRHOMAX 1 
#define MRHOMAX2 1
#define CLSMRHOMAX 1
#define CLSMRHOMAX2 1
#define CLSMAXN2 100
#define MAXTRI  10 /*WAS 48000*/
#define NSNAPSHPP 188 /*Preprocess all the frames now was NSNAPSH*/
#define CMAX 100 /*Maximum number of adjustable coefficients, once was 450 for meam*/
#define MAXATOM  300 /*Maximum number of atoms in set frames */
#define NSNAPSH 188 /*raise for Harris*/
#define CLSDATA FALSE /*Use special data*/
#define MAXRESID 16000 /*Needs to be increased with size of data set*/ 
#define MAXNEIGHB 40 

#define GAMMATYPE float /*Either float or double, data type for gammas in new powell routine*/
#define NSPMAX 20 /*Maximum number of spline points*/

#define CLSMAXATOM  100
#define CLSMAXNEIGHB  100

typedef struct {
  double x,y,z;
} position;

typedef struct {
  float a,b,c,d,e,f,g,h,m,n,nn;
  int i1,i2;
} ps_type;

/* OLD
typedef struct {
  int nbind[MAXATOM][MAXNEIGHB];
  int nb[MAXATOM];
  double r[MAXATOM][MAXNEIGHB];
  position grad[MAXATOM][MAXNEIGHB];
  ps_type phips[MAXATOM][MAXNEIGHB];
} framepp_type;
*/

typedef struct {
  /*For neighbor lists with reversed ordering: */
  int nbtrixx[MAXATOM];
  int nbindxx[MAXATOM][MAXNEIGHB];

  int nbind[MAXATOM][MAXNEIGHB];
  int nb[MAXATOM];
  int nbrho[MAXATOM];
  int nbtri[MAXATOM];
  double r[MAXATOM][MAXNEIGHB];
  position grad[MAXATOM][MAXNEIGHB];
  ps_type rhops[MAXATOM][MAXNEIGHB];
  ps_type phips[MAXATOM][MAXNEIGHB];
  ps_type fps[MAXATOM][MAXNEIGHB];
  int ntri; /*Number of triangles in frame*/
  double rtri[MAXTRI][3];
  position gradtri[MAXTRI][3];
  int indtri[MAXTRI][3];

/*Note as spline intervals are the same for every rho and f, don't
  need to store seperate ps records (just as for rhops): Also list
  of triangles is the same for every emb. function:  Thus rhops,
  fps, ntri, rtri, gradtri, indtri are all re-used*/
  double xbeta[MRHOMAX2][MAXATOM];
  double xfgrhomax[MRHOMAX2], xfgrhomin[MRHOMAX2];

/*Note also that fvalptr, cosval, and also gps and tusage can all be re-used*/

  /*These fields are not set by preprocess but rather by rhomaxmin:*/
  double beta[MRHOMAX][MAXATOM];
  double fgrhomax[MRHOMAX], fgrhomin[MRHOMAX];

  /*Values of f and its derivative precomputed during energy routine:*/
  double fval[MAXATOM][MAXNEIGHB][2+2*MRHOMAX2]; /*Fields store f and df/dr for sw & meam terms*/
  double *fvalptr[MAXTRI][3]; /*Pointer to fval[..][..][2] array for triangle*/
  double cosval[MAXTRI][3][4]; /*Precomputed cosines and derivatives*/
  ps_type gps[MAXTRI][3]; /*Precomputed splines for g*/
  int tusage[MAXTRI][3];
} framepp_type;


typedef struct {
  /*Note that total storage in bytes is 56AN + 16A+ 8AR
    where A= Number of atoms
          N = Number of neighbors max.
          R = Number of emb. potls max.

    Normally only the first term is important--the other
    terms being typically less than 100 kilobytes*/
  /*Added memory for MEAM terms:
      If Q = Number of emb. potls max.

      Then 8QA + 16ANQ = Added memory in bytes.
  */

  /*xx records are for neighbor lists with reversed ordering: */
  int nbtrixx[CLSMAXATOM];
  int nbindxx[CLSMAXATOM][CLSMAXNEIGHB];

  int nbind[CLSMAXATOM][CLSMAXNEIGHB];
  int nb[CLSMAXATOM];
  int nbrho[CLSMAXATOM];
  int nbtri[CLSMAXATOM];
  int iset[CLSMAXATOM];
  double r[CLSMAXATOM][CLSMAXNEIGHB];
  position grad[CLSMAXATOM][CLSMAXNEIGHB];

  /*These fields are not set by preprocess but rather by rhomaxmin:*/
  float beta[CLSMRHOMAX][CLSMAXATOM];

  /*Values of f and its derivative precomputed during energy routine:*/
  double fval[CLSMAXATOM][CLSMAXN2][2+2*CLSMRHOMAX2]; /*Fields store f and df/dr */
  double zval[CLSMRHOMAX2][CLSMAXATOM];
  double de_dz[CLSMRHOMAX2][CLSMAXATOM]; /* dxbeta/dZ */
  double xbeta[CLSMRHOMAX2][CLSMAXATOM];

  int cellptr[MAXCELLS][MAXCELLS][MAXCELLS];
  int linkptr[CLSMAXATOM];

} clsframepp_type;
typedef struct {
  position *p, *f, a;
  double lattforce[3][3];
  double en;
  int n,nsym,nsi;
  int id; /*identification number*/
  int *itype;
  int ntbcall;
  position *ldaforce;
  double ldaen;
  int stride;
  } frame_type;

typedef struct {
  double x[NSPMAX+1];
  double y[NSPMAX+1];
  double yp1;
  double ypn;
  double y2[NSPMAX+1];
  int npt;
  int iset1,isetn,ider1,idern;
} spline_type;

typedef struct {
  spline_type phi[3],f,g,u[MRHOMAX],rho[MRHOMAX],h[MHMAX];
  double eps[MEPSMAX];
  spline_type xf[MRHOMAX2], xg[MRHOMAX2],xu[MRHOMAX2],xrho[MRHOMAX2];
} potl_type;

typedef struct {
framepp_type *framepp;
frame_type *frame;
double fcomp, totenu, totenu2,totentri, totenphi,totentb,toten;
position resf[MAXATOM];
double rese;
} framerec_type;

typedef struct {
framerec_type data;
double eweight, fweight;
struct framereclist_type *next;
} framereclist_type;

#define NULLTEST(px) { if (px == NULL) printf("WARNING: Malloc failure...\n");  }
 
#define MIN(a,b) (a<b ? a : b)
#define MAX(a,b) (a>b ? a : b)

#define NINT(x) floor((x)+0.5)

#define VECT_SUBTRACT(a,b,c) \
  (c).x = (a).x - (b).x; \
  (c).y = (a).y - (b).y; \
  (c).z = (a).z - (b).z

#define VECT_ADD(a,b,c) \
  (c).x = (a).x + (b).x; \
  (c).y = (a).y + (b).y; \
  (c).z = (a).z + (b).z

#define CROSS(a,b,c) \
  (c).x = (a).y * (b).z - (a).z * (b).y; \
  (c).y = (a).z * (b).x - (a).x * (b).z; \
  (c).z = (a).x * (b).y - (a).y * (b).x

#define APPLY_PBC(a) \
  (a).x -= NINT((a).x); \
  (a).y -= NINT((a).y); \
  (a).z -= NINT((a).z)

#define CONVERT(s,r) \
  (r).x = gg[0][0]*(s).x + gg[0][1]*(s).y + gg[0][2]*(s).z; \
  (r).y = gg[1][0]*(s).x + gg[1][1]*(s).y + gg[1][2]*(s).z; \
  (r).z = gg[2][0]*(s).x + gg[2][1]*(s).y + gg[2][2]*(s).z

#define REVERT(r,s) \
  (s).x = gg_inv[0][0]*(r).x + gg_inv[0][1]*(r).y + gg_inv[0][2]*(r).z; \
  (s).y = gg_inv[1][0]*(r).x + gg_inv[1][1]*(r).y + gg_inv[1][2]*(r).z; \
  (s).z = gg_inv[2][0]*(r).x + gg_inv[2][1]*(r).y + gg_inv[2][2]*(r).z

#define CONVERT_TO_RECIPROCAL_LATTICE(r,s) \
		(s).x = gg_inv[0][0]*(r).x + gg_inv[1][0]*(r).y + gg_inv[2][0]*(r).z; \
		(s).y = gg_inv[0][1]*(r).x + gg_inv[1][1]*(r).y + gg_inv[2][1]*(r).z; \
		(s).z = gg_inv[0][2]*(r).x + gg_inv[1][2]*(r).y + gg_inv[2][2]*(r).z

#define LENGTH(a) sqrt(DOT(a,a))

#define DOT(a,b) ((a).x*(b).x + (a).y*(b).y + (a).z*(b).z)
#define EQUAL_ZERO(a,b) (fabs(a)<b)

#define MALLFR(fr,nat) {int kx; fr.ntbcall = 0; fr.p = (position *) malloc (nat * sizeof(position)); NULLTEST(fr.p); fr.f = (position *) malloc(nat*sizeof(position)); NULLTEST(fr.f); fr.itype = (int *) malloc (nat*sizeof(int)); NULLTEST(fr.itype); for(kx=0;kx<(nat);kx++) (fr).itype[kx] = 1; fr.nsym = nat;} 

#define SVAL(v,u) ((v).a * (u).y[(v).i1] + (v).b * (u).y[(v).i2] + (v).c * (u).y2[(v).i1] +  (v).d * (u).y2[(v).i2] + (v).m *((v).n*(u).yp1 + (v).nn*(u).ypn) )

#define SDER(a,b) ((a).e * (b).y[(a).i1] + (a).f * (b).y[(a).i2] + (a).g * (b).y2[(a).i1] + (a).h * (b).y2[(a).i2] + (a).n*(b).yp1 + (a).nn*(b).ypn)

#define SSETX(sp,low,high) {int j; for(j=1; j<=(sp).npt;j++) (sp).x[j] = low + ((high-low)*(j-1))/((double)(((sp).npt)-1)); }

#define SSETY(sp,px) {int j; double *q; q=px; for(j=2-(sp).iset1;j<=(sp).npt -1 + (sp).isetn;j++) {(sp).y[j] = q[0]; q++;}; if((sp).ider1==1) {(sp).yp1 = q[0]; q++; }; if((sp).idern==1) {(sp).ypn = q[0]; q++;}; spline(&(sp)); }

#define SENDP(sp,va,vb,vc,vd) {(sp).iset1 = va; (sp).isetn = vb; (sp).ider1=vc; (sp).idern=vd; }

/*HAS NO itype*/
#define COPYFRAME(xx,yy) {int i; yy=(frame_type *) malloc(sizeof(frame_type)); NULLTEST(yy); (*yy) = (*xx); (*yy).stride = (*xx).stride; (*yy).p =(position *) malloc(sizeof(position)*(*xx).n); NULLTEST((*yy).p); (*yy).f=(position *) malloc(sizeof(position)*(*xx).n); NULLTEST((*yy).f); (*yy).itype = (int *) malloc(sizeof(int)*(*xx).n); NULLTEST((*yy).itype); for(i=0;i<(*xx).n;i++) {(*yy).p[i]=(*xx).p[i]; (*yy).f[i] = (*xx).f[i]; (*yy).itype[i] = (*xx).itype[i];} (*yy).nsym = (*xx).nsym; }

#define CELLDIST(p1,p2,a,d) {position diff; VECT_SUBTRACT(p1,p2,diff); diff.x /= a.x; diff.y /=a.y; diff.z /= a.z; APPLY_PBC(diff); diff.x *= a.x; diff.y *= a.y; diff.z *= a.z; d=LENGTH(diff); }

#define CELLGRAD(p1,p2,a,g) {position diff; double r; VECT_SUBTRACT(p1,p2,diff); diff.x /= a.x; diff.y /= a.y; diff.z /= a.z; APPLY_PBC(diff); diff.x *= a.x; diff.y *=a.y; diff.z *=a.z; r=LENGTH(diff); g.x = diff.x / r; g.y = diff.y /r; g.z = diff.z/r; }

#define CELLDISTNO(p1,p2,a,d) {position diff,t; VECT_SUBTRACT(p1,p2,diff); diff.x /= a.x; diff.y /=a.y; diff.z /= a.z; APPLY_PBC(diff); diff.x *= a.x; diff.y *= a.y; diff.z *= a.z; CONVERT(diff,t); d=LENGTH(t); }

#define CELLGRADNO(p1,p2,a,g) {position diff,t; double r; VECT_SUBTRACT(p1,p2,diff); diff.x /= a.x; diff.y /= a.y; diff.z /= a.z; APPLY_PBC(diff); diff.x *= a.x; diff.y *=a.y; diff.z *=a.z; CONVERT(diff,t); r=LENGTH(t); g.x = t.x / r; g.y = t.y /r; g.z = t.z/r; }


#define ACPOTL(xv,d1,d2,d3,p1,p2,pg,psg) {double t0,t1,t2,t3,f1,f2,v1,v2,f3,v3; t0 = pg[0]; t1 = pg[1]; t2=pg[2]; t3=pg[3]; f1 = p1[0]; v1= p1[1];  f2 = p2[0]; v2=p2[1]; f3 = SVAL((psg),pp.g); v3 = SDER((psg),pp.g); d1 += v1*f2*f3 +f1*f2*v3*t1; d2 += f1*v2*f3+f1*f2*v3*t2; d3+=f1*f2*v3*t3; xv+= f1*f2*f3; }

#define GACPOTL(gxx,xv,d1,d2,d3,p1,p2,pg,psg) {double t0,t1,t2,t3,f1,f2,v1,v2,f3,v3; t0 = pg[0]; t1 = pg[1]; t2=pg[2]; t3=pg[3]; f1 = p1[0]; v1= p1[1];  f2 = p2[0]; v2=p2[1]; f3 = SVAL((psg),gxx); v3 = SDER((psg),gxx); d1 += v1*f2*f3 +f1*f2*v3*t1; d2 += f1*v2*f3+f1*f2*v3*t2; d3+=f1*f2*v3*t3; xv+= f1*f2*f3; }

#define NEWGACPOTL(gxx,xv,d1,d2,d3,p1,p2,pg,psg,rr) {double t0,t1,t2,t3,f1,f2,v1,v2,f3,v3,temp; t0 = pg[0]; t1 = pg[1]; t2=pg[2]; t3=pg[3]; f1 = p1[0]; v1= p1[1];  f2 = p2[0]; v2=p2[1]; temp = gxx.yp1; f3=exp(-temp*(rr+1)*(rr+1)); v3 = -2*temp*(rr+1)*f3; d1 += v1*f2*f3 +f1*f2*v3*t1; d2 += f1*v2*f3+f1*f2*v3*t2; d3+=f1*f2*v3*t3; xv+= f1*f2*f3; }

#define DGACPOTL(gxx,xv,p1,p2,pg,psg) {double f1,f2,f3; f1=p1[0]; f2 = p2[0]; f3 = SVAL((psg),gxx); xv+=f1*f2*f3; }

#define NEWDGACPOTL(gxx,xv,p1,p2,pg,psg,rr) {double f1,f2,f3,temp; f1=p1[0]; f2 = p2[0]; temp = gxx.yp1; f3 = exp(-temp*(rr+1)*(rr+1)); xv+=f1*f2*f3; }

#define CLSACPOTL(xv,d1,d2,d3,p1,p2,pg,f3,v3) {double t0,t1,t2,t3,f1,f2,v1,v2; t0=pg[0]; t1=pg[1]; t2=pg[2]; t3=pg[3]; f1=p1[0]; v1=p1[1]; f2= p2[0]; v2=p2[1]; d1 += v1*f2*f3 + f1*f2*v3*t1; d2 += f1*v2*f3+f1*f2*v3*t2; d3+= f1*f2*v3*t3; xv+= f1*f2*f3; }

#define CALCGRADR(r1,g1,r2,g2,r3,g3) {position sum; sum.x = r1*g1.x + r2*g2.x; sum.y =  r1*g1.y+r2*g2.y; sum.z = r1*g1.z + r2*g2.z; r3 = LENGTH(sum); g3.x = sum.x / r3; g3.y = sum.y / r3; g3.z =  sum.z/r3; }

