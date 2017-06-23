	#include <stdlib.h>
#include <stdio.h>
/*Note this package contains modified versions of numerical recipes 
  minimization routines (Numerical Recipes in C):

  see routines "brent", "linmin", "mnbrak", "f1dim", "frprmn" from which 
  these routines are adapted.

*/

#include <math.h>

extern double gmaxforce;

#define AHF 2.0 /* Ad hoc factor for ylinmin call to mnbrak */
#define USENL 0 /*Flag to use adaptive stepsize in ylinmin*/

#define USEFN 0
#define FNORM 1.0 /*Convergence parameter for CG routine */ 

#define USEFEXP 0
#define FDECAY 0.2 

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

	double brent(ax,bx,cx,f,tol,xmin,fxinit)
	double ax,bx,cx,tol,*xmin,fxinit;
	double (*f)();	/* ANSI: double (*f)(double); */
	{
		int iter;
		double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
		double e=0.0;
		void nrerror();

		a=((ax < cx) ? ax : cx);
		b=((ax > cx) ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=fxinit;
		for (iter=1;iter<=ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
			if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
				*xmin=x;
				return fx;
			}
			if (fabs(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=fabs(q);
				etemp=e;
				e=d;
				if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
			} else {
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fu=(*f)(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				SHFT(v,w,x,u)
				SHFT(fv,fw,fx,fu)
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		nrerror("Too many iterations in BRENT");
		*xmin=x;
		return fx;
	}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN

	/*VARIABLES FOR LINMIN:*/
	int ncom=0;	/* defining declarations */
	double *pcom=0,*xicom=0,(*nrfunc)();

	void ylinmin(p,xi,n,fret,func,lintol,oldstep)
	double p[],xi[],*fret,(*func)(), lintol,*oldstep;
	int n;
	/* Contains code for setting scale of initial guess using previous motion as a guide */
	{
		static double xdist = 1.0;
		int j;
		double xx,xmin,fx,fb,fa,bx,ax;
		double brent(),f1dim(),*dvector();
		void mnbrak(),free_dvector();
		ncom=n;
		pcom=dvector(1,n);
		xicom=dvector(1,n);
		nrfunc=func;
		for (j=1;j<=n;j++) {
			pcom[j]=p[j];
			xicom[j]=xi[j];
		}
#if (USENL)
		ax=0.0;
		xx=1.0*xdist/AHF;
		bx=2.0*xdist/AHF;
#endif
#if (!USENL)
		ax = 0.0;
		xx = 1.0;
		bx = 2.0;
#endif
		fa = *fret; /*setup for mnbrak*/ 
		mnbrak(&ax,&xx,&bx,&fa,&fx,&fb,f1dim);
		*fret=brent(ax,xx,bx,f1dim,lintol,&xmin,fx);
		xdist = sqrt(fabs(xmin*xdist)); /* New to prevent blowup while getting it ok*/
		for (j=1;j<=n;j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		free_dvector(xicom,1,n);
		free_dvector(pcom,1,n);
	}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak(ax,bx,cx,fa,fb,fc,func)
double *ax,*bx,*cx,*fa,*fb,*fc;
double (*func)();	/* ANSI: double (*func)(double); */
{
	double ulim,u,r,q,fu,dum;

	/* *fa=(*func)(*ax); REMOVED TO SPEED CODE*/
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT

extern int ncom;	/* defined in LINMIN */
extern double *pcom,*xicom,(*nrfunc)();

double f1dim(x)
double x;
{
	int j;
	double f,*xt,*dvector();
	void free_dvector();

	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_dvector(xt,1,ncom);
	return f;
}

#define EPS 1.0e-10 
#define FREEALL free_dvector(xi,1,n);free_dvector(h,1,n);free_dvector(g,1,n);

	void frprmn(p,n,ftol,iter,fret,itmax,lintol)
	double p[],ftol,*fret,lintol;

	int n,*iter;
	int itmax;
	{
		int j,its;
		double gg,gam,fp,dgg;
		double *g,*h,*xi,*dvector();
		void ylinminspecial(),nrerror(),free_dvector();
		void dylinmin();
                extern double dfunc();
                extern double func();

		g=dvector(1,n);
		h=dvector(1,n);
		xi=dvector(1,n);

		fp=dfunc(p,xi);
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j];
		}
		for (its=1;its<=itmax;its++) {

		*iter=its;
	        *fret = fp; /*Set up for new version of ylinmin*/	
// 		ylinminspecial(p,xi,n,fret,lintol);

               dylinmin(p,xi,n,fret,lintol);

		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
  
		/*fp=(*func)(p); Note this change speeds program*/
                fp = *fret;

		/*dfunc(p,xi);
*/
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
/*		  dgg += xi[j]*xi[j];	*/
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
}
#undef FREEALL
#undef EPS


#define DDEPS 0.00100
	void dylinmin(p,xi,n,fret,lintol)
	double p[],xi[],*fret,lintol;
	int n;
	/* Contains code for setting scale of initial guess using previous motion as a guide */
	{
		static double xdist = 0.01;
		int i,j;
		double xx,xmin,fx,fb,fa,bx,ax;
		double brentspecial(),f1dimspecial(),*dvector();
		void mnbrakspecial(),free_dvector();
 		void ylinminspecial();
                extern double func();
                extern double dfunc();
		static double force[300000],xiunit[300000];
                static double dtot = 0.0;
                static double norm = 0.0;
		static double aa,bb[2],cc[2],newval,v0,v1;
                static int qfail = 0;

/* *fret is value of starting point */
//if (qfail == 0 && 0==1) {
if (1==1) {
               aa = *fret;
               newval =0;
               norm = 0.0;
	       for(j=1;j<=n;j++)
                 norm+= xi[j]*xi[j];
               for(j=1;j<=n;j++)
                  xiunit[j] = xi[j] / sqrt(norm);

for(i=0;i<2;i++) {
               for(j=1;j<=n;j++)
		  p[j] += DDEPS*(1-2*i)*xiunit[j];
               v1 = dfunc(p,force);
               for(j=1;j<=n;j++)
                  p[j] -= DDEPS *(1-2*i)* xiunit[j];
               dtot = 0.0;
               for(j=1;j<=n;j++)
                  {
                  dtot -= force[j] * xiunit[j]; 
                  }

                cc[i] = (dtot * DDEPS*(1-2*i) + (aa - v1)) / (DDEPS * DDEPS);              
                bb[i] = -dtot + 2 * (v1 -aa) /(DDEPS*(1-2*i));
                    }
                cc[0] += cc[1];  cc[0] /= 2.0;
                bb[0] += bb[1];  bb[0] /= 2.0;
                for(j=1;j<=n;j++)
                  p[j] += -bb[0]/cc[0]/2.0 * xiunit[j];

                newval = dfunc(p,force);

                printf("STEPSIZE = %lf\n",-bb[0]/cc[0]/2.0);
                printf("aa, newval = %3.18le %3.18le\n",aa,newval);
		printf("v1 = %3.18le\n",v1);
                if (newval > aa )
                   {
                   printf("FAILURE OF QUAD STEP\n");
                   qfail = 1;
                   for(j=1;j<=n;j++)
                     p[j] -= -bb[0]/cc[0]/2.0 * xiunit[j];
                   ylinminspecial(p,xi,n,fret,lintol);
                   dfunc(p,xi);
                   }
                else
                  {
                  *fret = newval;
                  for(j=1;j<=n;j++)
                    xi[j] = force[j];
                  };
 }
else
 {ylinminspecial(p,xi,n,fret,lintol);
  dfunc(p,xi);
 }
	}

	void ylinminspecial(p,xi,n,fret,lintol)
	double p[],xi[],*fret,lintol;
	int n;
	/* Contains code for setting scale of initial guess using previous motion as a guide */
	{
		static double xdist = 0.1;;
		int j;
		double xx,xmin,fx,fb,fa,bx,ax;
		double brentspecial(),f1dimspecial(),*dvector();
		void mnbrakspecial(),free_dvector();
                extern double func();

xdist = 0.1;

		ncom=n;
		pcom=dvector(1,n);
		xicom=dvector(1,n);
		for (j=1;j<=n;j++) {
			pcom[j]=p[j];
			xicom[j]=xi[j];
		}
		ax = 0.0;
		xx = xdist;
		bx = 2*xdist;
		fa = *fret; /*setup for mnbrak*/ 
		mnbrakspecial(&ax,&xx,&bx,&fa,&fx,&fb,&f1dimspecial);
		*fret=brentspecial(ax,xx,bx,f1dimspecial,lintol,&xmin,fx);
		xdist = sqrt(fabs(xmin*xdist)); /* New to prevent blowup while getting it ok*/
		for (j=1;j<=n;j++) {
			xi[j] *= xmin;
			p[j] += xi[j];
		}
		free_dvector(xicom,1,n);
		free_dvector(pcom,1,n);
	}




/*  ****** */
  
#include <math.h>
#include <stdio.h>

#define ITMAX 100
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

	double brentspecial(ax,bx,cx,f,tol,xmin,fxinit)
	double ax,bx,cx,tol,*xmin,fxinit;
	double (*f)();	/* ANSI: double (*f)(double); */
	{
		int iter;
		double a,b,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
		double e=0.0;
                double d=0.0;
		void nrerror();

		a=((ax < cx) ? ax : cx);
		b=((ax > cx) ? ax : cx);
		x=w=v=bx;
		fw=fv=fx=fxinit;
		for (iter=1;iter<=ITMAX;iter++) {
			xm=0.5*(a+b);
			tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
			if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
				*xmin=x;
				return fx;
			}
			if (fabs(e) > tol1) {
				r=(x-w)*(fx-fv);
				q=(x-v)*(fx-fw);
				p=(x-v)*q-(x-w)*r;
				q=2.0*(q-r);
				if (q > 0.0) p = -p;
				q=fabs(q);
				etemp=e;
				e=d;
				if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
					d=CGOLD*(e=(x >= xm ? a-x : b-x));
				else {
					d=p/q;
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				}
			} else {
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			}
			u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
			fu=(*f)(u);
			if (fu <= fx) {
				if (u >= x) a=x; else b=x;
				SHFT(v,w,x,u)
				SHFT(fv,fw,fx,fu)
			} else {
				if (u < x) a=u; else b=u;
				if (fu <= fw || w == x) {
					v=w;
					w=u;
					fv=fw;
					fw=fu;
				} else if (fu <= fv || v == x || v == w) {
					v=u;
					fv=fu;
				}
			}
		}
		nrerror("Too many iterations in BRENT");
		*xmin=x;
		return fx;
	}

#undef ITMAX
#undef CGOLD
#undef ZEPS
#undef SIGN

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SIGN(a,b) ((b) > 0.0 ? fabs(a) : -fabs(a))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrakspecial(ax,bx,cx,fa,fb,fc,func)
double *ax,*bx,*cx,*fa,*fb,*fc;
double (*func)();	/* ANSI: double (*func)(double); */
{
	double ulim,u,r,q,fu,dum;

	/* *fa=(*func)(*ax); REMOVED TO SPEED CODE*/
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT

extern int ncom;	/* defined in LINMIN */
extern double *pcom,*xicom;

double f1dimspecial(x)
double x;
{
	int j;
	double f,*xt,*dvector();
	void free_dvector();
        extern double func();

	xt=dvector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=func(xt);
	free_dvector(xt,1,ncom);
	return f;
}

