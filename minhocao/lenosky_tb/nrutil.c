#include <stdlib.h>
#include <stdio.h>
#include <complex.h>

void nrerror(error_text)
char error_text[];
/* Numerical Recipes standard error handler */
{

	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
}

float *vector(nl,nh)
int nl,nh;
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((unsigned) (nh-nl+1)*sizeof(float));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl;
}

int *ivector(nl,nh)
int nl,nh;
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((unsigned) (nh-nl+1)*sizeof(int));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl;
}

double *dvector(nl,nh)
int nl,nh;
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl;
}

float **matrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((unsigned) (nrh-nrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure 1 in matrix()");
      m -= nrl;

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(float *) malloc((unsigned) (nch-ncl+1)*sizeof(float));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
            m[i] -= ncl;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
      m -= nrl;

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
            m[i] -= ncl;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

double complex **Complex_dmatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i;
	double complex **m;

	/* allocate pointers to rows */
	m=(double complex **) malloc((unsigned) (nrh-nrl+1)*sizeof(double complex*));
	if (!m) nrerror("allocation failure 1 in dmatrix()");
      m -= nrl;

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(double complex *) malloc((unsigned) (nch-ncl+1)*sizeof(double complex));
		if (!m[i]) nrerror("allocation failure 2 in dmatrix()");
            m[i] -= ncl;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(nrl,nrh,ncl,nch)
int nrl,nrh,ncl,nch;
/* allocate an int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	int i,**m;

	/* allocate pointers to rows */
	m=(int **)malloc((unsigned) (nrh-nrl+1)*sizeof(int*));
	if (!m) nrerror("allocation failure 1 in imatrix()");
      m -= nrl;

	/* allocate rows and set pointers to them */
	for(i=nrl;i<=nrh;i++) {
		m[i]=(int *)malloc((unsigned) (nch-ncl+1)*sizeof(int));
		if (!m[i]) nrerror("allocation failure 2 in imatrix()");
            m[i] -= ncl;
	}
	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(a,oldrl,oldrh,oldcl,oldch,newrl,newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
float **a;
int oldrl,oldrh,oldcl,oldch,newrl,newcl;
{
	int i,j;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((unsigned) (oldrh-oldrl+1)*sizeof(float*));
	if (!m) nrerror("allocation failure in submatrix()");
      m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+oldcl-newcl;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(a,nrl,nrh,ncl,nch)
float *a;
int nrl,nrh,ncl,nch;
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	int i,j,nrow,ncol;
	float **m;

	nrow=nrh-nrl+1;
	ncol=nch-ncl+1;
	/* allocate pointers to rows */
	if ((m=(float **) malloc((unsigned) (nrow)*sizeof(float*))) == NULL)
		nrerror("allocation failure in convert_matrix()");
	m -= nrl;

	/* set pointers to rows */
	for(i=0,j=nrl;i<=nrow-1;i++,j++) m[j]=a+ncol*i-ncl;
	/* return pointer to array of pointers to rows */
	return m;
}

void free_vector(v,nl,nh)
float *v;
int nl,nh;
/* free a float vector allocated with vector() */
{
	free((char*) (v+nl));
}

void free_ivector(v,nl,nh)
int *v,nl,nh;
/* free an int vector allocated with ivector() */
{
	free((char*) (v+nl));
}

void free_dvector(v,nl,nh)
double *v;
int nl,nh;
/* free a double vector allocated with dvector() */
{
	free((char*) (v+nl));
}

void free_matrix(m,nrl,nrh,ncl,nch)
float **m;
int nrl,nrh,ncl,nch;
/* free a float matrix allocated by matrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_dmatrix(m,nrl,nrh,ncl,nch)
double **m;
int nrl,nrh,ncl,nch;
/* free a double matrix allocated by dmatrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_imatrix(m,nrl,nrh,ncl,nch)
int **m;
int nrl,nrh,ncl,nch;
/* free an int matrix allocated by imatrix() */
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}

void free_submatrix(b,nrl,nrh,ncl,nch)
/* free a submatrix allocated by submatrix() */
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}

void free_convert_matrix(b,nrl,nrh,ncl,nch)
/* free a matrix allocated by convert_matrix() */
float **b;
int nrl,nrh,ncl,nch;
{
	free((char*) (b+nrl));
}
