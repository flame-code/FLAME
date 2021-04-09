//! @file
//! Various functions used during experiences (dgemmsy)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 
//! Eric Bainville, Mar 2010


#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <malloc.h>
#include <math.h>
#include <unistd.h>
#include <xmmintrin.h>

// Set to 1 to use BLAS axpy call to accumulate partial results.
#ifndef USE_AXPY
#define USE_AXPY 1
/*#define ACML_BLAS 1*/
#endif

#if USE_AXPY
#if defined(MKL_BLAS) || defined(WIN32)
#include <mkl_cblas.h>
#else
#if defined(ACML_BLAS)
#include <acml.h>
#define cblas_daxpy daxpy
#else
#include <cblas.h>
#endif
#endif
#endif

#include "utils.h"

// System-specific functions

#include <unistd.h>
void * allocPages(size_t sz)
{
  size_t page_size = sysconf(_SC_PAGESIZE);
  size_t extra = sz % page_size;
  if (extra != 0) sz = (sz - extra) + page_size; // round to next page
  void * address = 0;
  int status = posix_memalign(&address,page_size,sz);
  if (status != 0) return 0; // Failed
  return address;
}
void freePages(void * address)
{
  free(address);
}
double * allocDoubles(int n)
{
  void * address = 0;
  int status = posix_memalign(&address,16,n*sizeof(double));
  if (status != 0) return 0; // Failed
  return address;
}
float * allocFloats(int n)
{
  void * address = 0;
  int status = posix_memalign(&address,16,n*sizeof(float));
  if (status != 0) return 0; // Failed
  return address;
}
void freeDoubles(double * x)
{
  free(x);
}
void freeFloats(float * x)
{
  free(x);
}

// Reference implementations

void gemm_block_2x2_ref(const double * a,const double * b,long long int n,double * y,long long int ldy)
{
  int k;
  for (k=0;k<n;k+=2)
  {
    y[0] += a[2*k] * b[4*k];
    y[0] += a[2*k+1] * b[4*k+1];
    y[1] += a[2*k] * b[4*k+2];
    y[1] += a[2*k+1] * b[4*k+3];
    y[ldy] += a[2*k+2] * b[4*k];
    y[ldy] += a[2*k+3] * b[4*k+1];
    y[ldy+1] += a[2*k+2] * b[4*k+2];
    y[ldy+1] += a[2*k+3] * b[4*k+3];
  }
}

void gemm_block_2x4_ref(const double * a,const double * b,long long int n,double * y,long long int ldy)
{
  int k;
  double * work;
  for (k=0;k<n;k+=2)
  {
    work = y;
    work[0] += a[2*k]   * b[4*k];
    work[0] += a[2*k+1] * b[4*k+1];
    work += ldy;
    work[0] += a[2*k]   * b[4*k+2];
    work[0] += a[2*k+1] * b[4*k+3];
    work += ldy;
    work[0] += a[2*k]   * b[4*k+4];
    work[0] += a[2*k+1] * b[4*k+5];
    work += ldy;
    work[0] += a[2*k]   * b[4*k+6];
    work[0] += a[2*k+1] * b[4*k+7];

    work = y + 1;
    work[0] += a[2*k+2] * b[4*k];
    work[0] += a[2*k+3] * b[4*k+1];
    work += ldy;
    work[0] += a[2*k+2] * b[4*k+2];
    work[0] += a[2*k+3] * b[4*k+3];
    work += ldy;
    work[0] += a[2*k+2] * b[4*k+4];
    work[0] += a[2*k+3] * b[4*k+5];
    work += ldy;
    work[0] += a[2*k+2] * b[4*k+6];
    work[0] += a[2*k+3] * b[4*k+7];
  }
}

int dgemmsy_base(dgemmsyBaseArgs * args)
{
  int status = 0; // return value
  int row,col;
  size_t slice_size,ywork_size;
  int p2; // P rounded to next multiple of 4
  double *a_slice,*b_slice,*y_work;
  ComputeData cdata;
  TransposeData tdata;
  int slice_n,n,p;
  const double *a,*b;
  double *y;
  int lda,ldb,ldy;
  int transa,transb;
  BlockPattern_2x4_Proc pattern;
  void * pattern_arg;
  double alpha;

  if (args == 0) return -2;

  pattern = args->params.pattern;
  pattern_arg = args->params.pattern_arg;
  slice_n = args->params.slice_n;

  transa = args->transa;
  transb = args->transb;
  n = args->n;
  p = args->p;
  a = args->a;
  lda = args->lda;
  b = args->b;
  ldb = args->ldb;
  y = args->y;
  ldy = args->ldy;
  alpha = args->alpha;

  // Check dimensions
  if (slice_n < 8) return -1;
  if ( n <= 0 || p <= 0 ) return -1;
  p2 = (p+3) & 0x7FFFFFFC;

  // printf("N=%d P=%d LDA=%d LDB=%d LDY=%d  P2=%d\n",n,p,lda,ldb,ldy,p2);

  // Check other arguments
  if (a == 0 || b == 0 || y == 0 || pattern == 0) return -5;

  // Allocate memory
  slice_size = slice_n * p2 * sizeof(double); // Slice size in bytes
  ywork_size = p2 * p2 * sizeof(double); // Result size in bytes
  a_slice = (double *)allocPages(slice_size);
  b_slice = (double *)allocPages(slice_size);
  y_work = (double *)allocPages(ywork_size);
  if (a_slice == 0 || b_slice == 0 || y_work == 0) { status = -3; goto END; }
  memset(y_work,0,ywork_size);

  // Loop on all slices and accumulate products
  initComputeData(&cdata,p2,slice_n,a_slice,b_slice,y_work);
  for (row=0;row<n;row+=slice_n)
  {
    int s = slice_n;
    if (row+s > n) s = n-row; // Limit size of last slice if needed

    // 2-pack A slice
    if (transa) tpack_2(s,p,a+lda*row,lda, slice_n,p2,a_slice);
    else npack_2(s,p,a+row,lda, slice_n,p2,a_slice);

    // 4-pack B slice
    if (transb) tpack_4(s,p,b+ldb*row,ldb, slice_n,p2,b_slice);
    else npack_4(s,p,b+row,ldb, slice_n,p2,b_slice);

    pattern(pattern_arg,p2,Compute_visitor,&cdata);
  }
  cleanupComputeData(&cdata);

  // Complete result by symmetry
  initTransposeData(&tdata,p2,y_work);
  pattern(pattern_arg,p2,Transpose_visitor,&tdata);
  cleanupTransposeData(&tdata);

  // Combine and store (untransposed) result. If we are multithreading,
  // we must protect the update with a mutex, since all threads
  // will update the same Y.
  if (args->yMutex != 0)
    {
      int locked = pthread_mutex_lock(args->yMutex);
      if (locked != 0) status = -4;
    }
  for (col=0;col<p;col++)
    {
      double * yy = y+ldy*col;
      double * yy_work = y_work+p2*col;
#if USE_AXPY
      cblas_daxpy(p,alpha,yy_work,1,yy,1);
#else
      if (alpha == 1)
	{
	  for (row=0;row<p;row++) yy[row] += yy_work[row];
	}
      else
	{
	  for (row=0;row<p;row++) yy[row] += alpha * yy_work[row];
	}
#endif
    }
  if (args->yMutex != 0)
    {
      int unlocked = pthread_mutex_unlock(args->yMutex);
      if (unlocked != 0) status = -4;
    }

END:
  // Cleanup
  freePages(y_work);
  freePages(a_slice);
  freePages(b_slice);
  return status;
}

void npack_2(int n,int p,const double * a,int lda,int s,int q,double * slice)
{
  int row,col,icol;
  double * x = 0;
  int oddn = 0;
  if (n&1) { n--; oddn=1; }
  for (col=0;col<q;col++)
    {
      icol = col&1;
      x = slice + s*(col-icol) + (icol<<1);
      row = 0;
      if (col < p)
	{
	  for ( ;row<n;row+=2)
	    {
	      x[2*row] = a[row];
	      x[2*row+1] = a[row+1];
	    }
	  if (oddn) { x[2*row] = a[row]; x[2*row+1] = 0; row+=2; }
	}
      for ( ;row<s;row+=2)
	{
	  x[2*row] = 0;
	  x[2*row+1] = 0;
	}
      a += lda;
    }
}

void npack_4(int n,int p,const double * a,int lda,int s,int q,double * slice)
{
  int row,col,icol;
  double * x = 0;
  int oddn = 0;
  if (n&1) { n--; oddn=1; }
  for (col=0;col<q;col++)
    {
      icol = col&3;
      x = slice + s*(col-icol) + (icol<<1);
      row = 0;
      if (col < p)
	{
	  for ( ;row<n;row+=2)
	    {
	      x[4*row] = a[row];
	      x[4*row+1] = a[row+1];
	    }
	  if (oddn) { x[4*row] = a[row]; x[4*row+1] = 0; row+=2; }
	}
      for ( ;row<s;row+=2)
	{
	  x[4*row] = 0;
	  x[4*row+1] = 0;
	}
      a += lda;
    }
}

void tpack_2(int n,int p,const double * a,int lda,int s,int q,double * slice)
{
  int row,col,icol;
  double * x = 0;
  const double * aa;
  int da = lda<<1;
  int oddn = 0;
  if (n&1) { n--; oddn=1; }
  for (col=0;col<q;col++)
    {
      icol = col&1;
      x = slice + s*(col-icol) + (icol<<1);
      row = 0;
      if (col < p)
	{
	  aa = a;
	  for ( ;row<n;row+=2)
	    {
	      x[2*row] = aa[0];
	      x[2*row+1] = aa[lda];
	      aa += da;
	    }
	  if (oddn) { x[2*row] = aa[0]; x[2*row+1] = 0; row += 2; }
	}
      for ( ;row<s;row+=2)
	{
	  x[2*row] = 0;
	  x[2*row+1] = 0;
	}
      a++;
    }
}

void tpack_4(int n,int p,const double * a,int lda,int s,int q,double * slice)
{
  int row,col,icol;
  double * x = 0;
  const double * aa;
  int da = lda<<1;
  int oddn = 0;
  if (n&1) { n--; oddn=1; }
  for (col=0;col<q;col++)
    {
      icol = col&3;
      x = slice + s*(col-icol) + (icol<<1);
      row = 0;
      if (col < p)
	{
	  aa = a;
	  for ( ;row<n;row+=2)
	    {
	      x[4*row] = aa[0];
	      x[4*row+1] = aa[lda];
	      aa += da;
	    }
	  if (oddn) { x[4*row] = aa[0]; x[4*row+1] = 0; row += 2; }
	}
      for ( ;row<s;row+=2)
	{
	  x[4*row] = 0;
	  x[4*row+1] = 0;
	}
      a++;
    }
}
