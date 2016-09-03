//! @file
//! dgemmsy entry points
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 
//! Eric Bainville, Mar 2010

#define _GNU_SOURCE
#include <pthread.h>
#include <malloc.h>
#include <string.h>
#include "utils.h"
#include "patterns.h"
#include "dgemmsy.h"

static BlocksPatternData BEST_PATTERN_ARG_MT = { 32,128,1 }; // na,nb,order
static BlocksPatternData BEST_PATTERN_ARG_1T = { 32,512,1 }; // na,nb,order

static void setBestParams(int n_threads,dgemmsyParams * params)
{
  params->pattern = blocks_pattern;
  if (n_threads == 1)
    params->pattern_arg = &BEST_PATTERN_ARG_1T;
  else
    params->pattern_arg = &BEST_PATTERN_ARG_MT;
  params->slice_n = 512;
}

// Scale Y(P,P,LDY) by BETA
static void scaleOutput(int p,double beta,double * y,int ldy)
{
  // Special cases
  if (beta == 0.0) // Clear memory
    {
      size_t sz = p * sizeof(y[0]); // Column size (bytes)
      int c;
      for (c=0;c<p;c++) memset(y+c*ldy,0,sz);
      return;
    }
  if (beta == 1.0) return; // Nothing to do if BETA is exactly 1

  // Otherwise scale Y (not optimized)
  int c,r;
  for (c=0;c<p;c++)
    {
      double * yy = y + c*ldy;
      for (r=0;r<p;r++) yy[r] *= beta;
    }
}

void FC_FUNC(dgemmsy,DGEMMSY)(char *transa,char *transb,int *n,int *p,
	     double *alpha,const double * a,int *lda,const double * b,int *ldb,
	     double *beta,double * y,int *ldy)
{
  int status;
  dgemmsyBaseArgs U;

  // Scale Y
  scaleOutput(*p,*beta,y,*ldy);

  // Setup args
  setBestParams(1,&(U.params));
  U.transa = (*transa == 'T' || *transa == 't');
  U.transb = (*transb == 'T' || *transb == 't');
  U.n = *n;
  U.p = *p;
  U.alpha = *alpha;
  U.a = a;
  U.lda = *lda;
  U.b = b;
  U.ldb = *ldb;
  U.y = y;
  U.ldy = *ldy;
  U.yMutex = 0;

  status = dgemmsy_base(&U);

  // ZZZ: do something in case of error?
  if (status < 0)
    {
      fprintf(stderr,"dgemmsy_base failed: error %d\n",status);
    }
}
/*
void * dgemmsy_base_mt(void * arg)
{
  dgemmsyBaseArgs * x = (dgemmsyBaseArgs *)arg;
  dgemmsy_base(x);
  return 0;
}

void dgemmsy_mt(int n_threads,
		char transa,char transb,int n,int p,
		double alpha,const double * a,int lda,const double * b,int ldb,
		double beta,double * y,int ldy,void * params)
{
  pthread_mutex_t yMutex = PTHREAD_MUTEX_INITIALIZER;
  dgemmsyBaseArgs *U;
  pthread_t *T;
  int rows_per_thread,row0,t;

  // Single thread
  if (n_threads <= 1)
    {
      FC_FUNC(dgemmsy,DGEMMSY)(transa,transb,n,p,alpha,a,lda,b,ldb,beta,y,ldy,params);
      return;
    }

  rows_per_thread = (n+n_threads-1) / n_threads; // rounded up

  // Scale Y
  scaleOutput(p,beta,y,ldy);

  // Setup args
  U = (dgemmsyBaseArgs *)malloc(n_threads * sizeof(U[0]));
  row0 = 0;
  for (t=0;t<n_threads;t++)
    {
      // Rows to compute
      int rows = n - row0;
      if (rows > rows_per_thread) rows = rows_per_thread;

      if (params != 0) U[t].params = *((dgemmsyParams *)params);
      else setBestParams(n_threads,&(U[t].params));
      U[t].transa = (transa == 'T' || transa == 't');
      U[t].transb = (transb == 'T' || transb == 't');
      U[t].n = rows;
      U[t].p = p;
      U[t].alpha = alpha;
      U[t].a = a;
      U[t].lda = lda;
      U[t].b = b;
      U[t].ldb = ldb;
      U[t].y = y;
      U[t].ldy = ldy;
      U[t].yMutex = &yMutex;

      if (U[t].transa) U[t].a += row0 * lda; else U[t].a += row0;
      if (U[t].transb) U[t].b += row0 * ldb; else U[t].b += row0;

      row0 += rows;
    }

  // Run threads
  T = (pthread_t *)malloc(n_threads * sizeof(T[0]));
  for (t=0;t<n_threads;t++)
    {
      pthread_create(T+t,0,dgemmsy_base_mt,U+t);

#if 0 // disabled: no gain observed with thread affinity
      cpu_set_t cpu_set;
      CPU_ZERO(&cpu_set);
      CPU_SET(t&3,&cpu_set);
      pthread_setaffinity_np(T[t],sizeof(cpu_set_t),&cpu_set);
#endif
    }

  // Wait for all threads to terminate
  for (t=0;t<n_threads;t++)
    {
      pthread_join(T[t],0);
    }

  // Cleanup
  free(U);
  free(T);
}
*/
