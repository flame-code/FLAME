// Various functions used during experiences
//
// Copyright (c) 2010, Commissariat a l'Energie Atomique
// Eric Bainville, Mar 2010

#ifndef utils_h
#define utils_h

#include <pthread.h>

#include "visitors.h"
#include "patterns.h"

// Alloc/free system pages for SZ bytes.
// SZ is rounded up to the next multiple of the system page size.
// The returned address is multiple of the system page size.
// Return 0 if allocation failed.
void * allocPages(size_t sz);
void freePages(void * address);

// Alloc 16-byte aligned memory using the C runtime.
// Alloc/free aligned array of N double/float
double * allocDoubles(int n);
float * allocFloats(int n);
void freeDoubles(double * x);
void freeFloats(float * x);

// Parameters for high-level dgemmsy calls
typedef struct s_dgemmsyParams
{
  // Pattern to apply
  BlockPattern_2x4_Proc pattern;
  void * pattern_arg;
  // Slice size (rows)
  int slice_n;
} dgemmsyParams;

// Parameters of the base dgemmsy call
typedef struct s_dgemmsyBaseArgs
{
  // Pattern parameters
  dgemmsyParams params;

  // Memory representations of A and B (1=transpose, 0=normal)
  int transa,transb;
  // Dimensions of op(A) and op(B)
  int n,p;
  // A coefficient
  double alpha;
  // A and B
  const double * a;
  int lda;
  const double * b;
  int ldb;
  // The product is added to the passed Y.  Higher level functions
  // should handle the BETA coefficient.
  // Y
  double * y;
  int ldy;
  // Y access mutex
  pthread_mutex_t * yMutex;
} dgemmsyBaseArgs;

// dgemmsy low-level call, with a specific pattern and slice size.
// Return 0 if OK, and a <0 error code on failure:
// -1: invalid dimensions
// -2: invalid argument
// -3: memory allocation failed
// -4: mutex lock error
int dgemmsy_base(dgemmsyBaseArgs * args);

// 2-pack/4-pack normal matrix (N,P,A,LDA) N rows, P columns into SLICE, S rows, Q columns.
// S>=N multiple of 2, Q>=P multiple of 4.
void npack_2(int n,int p,const double * a,int lda,int s,int q,double * slice);
void npack_4(int n,int p,const double * a,int lda,int s,int q,double * slice);

// 2-pack/4-pack transpose matrix (N,P,A,LDA) N columns, P rows into SLICE, S rows, Q columns.
// S>=N multiple of 2, Q>=P multiple of 4.
void tpack_2(int n,int p,const double * a,int lda,int s,int q,double * slice);
void tpack_4(int n,int p,const double * a,int lda,int s,int q,double * slice);

// ASM functions.

// gemm_block_UxV(a,b,n,work) processes n coordinates from packed
// columns a[U*n] and b[V*n] into work[U*V].
// unsigned long long eRDTSC();
// void gemm_block_2x2(const double * a,const double * b,long long int n,double * y,long long int ldy);
// void gemm_block_2x4(const double * a,const double * b,long long int n,double * y,long long int ldy);

// Reference implementations in C (debug)
void gemm_block_2x2_ref(const double * a,const double * b,long long int n,double * y,long long int ldy);
void gemm_block_2x4_ref(const double * a,const double * b,long long int n,double * y,long long int ldy);

#endif // #ifndef utils_h
