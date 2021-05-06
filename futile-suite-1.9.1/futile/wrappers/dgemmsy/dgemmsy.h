// dgemmsy main include file
//
// Copyright (c) 2010, Commissariat a l'Energie Atomique
// Eric Bainville, Mar 2010

#ifndef dgemmsy_h
#define dgemmsy_h
#include <config.h>

// Compute Y' = ALPHA * op(A)^T*op(B) + BETA * Y
//
// op(A) and op(B) both have N rows and P columns, representing
// P vectors of dimension N.
//
// If TRANSA='N' or 'n', then op(A)=A, the increment in
// memory between consecutive components of a vector is 1,
// and between two vectors is LDA.
//
// If TRANSA='T' or 't', then op(A)=A^T, the increment in
// memory between consecutive components of a vector is LDA,
// and between two vectors is 1.
//
// The same storage convention applies to TRANSB and B.
// 
// Y is a symmetric (in input and in output) P*P matrix stored
// in column major order: Y(I,J)=Y[I-1+(J-1)*LDY].
//
// PARAMS is an internal parameter, must be set to 0 for general use.
//
// All coefficients of Y are updated.
void FC_FUNC(dgemmsy,DGEMMSY)(char *transa,char *transb,int *n,int *p,
	     double *alpha,const double * a,int *lda,const double * b,int *ldb,
	     double *beta,double * y,int *ldy);

// Same call, but dispatched on N_THREADS threads.
// On Windows, N_THREADS is ignored, and the single thread function will be called.
void dgemmsy_mt(int n_threads,
		char transa,char transb,int n,int p,
		double alpha,const double * a,int lda,const double * b,int ldb,
		double beta,double * y,int ldy,
		void * params);

#endif // #ifndef dgemmsy_h
