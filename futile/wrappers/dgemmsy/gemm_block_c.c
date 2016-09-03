//! @file
//!  Matrix-matrix multiplication (block 2x2)
//!
//! @author
//!    Copyright (C) 2009-2011 BigDFT group 
//!    This file is distributed under the terms of the
//!    GNU General Public License, see ~/COPYING file
//!    or http://www.gnu.org/copyleft/gpl.txt .
//!    For the list of contributors, see ~/AUTHORS 


#include <emmintrin.h>
#include <pmmintrin.h>
#include "gemm_block_c.h"

/*void gemm_block_2x2_c(const double * a,const double * b,long long int n,double * y,long long int ldy)
{
  __m128d A0,A1,B0,B1;
  __m128d S00,S01,S10,S11;
  S00 = _mm_setzero_pd();
  S01 = _mm_setzero_pd();
  S10 = _mm_setzero_pd();
  S11 = _mm_setzero_pd();
  unsigned int k;
  for (k=0;k<n;k+=2)
  {
    A0 = _mm_load_pd(a);
    a += 2;
    B0 = _mm_load_pd(b);
    b += 2;
    S00 = _mm_add_pd(S00,_mm_mul_pd(A0,B0));
    B1 = _mm_load_pd(b);
    b += 2;
    S01 = _mm_add_pd(S01,_mm_mul_pd(A0,B1));
    A1 = _mm_load_pd(a);
    a += 2;
    S10 = _mm_add_pd(S10,_mm_mul_pd(A1,B0));
    S11 = _mm_add_pd(S11,_mm_mul_pd(A1,B1));
//    y[0] += a[2*k] * b[4*k];
//    y[0] += a[2*k+1] * b[4*k+1];
//    y[1] += a[2*k] * b[4*k+2];
//    y[1] += a[2*k+1] * b[4*k+3];
//    y[ldy] += a[2*k+2] * b[4*k];
//    y[ldy] += a[2*k+3] * b[4*k+1];
//    y[ldy+1] += a[2*k+2] * b[4*k+2];
//    y[ldy+1] += a[2*k+3] * b[4*k+3];
  }
  A0 = _mm_load_pd(y);
  _mm_store_pd(y, _mm_add_pd(A0, _mm_hadd_pd(S00,S01)));
  A0 = _mm_load_pd(y+ldy);
  _mm_store_pd(y+ldy, _mm_add_pd(A0, _mm_hadd_pd(S10,S11)));
}*/

void gemm_block_2x4_c(const double * a,const double * b,long long int n,double * y,long long int ldy)
{
  __m128d A0,A1,B0,B1,B2,B3;
  __m128d S00,S01,S02,S03,S10,S11,S12,S13;
  S00 = _mm_setzero_pd();
  S01 = _mm_setzero_pd();
  S02 = _mm_setzero_pd();
  S03 = _mm_setzero_pd();
  S10 = _mm_setzero_pd();
  S11 = _mm_setzero_pd();
  S12 = _mm_setzero_pd();
  S13 = _mm_setzero_pd();
  unsigned long long int k=n>>3;
//  double * work;
  do {
    A0 = _mm_load_pd(a);
    a += 4;
    B0 = _mm_load_pd(b);
    S00 = _mm_add_pd(S00,_mm_mul_pd(A0,B0));
    B1 = _mm_load_pd(b+2);
    b += 8;
    S01 = _mm_add_pd(S01,_mm_mul_pd(A0,B1));
    B2 = _mm_load_pd(b-4);
    S02 = _mm_add_pd(S02,_mm_mul_pd(A0,B2));
    B3 = _mm_load_pd(b-2);
    S03 = _mm_add_pd(S03,_mm_mul_pd(A0,B3));
    A1 = _mm_load_pd(a-2);
    S10 = _mm_add_pd(S10,_mm_mul_pd(A1,B0));
    S11 = _mm_add_pd(S11,_mm_mul_pd(A1,B1));
    S12 = _mm_add_pd(S12,_mm_mul_pd(A1,B2));
    S13 = _mm_add_pd(S13,_mm_mul_pd(A1,B3));
    A0 = _mm_load_pd(a);
    a += 4;
    B0 = _mm_load_pd(b);
    S00 = _mm_add_pd(S00,_mm_mul_pd(A0,B0));
    B1 = _mm_load_pd(b+2);
    b += 8;
    S01 = _mm_add_pd(S01,_mm_mul_pd(A0,B1));
    B2 = _mm_load_pd(b-4);
    S02 = _mm_add_pd(S02,_mm_mul_pd(A0,B2));
    B3 = _mm_load_pd(b-2);
    S03 = _mm_add_pd(S03,_mm_mul_pd(A0,B3));
    A1 = _mm_load_pd(a-2);
    S10 = _mm_add_pd(S10,_mm_mul_pd(A1,B0));
    S11 = _mm_add_pd(S11,_mm_mul_pd(A1,B1));
    S12 = _mm_add_pd(S12,_mm_mul_pd(A1,B2));
    S13 = _mm_add_pd(S13,_mm_mul_pd(A1,B3));
    A0 = _mm_load_pd(a);
    a += 4;
    B0 = _mm_load_pd(b);
    S00 = _mm_add_pd(S00,_mm_mul_pd(A0,B0));
    B1 = _mm_load_pd(b+2);
    b += 8;
    S01 = _mm_add_pd(S01,_mm_mul_pd(A0,B1));
    B2 = _mm_load_pd(b-4);
    S02 = _mm_add_pd(S02,_mm_mul_pd(A0,B2));
    B3 = _mm_load_pd(b-2);
    S03 = _mm_add_pd(S03,_mm_mul_pd(A0,B3));
    A1 = _mm_load_pd(a-2);
    S10 = _mm_add_pd(S10,_mm_mul_pd(A1,B0));
    S11 = _mm_add_pd(S11,_mm_mul_pd(A1,B1));
    S12 = _mm_add_pd(S12,_mm_mul_pd(A1,B2));
    S13 = _mm_add_pd(S13,_mm_mul_pd(A1,B3));
    A0 = _mm_load_pd(a);
    a += 4;
    B0 = _mm_load_pd(b);
    S00 = _mm_add_pd(S00,_mm_mul_pd(A0,B0));
    B1 = _mm_load_pd(b+2);
    b += 8;
    S01 = _mm_add_pd(S01,_mm_mul_pd(A0,B1));
    B2 = _mm_load_pd(b-4);
    S02 = _mm_add_pd(S02,_mm_mul_pd(A0,B2));
    B3 = _mm_load_pd(b-2);
    S03 = _mm_add_pd(S03,_mm_mul_pd(A0,B3));
    A1 = _mm_load_pd(a-2);
    S10 = _mm_add_pd(S10,_mm_mul_pd(A1,B0));
    S11 = _mm_add_pd(S11,_mm_mul_pd(A1,B1));
    S12 = _mm_add_pd(S12,_mm_mul_pd(A1,B2));
    S13 = _mm_add_pd(S13,_mm_mul_pd(A1,B3));

  } while (--k>0);
  A0 = _mm_load_pd(y);
  _mm_store_pd(y, _mm_add_pd(A0, _mm_hadd_pd(S00,S10)));
  A0 = _mm_load_pd(y+ldy);
  _mm_store_pd(y+ldy, _mm_add_pd(A0, _mm_hadd_pd(S01,S11)));
  A0 = _mm_load_pd(y+2*ldy);
  _mm_store_pd(y+2*ldy, _mm_add_pd(A0, _mm_hadd_pd(S02,S12)));
  A0 = _mm_load_pd(y+3*ldy);
  _mm_store_pd(y+3*ldy, _mm_add_pd(A0, _mm_hadd_pd(S03,S13)));
}

