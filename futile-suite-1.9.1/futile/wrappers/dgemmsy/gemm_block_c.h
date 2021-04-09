#include <emmintrin.h>
#include <pmmintrin.h>

void gemm_block_2x4_c(const double * a,const double * b,long long int n,double * y,long long int ldy);
/*void gemm_block_2x2_c(const double * a,const double * b,long long int n,double * y,long long int ldy);*/
