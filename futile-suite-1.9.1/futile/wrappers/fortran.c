/*!
 * Change T.D. 2011-11-18: pb of dereferencing type-punned
 *
 * Copyright 1993-2008 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO USER:   
 *
 * This source code is subject to NVIDIA ownership rights under U.S. and 
 * international Copyright laws.  Users and possessors of this source code 
 * are hereby granted a nonexclusive, royalty-free license to use this code 
 * in individual and commercial software.
 *
 * NVIDIA MAKES NO REPRESENTATION ABOUT THE SUITABILITY OF THIS SOURCE 
 * CODE FOR ANY PURPOSE.  IT IS PROVIDED "AS IS" WITHOUT EXPRESS OR 
 * IMPLIED WARRANTY OF ANY KIND.  NVIDIA DISCLAIMS ALL WARRANTIES WITH 
 * REGARD TO THIS SOURCE CODE, INCLUDING ALL IMPLIED WARRANTIES OF 
 * MERCHANTABILITY, NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY SPECIAL, INDIRECT, INCIDENTAL, 
 * OR CONSEQUENTIAL DAMAGES, OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS 
 * OF USE, DATA OR PROFITS,  WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE 
 * OR OTHER TORTIOUS ACTION,  ARISING OUT OF OR IN CONNECTION WITH THE USE 
 * OR PERFORMANCE OF THIS SOURCE CODE.  
 *
 * U.S. Government End Users.   This source code is a "commercial item" as 
 * that term is defined at  48 C.F.R. 2.101 (OCT 1995), consisting  of 
 * "commercial computer  software"  and "commercial computer software 
 * documentation" as such terms are  used in 48 C.F.R. 12.212 (SEPT 1995) 
 * and is provided to the U.S. Government only as a commercial end item.  
 * Consistent with 48 C.F.R.12.212 and 48 C.F.R. 227.7202-1 through 
 * 227.7202-4 (JUNE 1995), all U.S. Government End Users acquire the 
 * source code with only those rights set forth herein. 
 *
 * Any use of this source code in individual and commercial software must 
 * include, in the user documentation and internal comments to the code,
 * the above Disclaimer and U.S. Government End Users Notice.
 */

/*
 * This file contains example Fortran bindings for the CUBLAS library, These
 * bindings have been tested with Intel Fortran 9.0 on 32-bit and 64-bit 
 * Windows, and with g77 3.4.5 on 32-bit and 64-bit Linux. They will likely
 * have to be adjusted for other Fortran compilers and platforms.
 */

#include <ctype.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <stdlib.h>
#if defined(__GNUC__)
#include <stdint.h>
#endif /* __GNUC__ */
#include "cublas.h"   /* CUBLAS public header file  */
#include <config.h>

#define imin(a,b) (((a)<(b))?(a):(b))
#define imax(a,b) (((a)<(b))?(b):(a))

#define CUBLAS_G77              1
#define CUBLAS_INTEL_FORTRAN    2
#define CUBLAS_ONE_UNDERSCORE   3
/* Default to g77 on Linux, and Intel Fortran on Win32 */
#if defined(_WIN32)
#define CUBLAS_FORTRAN_COMPILER CUBLAS_INTEL_FORTRAN
#elif defined(__linux)
#define CUBLAS_FORTRAN_COMPILER CUBLAS_G77
#elif defined(__APPLE__)
#define CUBLAS_FORTRAN_COMPILER CUBLAS_G77
#else
#error unsupported platform
#endif


/*#define CUBLAS_FORTRAN_COMPILER CUBLAS_ONE_UNDERSCORE*/


#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
/* NOTE: Must use -fno-second-underscore when building Fortran source with g77
 *       g77 invocation may not use -fno-f2c, which forces different return 
 *       type conventions than the one used below
 */
#define CUBLAS_INIT             cublas_init__
#define CUBLAS_SHUTDOWN         cublas_shutdown__
#define CUBLAS_ALLOC            cublas_alloc__
#define CUBLAS_FREE             cublas_free__
#define CUBLAS_SET_VECTOR       cublas_set_vector__
#define CUBLAS_GET_VECTOR       cublas_get_vector__
#define CUBLAS_SET_MATRIX       cublas_set_matrix__
#define CUBLAS_GET_MATRIX       cublas_get_matrix__
#define CUBLAS_GET_ERROR        cublas_get_error__
#define CUBLAS_XERBLA           cublas_xerbla__
#define CUBLAS_ISAMAX           cublas_isamax__
#define CUBLAS_ISAMIN           cublas_isamin__
#define CUBLAS_SASUM            cublas_sasum__
#define CUBLAS_SAXPY            cublas_saxpy__
#define CUBLAS_SCOPY            cublas_scopy__
#define CUBLAS_SDOT             cublas_sdot__
#define CUBLAS_SNRM2            cublas_snrm2__
#define CUBLAS_SROT             cublas_srot__
#define CUBLAS_SROTG            cublas_srotg__
#define CUBLAS_SROTM            cublas_srotm__
#define CUBLAS_SROTMG           cublas_srotmg__
#define CUBLAS_SSCAL            cublas_sscal__
#define CUBLAS_SSWAP            cublas_sswap__
#define CUBLAS_CAXPY            cublas_caxpy__
#define CUBLAS_CCOPY            cublas_ccopy__
#define CUBLAS_CROT             cublas_crot__
#define CUBLAS_CROTG            cublas_crotg__
#define CUBLAS_CSCAL            cublas_cscal__
#define CUBLAS_CSROT            cublas_csrot__
#define CUBLAS_CSSCAL           cublas_csscal__
#define CUBLAS_CSWAP            cublas_cswap__
#define CUBLAS_CDOTU            cublas_cdotu__
#define CUBLAS_CDOTC            cublas_cdotc__
#define CUBLAS_ICAMAX           cublas_icamax__
#define CUBLAS_SCASUM           cublas_scasum__
#define CUBLAS_SCNRM2           cublas_scnrm2__
#define CUBLAS_SGBMV            cublas_sgbmv__
#define CUBLAS_SGEMV            cublas_sgemv__
#define CUBLAS_SGER             cublas_sger__
#define CUBLAS_SSBMV            cublas_ssbmv__
#define CUBLAS_SSPMV            cublas_sspmv__
#define CUBLAS_SSPR             cublas_sspr__
#define CUBLAS_SSPR2            cublas_sspr2__
#define CUBLAS_SSYMV            cublas_ssymv__
#define CUBLAS_SSYR             cublas_ssyr__
#define CUBLAS_SSYR2            cublas_ssyr2__
#define CUBLAS_STBMV            cublas_stbmv__
#define CUBLAS_STBSV            cublas_stbsv__
#define CUBLAS_STPMV            cublas_stpmv__
#define CUBLAS_STPSV            cublas_stpsv__
#define CUBLAS_STRMV            cublas_strmv__
#define CUBLAS_STRSV            cublas_strsv__
#define CUBLAS_SGEMM            cublas_sgemm__
#define CUBLAS_SSYMM            cublas_ssymm__
#define CUBLAS_SSYR2K           cublas_ssyr2k__
#define CUBLAS_SSYRK            cublas_ssyrk__
#define CUBLAS_STRMM            cublas_strmm__
#define CUBLAS_STRSM            cublas_strsm__
#define CUBLAS_CGEMM            cublas_cgemm__
#define CUBLAS_CHEMM            cublas_chemm__
#define CUBLAS_CSYMM            cublas_csymm__
#define CUBLAS_CTRMM            cublas_ctrmm__
#define CUBLAS_CTRSM            cublas_ctrsm__
#define CUBLAS_CHERK            cublas_cherk__
#define CUBLAS_CSYRK            cublas_csyrk__
#define CUBLAS_CHER2K           cublas_cher2k__
#define CUBLAS_CSYR2K           cublas_csyr2k__
#define CUBLAS_IDAMAX           cublas_idamax__
#define CUBLAS_IDAMIN           cublas_idamin__
#define CUBLAS_DASUM            cublas_dasum__
#define CUBLAS_DAXPY            cublas_daxpy__
#define CUBLAS_DCOPY            cublas_dcopy__
#define CUBLAS_DDOT             cublas_ddot__
#define CUBLAS_DNRM2            cublas_dnrm2__
#define CUBLAS_DROT             cublas_drot__
#define CUBLAS_DROTG            cublas_drotg__
#define CUBLAS_DROTM            cublas_drotm__
#define CUBLAS_DROTMG           cublas_drotmg__
#define CUBLAS_DSCAL            cublas_dscal__
#define CUBLAS_DSWAP            cublas_dswap__
#define CUBLAS_ZAXPY            cublas_zaxpy__
#define CUBLAS_ZCOPY            cublas_zcopy__
#define CUBLAS_ZROT             cublas_zrot__
#define CUBLAS_ZROTG            cublas_zrotg__
#define CUBLAS_ZSCAL            cublas_zscal__
#define CUBLAS_ZDROT            cublas_zdrot__
#define CUBLAS_ZDSCAL           cublas_zdscal__
#define CUBLAS_ZSWAP            cublas_zswap__
#define CUBLAS_ZDOTU            cublas_zdotu__
#define CUBLAS_ZDOTC            cublas_zdotc__
#define CUBLAS_IZAMAX           cublas_izamax__
#define CUBLAS_DZASUM           cublas_dzasum__
#define CUBLAS_DZNRM2           cublas_dznrm2__
#define CUBLAS_DGBMV            cublas_dgbmv__
#define CUBLAS_DGEMV            cublas_dgemv__
#define CUBLAS_DGER             cublas_dger__
#define CUBLAS_DSBMV            cublas_dsbmv__
#define CUBLAS_DSPMV            cublas_dspmv__
#define CUBLAS_DSPR             cublas_dspr__
#define CUBLAS_DSPR2            cublas_dspr2__
#define CUBLAS_DSYMV            cublas_dsymv__
#define CUBLAS_DSYR             cublas_dsyr__
#define CUBLAS_DSYR2            cublas_dsyr2__
#define CUBLAS_DTBMV            cublas_dtbmv__
#define CUBLAS_DTBSV            cublas_dtbsv__
#define CUBLAS_DTPMV            cublas_dtpmv__
#define CUBLAS_DTPSV            cublas_dtpsv__
#define CUBLAS_DTRMV            cublas_dtrmv__
#define CUBLAS_DTRSV            cublas_dtrsv__
#define CUBLAS_DGEMM            cublas_dgemm__
#define CUBLAS_DSYMM            cublas_dsymm__
#define CUBLAS_DSYR2K           cublas_dsyr2k__
#define CUBLAS_DSYRK            cublas_dsyrk__
#define CUBLAS_DTRMM            cublas_dtrmm__
#define CUBLAS_DTRSM            cublas_dtrsm__
#define CUBLAS_ZGEMM            cublas_zgemm__
#define CUBLAS_ZHEMM            cublas_zhemm__
#define CUBLAS_ZSYMM            cublas_zsymm__
#define CUBLAS_ZTRMM            cublas_ztrmm__
#define CUBLAS_ZTRSM            cublas_ztrsm__
#define CUBLAS_ZHERK            cublas_zherk__
#define CUBLAS_ZSYRK            cublas_zsyrk__
#define CUBLAS_ZHER2K           cublas_zher2k__
#define CUBLAS_ZSYR2K           cublas_zsyr2k__

#elif CUBLAS_FORTRAN_COMPILER==CUBLAS_ONE_UNDERSCORE

#define CUBLAS_INIT             cublas_init_
#define CUBLAS_SHUTDOWN         cublas_shutdown_
#define CUBLAS_ALLOC            cublas_alloc_
#define CUBLAS_FREE             cublas_free_
#define CUBLAS_SET_VECTOR       cublas_set_vector_
#define CUBLAS_GET_VECTOR       cublas_get_vector_
#define CUBLAS_SET_MATRIX       cublas_set_matrix_
#define CUBLAS_GET_MATRIX       cublas_get_matrix_
#define CUBLAS_GET_ERROR        cublas_get_error_
#define CUBLAS_XERBLA           cublas_xerbla_
#define CUBLAS_ISAMAX           cublas_isamax_
#define CUBLAS_ISAMIN           cublas_isamin_
#define CUBLAS_SASUM            cublas_sasum_
#define CUBLAS_SAXPY            cublas_saxpy_
#define CUBLAS_SCOPY            cublas_scopy_
#define CUBLAS_SDOT             cublas_sdot_
#define CUBLAS_SNRM2            cublas_snrm2_
#define CUBLAS_SROT             cublas_srot_
#define CUBLAS_SROTG            cublas_srotg_
#define CUBLAS_SROTM            cublas_srotm_
#define CUBLAS_SROTMG           cublas_srotmg_
#define CUBLAS_SSCAL            cublas_sscal_
#define CUBLAS_SSWAP            cublas_sswap_
#define CUBLAS_CAXPY            cublas_caxpy_
#define CUBLAS_CCOPY            cublas_ccopy_
#define CUBLAS_CROT             cublas_crot_
#define CUBLAS_CROTG            cublas_crotg_
#define CUBLAS_CSCAL            cublas_cscal_
#define CUBLAS_CSROT            cublas_csrot_
#define CUBLAS_CSSCAL           cublas_csscal_
#define CUBLAS_CSWAP            cublas_cswap_
#define CUBLAS_CDOTU            cublas_cdotu_
#define CUBLAS_CDOTC            cublas_cdotc_
#define CUBLAS_ICAMAX           cublas_icamax_
#define CUBLAS_SCASUM           cublas_scasum_
#define CUBLAS_SCNRM2           cublas_scnrm2_
#define CUBLAS_SGBMV            cublas_sgbmv_
#define CUBLAS_SGEMV            cublas_sgemv_
#define CUBLAS_SGER             cublas_sger_
#define CUBLAS_SSBMV            cublas_ssbmv_
#define CUBLAS_SSPMV            cublas_sspmv_
#define CUBLAS_SSPR             cublas_sspr_
#define CUBLAS_SSPR2            cublas_sspr2_
#define CUBLAS_SSYMV            cublas_ssymv_
#define CUBLAS_SSYR             cublas_ssyr_
#define CUBLAS_SSYR2            cublas_ssyr2_
#define CUBLAS_STBMV            cublas_stbmv_
#define CUBLAS_STBSV            cublas_stbsv_
#define CUBLAS_STPMV            cublas_stpmv_
#define CUBLAS_STPSV            cublas_stpsv_
#define CUBLAS_STRMV            cublas_strmv_
#define CUBLAS_STRSV            cublas_strsv_
#define CUBLAS_SGEMM            cublas_sgemm_
#define CUBLAS_SSYMM            cublas_ssymm_
#define CUBLAS_SSYR2K           cublas_ssyr2k_
#define CUBLAS_SSYRK            cublas_ssyrk_
#define CUBLAS_STRMM            cublas_strmm_
#define CUBLAS_STRSM            cublas_strsm_
#define CUBLAS_CGEMM            cublas_cgemm_
#define CUBLAS_CHEMM            cublas_chemm_
#define CUBLAS_CSYMM            cublas_csymm_
#define CUBLAS_CTRMM            cublas_ctrmm_
#define CUBLAS_CTRSM            cublas_ctrsm_
#define CUBLAS_CHERK            cublas_cherk_
#define CUBLAS_CSYRK            cublas_csyrk_
#define CUBLAS_CHER2K           cublas_cher2k_
#define CUBLAS_CSYR2K           cublas_csyr2k_
#define CUBLAS_IDAMAX           cublas_idamax_
#define CUBLAS_IDAMIN           cublas_idamin_
#define CUBLAS_DASUM            cublas_dasum_
#define CUBLAS_DAXPY            cublas_daxpy_
#define CUBLAS_DCOPY            cublas_dcopy_
#define CUBLAS_DDOT             cublas_ddot_
#define CUBLAS_DNRM2            cublas_dnrm2_
#define CUBLAS_DROT             cublas_drot_
#define CUBLAS_DROTG            cublas_drotg_
#define CUBLAS_DROTM            cublas_drotm_
#define CUBLAS_DROTMG           cublas_drotmg_
#define CUBLAS_DSCAL            cublas_dscal_
#define CUBLAS_DSWAP            cublas_dswap_
#define CUBLAS_ZAXPY            cublas_zaxpy_
#define CUBLAS_ZCOPY            cublas_zcopy_
#define CUBLAS_ZROT             cublas_zrot_
#define CUBLAS_ZROTG            cublas_zrotg_
#define CUBLAS_ZSCAL            cublas_zscal_
#define CUBLAS_ZDROT            cublas_zdrot_
#define CUBLAS_ZDSCAL           cublas_zdscal_
#define CUBLAS_ZSWAP            cublas_zswap_
#define CUBLAS_ZDOTU            cublas_zdotu_
#define CUBLAS_ZDOTC            cublas_zdotc_
#define CUBLAS_IZAMAX           cublas_izamax_
#define CUBLAS_DZASUM           cublas_dzasum_
#define CUBLAS_DZNRM2           cublas_dznrm2_
#define CUBLAS_DGBMV            cublas_dgbmv_
#define CUBLAS_DGEMV            cublas_dgemv_
#define CUBLAS_DGER             cublas_dger_
#define CUBLAS_DSBMV            cublas_dsbmv_
#define CUBLAS_DSPMV            cublas_dspmv_
#define CUBLAS_DSPR             cublas_dspr_
#define CUBLAS_DSPR2            cublas_dspr2_
#define CUBLAS_DSYMV            cublas_dsymv_
#define CUBLAS_DSYR             cublas_dsyr_
#define CUBLAS_DSYR2            cublas_dsyr2_
#define CUBLAS_DTBMV            cublas_dtbmv_
#define CUBLAS_DTBSV            cublas_dtbsv_
#define CUBLAS_DTPMV            cublas_dtpmv_
#define CUBLAS_DTPSV            cublas_dtpsv_
#define CUBLAS_DTRMV            cublas_dtrmv_
#define CUBLAS_DTRSV            cublas_dtrsv_
#define CUBLAS_DGEMM            cublas_dgemm_
#define CUBLAS_DSYMM            cublas_dsymm_
#define CUBLAS_DSYR2K           cublas_dsyr2k_
#define CUBLAS_DSYRK            cublas_dsyrk_
#define CUBLAS_DTRMM            cublas_dtrmm_
#define CUBLAS_DTRSM            cublas_dtrsm_
#define CUBLAS_ZGEMM            cublas_zgemm_
#define CUBLAS_ZHEMM            cublas_zhemm_
#define CUBLAS_ZSYMM            cublas_zsymm_
#define CUBLAS_ZTRMM            cublas_ztrmm_
#define CUBLAS_ZTRSM            cublas_ztrsm_
#define CUBLAS_ZHERK            cublas_zherk_
#define CUBLAS_ZSYRK            cublas_zsyrk_
#define CUBLAS_ZHER2K           cublas_zher2k_
#define CUBLAS_ZSYR2K           cublas_zsyr2k_

#elif CUBLAS_FORTRAN_COMPILER==CUBLAS_INTEL_FORTRAN

#define CUBLAS_INIT             CUBLAS_INIT 
#define CUBLAS_SHUTDOWN         CUBLAS_SHUTDOWN
#define CUBLAS_ALLOC            CUBLAS_ALLOC
#define CUBLAS_FREE             CUBLAS_FREE
#define CUBLAS_SET_VECTOR       CUBLAS_SET_VECTOR
#define CUBLAS_GET_VECTOR       CUBLAS_GET_VECTOR
#define CUBLAS_SET_MATRIX       CUBLAS_SET_MATRIX
#define CUBLAS_GET_MATRIX       CUBLAS_GET_MATRIX
#define CUBLAS_GET_ERROR        CUBLAS_GET_ERROR
#define CUBLAS_XERBLA           CUBLAS_XERBLA
#define CUBLAS_ISAMAX           CUBLAS_ISAMAX
#define CUBLAS_ISAMIN           CUBLAS_ISAMIN
#define CUBLAS_SASUM            CUBLAS_SASUM
#define CUBLAS_SAXPY            CUBLAS_SAXPY
#define CUBLAS_SCOPY            CUBLAS_SCOPY
#define CUBLAS_SDOT             CUBLAS_SDOT
#define CUBLAS_SNRM2            CUBLAS_SNRM2
#define CUBLAS_SROT             CUBLAS_SROT
#define CUBLAS_SROTG            CUBLAS_SROTG
#define CUBLAS_SROTM            CUBLAS_SROTM
#define CUBLAS_SROTMG           CUBLAS_SROTMG
#define CUBLAS_SSCAL            CUBLAS_SSCAL
#define CUBLAS_SSWAP            CUBLAS_SSWAP
#define CUBLAS_CAXPY            CUBLAS_CAXPY
#define CUBLAS_CCOPY            CUBLAS_CCOPY
#define CUBLAS_CROT             CUBLAS_CROT
#define CUBLAS_CROTG            CUBLAS_CROTG
#define CUBLAS_CSCAL            CUBLAS_CSCAL
#define CUBLAS_CSROT            CUBLAS_CSROT
#define CUBLAS_CSSCAL           CUBLAS_CSSCAL
#define CUBLAS_CSWAP            CUBLAS_CSWAP 
#define CUBLAS_CDOTU            CUBLAS_CDOTU
#define CUBLAS_CDOTC            CUBLAS_CDOTC
#define CUBLAS_ICAMAX           CUBLAS_ICAMAX
#define CUBLAS_SCASUM           CUBLAS_SCASUM
#define CUBLAS_SCNRM2           CUBLAS_SCNRM2
#define CUBLAS_SGBMV            CUBLAS_SGBMV
#define CUBLAS_SGEMV            CUBLAS_SGEMV
#define CUBLAS_SGER             CUBLAS_SGER
#define CUBLAS_SSBMV            CUBLAS_SSBMV
#define CUBLAS_SSPMV            CUBLAS_SSPMV
#define CUBLAS_SSPR             CUBLAS_SSPR
#define CUBLAS_SSPR2            CUBLAS_SSPR2
#define CUBLAS_SSYMV            CUBLAS_SSYMV
#define CUBLAS_SSYR             CUBLAS_SSYR
#define CUBLAS_SSYR2            CUBLAS_SSYR2
#define CUBLAS_STBMV            CUBLAS_STBMV
#define CUBLAS_STBSV            CUBLAS_STBSV
#define CUBLAS_STPMV            CUBLAS_STPMV
#define CUBLAS_STPSV            CUBLAS_STPSV
#define CUBLAS_STRMV            CUBLAS_STRMV
#define CUBLAS_STRSV            CUBLAS_STRSV
#define CUBLAS_SGEMM            CUBLAS_SGEMM
#define CUBLAS_SSYMM            CUBLAS_SSYMM
#define CUBLAS_SSYR2K           CUBLAS_SSYR2K
#define CUBLAS_SSYRK            CUBLAS_SSYRK
#define CUBLAS_STRMM            CUBLAS_STRMM
#define CUBLAS_STRSM            CUBLAS_STRSM
#define CUBLAS_CGEMM            CUBLAS_CGEMM
#define CUBLAS_CHEMM            CUBLAS_CHEMM
#define CUBLAS_CSYMM            CUBLAS_CSYMM
#define CUBLAS_CTRMM            CUBLAS_CTRMM
#define CUBLAS_CTRSM            CUBLAS_CTRSM
#define CUBLAS_CHERK            CUBLAS_CHERK
#define CUBLAS_CSYRK            CUBLAS_CSYRK
#define CUBLAS_CHER2K           CUBLAS_CHER2K
#define CUBLAS_CSYR2K           CUBLAS_CSYR2K
#define CUBLAS_IDAMAX           CUBLAS_IDAMAX
#define CUBLAS_IDAMIN           CUBLAS_IDAMIN
#define CUBLAS_DASUM            CUBLAS_DASUM
#define CUBLAS_DAXPY            CUBLAS_DAXPY
#define CUBLAS_DCOPY            CUBLAS_DCOPY
#define CUBLAS_DDOT             CUBLAS_DDOT
#define CUBLAS_DNRM2            CUBLAS_DNRM2
#define CUBLAS_DROT             CUBLAS_DROT
#define CUBLAS_DROTG            CUBLAS_DROTG
#define CUBLAS_DROTM            CUBLAS_DROTM
#define CUBLAS_DROTMG           CUBLAS_DROTMG
#define CUBLAS_DSCAL            CUBLAS_DSCAL
#define CUBLAS_DSWAP            CUBLAS_DSWAP
#define CUBLAS_ZAXPY            CUBLAS_ZAXPY
#define CUBLAS_ZCOPY            CUBLAS_ZCOPY
#define CUBLAS_ZROT             CUBLAS_ZROT
#define CUBLAS_ZROTG            CUBLAS_ZROTG
#define CUBLAS_ZSCAL            CUBLAS_ZSCAL
#define CUBLAS_ZDROT            CUBLAS_ZDROT
#define CUBLAS_ZDSCAL           CUBLAS_ZDSCAL
#define CUBLAS_ZSWAP            CUBLAS_ZSWAP 
#define CUBLAS_ZDOTU            CUBLAS_ZDOTU
#define CUBLAS_ZDOTC            CUBLAS_ZDOTC
#define CUBLAS_IZAMAX           CUBLAS_IZAMAX
#define CUBLAS_DZASUM           CUBLAS_DZASUM
#define CUBLAS_DZNRM2           CUBLAS_DZNRM2
#define CUBLAS_DGBMV            CUBLAS_DGBMV
#define CUBLAS_DGEMV            CUBLAS_DGEMV
#define CUBLAS_DGER             CUBLAS_DGER
#define CUBLAS_DSBMV            CUBLAS_DSBMV
#define CUBLAS_DSPMV            CUBLAS_DSPMV
#define CUBLAS_DSPR             CUBLAS_DSPR
#define CUBLAS_DSPR2            CUBLAS_DSPR2
#define CUBLAS_DSYMV            CUBLAS_DSYMV
#define CUBLAS_DSYR             CUBLAS_DSYR
#define CUBLAS_DSYR2            CUBLAS_DSYR2
#define CUBLAS_DTBMV            CUBLAS_DTBMV
#define CUBLAS_DTBSV            CUBLAS_DTBSV
#define CUBLAS_DTPMV            CUBLAS_DTPMV
#define CUBLAS_DTPSV            CUBLAS_DTPSV
#define CUBLAS_DTRMV            CUBLAS_DTRMV
#define CUBLAS_DTRSV            CUBLAS_DTRSV
#define CUBLAS_DGEMM            CUBLAS_DGEMM
#define CUBLAS_DSYMM            CUBLAS_DSYMM
#define CUBLAS_DSYR2K           CUBLAS_DSYR2K
#define CUBLAS_DSYRK            CUBLAS_DSYRK
#define CUBLAS_DTRMM            CUBLAS_DTRMM
#define CUBLAS_DTRSM            CUBLAS_DTRSM
#define CUBLAS_ZGEMM            CUBLAS_ZGEMM
#define CUBLAS_ZHEMM            CUBLAS_ZHEMM
#define CUBLAS_ZSYMM            CUBLAS_ZSYMM
#define CUBLAS_ZTRMM            CUBLAS_ZTRMM
#define CUBLAS_ZTRSM            CUBLAS_ZTRSM
#define CUBLAS_ZHERK            CUBLAS_ZHERK
#define CUBLAS_ZSYRK            CUBLAS_ZSYRK
#define CUBLAS_ZHER2K           CUBLAS_ZHER2K
#define CUBLAS_ZSYR2K           CUBLAS_ZSYR2K

#else
#error unsupported Fortran compiler
#endif

/* For now, the GPU only supports a 32-bit address space, so device pointers
   can be represented as INTEGER*4 in Fortran. In the future, device pointers
   may become 64-bit pointers, and will have to be represented as INTEGER*8 in
   Fortran, at which point devptr_t needs to be typedef'ed as long long.
*/
typedef int devptr_t;

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
int FC_FUNC_(cublas_init,CUBLAS_INIT) (void);
int FC_FUNC_(cublas_shutdown,CUBLAS_SHUTDOWN) (void);
int FC_FUNC_(cublas_alloc,CUBLAS_ALLOC) (const int *n, const int *elemSize, devptr_t *devicePtr);
int FC_FUNC_(cublas_free,CUBLAS_FREE) (const devptr_t *devicePtr);
int FC_FUNC_(cublas_set_vector,CUBLAS_SET_VECTOR) (const int *n, const int *elemSize, const void *x,
                       const int *incx, const devptr_t *y, const int *incy);
int FC_FUNC_(cublas_get_vector,CUBLAS_GET_VECTOR) (const int *n, const int *elemSize, const devptr_t *x,
                       const int *incx, void *y, const int *incy);
int FC_FUNC_(cublas_set_matrix,CUBLAS_SET_MATRIX) (const int *rows, const int *cols, const int *elemSize,
                       const void *A, const int *lda, const int *B, 
                       const int *ldb);
int FC_FUNC_(cublas_get_matrix,CUBLAS_GET_MATRIX) (const int *rows, const int *cols, const int *elemSize,
                       const int *A, const int *lda, void *B, const int *ldb);

/* BLAS util */
void FC_FUNC_(cublas_xerbla,CUBLAS_XERBLA) (const char *srName, int *info);
#if defined(__cplusplus)
}
#endif /* __cplusplus */

int FC_FUNC_(cublas_init,CUBLAS_INIT) (void) 
{
    return (int)cublasInit ();
}

int FC_FUNC_(cublas_shutdown,CUBLAS_SHUTDOWN) (void) 
{
    return (int)cublasShutdown ();
}

int FC_FUNC_(cublas_alloc,CUBLAS_ALLOC) (const int *n, const int *elemSize, devptr_t *devicePtr)
{    
    void *tPtr;
    int retVal;
    retVal = (int)cublasAlloc (*n, *elemSize, &tPtr);
    *devicePtr = (devptr_t)(uintptr_t)tPtr;
    return retVal;
}

int FC_FUNC_(cublas_free,CUBLAS_FREE) (const devptr_t *devicePtr)
{
    void *tPtr;
    tPtr = (void *)(uintptr_t)(*devicePtr);
    return (int)cublasFree (tPtr);
}

int FC_FUNC_(cublas_set_vector,CUBLAS_SET_VECTOR) (const int *n, const int *elemSize, const void *x,
                       const int *incx, const devptr_t *y, const int *incy)
{
    void *tPtr = (void *)(uintptr_t)(*y);
    return (int)cublasSetVector (*n, *elemSize, x, *incx, tPtr, *incy);
}

int FC_FUNC_(cublas_get_vector,CUBLAS_GET_VECTOR) (const int *n, const int *elemSize, const devptr_t *x,
                       const int *incx, void *y, const int *incy)
{
    const void *tPtr = (const void *)(uintptr_t)(*x);
    return (int)cublasGetVector (*n, *elemSize, tPtr, *incx, y, *incy);
}

int FC_FUNC_(cublas_set_matrix,CUBLAS_SET_MATRIX) (const int *rows, const int *cols, const int *elemSize,
                       const void *A, const int *lda, const devptr_t *B, 
                       const int *ldb)
{
    void *tPtr = (void *)(uintptr_t)(*B);
    return (int)cublasSetMatrix (*rows, *cols, *elemSize, A, *lda, tPtr,*ldb);
}

int FC_FUNC_(cublas_get_matrix,CUBLAS_GET_MATRIX) (const int *rows, const int *cols, const int *elemSize,
                       const devptr_t *A, const int *lda, void *B, 
                       const int *ldb)
{
    const void *tPtr = (const void *)(uintptr_t)(*A);
    return (int)cublasGetMatrix (*rows, *cols, *elemSize, tPtr, *lda, B, *ldb);
}

int FC_FUNC_(cublas_get_error,CUBLAS_GET_ERROR) (void)
{
    return (int)cublasGetError();
}

void FC_FUNC_(cublas_xerbla,CUBLAS_XERBLA) (const char *srName, int *info)
{
    cublasXerbla (srName, *info);
}

/*
 *  Fortran callable BLAS functions that include GPU memory allocation and
 *  copy-up and copy-down code. These can be called from unmodified Fortran 
 *  code, but they are inefficient due to the data constantly bouncing back 
 *  and forth between CPU and GPU.
 */
#if defined(CUBLAS_USE_THUNKING)

#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
/* BLAS1 */
#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const float *x, const int *incx, float *y, 
                    const int *incy);
double FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const float *x, const int *incx);
double FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const float *x, const int *incx);
double FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const cuComplex *x, const int *incx);
double FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const cuComplex *x, const int *incx);
#else
float FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const float *x, const int *incx, float *y, 
                   const int *incy);
float FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const float *x, const int *incx);
float FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const float *x, const int *incx);
float FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const cuComplex *x, const int *incx);
float FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const cuComplex *x, const int *incx);
#endif
int FC_FUNC_(cublas_isamax,CUBLAS_ISAMAX) (const int *n, const float *x, const int *incx);
int FC_FUNC_(cublas_isamin,CUBLAS_ISAMIN) (const int *n, const float *x, const int *incx);
void FC_FUNC_(cublas_saxpy,CUBLAS_SAXPY) (const int *n, const float *alpha, const float *x, 
                   const int *incx, float *y, const int *incy);
void FC_FUNC_(cublas_scopy,CUBLAS_SCOPY) (const int *n, const float *x, const int *incx, float *y, 
                   const int *incy);
void FC_FUNC_(cublas_srot,CUBLAS_SROT) (const int *n, float *x, const int *incx, float *y, 
                  const int *incy, const float *sc, const float *ss);
void FC_FUNC_(cublas_srotg,CUBLAS_SROTG) (float *sa, float *sb, float *sc, float *ss);
void FC_FUNC_(cublas_srotm,CUBLAS_SROTM) (const int *n, float *x, const int *incx, float *y,
                   const int *incy, const float* sparam);
void FC_FUNC_(cublas_srotmg,CUBLAS_SROTMG) (float *sd1, float *sd2, float *sx1, const float *sy1, 
                    float* sparam);
void FC_FUNC_(cublas_sscal,CUBLAS_SSCAL) (const int *n, const float *alpha, float *x,
                   const int *incx);
void FC_FUNC_(cublas_sswap,CUBLAS_SSWAP)(const int *n, float *x, const int *incx, float *y,
                  const int *incy);

void FC_FUNC_(cublas_caxpy,CUBLAS_CAXPY) (const int *n, const cuComplex *alpha, const cuComplex *x, 
                   const int *incx, cuComplex *y, const int *incy);
void FC_FUNC_(cublas_ccopy,CUBLAS_CCOPY) (const int *n, const cuComplex *x, const int *incx, 
                   cuComplex *y, const int *incy);
void FC_FUNC_(cublas_crot,CUBLAS_CROT) (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                  const int *incy, const float *sc, const cuComplex *cs);
void FC_FUNC_(cublas_crotg,CUBLAS_CROTG) (cuComplex *ca, const cuComplex *cb, float *sc, 
                   cuComplex *cs);
void FC_FUNC_(cublas_cscal,CUBLAS_CSCAL) (const int *n, const cuComplex *alpha, cuComplex *x, 
                   const int *incx);
void FC_FUNC_(cublas_csrot,CUBLAS_CSROT) (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                   const int *incy, const float *sc, const float *ss);
void FC_FUNC_(cublas_csscal,CUBLAS_CSSCAL) (const int *n, const float *alpha, cuComplex *x,
                    const int *incx);
void FC_FUNC_(cublas_cswap,CUBLAS_CSWAP) (const int *n, cuComplex *x, const int *incx, cuComplex *y,
                   const int *incy);
void FC_FUNC_(cublas_cdotu,CUBLAS_CDOTU) (cuComplex *retVal, const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy);
void FC_FUNC_(cublas_cdotc,CUBLAS_CDOTC) (cuComplex *retVal,const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy);
int FC_FUNC_(cublas_icamax,CUBLAS_ICAMAX) (const int *n, const cuComplex *x, const int *incx);
int FC_FUNC_(cublas_icamin,CUBLAS_ICAMIN) (const int *n, const cuComplex *x, const int *incx);

/* BLAS2 */
void FC_FUNC_(cublas_sgbmv,CUBLAS_SGBMV) (const char *trans, const int *m, const int *n, 
                   const int *kl, const int *ku, const float *alpha, 
                   const float *A, const int *lda, const float *x,
                   const int *incx, const float *beta, float *y, 
                   const int *incy);
void FC_FUNC_(cublas_sgemv,CUBLAS_SGEMV) (const char *trans, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta, 
                   float *y, const int *incy);
void FC_FUNC_(cublas_sger,CUBLAS_SGER) (const int *m, const int *n, const float *alpha,
                  const float *x, const int *incx, const float *y,
                  const int *incy, float *A, const int *lda);
void FC_FUNC_(cublas_ssbmv,CUBLAS_SSBMV) (const char *uplo, const int *n, const int *k, 
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta,
                   float *y, const int *incy);
void FC_FUNC_(cublas_sspmv,CUBLAS_SSPMV) (const char *uplo, const int *n, const float *alpha, 
                   const float *AP, const float *x, const int *incx, 
                   const float *beta, float *y, const int *incy);
void FC_FUNC_(cublas_sspr,CUBLAS_SSPR) (const char *uplo, const int *n, const float *alpha,
                  const float *x, const int *incx, float *AP);
void FC_FUNC_(cublas_sspr2,CUBLAS_SSPR2) (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y,
                   const int *incy, float *AP);
void FC_FUNC_(cublas_ssymv,CUBLAS_SSYMV) (const char *uplo, const int *n, const float *alpha,
                   const float *A, const int *lda, const float *x,
                   const int *incx, const float *beta, float *y, 
                   const int *incy);
void FC_FUNC_(cublas_ssyr,CUBLAS_SSYR) (const char *uplo, const int *n, const float *alpha,
                  const float *x, const int *incx, float *A, const int *lda);
void FC_FUNC_(cublas_ssyr2,CUBLAS_SSYR2) (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y,
                   const int *incy, float *A, const int *lda);
void FC_FUNC_(cublas_stbmv,CUBLAS_STBMV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const float *A, const int *lda, 
                   float *x, const int *incx);
void FC_FUNC_(cublas_stbsv,CUBLAS_STBSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const float *A, const int *lda, 
                   float *x, const int *incx);
void FC_FUNC_(cublas_stpmv,CUBLAS_STPMV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *AP, float *x, const int *incx);
void FC_FUNC_(cublas_stpsv,CUBLAS_STPSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *AP, float *x, const int *incx);
void FC_FUNC_(cublas_strmv,CUBLAS_STRMV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *A, const int *lda, float *x, 
                   const int *incx);
void FC_FUNC_(cublas_strsv,CUBLAS_STRSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *A, const int *lda, float *x, 
                   const int *incx);

/* BLAS3 */
void FC_FUNC_(cublas_sgemm,CUBLAS_SGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha, 
                   const float *A, const int *lda, const float *B, 
                   const int *ldb, const float *beta, float *C, 
                   const int *ldc);
void FC_FUNC_(cublas_ssymm,CUBLAS_SSYMM) (const char *side, const char *uplo, const int *m,
                   const int *n, const float *alpha, const float *A,
                   const int *lda, const float *B, const int *ldb, 
                   const float *beta, float *C, const int *ldc);
void FC_FUNC_(cublas_ssyr2k,CUBLAS_SSYR2K) (const char *uplo, const char *trans, const int *n, 
                    const int *k, const float *alpha, const float *A,
                    const int *lda,const float *B, const int *ldb, 
                    const float *beta, float *C, const int *ldc);
void FC_FUNC_(cublas_ssyrk,CUBLAS_SSYRK) (const char *uplo, const char *trans, const int *n,
                   const int *k, const float *alpha, const float *A,
                   const int *lda, const float *beta, float *C,
                   const int *ldc);
void FC_FUNC_(cublas_strmm,CUBLAS_STRMM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb);
void FC_FUNC_(cublas_strsm,CUBLAS_STRSM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb);

void FC_FUNC_(cublas_cgemm,CUBLAS_CGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const cuComplex *A, const int *lda, const cuComplex *B,
                   const int *ldb, const cuComplex *beta, cuComplex *C, 
                   const int *ldc);

/* DP BLAS1 */
double FC_FUNC_(cublas_ddot,CUBLAS_DDOT) (const int *n, const double *x, const int *incx, double *y, 
                   const int *incy);
double FC_FUNC_(cublas_dasum,CUBLAS_DASUM) (const int *n, const double *x, const int *incx);
double FC_FUNC_(cublas_dnrm2,CUBLAS_DNRM2) (const int *n, const double *x, const int *incx);
int FC_FUNC_(cublas_idamax,CUBLAS_IDAMAX) (const int *n, const double *x, const int *incx);
int FC_FUNC_(cublas_idamin,CUBLAS_IDAMIN) (const int *n, const double *x, const int *incx);
void FC_FUNC_(cublas_daxpy,CUBLAS_DAXPY) (const int *n, const double *alpha, const double *x, 
                   const int *incx, double *y, const int *incy);
void FC_FUNC_(cublas_dcopy,CUBLAS_DCOPY) (const int *n, const double *x, const int *incx, double *y, 
                   const int *incy);
void FC_FUNC_(cublas_drot,CUBLAS_DROT) (const int *n, double *x, const int *incx, double *y, 
                  const int *incy, const double *sc, const double *ss);
void FC_FUNC_(cublas_drotg,CUBLAS_DROTG) (double *sa, double *sb, double *sc, double *ss);
void FC_FUNC_(cublas_drotm,CUBLAS_DROTM) (const int *n, double *x, const int *incx, double *y,
                   const int *incy, const double* sparam);
void FC_FUNC_(cublas_drotmg,CUBLAS_DROTMG) (double *sd1, double *sd2, double *sx1, const double *sy1, 
                    double* sparam);
void FC_FUNC_(cublas_dscal,CUBLAS_DSCAL) (const int *n, const double *alpha, double *x,
                   const int *incx);
void FC_FUNC_(cublas_dswap,CUBLAS_DSWAP)(const int *n, double *x, const int *incx, double *y,
                  const int *incy);

/* DP BLAS2 */
void FC_FUNC_(cublas_dgemv,CUBLAS_DGEMV) (const char *trans, const int *m, const int *n,
                   const double *alpha, const double *A, const int *lda,
                   const double *x, const int *incx, const double *beta, 
                   double *y, const int *incy);
void FC_FUNC_(cublas_dger,CUBLAS_DGER) (const int *m, const int *n, const double *alpha,
                  const double *x, const int *incx, const double *y,
                  const int *incy, double *A, const int *lda);
void FC_FUNC_(cublas_dsyr,CUBLAS_DSYR) (const char *uplo, const int *n, const double *alpha,
                  const double *x, const int *incx, double *A, const int *lda);
void FC_FUNC_(cublas_dtrsv,CUBLAS_DTRSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const double *A, const int *lda, double *x, 
                   const int *incx);

/* DP BLAS3 */
void FC_FUNC_(cublas_dgemm,CUBLAS_DGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const double *alpha, 
                   const double *A, const int *lda, const double *B, 
                   const int *ldb, const double *beta, double *C, 
                   const int *ldc);
void FC_FUNC_(cublas_dsymm,CUBLAS_DSYMM) (const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const double *A,
                   const int *lda, const double *B, const int *ldb, 
                   const double *beta, double *C, const int *ldc);
void FC_FUNC_(cublas_dsyr2k,CUBLAS_DSYR2K) (const char *uplo, const char *trans, const int *n, 
                    const int *k, const double *alpha, const double *A,
                    const int *lda,const double *B, const int *ldb, 
                    const double *beta, double *C, const int *ldc);
void FC_FUNC_(cublas_dsyrk,CUBLAS_DSYRK) (const char *uplo, const char *trans, const int *n,
                   const int *k, const double *alpha, const double *A,
                   const int *lda, const double *beta, double *C,
                   const int *ldc);
void FC_FUNC_(cublas_dtrmm,CUBLAS_DTRMM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const double *alpha, const double *A, const int *lda,
                   double *B, const int *ldb);
void FC_FUNC_(cublas_dtrsm,CUBLAS_DTRSM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const double *alpha, const double *A, const int *lda,
                   double *B, const int *ldb);

void FC_FUNC_(cublas_zgemm,CUBLAS_ZGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuDoubleComplex *alpha,
                   const cuDoubleComplex *A, const int *lda, 
                   const cuDoubleComplex *B, const int *ldb, 
                   const cuDoubleComplex *beta, cuDoubleComplex *C, 
                   const int *ldc);

#if defined(__cplusplus)
}
#endif /* __cplusplus */

#define CUBLAS_WRAPPER_ERROR_NOERR      0
#define CUBLAS_WRAPPER_ERROR_ALLOC      1
#define CUBLAS_WRAPPER_ERROR_SET        2
#define CUBLAS_WRAPPER_ERROR_GET        3
#define CUBLAS_WRAPPER_ERROR_STUB       4

static char *errMsg[5] = 
{
    "no error",
    "allocation error",
    "setVector/setMatrix error",
    "getVector/getMatrix error",
    "not implemented"
};

static void wrapperError (const char *funcName, int error)
{
    printf ("cublas%s wrapper: %s\n", funcName, errMsg[error]);
    fflush (stdout);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS1 ----------------------------------*/
/*---------------------------------------------------------------------------*/

int FC_FUNC_(cublas_isamax,CUBLAS_ISAMAX) (const int *n, const float *x, const int *incx)
{
    void *devPtrx = 0;
    int retVal = 0;
    cublasStatus stat;

    if (*n <= 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx), sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Isamax", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Isamax", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasIsamax (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

int FC_FUNC_(cublas_isamin,CUBLAS_ISAMIN) (const int *n, const float *x, const int *incx)
{
    void *devPtrx = 0;
    int retVal = 0;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx), sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Isamin", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Isamin", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasIsamin (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const float *x, const int *incx)
#else
float FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const float *x, const int *incx)
#endif
{
    void *devPtrx = 0;
    float retVal = 0.0f;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx), sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sasum", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sasum", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
   
    retVal = cublasSasum (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

void FC_FUNC_(cublas_saxpy,CUBLAS_SAXPY) (const int *n, const float *alpha, const float *x, 
                   const int *incx, float *y, const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Saxpy", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Saxpy", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasSaxpy (*n, *alpha, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Saxpy", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_scopy,CUBLAS_SCOPY) (const int *n, const float *x, const int *incx, float *y,
                   const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Scopy", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Scopy", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasScopy (*n, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Scopy", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const float *x, const int *incx, float *y,
                    const int *incy)
#else
float FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const float *x, const int *incx, float *y,
                   const int *incy)
#endif
{
    void *devPtrx = 0, *devPtry = 0;
    float retVal = 0.0f;
    cublasStatus stat1, stat2;

    if (*n == 0) return retVal;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sdot", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return retVal;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sdot", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return retVal;
    }
retVal = cublasSdot (*n, devPtrx, *incx, devPtry, *incy);
    cublasFree (devPtrx);
    cublasFree (devPtry);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const float *x, const int *incx)
#else
float FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const float *x, const int *incx)
#endif
{
    void *devPtrx = 0;
    float retVal = 0.0f;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Snrm2", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Snrm2", CUBLAS_WRAPPER_ERROR_SET);
        return retVal;
    }
    retVal = cublasSnrm2 (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

void FC_FUNC_(cublas_srot,CUBLAS_SROT) (const int *n, float *x, const int *incx, float *y, 
                  const int *incy, const float *sc, const float *ss)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Srot", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Srot", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasSrot (*n, devPtrx, *incx, devPtry, *incy, *sc, *ss);
    stat1 = cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Srot", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_srotg,CUBLAS_SROTG) (float *sa, float *sb, float *sc, float *ss)
{
    cublasSrotg (sa, sb, sc, ss);
}

void FC_FUNC_(cublas_srotm,CUBLAS_SROTM) (const int *n, float *x, const int *incx, float *y, 
                   const int *incy, const float* sparam)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Srotm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Srotm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasSrotm (*n, devPtrx, *incx, devPtry, *incy, sparam);
    stat1 = cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Srotm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_srotmg,CUBLAS_SROTMG) (float *sd1, float *sd2, float *sx1, const float *sy1,
                    float* sparam)
{
    cublasSrotmg (sd1, sd2, sx1, sy1, sparam);
}

void FC_FUNC_(cublas_sscal,CUBLAS_SSCAL) (const int *n, const float *alpha, float *x, const int *incx)
{
    void *devPtrx = 0;
    cublasStatus stat;

    if (*n == 0) return;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sscal", CUBLAS_WRAPPER_ERROR_ALLOC);
        return;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sscal", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return;
    }
    cublasSscal (*n, *alpha, devPtrx, *incx);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, *incx, x, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sscal", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx); 
}

void FC_FUNC_(cublas_sswap,CUBLAS_SSWAP) (const int *n, float *x, const int *incx, float *y, 
                   const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sswap", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sswap", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasSswap (*n, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sswap", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_caxpy,CUBLAS_CAXPY) (const int *n, const cuComplex *alpha, const cuComplex *x, 
                   const int *incx, cuComplex *y, const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Caxpy", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Caxpy", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasCaxpy (*n, *alpha, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Caxpy", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_ccopy,CUBLAS_CCOPY) (const int *n, const cuComplex *x, const int *incx, 
                   cuComplex *y, const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ccopy", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ccopy", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasCcopy (*n, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ccopy", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_crot,CUBLAS_CROT) (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                  const int *incy, const float *sc, const cuComplex *cs)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Crot", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Crot", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasCrot (*n, devPtrx, *incx, devPtry, *incy, *sc, *cs);
    stat1 = cublasGetVector (*n, sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Crot", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_crotg,CUBLAS_CROTG) (cuComplex *ca, const cuComplex *cb, float *sc,
                   cuComplex *cs)
{
    cublasCrotg (ca, *cb, sc, cs);
}

void FC_FUNC_(cublas_cscal,CUBLAS_CSCAL) (const int *n, const cuComplex *alpha, cuComplex *x, 
                   const int *incx)
{
    void *devPtrx = 0;
    cublasStatus stat;
    
    if (*n == 0) return;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Cscal", CUBLAS_WRAPPER_ERROR_ALLOC);
        return;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Cscal", CUBLAS_WRAPPER_ERROR_SET);
        return;
    }
    cublasCscal (*n, *alpha, devPtrx, *incx);
    stat = cublasGetVector (*n, sizeof(x[0]), devPtrx, *incx, x, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Cscal", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx); 
}

void FC_FUNC_(cublas_csrot,CUBLAS_CSROT) (const int *n, cuComplex *x, const int *incx, cuComplex *y, 
                   const int *incy, const float *sc, const float *ss)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Csrot", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Csrot", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasCsrot (*n, devPtrx, *incx, devPtry, *incy, *sc, *ss);
    stat1 = cublasGetVector (*n, sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Csrot", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_csscal,CUBLAS_CSSCAL) (const int *n, const float *alpha, cuComplex *x, 
                    const int *incx)
{
    void *devPtrx = 0;
    cublasStatus stat;

    if (*n == 0) return;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Csscal", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        return;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Csscal", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return;
    }
    cublasCsscal (*n, *alpha, devPtrx, *incx);
    cublasGetVector (*n, sizeof(x[0]), devPtrx, *incx, x, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Csscal", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx); 
}

void FC_FUNC_(cublas_cswap,CUBLAS_CSWAP) (const int *n, cuComplex *x, const int *incx, cuComplex *y,
                   const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (imax(1,*n *abs(*incy)),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cswap", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cswap", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasCswap (*n, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n, sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cswap", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_cdotu,CUBLAS_CDOTU) (cuComplex *retVal, const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    *retVal = make_cuComplex (0.0f, 0.0f);
    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cdotu", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cdotu", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    *retVal = cublasCdotu (*n, devPtrx, *incx, devPtry, *incy);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_cdotc,CUBLAS_CDOTC) (cuComplex *retVal, const int *n, const cuComplex *x, 
                   const int *incx, const cuComplex *y, const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    *retVal = make_cuComplex (0.0f, 0.0f);
    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cdotc", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cdotc", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    *retVal = cublasCdotc (*n, devPtrx, *incx, devPtry, *incy);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

int FC_FUNC_(cublas_icamax,CUBLAS_ICAMAX) (const int *n, const cuComplex *x, const int *incx)
{
    void *devPtrx = 0;
    int retVal = 0;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Icamax", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Icamax", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasIcamax (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

int FC_FUNC_(cublas_icamin,CUBLAS_ICAMIN) (const int *n, const cuComplex *x, const int *incx)
{
    void *devPtrx = 0;
    int retVal = 0;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Icamin", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Icamin", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasIcamin (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const cuComplex *x, const int *incx)
#else
float FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const cuComplex *x, const int *incx)
#endif
{
    void *devPtrx = 0;
    float retVal = 0.0f;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Scasum", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Scasum", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasScasum (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const cuComplex *x, const int *incx)
#else
float FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const cuComplex *x, const int *incx)
#endif
{
    void *devPtrx = 0;
    float retVal = 0.0f;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Scnrm2", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Scnrm2", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasScnrm2 (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

int FC_FUNC_(cublas_idamax,CUBLAS_IDAMAX) (const int *n, const double *x, const int *incx)
{
    void *devPtrx = 0;
    int retVal = 0;
    cublasStatus stat;

    if (*n == 0) return retVal;;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Idamax", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Idamax", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasIdamax (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

int FC_FUNC_(cublas_idamin,CUBLAS_IDAMIN) (const int *n, const double *x, const int *incx)
{
    void *devPtrx = 0;
    int retVal = 0;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Idamin", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Idamin", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasIdamin (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

double FC_FUNC_(cublas_dasum,CUBLAS_DASUM) (const int *n, const double *x, const int *incx)
{
    void *devPtrx = 0;
    double retVal = 0;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dasum", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dasum", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        return retVal;
    }
    retVal = cublasDasum (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

void FC_FUNC_(cublas_daxpy,CUBLAS_DAXPY) (const int *n, const double *alpha, const double *x, 
                   const int *incx, double *y, const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Daxpy", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Daxpy", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasDaxpy (*n, *alpha, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Daxpy", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_dcopy,CUBLAS_DCOPY) (const int *n, const double *x, const int *incx, double *y,
                   const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dcopy", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dcopy", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasDcopy (*n, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dcopy", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

double FC_FUNC_(cublas_ddot,CUBLAS_DDOT) (const int *n, const double *x, const int *incx, double *y,
                    const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    double retVal = 0.0;
    cublasStatus stat1, stat2;

    if (*n == 0) return retVal;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ddot", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return retVal;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ddot", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return retVal;
    }
    retVal = cublasDdot (*n, devPtrx, *incx, devPtry, *incy);
    cublasFree (devPtrx);
    cublasFree (devPtry);
    return retVal;
}

double FC_FUNC_(cublas_dnrm2,CUBLAS_DNRM2) (const int *n, const double *x, const int *incx)
{
    void *devPtrx = 0;
    double retVal = 0.0;
    cublasStatus stat;

    if (*n == 0) return retVal;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dnrm2", CUBLAS_WRAPPER_ERROR_ALLOC);
        return retVal;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dnrm2", CUBLAS_WRAPPER_ERROR_SET);
        return retVal;
    }
    retVal = cublasDnrm2 (*n, devPtrx, *incx);
    cublasFree (devPtrx);
    return retVal;
}

void FC_FUNC_(cublas_drot,CUBLAS_DROT) (const int *n, double *x, const int *incx, double *y, 
                  const int *incy, const double *sc, const double *ss)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Drot", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat1 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Drot", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasDrot (*n, devPtrx, *incx, devPtry, *incy, *sc, *ss);
    stat1 = cublasGetVector (*n, sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Drot", CUBLAS_WRAPPER_ERROR_GET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_drotg,CUBLAS_DROTG) (double *sa, double *sb, double *sc, double *ss)
{
    cublasDrotg (sa, sb, sc, ss);
}

void FC_FUNC_(cublas_drotm,CUBLAS_DROTM) (const int *n, double *x, const int *incx, double *y, 
                   const int *incy, const double* sparam)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Drotm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Drotm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasDrotm (*n, devPtrx, *incx, devPtry, *incy, sparam);
    stat1 = cublasGetVector (*n, sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Drotm", CUBLAS_WRAPPER_ERROR_GET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_drotmg,CUBLAS_DROTMG) (double *sd1, double *sd2, double *sx1, const double *sy1,
                    double* sparam)
{
    cublasDrotmg (sd1, sd2, sx1, sy1, sparam);
}

void FC_FUNC_(cublas_dscal,CUBLAS_DSCAL) (const int *n, const double *alpha, double *x, 
                   const int *incx)
{
    void *devPtrx = 0;
    cublasStatus stat;

    if (*n == 0) return;
    stat = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dscal", CUBLAS_WRAPPER_ERROR_ALLOC);
        return;
    }
    stat = cublasSetVector (*n, sizeof(x[0]), x, *incx, devPtrx, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dscal", CUBLAS_WRAPPER_ERROR_SET);
        return;
    }
    cublasDscal (*n, *alpha, devPtrx, *incx);
    stat = cublasGetVector (*n, sizeof(x[0]), devPtrx, *incx, x, *incx);
    if (stat != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dscal", CUBLAS_WRAPPER_ERROR_GET);
        return;
    }
    cublasFree (devPtrx); 
}

void FC_FUNC_(cublas_dswap,CUBLAS_DSWAP) (const int *n, double *x, const int *incx, double *y, 
                   const int *incy)
{
    void *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    stat1 = cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2 = cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dswap", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    stat1 = cublasSetVector (*n, sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n, sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dswap", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        return;
    }
    cublasDswap (*n, devPtrx, *incx, devPtry, *incy);
    stat1 = cublasGetVector (*n, sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    stat2 = cublasGetVector (*n, sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dswap", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS2 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void FC_FUNC_(cublas_sgbmv,CUBLAS_SGBMV) (const char *trans, const int *m, const int *n, 
                   const int *kl, const int *ku, const float *alpha, 
                   const float *A, const int *lda, const float *x, 
                   const int *incx, const float *beta, float *y, 
                   const int *incy)
{
    void *devPtrx = 0, *devPtry = 0, *devPtrA = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return;

    /*  X      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
     *  Y      - REAL             array of DIMENSION at least
     *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
     *  A      - REAL             array of DIMENSION ( LDA, n ). 
     * Before entry, the leading ( kl + ku + 1 ) by n part of the
     * array A must contain the matrix of coefficients
     */
    if (toupper(trans[0]) == 'N') {
        stat1 = cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
        stat2 = cublasAlloc(1+(*m-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    } else {
        stat1 = cublasAlloc(1+(*m-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
        stat2 = cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    }
    stat3 = cublasAlloc ((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sgbmv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    if (toupper(trans[0]) == 'N') {
        stat1=cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
        stat2=cublasSetVector(*m,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    } else {
        stat1=cublasSetVector(*m,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
        stat2=cublasSetVector(*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    }
    stat3 = cublasSetMatrix (imin(*kl+*ku+1,*lda), *n, sizeof(A[0]), A, *lda, 
                             devPtrA, *lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sgbmv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasSgbmv (trans[0], *m, *n, *kl, *ku, *alpha, devPtrA, *lda, devPtrx, 
                 *incx, *beta, devPtry, *incy);
    if (toupper(trans[0]) == 'N') {
        stat1=cublasGetVector(*m,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    } else {
        stat1=cublasGetVector(*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    }
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sgbmv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_sgemv,CUBLAS_SGEMV) (const char *trans, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta,
                   float *y, const int *incy)
{
    void *devPtrA = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return;
    
    /*  X      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
     *  Y      - REAL             array of DIMENSION at least
     *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry, the leading m by n part of the array A must
     *           contain the matrix of coefficients.
     */
    if (toupper(trans[0]) == 'N') {
        stat1=cublasAlloc (1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
        stat2=cublasAlloc (1+(*m-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    } else {
        stat1=cublasAlloc (1+(*m-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
        stat2=cublasAlloc (1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    }
    stat3=cublasAlloc ((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sgemv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    if (toupper(trans[0]) == 'N') {
        stat1=cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
        stat2=cublasSetVector(*m,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    } else {
        stat1=cublasSetVector(*m,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
        stat2=cublasSetVector(*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    }       
    stat3=cublasSetMatrix (imin(*m,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sgemv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasSgemv (trans[0], *m, *n, *alpha, devPtrA, *lda, devPtrx, *incx,
                 *beta, devPtry, *incy);
    if (toupper(trans[0]) == 'N') {
        stat1=cublasGetVector(*m,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    } else {
        stat1=cublasGetVector(*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    }       
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sgemv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_sger,CUBLAS_SGER) (const int *m, const int *n, const float *alpha, 
                  const float *x, const int *incx, const float *y,
                  const int *incy, float *A, const int *lda)
{
    void *devPtrA = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ).
     *  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *  A      - REAL array of DIMENSION ( LDA, n ).
     *           Before entry, the leading m by n part of the array A must
     *           contain the matrix of coefficients. On exit, A is
     */
    stat1=cublasAlloc(1+(*m-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    stat3=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sger", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector(*m,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetVector(*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    stat3=cublasSetMatrix(imin(*m,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sger", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasSger (*m, *n, *alpha, devPtrx, *incx, devPtry, *incy, devPtrA, *lda);
    stat1 = cublasGetMatrix(imin(*m,*lda),*n,sizeof(A[0]),devPtrA,*lda,A,*lda);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sger", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_ssbmv,CUBLAS_SSBMV) (const char *uplo, const int *n, const int *k, 
                   const float *alpha, const float *A, const int *lda,
                   const float *x, const int *incx, const float *beta, 
                   float *y, const int *incy)
{
    void *devPtrA = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if (*n == 0) return;
    
    /*  X      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *
     *  Y      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    stat3=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssbmv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    stat1 = cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector(*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    stat3 = cublasSetMatrix (imin(*k+1,*lda), *n, sizeof(A[0]), A, *lda, 
                             devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssbmv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasSsbmv (uplo[0], *n, *k, *alpha, devPtrA, *lda, devPtrx, *incx, *beta,
                 devPtry, *incy);
    stat1 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ssbmv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_sspmv,CUBLAS_SSPMV) (const char *uplo, const int *n, const float *alpha,
                   const float *AP, const float *x, const int *incx, 
                   const float *beta, float *y, const int *incy)
{
    void *devPtrAP = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular part of the symmetric matrix
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    stat3=cublasAlloc(((*n)*(*n+1))/2, sizeof(AP[0]), (void**)&devPtrAP);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sspmv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrAP);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    stat3 = cublasSetVector (((*n)*(*n+1))/2,sizeof(AP[0]),AP,1,devPtrAP,1);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sspmv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrAP);
        return;
    }
    cublasSspmv (*uplo, *n, *alpha, devPtrAP, devPtrx, *incx, *beta, devPtry,
                 *incy);
    stat1 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sspmv", CUBLAS_WRAPPER_ERROR_GET); 
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrAP);
}

void FC_FUNC_(cublas_sspr,CUBLAS_SSPR) (const char *uplo, const int *n, const float *alpha, 
                  const float *x, const int *incx, float *AP)
{
    void *devPtrAP = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(((*n)*(*n+1))/2, sizeof(AP[0]), (void**)&devPtrAP);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sspr", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrAP);
        return;
    }
    stat1=cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetVector (((*n)*(*n+1))/2,sizeof(AP[0]),AP,1,devPtrAP,1);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sspr", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrAP);
        return;
    }
    cublasSspr (uplo[0], *n, *alpha, devPtrx, *incx, devPtrAP);
    stat1=cublasGetVector (((*n)*(*n+1))/2,sizeof(AP[0]),devPtrAP,1,AP,1);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sspr", CUBLAS_WRAPPER_ERROR_GET); 
    }
    cublasFree (devPtrx);
    cublasFree (devPtrAP);
}

void FC_FUNC_(cublas_sspr2,CUBLAS_SSPR2) (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y, 
                   const int *incy, float *AP)
{
    void *devPtrAP = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     *           Before entry with  UPLO = 'U' or 'u', the array AP must
     *           contain the upper triangular part of the symmetric matrix
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    stat3=cublasAlloc(((*n)*(*n+1))/2, sizeof(AP[0]), (void**)&devPtrAP);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sspr2", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrAP);
        return;
    }
    stat1 = cublasSetVector (1+(*n-1)*abs(*incx),sizeof(x[0]),x,1,devPtrx,1);
    stat2 = cublasSetVector (1+(*n-1)*abs(*incy),sizeof(y[0]),y,1,devPtry,1);
    stat3 = cublasSetVector (((*n)*(*n+1))/2,sizeof(AP[0]),AP,1,devPtrAP,1);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sspr2", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrAP);
        return;
    }
    cublasSspr2 (uplo[0], *n, *alpha, devPtrx, *incx, devPtry, *incy,devPtrAP);
    stat1 = cublasGetVector (((*n)*(*n+1))/2,sizeof(AP[0]),devPtrAP,1,AP,1);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sspr2", CUBLAS_WRAPPER_ERROR_GET); 
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrAP);
}

void FC_FUNC_(cublas_ssymv,CUBLAS_SSYMV) (const char *uplo, const int *n, const float *alpha,
                   const float *A, const int *lda, const float *x, 
                   const int *incx, const float *beta, float *y, 
                   const int *incy)
{
    void *devPtrA = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    stat3=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssymv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    stat3 = cublasSetMatrix (imin(*n,*lda), *n, sizeof(A[0]), A, *lda, devPtrA,
                             *lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssymv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasSsymv (uplo[0], *n, *alpha, devPtrA, *lda, devPtrx, *incx, *beta,
                 devPtry, *incy);
    stat1 = cublasGetVector (*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ssymv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_ssyr,CUBLAS_SSYR) (const char *uplo, const int *n, const float *alpha, 
                  const float *x, const int *incx, float *A, const int *lda)
{
    void *devPtrA = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc((*lda)*(*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyr", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    stat1 = cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetMatrix(imin(*n,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyr", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    cublasSsyr (uplo[0], *n, *alpha, devPtrx, *incx, devPtrA, *lda);
    stat1 = cublasGetMatrix(imin(*n,*lda),*n,sizeof(A[0]),devPtrA,*lda,A,*lda);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ssyr", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_ssyr2,CUBLAS_SSYR2) (const char *uplo, const int *n, const float *alpha,
                   const float *x, const int *incx, const float *y,
                   const int *incy, float *A, const int *lda)
{
    void *devPtrA = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    stat3=cublasAlloc((*lda)*(*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyr2", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetVector (*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    stat3=cublasSetMatrix (imin(*n,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyr2", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasSsyr2 (uplo[0], *n, *alpha, devPtrx, *incx, devPtry, *incy, devPtrA,
                 *lda);
    stat1=cublasGetMatrix (imin(*n,*lda),*n,sizeof(A[0]),devPtrA,*lda,A,*lda);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ssyr2", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_stbmv,CUBLAS_STBMV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const float *A, const int *lda,
                   float *x, const int *incx)
{
    void *devPtrA = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    
    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
     *           by n part of the array A must contain the lower triangular
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stbmv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetMatrix(imin(*k+1,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stbmv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    cublasStbmv (uplo[0],trans[0],diag[0],*n,*k,devPtrA,*lda,devPtrx,*incx);
    stat1 = cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Stbmv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_stbsv,CUBLAS_STBSV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const float *A, const int *lda,
                   float *x, const int *incx)
{
    void *devPtrA = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
     *           by n part of the array A must contain the upper triangular
     *           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
     *           by n part of the array A must contain the lower triangular
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stbsv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetMatrix(imin(*k+1,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stbsv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    cublasStbsv (uplo[0],trans[0],diag[0],*n,*k,devPtrA,*lda,devPtrx,*incx);
    stat1 = cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Stbsv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_stpmv,CUBLAS_STPMV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const float *AP, float *x, const int *incx)
{
    void *devPtrAP = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(((*n)*(*n+1))/2, sizeof(AP[0]), (void**)&devPtrAP);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stpmv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrAP);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (((*n)*(*n+1))/2,sizeof(AP[0]),AP,1,devPtrAP,1);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stpmv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrAP);
        return;
    }
    cublasStpmv (uplo[0], trans[0], diag[0], *n, devPtrAP, devPtrx, *incx);
    stat1 = cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Stpmv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrAP);
}

void FC_FUNC_(cublas_stpsv,CUBLAS_STPSV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const float *AP, float *x, const int *incx)
{
    void *devPtrAP = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  AP     - REAL             array of DIMENSION at least
     *           ( ( n*( n + 1 ) )/2 ).
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(((*n)*(*n+1))/2, sizeof(AP[0]), (void**)&devPtrAP);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stpsv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrAP);
        return;
    }
    stat1 = cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetVector (((*n)*(*n+1))/2,sizeof(AP[0]),AP,1,devPtrAP,1);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Stpsv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrAP);
        return;
    }
    cublasStpsv (uplo[0], trans[0], diag[0], *n, devPtrAP, devPtrx, *incx);
    stat1 = cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Stpsv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrAP);
}

void FC_FUNC_(cublas_strmv,CUBLAS_STRMV) (const char *uplo, const char *trans,
                            const char *diag, const int *n, const float *A,
                            const int *lda, float *x, const int *incx)
{
    void *devPtrA = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;
    
    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strmv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetMatrix (imin(*n,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strmv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    cublasStrmv (uplo[0], trans[0], diag[0], *n, devPtrA, *lda, devPtrx,*incx);
    stat1=cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Strmv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrx);
}

void FC_FUNC_(cublas_strsv,CUBLAS_STRSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const float *A, const int *lda, float *x, 
                   const int *incx)
{
    void *devPtrA = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strsv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetMatrix (imin(*n,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strsv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    cublasStrsv (uplo[0], trans[0], diag[0], *n, devPtrA, *lda, devPtrx,*incx);
    stat1=cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Strsv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_dgemv,CUBLAS_DGEMV) (const char *trans, const int *m, const int *n,
                   const double *alpha, const double *A, const int *lda,
                   const double *x, const int *incx, const double *beta,
                   double *y, const int *incy)
{
    void *devPtrA = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;
    
    if ((*m == 0) || (*n == 0)) return;

    /*  X      - REAL             array of DIMENSION at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
     *  Y      - REAL             array of DIMENSION at least
     *           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
     *           and at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry, the leading m by n part of the array A must
     *           contain the matrix of coefficients.
     */
    if (toupper(trans[0]) == 'N') {
        stat1 = cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
        stat2 = cublasAlloc(1+(*m-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    } else {
        stat1 = cublasAlloc(1+(*m-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
        stat2 = cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    }
    stat3 = cublasAlloc ((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dgemv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    if (toupper(trans[0]) == 'N') {
        stat1=cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
        stat2=cublasSetVector(*m,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    } else {
        stat1=cublasSetVector(*m,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
        stat2=cublasSetVector(*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    }       
    stat3=cublasSetMatrix (imin(*m,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    stat3=cublasSetMatrix (imin(*m,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dgemv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasDgemv (trans[0], *m, *n, *alpha, devPtrA, *lda, devPtrx, *incx,
                 *beta, devPtry, *incy);
    if (toupper(trans[0]) == 'N') {
        stat1=cublasGetVector(*m,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    } else {
        stat1=cublasGetVector(*n,sizeof(y[0]),devPtry,abs(*incy),y,abs(*incy));
    }       
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dgemv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrx);
    cublasFree (devPtry);
}

void FC_FUNC_(cublas_dger,CUBLAS_DGER) (const int *m, const int *n, const double *alpha, 
                  const double *x, const int *incx, const double *y,
                  const int *incy, double *A, const int *lda)
{
    void *devPtrA = 0, *devPtrx = 0, *devPtry = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( m - 1 )*abs( INCX ) ).
     *
     *  Y      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCY ) ).
     *
     * A       - REAL array of DIMENSION ( LDA, n ).
     *           Before entry, the leading m by n part of the array A must
     *           contain the matrix of coefficients. On exit, A is
     */
    stat1=cublasAlloc(1+(*m-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc(1+(*n-1)*abs(*incy),sizeof(y[0]),(void**)&devPtry);
    stat3=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dger", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector(*m,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetVector(*n,sizeof(y[0]),y,abs(*incy),devPtry,abs(*incy));
    stat3=cublasSetMatrix(imin(*m,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dger", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtry);
        cublasFree (devPtrA);
        return;
    }
    cublasDger (*m, *n, *alpha, devPtrx, *incx, devPtry, *incy, devPtrA, *lda);
    stat1 = cublasGetMatrix(imin(*m,*lda),*n,sizeof(A[0]),devPtrA,*lda,A,*lda);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dger", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtry);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_dsyr,CUBLAS_DSYR) (const char *uplo, const int *n, const double *alpha, 
                  const double *x, const int *incx, double *A, const int *lda)
{
    void *devPtrA = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsyr", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    stat1 = cublasSetVector(*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2 = cublasSetMatrix(imin(*n,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsyr", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    cublasDsyr (uplo[0], *n, *alpha, devPtrx, *incx, devPtrA, *lda);
    stat1 = cublasGetMatrix(imin(*n,*lda),*n,sizeof(A[0]),devPtrA,*lda,A,*lda);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dsyr", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_dtrsv,CUBLAS_DTRSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const double *A, const int *lda, double *x, 
                   const int *incx)
{
    void *devPtrA = 0, *devPtrx = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /*  X      - REAL             array of dimension at least
     *           ( 1 + ( n - 1 )*abs( INCX ) ).
     *  A      - REAL             array of DIMENSION ( LDA, n ).
     *           Before entry with  UPLO = 'U' or 'u', the leading n by n
     *           upper triangular part of the array A must contain the upper
     *           Before entry with UPLO = 'L' or 'l', the leading n by n
     *           lower triangular part of the array A must contain the lower
     */
    stat1=cublasAlloc(1+(*n-1)*abs(*incx),sizeof(x[0]),(void**)&devPtrx);
    stat2=cublasAlloc((*lda) * (*n), sizeof(A[0]), (void**)&devPtrA);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dtrsv", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    stat1=cublasSetVector (*n,sizeof(x[0]),x,abs(*incx),devPtrx,abs(*incx));
    stat2=cublasSetMatrix (imin(*n,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) ||
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dtrsv", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrx);
        cublasFree (devPtrA);
        return;
    }
    cublasDtrsv (uplo[0], trans[0], diag[0], *n, devPtrA, *lda, devPtrx,*incx);
    stat1=cublasGetVector (*n,sizeof(x[0]),devPtrx,abs(*incx),x,abs(*incx));
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dtrsv", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrx);
    cublasFree (devPtrA);
}

void FC_FUNC_(cublas_dspr2,CUBLAS_DSPR2) (void)
{
    wrapperError ("Dspr2", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dsyr2,CUBLAS_DSYR2) (void)
{
    wrapperError ("Dsyr2", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dspr,CUBLAS_DSPR) (void)
{
    wrapperError ("Dspr", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dtpsv,CUBLAS_DTPSV) (void)
{
    wrapperError ("Dtpsv", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dtbsv,CUBLAS_DTBSV) (void)
{
    wrapperError ("Dtbsv", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dtpmv,CUBLAS_DTPMV) (void)
{
    wrapperError ("Dtpmv", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dtbmv,CUBLAS_DTBMV) (void)
{
    wrapperError ("Dtbmv", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dtrmv,CUBLAS_DTRMV) (void)
{
    wrapperError ("Dtrmv", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dspmv,CUBLAS_DSPMV) (void)
{
    wrapperError ("Dspmv", CUBLAS_WRAPPER_ERROR_STUB); 
}

void FC_FUNC_(cublas_dsbmv,CUBLAS_DSBMV) (void)
{
    wrapperError ("Dsbmv", CUBLAS_WRAPPER_ERROR_STUB); 
} 

void FC_FUNC_(cublas_dsymv,CUBLAS_DSYMV) (void)
{
    wrapperError ("Dsymv", CUBLAS_WRAPPER_ERROR_STUB); 
}

void FC_FUNC_(cublas_dgbmv,CUBLAS_DGBMV) (void)
{
    wrapperError ("Dgbmv", CUBLAS_WRAPPER_ERROR_STUB);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS3 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void FC_FUNC_(cublas_sgemm,CUBLAS_SGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha,
                   const float *A, const int *lda, const float *B,
                   const int *ldb, const float *beta, float *C, const int *ldc)
{
    int ka, kb;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return; 

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     *  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B.
     *  C      - REAL             array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     */
    ka = (toupper(transa[0]) == 'N') ? *k : *m;
    kb = (toupper(transb[0]) == 'N') ? *n : *k;
    stat1 = cublasAlloc (imax(1,*lda*ka),sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc (imax(1,*ldb*kb),sizeof(B[0]), (void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n),  sizeof(C[0]), (void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sgemm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(transa[0]) == 'N') {
        stat1=cublasSetMatrix(imin(*m,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    } else {
        stat1=cublasSetMatrix(imin(*k,*lda),*m,sizeof(A[0]),A,*lda,devPtrA,*lda);
    }
    if (toupper(transb[0]) == 'N') {
        stat2=cublasSetMatrix(imin(*k,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    } else {
        stat2=cublasSetMatrix(imin(*n,*ldb),*k,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    }
    stat3=cublasSetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Sgemm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasSgemm (transa[0], transb[0], *m, *n, *k, *alpha, devPtrA, *lda, 
                 devPtrB, *ldb, *beta, devPtrC, *ldc);
    stat1=cublasGetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Sgemm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_ssymm,CUBLAS_SSYMM) (const char *side, const char *uplo, const int *m, 
                   const int *n, const float *alpha, const float *A, 
                   const int *lda, const float *B, const int *ldb, 
                   const float *beta, float *C, const int *ldc)
{
    int ka;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return;
    
    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
     *           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
     *           the array  A  must contain the  symmetric matrix, [..]
     *           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
     *           the array  A  must contain the  symmetric matrix, [..]
     *  B      - REAL             array of DIMENSION ( LDB, n ).
     *           Before entry, the leading  m by n part of the array  B  must
     *           contain the matrix B.
     *  C      - REAL             array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     */
    ka = (toupper(side[0]) == 'L') ? *m : *n;
    stat1 = cublasAlloc ((*lda) * ka, sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc ((*ldb) * (*n), sizeof(B[0]), (void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n), sizeof(C[0]), (void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssymm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    stat1 = cublasSetMatrix(imin(ka,*lda),ka,sizeof(A[0]),A,*lda,devPtrA,*lda);
    stat2 = cublasSetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    stat3 = cublasSetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssymm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasSsymm (side[0], uplo[0], *m, *n, *alpha, devPtrA, *lda, devPtrB,
                 *ldb, *beta, devPtrC, *ldc);
    stat1 = cublasGetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ssymm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_ssyr2k,CUBLAS_SSYR2K) (const char *uplo, const char *trans, const int *n,
                    const int *k, const float *alpha, const float *A, 
                    const int *lda, const float *B, const int *ldb, 
                    const float *beta, float *C, const int *ldc)
{
    int ka, kb;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;

    if (*n == 0) return;

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the leading  n by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by n  part of the array  A  must contain  the
     *           matrix A.
     *  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the leading  n by k
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  k by n  part of the array  B  must contain  the
     *           matrix B.
     * C       - single precision array of dimensions (ldc, n). If uplo == 'U' 
     *           or 'u', the leading n x n triangular part of the array C must 
     */
    ka = (toupper(trans[0]) == 'N') ? *k : *n;
    kb = (toupper(trans[0]) == 'N') ? *k : *n;
    stat1 = cublasAlloc (imax(1,*lda*ka),sizeof(A[0]),(void**)&devPtrA);
    stat2 = cublasAlloc (imax(1,*ldb*kb),sizeof(B[0]),(void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n),  sizeof(C[0]),(void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyr2k", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(trans[0]) == 'N') {
      stat1=cublasSetMatrix(imin(*n,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
      stat2=cublasSetMatrix(imin(*n,*ldb),*k,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    } else {
      stat1=cublasSetMatrix(imin(*k,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
      stat2=cublasSetMatrix(imin(*k,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    }
    stat3 = cublasSetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyr2k", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasSsyr2k (uplo[0], trans[0], *n, *k, *alpha, devPtrA, *lda, devPtrB, 
                  *ldb, *beta, devPtrC, *ldc);
    stat1 = cublasGetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ssyr2k", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_ssyrk,CUBLAS_SSYRK) (const char *uplo, const char *trans, const int *n, 
                   const int *k, const float *alpha, const float *A, 
                   const int *lda, const float *beta, float *C, const int *ldc)
{
    int ka;
    void *devPtrA = 0, *devPtrC = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /* A      single precision array of dimensions (lda, ka), where ka is k 
     *        when trans == 'N' or 'n', and is n otherwise. When trans == 'N' 
     *        or 'n', the leading n x k part of array A must contain the matrix
     *        A, otherwise the leading k x n part of the array must contain the
     *        matrix A.
     * C      single precision array of dimensions (ldc, n). If uplo='U'or'u',
     *        the leading n x n triangular part of the array C must contain the
     *        upper triangular part of the symmetric matrix C and the strictly 
     */
    ka = (toupper(trans[0]) == 'N') ? *k : *n;
    stat1 = cublasAlloc (imax(1,*lda*ka),sizeof(A[0]),(void**)&devPtrA);
    stat2 = cublasAlloc ((*ldc)*(*n), sizeof(C[0]), (void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyrk", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(trans[0]) == 'N') {
      stat1=cublasSetMatrix(imin(*n,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    } else {
      stat1=cublasSetMatrix(imin(*k,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    }
    stat2 = cublasSetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Ssyrk", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrC);
        return;
    }
    cublasSsyrk (uplo[0], trans[0], *n, *k, *alpha, devPtrA, *lda, *beta,
                 devPtrC, *ldc);
    stat1 = cublasGetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Ssyrk", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_strmm,CUBLAS_STRMM) (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb)
{
    int k;
    void *devPtrA = 0, *devPtrB = 0;
    cublasStatus stat1, stat2;
    
    if ((*m == 0) || (*n == 0)) return;

    /* A      single precision array of dimensions (lda, k). k = m if side =
     *        'L' or 'l', k = n if side = 'R' or 'r'. If uplo = 'U' or 'u'
     *        the leading k x k upper triangular part of the array A must
     *        contain the upper triangular matrix, and the strictly lower
     *        triangular part of A is not referenced. If uplo = 'L' or 'l'
     *        the leading k x k lower triangular part of the array A must
     *        contain the lower triangular matrix, and the strictly upper
     * B      single precision array of dimensions (ldb, n). On entry, the 
     *        leading m x n part of the array contains the matrix B. It is
     *        overwritten with the transformed matrix on exit.
     */
    k = (toupper(side[0]) == 'L') ? *m : *n;
    stat1 = cublasAlloc (*lda * k, sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc (*ldb * (*n), sizeof(B[0]), (void**)&devPtrB);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strmm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    stat1 = cublasSetMatrix(imin(k,*lda),k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    stat2 = cublasSetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strmm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    cublasStrmm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha, devPtrA,
                 *lda, devPtrB, *ldb);
    stat1 = cublasGetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),devPtrB,*ldb,B,*ldb);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Strmm", CUBLAS_WRAPPER_ERROR_GET);
    }        
    cublasFree (devPtrA);
    cublasFree (devPtrB);
}

void FC_FUNC_(cublas_strsm,CUBLAS_STRSM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n, 
                   const float *alpha, const float *A, const int *lda,
                   float *B, const int *ldb)
{
    void *devPtrA = 0, *devPtrB = 0;
    int k;
    cublasStatus stat1, stat2;

    if ((*m == 0) || (*n == 0)) return;

    /*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
     *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
     *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
     *           upper triangular part of the array  A must contain the upper
     *  B      - REAL             array of DIMENSION ( LDB, n ).
     *           Before entry,  the leading  m by n part of the array  B must
     *           contain  the  right-hand  side  matrix  B,  and  on exit  is
     */
    k = (toupper(side[0]) == 'L') ? *m : *n;
    stat1 = cublasAlloc (*lda * k, sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc (*ldb * (*n), sizeof(B[0]), (void**)&devPtrB);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strsm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    stat1 = cublasSetMatrix(imin(k,*lda),k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    stat2 = cublasSetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Strsm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    cublasStrsm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha, devPtrA,
                 *lda, devPtrB, *ldb);
    stat1 = cublasGetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),devPtrB,*ldb,B,*ldb);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Strsm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
}

void FC_FUNC_(cublas_cgemm,CUBLAS_CGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const cuComplex *A, const int *lda, const cuComplex *B,
                   const int *ldb, const cuComplex *beta, cuComplex *C, 
                   const int *ldc)
{
    int ka, kb;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return; 

    /*  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     *  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B.
     *  C      - COMPLEX          array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     */
    ka = (toupper(transa[0]) == 'N') ? *k : *m;
    kb = (toupper(transb[0]) == 'N') ? *n : *k;
    stat1 = cublasAlloc (imax(1,*lda*ka),sizeof(A[0]),(void**)&devPtrA);
    stat2 = cublasAlloc (imax(1,*ldb*kb),sizeof(B[0]),(void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n),  sizeof(C[0]),(void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cgemm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(transa[0]) == 'N') {
      stat1=cublasSetMatrix(imin(*m,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    } else {
      stat1=cublasSetMatrix(imin(*k,*lda),*m,sizeof(A[0]),A,*lda,devPtrA,*lda);
    }
    if (toupper(transb[0]) == 'N') {
      stat2=cublasSetMatrix(imin(*k,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    } else {
      stat2=cublasSetMatrix(imin(*n,*ldb),*k,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    }
    stat3=cublasSetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Cgemm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasCgemm (transa[0], transb[0], *m, *n, *k, *alpha, devPtrA, *lda, 
                 devPtrB, *ldb, *beta, devPtrC, *ldc);
    stat1=cublasGetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Cgemm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_chemm,CUBLAS_CHEMM) (void)
{
    wrapperError ("Chemm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_csymm,CUBLAS_CSYMM) (void)
{
    wrapperError ("Csymm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_ctrmm,CUBLAS_CTRMM) (void)
{
    wrapperError ("Ctrmm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_ctrsm,CUBLAS_CTRSM) (void)
{
    wrapperError ("Ctrsm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_cherk,CUBLAS_CHERK) (void)
{
    wrapperError ("CHerk", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_csyrk,CUBLAS_CSYRK) (void)
{
    wrapperError ("CSyrk", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_cher2k,CUBLAS_CHER2K) (void)
{
    wrapperError ("CHer2k", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_csyr2k,CUBLAS_CSYR2K) (void)
{
    wrapperError ("CSyr2k", CUBLAS_WRAPPER_ERROR_STUB);
}

void FC_FUNC_(cublas_dgemm,CUBLAS_DGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const double *alpha,
                   const double *A, const int *lda, const double *B,
                   const int *ldb, const double *beta, double *C,
                   const int *ldc)
{
    int ka, kb;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return;

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     *  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B. 
     *  C      - REAL             array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     */
    ka = (toupper(transa[0]) == 'N') ? *k : *m;
    kb = (toupper(transb[0]) == 'N') ? *n : *k;
    stat1 = cublasAlloc (imax(1,*lda*ka),sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc (imax(1,*ldb*kb),sizeof(B[0]), (void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n),  sizeof(C[0]), (void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dgemm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(transa[0]) == 'N') {
      stat1=cublasSetMatrix(imin(*m,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    } else {
      stat1=cublasSetMatrix(imin(*k,*lda),*m,sizeof(A[0]),A,*lda,devPtrA,*lda);
    }
    if (toupper(transb[0]) == 'N') {
      stat2=cublasSetMatrix(imin(*k,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    } else {
      stat2=cublasSetMatrix(imin(*n,*ldb),*k,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    }
    stat3=cublasSetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dgemm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasDgemm (transa[0], transb[0], *m, *n, *k, *alpha, devPtrA, *lda, 
                 devPtrB, *ldb, *beta, devPtrC, *ldc);
    stat1=cublasGetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dgemm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_dsymm,CUBLAS_DSYMM) (const char *side, const char *uplo, const int *m, 
                   const int *n, const double *alpha, const double *A, 
                   const int *lda, const double *B, const int *ldb, 
                   const double *beta, double *C, const int *ldc)
{
    int ka;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;
    
    if ((*m == 0) || (*n == 0)) return;

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
     *           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
     *           the array  A  must contain the  symmetric matrix, [..]
     *           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
     *           the array  A  must contain the  symmetric matrix, [..]
     *  B      - REAL             array of DIMENSION ( LDB, n ).
     *           Before entry, the leading  m by n part of the array  B  must
     *           contain the matrix B.
     *  C      - REAL             array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     */
    ka = (toupper(side[0]) == 'L') ? *m : *n;
    stat1 = cublasAlloc ((*lda) * ka, sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc ((*ldb) * (*n), sizeof(B[0]), (void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n), sizeof(C[0]), (void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsymm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    stat1 = cublasSetMatrix(imin(ka,*lda),ka,sizeof(A[0]),A,*lda,devPtrA,*lda);
    stat2 = cublasSetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    stat3 = cublasSetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsymm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasDsymm (side[0], uplo[0], *m, *n, *alpha, devPtrA, *lda, devPtrB,
                 *ldb, *beta, devPtrC, *ldc);
    stat1 = cublasGetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dsymm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_dsyr2k,CUBLAS_DSYR2K) (const char *uplo, const char *trans, const int *n,
                    const int *k, const double *alpha, const double *A, 
                    const int *lda, const double *B, const int *ldb, 
                    const double *beta, double *C, const int *ldc)
{
    int ka, kb;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;

    if (*n == 0) return;

    /*  A      - REAL             array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the leading  n by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by n  part of the array  A  must contain  the
     *           matrix A.
     *  B      - REAL             array of DIMENSION ( LDB, kb ), where kb is
     *           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
     *           Before entry with  TRANS = 'N' or 'n',  the leading  n by k
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  k by n  part of the array  B  must contain  the
     *           matrix B.
     *  C      - single precision array of dimensions (ldc, n). If uplo == 'U' 
     *           or  'u', the leading n x n triangular part of the array C must
     */
    ka = (toupper(trans[0]) == 'N') ? *k : *n;
    kb = (toupper(trans[0]) == 'N') ? *k : *n;
    stat1 = cublasAlloc (imax(1,*lda*ka),sizeof(A[0]),(void**)&devPtrA);
    stat2 = cublasAlloc (imax(1,*ldb*kb),sizeof(B[0]),(void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n),  sizeof(C[0]),(void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsyr2k", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(trans[0]) == 'N') {
      stat1=cublasSetMatrix(imin(*n,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
      stat2=cublasSetMatrix(imin(*n,*ldb),*k,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    } else {
      stat1=cublasSetMatrix(imin(*k,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
      stat2=cublasSetMatrix(imin(*k,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    }
    stat3=cublasSetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsyr2k", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasDsyr2k (uplo[0], trans[0], *n, *k, *alpha, devPtrA, *lda, devPtrB, 
                  *ldb, *beta, devPtrC, *ldc);
    stat1=cublasGetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dsyr2k", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_dsyrk,CUBLAS_DSYRK) (const char *uplo, const char *trans, const int *n, 
                   const int *k, const double *alpha, const double *A, 
                   const int *lda, const double *beta, double *C, 
                   const int *ldc)
{
    int ka;
    void *devPtrA = 0, *devPtrC = 0;
    cublasStatus stat1, stat2;

    if (*n == 0) return;

    /* A      single precision array of dimensions (lda, ka), where ka is k 
     *        when trans == 'N' or 'n', and is n otherwise. When trans == 'N' 
     *        or 'n', the leading n x k part of array A must contain the matrix
     *        A, otherwise the leading k x n part of the array must contain the
     *        matrix A.
     * C      single precision array of dimensions (ldc, n). If uplo='U'or'u',
     *        the leading n x n triangular part of the array C must contain the
     *        upper triangular part of the symmetric matrix C and the strictly 
     */
    ka = (toupper(trans[0]) == 'N') ? *k : *n;
    stat1 = cublasAlloc(imax(1,*lda*ka),sizeof(A[0]),(void**)&devPtrA);
    stat2 = cublasAlloc(imax(1,*ldc*(*n)),sizeof(C[0]),(void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsyrk", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(trans[0]) == 'N') {
      stat1=cublasSetMatrix(imin(*n,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    } else {
      stat1=cublasSetMatrix(imin(*k,*lda),*n,sizeof(A[0]),A,*lda,devPtrA,*lda);
    }
    stat2 = cublasSetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dsyrk", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrC);
        return;
    }
    cublasDsyrk (uplo[0], trans[0], *n, *k, *alpha, devPtrA, *lda, *beta,
                 devPtrC, *ldc);
    stat1 = cublasGetMatrix(imin(*n,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dsyrk", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_dtrmm,CUBLAS_DTRMM) (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n,
                   const double *alpha, const double *A, const int *lda,
                   double *B, const int *ldb)
{
    int k;
    void *devPtrA = 0, *devPtrB = 0;
    cublasStatus stat1, stat2;

    if ((*m == 0) || (*n == 0)) return;

    /* A      single precision array of dimensions (lda, k). k = m if side =
     *        'L' or 'l', k = n if side = 'R' or 'r'. If uplo = 'U' or 'u'
     *        the leading k x k upper triangular part of the array A must
     *        contain the upper triangular matrix, and the strictly lower
     *        triangular part of A is not referenced. If uplo = 'L' or 'l'
     *        the leading k x k lower triangular part of the array A must
     *        contain the lower triangular matrix, and the strictly upper
     * B      single precision array of dimensions (ldb, n). On entry, the 
     *        leading m x n part of the array contains the matrix B. It is
     *        overwritten with the transformed matrix on exit.
     */
    k = (toupper(side[0]) == 'L') ? *m : *n;
    stat1 = cublasAlloc (*lda * k, sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc (*ldb * (*n), sizeof(B[0]), (void**)&devPtrB);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dtrmm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    stat1 = cublasSetMatrix(imin(k,*lda),k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    stat2 = cublasSetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dtrmm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    cublasDtrmm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha, devPtrA,
                 *lda, devPtrB, *ldb);
    stat1 = cublasGetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),devPtrB,*ldb,B,*ldb);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dtrmm", CUBLAS_WRAPPER_ERROR_GET);
    }        
    cublasFree (devPtrA);
    cublasFree (devPtrB);
}

void FC_FUNC_(cublas_dtrsm,CUBLAS_DTRSM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n, 
                   const double *alpha, const double *A, const int *lda,
                   double *B, const int *ldb)
{
    void *devPtrA = 0, *devPtrB = 0;
    int k;
    cublasStatus stat1, stat2;

    if ((*m == 0) || (*n == 0)) return;

    /*  A      - REAL             array of DIMENSION ( LDA, k ), where k is m
     *           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
     *           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
     *           upper triangular part of the array  A must contain the upper
     *  B      - REAL             array of DIMENSION ( LDB, n ).
     *           Before entry,  the leading  m by n part of the array  B must
     *           contain  the  right-hand  side  matrix  B,  and  on exit  is
     */
    k = (toupper(side[0]) == 'L') ? *m : *n;
    stat1 = cublasAlloc (*lda * k, sizeof(A[0]), (void**)&devPtrA);
    stat2 = cublasAlloc (*ldb * (*n), sizeof(B[0]), (void**)&devPtrB);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dtrsm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    stat1 = cublasSetMatrix(imin(k,*lda),k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    stat2 = cublasSetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Dtrsm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        return;
    }
    cublasDtrsm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha, devPtrA,
                 *lda, devPtrB, *ldb);
    stat1=cublasGetMatrix(imin(*m,*ldb),*n,sizeof(B[0]),devPtrB,*ldb,B,*ldb);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Dtrsm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
}

void FC_FUNC_(cublas_zgemm,CUBLAS_ZGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuDoubleComplex *alpha,
                   const cuDoubleComplex *A, const int *lda, 
                   const cuDoubleComplex *B, const int *ldb, 
                   const cuDoubleComplex *beta, cuDoubleComplex *C, 
                   const int *ldc)
{
    int ka, kb;
    void *devPtrA = 0, *devPtrB = 0, *devPtrC = 0;
    cublasStatus stat1, stat2, stat3;

    if ((*m == 0) || (*n == 0)) return; 
    
    /*  A      - COMPLEX          array of DIMENSION ( LDA, ka ), where ka is
     *           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
     *           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
     *           part of the array  A  must contain the matrix  A,  otherwise
     *           the leading  k by m  part of the array  A  must contain  the
     *           matrix A.
     *  B      - COMPLEX          array of DIMENSION ( LDB, kb ), where kb is
     *           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
     *           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
     *           part of the array  B  must contain the matrix  B,  otherwise
     *           the leading  n by k  part of the array  B  must contain  the
     *           matrix B.
     *  C      - COMPLEX          array of DIMENSION ( LDC, n ).
     *           Before entry, the leading  m by n  part of the array  C must
     *           contain the matrix  C,  except when  beta  is zero, in which
     *           case C need not be set on entry.
     *           On exit, the array  C  is overwritten by the  m by n  matrix
     */
    ka = (toupper(transa[0]) == 'N') ? *k : *m;
    kb = (toupper(transb[0]) == 'N') ? *n : *k;
    stat1 = cublasAlloc (imax(1,*lda*ka),sizeof(A[0]),(void**)&devPtrA);
    stat2 = cublasAlloc (imax(1,*ldb*kb),sizeof(B[0]),(void**)&devPtrB);
    stat3 = cublasAlloc ((*ldc) * (*n),  sizeof(C[0]),(void**)&devPtrC);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Zgemm", CUBLAS_WRAPPER_ERROR_ALLOC);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    if (toupper(transa[0]) == 'N') {
      stat1=cublasSetMatrix(imin(*m,*lda),*k,sizeof(A[0]),A,*lda,devPtrA,*lda);
    } else {
      stat1=cublasSetMatrix(imin(*k,*lda),*m,sizeof(A[0]),A,*lda,devPtrA,*lda);
    }
    if (toupper(transb[0]) == 'N') {
      stat2=cublasSetMatrix(imin(*k,*ldb),*n,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    } else {
      stat2=cublasSetMatrix(imin(*n,*ldb),*k,sizeof(B[0]),B,*ldb,devPtrB,*ldb);
    }
    stat3=cublasSetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),C,*ldc,devPtrC,*ldc);
    if ((stat1 != CUBLAS_STATUS_SUCCESS) || 
        (stat2 != CUBLAS_STATUS_SUCCESS) ||
        (stat3 != CUBLAS_STATUS_SUCCESS)) {
        wrapperError ("Zgemm", CUBLAS_WRAPPER_ERROR_SET);
        cublasFree (devPtrA);
        cublasFree (devPtrB);
        cublasFree (devPtrC);
        return;
    }
    cublasZgemm (transa[0], transb[0], *m, *n, *k, *alpha, devPtrA, *lda, 
                 devPtrB, *ldb, *beta, devPtrC, *ldc);
    stat1=cublasGetMatrix(imin(*m,*ldc),*n,sizeof(C[0]),devPtrC,*ldc,C,*ldc);
    if (stat1 != CUBLAS_STATUS_SUCCESS) {
        wrapperError ("Zgemm", CUBLAS_WRAPPER_ERROR_GET);
    }
    cublasFree (devPtrA);
    cublasFree (devPtrB);
    cublasFree (devPtrC);
}

void FC_FUNC_(cublas_zhemm,CUBLAS_ZHEMM) (void)
{
    wrapperError ("ZHemm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_zsymm,CUBLAS_ZSYMM) (void)
{
    wrapperError ("ZSymm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_ztrmm,CUBLAS_ZTRMM) (void)
{
    wrapperError ("ZTrmm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_ztrsm,CUBLAS_ZTRSM) (void)
{
    wrapperError ("ZTrsm", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_zherk,CUBLAS_ZHERK) (void)
{
    wrapperError ("ZHerk", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_zsyrk,CUBLAS_ZSYRK) (void)
{
    wrapperError ("ZSyrk", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_zher2k,CUBLAS_ZHER2K) (void)
{
    wrapperError ("ZHer2k", CUBLAS_WRAPPER_ERROR_STUB);
}
void FC_FUNC_(cublas_zsyr2k,CUBLAS_ZSYR2K) (void)
{
    wrapperError ("ZSyr2k", CUBLAS_WRAPPER_ERROR_STUB);
}

#else /* defined(CUBLAS_USE_THUNKING) */

/*
 * Fortran callable thin wrappers. Fortran application must allocate and
 * deallocate GPU memory, and copy data up and down.
 */
#if defined(__cplusplus)
extern "C" {
#endif /* __cplusplus */
#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                    const devptr_t *devPtry, const int *incy);
double FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const devptr_t *devPtrx, const int *incx);
double FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const devptr_t *devPtrx, const int *incx);
double FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const devptr_t *devPtrx, const int *incx);
double FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const devptr_t *devPtrx, const int *incx);
#else
float FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
float FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const devptr_t *devPtrx, const int *incx);
float FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const devptr_t *devPtrx, const int *incx);
float FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const devptr_t *devPtrx, const int *incx);
float FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const devptr_t *devPtrx, const int *incx);
#endif

int FC_FUNC_(cublas_isamax,CUBLAS_ISAMAX) (const int *n, const devptr_t *devPtrx, const int *incx);
int FC_FUNC_(cublas_isamin,CUBLAS_ISAMIN) (const int *n, const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_saxpy,CUBLAS_SAXPY) (const int *n, const float *alpha, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_scopy,CUBLAS_SCOPY) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_srot,CUBLAS_SROT) (const int *n, const devptr_t *devPtrX, const int *incx, 
                  const devptr_t *devPtrY, const int *incy, const float *sc, 
                  const float *ss);
void FC_FUNC_(cublas_srotg,CUBLAS_SROTG) (float *sa, float *sb, float *sc, float *ss);
void FC_FUNC_(cublas_srotm,CUBLAS_SROTM) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const float* sparam);
void FC_FUNC_(cublas_srotmg,CUBLAS_SROTMG) (float *sd1, float *sd2, float *sx1, const float *sy1, 
                    float* sparam);
void FC_FUNC_(cublas_sscal,CUBLAS_SSCAL) (const int *n, const float *alpha, const devptr_t *devPtrx,
                   const int *incx);
void FC_FUNC_(cublas_sswap,CUBLAS_SSWAP) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);

void FC_FUNC_(cublas_caxpy,CUBLAS_CAXPY) (const int *n, const cuComplex *alpha,
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_ccopy,CUBLAS_CCOPY) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_crot,CUBLAS_CROT) (const int *n, const devptr_t *devPtrX, const int *incx, 
                  const devptr_t *devPtrY, const int *incy, const float *sc, 
                  const cuComplex *cs);
void FC_FUNC_(cublas_crotg,CUBLAS_CROTG) (cuComplex *ca, const cuComplex *cb, float *sc,
                   cuComplex *cs);
void FC_FUNC_(cublas_cscal,CUBLAS_CSCAL) (const int *n, const cuComplex *alpha, 
                   const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_csrot,CUBLAS_CSROT) (const int *n, const devptr_t *devPtrX, const int *incx, 
                   const devptr_t *devPtrY, const int *incy, const float *sc, 
                   const float *ss);
void FC_FUNC_(cublas_csscal,CUBLAS_CSSCAL) (const int *n, const float *alpha, const devptr_t *devPtrx, 
                    const int *incx);
void FC_FUNC_(cublas_cswap,CUBLAS_CSWAP) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_cdotu,CUBLAS_CDOTU) (cuComplex *retVal, const int *n, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_cdotc,CUBLAS_CDOTC) (cuComplex *retVal, const int *n, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy);
int FC_FUNC_(cublas_icamax,CUBLAS_ICAMAX) (const int *n, const devptr_t *devPtrx, const int *incx);
int FC_FUNC_(cublas_icamin,CUBLAS_ICAMIN) (const int *n, const devptr_t *devPtrx, const int *incx);
double FC_FUNC_(cublas_ddot,CUBLAS_DDOT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
double FC_FUNC_(cublas_dasum,CUBLAS_DASUM) (const int *n, const devptr_t *devPtrx, const int *incx);
double FC_FUNC_(cublas_dnrm2,CUBLAS_DNRM2) (const int *n, const devptr_t *devPtrx, const int *incx);
int FC_FUNC_(cublas_idamax,CUBLAS_IDAMAX) (const int *n, const devptr_t *devPtrx, const int *incx);
int FC_FUNC_(cublas_idamin,CUBLAS_IDAMIN) (const int *n, const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_daxpy,CUBLAS_DAXPY) (const int *n, const double *alpha, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_dcopy,CUBLAS_DCOPY) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_drot,CUBLAS_DROT) (const int *n, const devptr_t *devPtrX, const int *incx, 
                  const devptr_t *devPtrY, const int *incy, const double *sc, 
                  const double *ss);
void FC_FUNC_(cublas_drotg,CUBLAS_DROTG) (double *sa, double *sb, double *sc, double *ss);
void FC_FUNC_(cublas_drotm,CUBLAS_DROTM) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const double* sparam);
void FC_FUNC_(cublas_drotmg,CUBLAS_DROTMG) (double *sd1, double *sd2, double *sx1, const double *sy1, 
                    double* sparam);
void FC_FUNC_(cublas_dscal,CUBLAS_DSCAL) (const int *n, const double *alpha, const devptr_t *devPtrx,
                   const int *incx);
void FC_FUNC_(cublas_dswap,CUBLAS_DSWAP) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy);


/* BLAS2 */
void FC_FUNC_(cublas_sgbmv,CUBLAS_SGBMV) (const char *trans, const int *m, const int *n,
                   const int *kl, const int *ku, const float *alpha, 
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_sgemv,CUBLAS_SGEMV) (const char *trans, const int *m, const int *n,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_sger,CUBLAS_SGER) (const int *m, const int *n, const float *alpha, 
                  const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, 
                  const devptr_t *devPtrA, const int *lda);
void FC_FUNC_(cublas_ssbmv,CUBLAS_SSBMV) (const char *uplo, const int *n, const int *k, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_sspmv,CUBLAS_SSPMV) (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrAP, const devptr_t *devPtrx, 
                   const int *incx, const float *beta, const devptr_t *devPtry,
                   const int *incy);
void FC_FUNC_(cublas_sspr,CUBLAS_SSPR) (const char *uplo, const int *n, const float *alpha,
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtrAP);
void FC_FUNC_(cublas_sspr2,CUBLAS_SSPR2) (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const devptr_t *devPtrAP);
void FC_FUNC_(cublas_ssymv,CUBLAS_SSYMV) (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy);
void FC_FUNC_(cublas_ssyr,CUBLAS_SSYR) (const char *uplo, const int *n, const float *alpha,
                  const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtrA, const int *lda);
void FC_FUNC_(cublas_ssyr2,CUBLAS_SSYR2) (const char *uplo, const int *n, const float *alpha, 
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const devptr_t *devPtrA, const int *lda);
void FC_FUNC_(cublas_stbmv,CUBLAS_STBMV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_stbsv,CUBLAS_STBSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const int *k, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_stpmv,CUBLAS_STPMV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_stpsv,CUBLAS_STPSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_strmv,CUBLAS_STRMV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_strsv,CUBLAS_STRSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx);
void FC_FUNC_(cublas_dgemv,CUBLAS_DGEMV) (const char *trans, const int *m, const int *n,
                   const double *alpha, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrx, const int *incx, 
                   const double *beta, const devptr_t *devPtry, 
                   const int *incy);
void FC_FUNC_(cublas_dger,CUBLAS_DGER) (const int *m, const int *n, const double *alpha, 
                  const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, 
                  const devptr_t *devPtrA, const int *lda);
void FC_FUNC_(cublas_dsyr,CUBLAS_DSYR) (const char *uplo, const int *n, const double *alpha,
                  const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtrA, const int *lda);
void FC_FUNC_(cublas_dtrsv,CUBLAS_DTRSV) (const char *uplo, const char *trans, const char *diag, 
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx);

/* BLAS 3 */
void FC_FUNC_(cublas_sgemm,CUBLAS_SGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha, 
                   const devptr_t *A, const int *lda, const devptr_t *B, 
                   const int *ldb, const float *beta, const devptr_t *C, 
                   const int *ldc);
void FC_FUNC_(cublas_ssymm,CUBLAS_SSYMM) (const char *side, const char *uplo, const int *m,
                   const int *n, const float *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrB, const int *ldb,
                   const float *beta, const devptr_t *devPtrC, const int *ldc);
void FC_FUNC_(cublas_ssyr2k,CUBLAS_SSYR2K) (const char *uplo, const char *trans, const int *n,
                    const int *k, const float *alpha, const devptr_t *devPtrA,
                    const int *lda, const devptr_t *devPtrB, const int *ldb,
                    const float *beta, const devptr_t *devPtrC, const int *ldc);
void FC_FUNC_(cublas_ssyrk,CUBLAS_SSYRK) (const char *uplo, const char *trans, const int *n,
                   const int *k, const float *alpha, const devptr_t *devPtrA,
                   const int *lda, const float *beta, const devptr_t *devPtrC,
                   const int *ldc);
void FC_FUNC_(cublas_strmm,CUBLAS_STRMM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb);
void FC_FUNC_(cublas_strsm,CUBLAS_STRSM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrB, const int *ldb);

void FC_FUNC_(cublas_cgemm,CUBLAS_CGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb, 
                   const cuComplex *beta, const devptr_t *devPtrC,
                   const int *ldc);

void FC_FUNC_(cublas_dgemm,CUBLAS_DGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const double *alpha, 
                   const devptr_t *A, const int *lda, const devptr_t *B, 
                   const int *ldb, const double *beta, const devptr_t *C, 
                   const int *ldc);
void FC_FUNC_(cublas_dsymm,CUBLAS_DSYMM) (const char *side, const char *uplo, const int *m,
                   const int *n, const double *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrB, const int *ldb,
                   const double *beta, const devptr_t *devPtrC,
                   const int *ldc);
void FC_FUNC_(cublas_dsyr2k,CUBLAS_DSYR2K) (const char *uplo, const char *trans, const int *n,
                    const int *k, const double *alpha, const devptr_t *devPtrA,
                    const int *lda, const devptr_t *devPtrB, const int *ldb,
                    const double *beta, const devptr_t *devPtrC, 
                    const int *ldc);
void FC_FUNC_(cublas_dsyrk,CUBLAS_DSYRK) (const char *uplo, const char *trans, const int *n,
                   const int *k, const double *alpha, const devptr_t *devPtrA,
                   const int *lda, const double *beta, const devptr_t *devPtrC,
                   const int *ldc);
void FC_FUNC_(cublas_dtrmm,CUBLAS_DTRMM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n,
                   const double *alpha, const devptr_t *devPtrA, 
                   const int *lda,
                   const devptr_t *devPtrB, const int *ldb);
void FC_FUNC_(cublas_dtrsm,CUBLAS_DTRSM) (const char *side, const char *uplo, const char *transa, 
                   const char *diag, const int *m, const int *n, 
                   const double *alpha, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrB, const int *ldb);

void FC_FUNC_(cublas_zgemm,CUBLAS_ZGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuDoubleComplex *alpha,
                   const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb, 
                   const cuDoubleComplex *beta, const devptr_t *devPtrC,
                   const int *ldc);

#if defined(__cplusplus)
}
#endif /* __cplusplus */

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS1 ----------------------------------*/
/*---------------------------------------------------------------------------*/

int FC_FUNC_(cublas_isamax,CUBLAS_ISAMAX) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    int retVal;
    retVal = cublasIsamax (*n, x, *incx);
    return retVal;
}

int FC_FUNC_(cublas_isamin,CUBLAS_ISAMIN) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    int retVal;
    retVal = cublasIsamin (*n, x, *incx);
    return retVal;
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float FC_FUNC_(cublas_sasum,CUBLAS_SASUM) (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float retVal;
    retVal = cublasSasum (*n, x, *incx);
    return retVal;
}

void FC_FUNC_(cublas_saxpy,CUBLAS_SAXPY) (const int *n, const float *alpha, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSaxpy (*n, *alpha, x, *incx, y, *incy);
}

void FC_FUNC_(cublas_scopy,CUBLAS_SCOPY) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasScopy (*n, x, *incx, y, *incy);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                    const devptr_t *devPtry, const int *incy)
#else
float FC_FUNC_(cublas_sdot,CUBLAS_SDOT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
#endif
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    return cublasSdot (*n, x, *incx, y, *incy);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float FC_FUNC_(cublas_snrm2,CUBLAS_SNRM2) (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    return cublasSnrm2 (*n, x, *incx);
}

void FC_FUNC_(cublas_srot,CUBLAS_SROT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, const float *sc, 
                  const float *ss)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSrot (*n, x, *incx, y, *incy, *sc, *ss);
}

void FC_FUNC_(cublas_srotg,CUBLAS_SROTG) (float *sa, float *sb, float *sc, float *ss)
{
    cublasSrotg (sa, sb, sc, ss);
}

void FC_FUNC_(cublas_srotm,CUBLAS_SROTM) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const float* sparam) 
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSrotm (*n, x, *incx, y, *incy, sparam);
}

void FC_FUNC_(cublas_srotmg,CUBLAS_SROTMG) (float *sd1, float *sd2, float *sx1, const float *sy1,
                    float* sparam)
{
    cublasSrotmg (sd1, sd2, sx1, sy1, sparam);
}

void FC_FUNC_(cublas_sscal,CUBLAS_SSCAL) (const int *n, const float *alpha, const devptr_t *devPtrx,
                   const int *incx)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    cublasSscal (*n, *alpha, x, *incx);
}

void FC_FUNC_(cublas_sswap,CUBLAS_SSWAP) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSswap (*n, x, *incx, y, *incy);
}

void FC_FUNC_(cublas_caxpy,CUBLAS_CAXPY) (const int *n, const cuComplex *alpha, 
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCaxpy (*n, *alpha, x, *incx, y, *incy);
}

void FC_FUNC_(cublas_ccopy,CUBLAS_CCOPY) (const int *n, const devptr_t *devPtrx, const int *incx,
                   const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCcopy (*n, x, *incx, y, *incy);
}

void FC_FUNC_(cublas_crot,CUBLAS_CROT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, const float *sc, 
                  const cuComplex *cs)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCrot (*n, x, *incx, y, *incy, *sc, *cs);
}

void FC_FUNC_(cublas_crotg,CUBLAS_CROTG) (cuComplex *ca, const cuComplex *cb, float *sc,
                   cuComplex *cs)
{
    cublasCrotg (ca, *cb, sc, cs);
}

void FC_FUNC_(cublas_cscal,CUBLAS_CSCAL) (const int *n, const cuComplex *alpha, 
                   const devptr_t *devPtrx, const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cublasCscal (*n, *alpha, x, *incx);
}

void FC_FUNC_(cublas_csrot,CUBLAS_CSROT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, const float *sc, 
                   const float *ss)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCsrot (*n, x, *incx, y, *incy, *sc, *ss);
}

void FC_FUNC_(cublas_csscal,CUBLAS_CSSCAL) (const int *n, const float *alpha, const devptr_t *devPtrx,
                    const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cublasCsscal (*n, *alpha, x, *incx);
}

void FC_FUNC_(cublas_cswap,CUBLAS_CSWAP) (const int *n, const devptr_t *devPtrx, const int *incx,
                   const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    cublasCswap (*n, x, *incx, y, *incy);
}

void FC_FUNC_(cublas_cdotu,CUBLAS_CDOTU) (cuComplex *retVal, const int *n, const devptr_t *devPtrx,
                   const int *incx, const devptr_t *devPtry,const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    *retVal = cublasCdotu (*n, x, *incx, y, *incy);
}

void FC_FUNC_(cublas_cdotc,CUBLAS_CDOTC) (cuComplex *retVal, const int *n, const devptr_t *devPtrx,
                   const int *incx, const devptr_t *devPtry, const int *incy)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    cuComplex *y = (cuComplex *)(uintptr_t)(*devPtry);
    *retVal = cublasCdotc (*n, x, *incx, y, *incy);
}

int FC_FUNC_(cublas_icamax,CUBLAS_ICAMAX) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasIcamax (*n, x, *incx);
}

int FC_FUNC_(cublas_icamin,CUBLAS_ICAMIN) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasIcamin (*n, x, *incx);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float FC_FUNC_(cublas_scasum,CUBLAS_SCASUM) (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasScasum (*n, x, *incx);
}

#if CUBLAS_FORTRAN_COMPILER==CUBLAS_G77
double FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const devptr_t *devPtrx, const int *incx)
#else
float FC_FUNC_(cublas_scnrm2,CUBLAS_SCNRM2) (const int *n, const devptr_t *devPtrx, const int *incx)
#endif
{
    cuComplex *x = (cuComplex *)(uintptr_t)(*devPtrx);
    return cublasScnrm2 (*n, x, *incx);
}

int FC_FUNC_(cublas_idamax,CUBLAS_IDAMAX) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    int retVal;
    retVal = cublasIdamax (*n, x, *incx);
    return retVal;
}

int FC_FUNC_(cublas_idamin,CUBLAS_IDAMIN) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    int retVal;
    retVal = cublasIdamin (*n, x, *incx);
    return retVal;
}

double FC_FUNC_(cublas_dasum,CUBLAS_DASUM) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    double retVal;
    retVal = cublasDasum (*n, x, *incx);
    return retVal;
}

void FC_FUNC_(cublas_daxpy,CUBLAS_DAXPY) (const int *n, const double *alpha, const devptr_t *devPtrx, 
                   const int *incx, const devptr_t *devPtry, const int *incy)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);
    cublasDaxpy (*n, *alpha, x, *incx, y, *incy);
}

void FC_FUNC_(cublas_dcopy,CUBLAS_DCOPY) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);
    cublasDcopy (*n, x, *incx, y, *incy);
}

double FC_FUNC_(cublas_ddot,CUBLAS_DDOT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                    const devptr_t *devPtry, const int *incy)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);
    return cublasDdot (*n, x, *incx, y, *incy);
}

double FC_FUNC_(cublas_dnrm2,CUBLAS_DNRM2) (const int *n, const devptr_t *devPtrx, const int *incx)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    return cublasDnrm2 (*n, x, *incx);
}

void FC_FUNC_(cublas_drot,CUBLAS_DROT) (const int *n, const devptr_t *devPtrx, const int *incx, 
                  const devptr_t *devPtry, const int *incy, const double *sc, 
                  const double *ss)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);
    cublasDrot (*n, x, *incx, y, *incy, *sc, *ss);
}

void FC_FUNC_(cublas_drotg,CUBLAS_DROTG) (double *sa, double *sb, double *sc, double *ss)
{
    cublasDrotg (sa, sb, sc, ss);
}

void FC_FUNC_(cublas_drotm,CUBLAS_DROTM) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy, 
                   const double* sparam) 
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);
    cublasDrotm (*n, x, *incx, y, *incy, sparam);
}

void FC_FUNC_(cublas_drotmg,CUBLAS_DROTMG) (double *sd1, double *sd2, double *sx1, const double *sy1,
                    double* sparam)
{
    cublasDrotmg (sd1, sd2, sx1, sy1, sparam);
}

void FC_FUNC_(cublas_dscal,CUBLAS_DSCAL) (const int *n, const double *alpha, const devptr_t *devPtrx,
                   const int *incx)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    cublasDscal (*n, *alpha, x, *incx);
}

void FC_FUNC_(cublas_dswap,CUBLAS_DSWAP) (const int *n, const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy)
{
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);
    cublasDswap (*n, x, *incx, y, *incy);
}


/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS2 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void FC_FUNC_(cublas_sgbmv,CUBLAS_SGBMV) (const char *trans, const int *m, const int *n, 
                   const int *kl, const int *ku, const float *alpha,
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSgbmv (trans[0], *m, *n, *kl, *ku, *alpha, A, *lda, x, *incx, *beta,
                 y, *incy);
}

void FC_FUNC_(cublas_sgemv,CUBLAS_SGEMV) (const char *trans, const int *m, const int *n, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);
    cublasSgemv (trans[0], *m, *n, *alpha, A, *lda, x, *incx, *beta, y, *incy);
}

void FC_FUNC_(cublas_sger,CUBLAS_SGER) (const int *m, const int *n, const float *alpha, 
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtry, const int *incy,
                  const devptr_t *devPtrA, const int *lda)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSger (*m, *n, *alpha, x, *incx, y, *incy, A, *lda);
}

void FC_FUNC_(cublas_ssbmv,CUBLAS_SSBMV) (const char *uplo, const int *n, const int *k,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry, const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSsbmv (uplo[0], *n, *k, *alpha, A, *lda, x, *incx, *beta, y, *incy);
}

void FC_FUNC_(cublas_sspmv,CUBLAS_SSPMV) (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrAP, const devptr_t *devPtrx,
                   const int *incx, const float *beta, const devptr_t *devPtry,
                   const int *incy)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSspmv (uplo[0], *n, *alpha, AP, x, *incx, *beta, y, *incy);
}

void FC_FUNC_(cublas_sspr,CUBLAS_SSPR) (const char *uplo, const int *n, const float *alpha, 
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtrAP)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);
    cublasSspr (uplo[0], *n, *alpha, x, *incx, AP);
}

void FC_FUNC_(cublas_sspr2,CUBLAS_SSPR2) (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrx, const int *incx, 
                   const devptr_t *devPtry, const int *incy,
                   const devptr_t *devPtrAP)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSspr2 (uplo[0], *n, *alpha, x, *incx, y, *incy, AP);
}

void FC_FUNC_(cublas_ssymv,CUBLAS_SSYMV) (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrx, const int *incx, const float *beta,
                   const devptr_t *devPtry,
                   const int *incy)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSsymv (uplo[0], *n, *alpha, A, *lda, x, *incx, *beta, y, *incy);
}

void FC_FUNC_(cublas_ssyr,CUBLAS_SSYR) (const char *uplo, const int *n, const float *alpha,
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtrA, const int *lda)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);    
    cublasSsyr (uplo[0], *n, *alpha, x, *incx, A, *lda);
}

void FC_FUNC_(cublas_ssyr2,CUBLAS_SSYR2) (const char *uplo, const int *n, const float *alpha,
                   const devptr_t *devPtrx, const int *incx,
                   const devptr_t *devPtry, const int *incy, 
                   const devptr_t *devPtrA, const int *lda)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);
    float *y = (float *)(uintptr_t)(*devPtry);    
    cublasSsyr2 (uplo[0], *n, *alpha, x, *incx, y, *incy, A, *lda);
}

void FC_FUNC_(cublas_stbmv,CUBLAS_STBMV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);    
    cublasStbmv (uplo[0], trans[0], diag[0], *n, *k, A, *lda, x, *incx);
}

void FC_FUNC_(cublas_stbsv,CUBLAS_STBSV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const int *k, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStbsv (uplo[0], trans[0], diag[0], *n, *k, A, *lda, x, *incx);
}

void FC_FUNC_(cublas_stpmv,CUBLAS_STPMV) (const char *uplo, const char *trans, const char *diag,
                   const int *n,  const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStpmv (uplo[0], trans[0], diag[0], *n, AP, x, *incx);
}

void FC_FUNC_(cublas_stpsv,CUBLAS_STPSV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const devptr_t *devPtrAP, 
                   const devptr_t *devPtrx, const int *incx)
{
    float *AP = (float *)(uintptr_t)(*devPtrAP);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStpsv (uplo[0], trans[0], diag[0], *n, AP, x, *incx);
}

void FC_FUNC_(cublas_strmv,CUBLAS_STRMV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStrmv (uplo[0], trans[0], diag[0], *n, A, *lda, x, *incx);
}

void FC_FUNC_(cublas_strsv,CUBLAS_STRSV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *x = (float *)(uintptr_t)(*devPtrx);       
    cublasStrsv (uplo[0], trans[0], diag[0], *n, A, *lda, x, *incx);
}

void FC_FUNC_(cublas_dgemv,CUBLAS_DGEMV) (const char *trans, const int *m, const int *n, 
                   const double *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrx, const int *incx,
                   const double *beta, const devptr_t *devPtry,
                   const int *incy)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);
    cublasDgemv (trans[0], *m, *n, *alpha, A, *lda, x, *incx, *beta, y, *incy);
}

void FC_FUNC_(cublas_dger,CUBLAS_DGER) (const int *m, const int *n, const double *alpha, 
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtry, const int *incy,
                  const devptr_t *devPtrA, const int *lda)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *x = (double *)(uintptr_t)(*devPtrx);
    double *y = (double *)(uintptr_t)(*devPtry);    
    cublasDger (*m, *n, *alpha, x, *incx, y, *incy, A, *lda);
}

void FC_FUNC_(cublas_dsyr,CUBLAS_DSYR) (const char *uplo, const int *n, const double *alpha,
                  const devptr_t *devPtrx, const int *incx,
                  const devptr_t *devPtrA, const int *lda)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *x = (double *)(uintptr_t)(*devPtrx);    
    cublasDsyr (uplo[0], *n, *alpha, x, *incx, A, *lda);
}

void FC_FUNC_(cublas_dtrsv,CUBLAS_DTRSV) (const char *uplo, const char *trans, const char *diag,
                   const int *n, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrx, const int *incx)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *x = (double *)(uintptr_t)(*devPtrx);       
    cublasDtrsv (uplo[0], trans[0], diag[0], *n, A, *lda, x, *incx);
}

/*---------------------------------------------------------------------------*/
/*---------------------------------- BLAS3 ----------------------------------*/
/*---------------------------------------------------------------------------*/

void FC_FUNC_(cublas_sgemm,CUBLAS_SGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const float *alpha,
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrB, const int *ldb, const float *beta,
                   const devptr_t *devPtrC, const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSgemm (transa[0], transb[0], *m, *n, *k, *alpha, A, *lda, 
                 B, *ldb, *beta, C, *ldc);
}

void FC_FUNC_(cublas_ssymm,CUBLAS_SSYMM) (const char *side, const char *uplo, const int *m, 
                   const int *n, const float *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrB, const int *ldb, 
                   const float *beta, const devptr_t *devPtrC, const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSsymm (*side, *uplo, *m, *m, *alpha, A, *lda, B, *ldb, *beta, C,
                 *ldc);
}

void FC_FUNC_(cublas_ssyr2k,CUBLAS_SSYR2K) (const char *uplo, const char *trans, const int *n,
                    const int *k, const float *alpha, const devptr_t *devPtrA,
                    const int *lda, const devptr_t *devPtrB, const int *ldb, 
                    const float *beta, const devptr_t *devPtrC, const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSsyr2k (*uplo, *trans, *n, *k, *alpha, A, *lda, B, *ldb, *beta, 
                  C, *ldc);
}

void FC_FUNC_(cublas_ssyrk,CUBLAS_SSYRK) (const char *uplo, const char *trans, const int *n, 
                   const int *k, const float *alpha, const devptr_t *devPtrA, 
                   const int *lda, const float *beta, const devptr_t *devPtrC,
                   const int *ldc)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *C = (float *)(uintptr_t)(*devPtrC);
    cublasSsyrk (*uplo, *trans, *n, *k, *alpha, A, *lda, *beta, C, *ldc);
}

void FC_FUNC_(cublas_strmm,CUBLAS_STRMM) (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n,
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb)
{
    float *A = (float *)(uintptr_t)(*devPtrA);
    float *B = (float *)(uintptr_t)(*devPtrB);
    cublasStrmm (*side, *uplo, *transa, *diag, *m, *n, *alpha, A, *lda, B,
                 *ldb);
}

void FC_FUNC_(cublas_strsm,CUBLAS_STRSM) (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n, 
                   const float *alpha, const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb)
{
    float *A = (float *)(uintptr_t)*devPtrA;
    float *B = (float *)(uintptr_t)*devPtrB;
    cublasStrsm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha,
                 A, *lda, B, *ldb);
}

void FC_FUNC_(cublas_cgemm,CUBLAS_CGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuComplex *alpha,
                   const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb, 
                   const cuComplex *beta, const devptr_t *devPtrC,
                   const int *ldc)
{
    cuComplex *A = (cuComplex *)(uintptr_t)*devPtrA;
    cuComplex *B = (cuComplex *)(uintptr_t)*devPtrB;
    cuComplex *C = (cuComplex *)(uintptr_t)*devPtrC;    
    cublasCgemm (transa[0], transb[0], *m, *n, *k, *alpha, A, *lda, B, *ldb, 
                 *beta, C, *ldc);
}

void FC_FUNC_(cublas_chemm,CUBLAS_CHEMM) (void)
{
    printf ("CUBLAS_CHEMM stub\n");
}
void FC_FUNC_(cublas_csymm,CUBLAS_CSYMM) (void)
{
    printf ("CUBLAS_CSYMM stub\n");
}
void FC_FUNC_(cublas_ctrmm,CUBLAS_CTRMM) (void)
{
    printf ("CUBLAS_CTRMM stub\n");
}
void FC_FUNC_(cublas_ctrsm,CUBLAS_CTRSM) (void)
{
    printf ("CUBLAS_CTRSM stub\n");
}
void FC_FUNC_(cublas_cherk,CUBLAS_CHERK) (void)
{
    printf ("CUBLAS_CHERK stub\n");
}
void FC_FUNC_(cublas_csyrk,CUBLAS_CSYRK) (void)
{
    printf ("CUBLAS_CSYRK stub\n");
}
void FC_FUNC_(cublas_cher2k,CUBLAS_CHER2K) (void)
{
    printf ("CUBLAS_CHER2K stub\n");
}
void FC_FUNC_(cublas_csyr2k,CUBLAS_CSYR2K) (void)
{
    printf ("CUBLAS_CSYR2K stub\n");
}

void FC_FUNC_(cublas_dgemm,CUBLAS_DGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const double *alpha,
                   const devptr_t *devPtrA, const int *lda, 
                   const devptr_t *devPtrB, const int *ldb, const double *beta,
                   const devptr_t *devPtrC, const int *ldc)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *B = (double *)(uintptr_t)(*devPtrB);
    double *C = (double *)(uintptr_t)(*devPtrC);
    cublasDgemm (transa[0], transb[0], *m, *n, *k, *alpha, A, *lda, 
                 B, *ldb, *beta, C, *ldc);
}

void FC_FUNC_(cublas_dsymm,CUBLAS_DSYMM) (const char *side, const char *uplo, const int *m, 
                   const int *n, const double *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrB, const int *ldb, 
                   const double *beta, const devptr_t *devPtrC, const int *ldc)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *B = (double *)(uintptr_t)(*devPtrB);
    double *C = (double *)(uintptr_t)(*devPtrC);
    cublasDsymm (*side, *uplo, *m, *m, *alpha, A, *lda, B, *ldb, *beta, C,
                 *ldc);
}

void FC_FUNC_(cublas_dsyr2k,CUBLAS_DSYR2K) (const char *uplo, const char *trans, const int *n,
                    const int *k, const double *alpha, const devptr_t *devPtrA,
                    const int *lda, const devptr_t *devPtrB, const int *ldb, 
                    const double *beta, const devptr_t *devPtrC,
                    const int *ldc)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *B = (double *)(uintptr_t)(*devPtrB);
    double *C = (double *)(uintptr_t)(*devPtrC);
    cublasDsyr2k (*uplo, *trans, *n, *k, *alpha, A, *lda, B, *ldb, *beta, 
                  C, *ldc);
}

void FC_FUNC_(cublas_dsyrk,CUBLAS_DSYRK) (const char *uplo, const char *trans, const int *n, 
                   const int *k, const double *alpha, const devptr_t *devPtrA, 
                   const int *lda, const double *beta, const devptr_t *devPtrC,
                   const int *ldc)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *C = (double *)(uintptr_t)(*devPtrC);
    cublasDsyrk (*uplo, *trans, *n, *k, *alpha, A, *lda, *beta, C, *ldc);
}

void FC_FUNC_(cublas_dtrmm,CUBLAS_DTRMM) (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n,
                   const double *alpha, const devptr_t *devPtrA, 
                   const int *lda, const devptr_t *devPtrB, const int *ldb)
{
    double *A = (double *)(uintptr_t)(*devPtrA);
    double *B = (double *)(uintptr_t)(*devPtrB);
    cublasDtrmm (*side, *uplo, *transa, *diag, *m, *n, *alpha, A, *lda, B,
                 *ldb);
}

void FC_FUNC_(cublas_dtrsm,CUBLAS_DTRSM) (const char *side, const char *uplo, const char *transa,
                   const char *diag, const int *m, const int *n, 
                   const double *alpha, const devptr_t *devPtrA,
                   const int *lda, const devptr_t *devPtrB, const int *ldb)
{
    double *A = (double *)(uintptr_t)*devPtrA;
    double *B = (double *)(uintptr_t)*devPtrB;
    cublasDtrsm (side[0], uplo[0], transa[0], diag[0], *m, *n, *alpha,
                 A, *lda, B, *ldb);
}

void FC_FUNC_(cublas_zgemm,CUBLAS_ZGEMM) (const char *transa, const char *transb, const int *m,
                   const int *n, const int *k, const cuDoubleComplex *alpha,
                   const devptr_t *devPtrA, const int *lda,
                   const devptr_t *devPtrB, const int *ldb, 
                   const cuDoubleComplex *beta, const devptr_t *devPtrC,
                   const int *ldc)
{
    cuDoubleComplex *A = (cuDoubleComplex *)(uintptr_t)*devPtrA;
    cuDoubleComplex *B = (cuDoubleComplex *)(uintptr_t)*devPtrB;
    cuDoubleComplex *C = (cuDoubleComplex *)(uintptr_t)*devPtrC;    
    cublasZgemm (transa[0], transb[0], *m, *n, *k, *alpha, A, *lda, B, *ldb, 
                 *beta, C, *ldc);
}

#endif  /* defined(CUBLAS_USE_THUNKING) */

