!> @file
!!  Define a module to wrap the linear algebra routines
!! @author
!!    Copyright (C) 2013-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Modules which defines wrappers for the linear alegra.
module wrapper_linalg
  use time_profiling, only: TIMING_UNINITIALIZED
  use f_precisions
  implicit none


!!$  !> central type to handle the behaviour of the linear
!!$  !! problem
!!$  type, private :: f_linalg_data
!!$     !> all the typical arguments of ev and sv routines
!!$     integer :: n
!!$     [etc,etc..]
!!$  end type f_linalg_data
!!$
!!$  !exemple of API, inspired from numpy
!!$  eig=f_linalg_eig(A=mat,n=n,symmetric=.true.,hermitian=.false.,complex=.true.,work=work,lwork=lwork,eigenvalues=eval)
!!$
!!$  eig=f_linalg_solve(A=mat,n=n,symmetric=.true.,hermitian=.false.,complex=.true.,work=work,lwork=lwork,eigenvalues=eval)
!!$
!!$  !then retrieve the information, if not specified in the caller
!!$  call f_get_evals(eig,eval)
!!$  !or maybe by 
!!$  call f_memcpy(eig%eval_ptr,eval)
!!$
!!$  !all the thinkg which have not been specified are then destroyed by
!!$  call f_linalg_free(eig)


  !> Flag for GPU computing, if CUDA libraries are present
  !! in that case if a GPU is present a given MPI processor may or not perform a GPU calculation
  !! this value can be changed in the read_input_variables routine
  logical :: GPUblas=.false.

  !>timing categories of the module
  integer, save, private :: TCAT_COPY_ARRAYS=TIMING_UNINITIALIZED
  integer, save, private :: TCAT_BLAS_GEMM  =TIMING_UNINITIALIZED
  integer, save, private :: TCAT_LAPACK_EV  =TIMING_UNINITIALIZED

  !> interfaces for LAPACK routines
  interface potrf
     module procedure potrf_simple,potrf_double
  end interface
  interface c_potrf
     module procedure c_potrf_simple,c_potrf_double
  end interface
  interface trtri
     module procedure trtri_simple,trtri_double
  end interface 
  interface c_trtri
     module procedure c_trtri_simple,c_trtri_double
  end interface
  interface syev
     module procedure syev_simple,syev_double
  end interface
  interface heev
     module procedure heev_simple,heev_double
  end interface
  interface sygv
     module procedure sygv_simple,sygv_double
  end interface
  interface hegv
     module procedure hegv_simple,hegv_double
  end interface
  interface gesv
     module procedure gesv_simple,gesv_double
  end interface
  interface c_gesv
     module procedure c_gesv_simple,c_gesv_double
  end interface


  !> interfaces for BLAS routines
  interface gemm
     module procedure gemm_simple,gemm_double
  end interface
  interface gemmsy
     module procedure gemm_simple,gemmsy_double_wrap
  end interface
  interface c_gemm
     module procedure c_gemm_simple,c_gemm_double
  end interface
  interface dot
     module procedure dot_simple,dot_double
  end interface
  interface dotc
     module procedure dotc_simple,dotc_double
  end interface
  interface nrm2
     module procedure nrm2_simple,nrm2_double
  end interface
  interface vscal
     module procedure scal_simple,scal_double
  end interface
  interface vcopy
     module procedure copy_integer,copy_simple,copy_double,copy_double_to_simple,&
          copy_complex_real_simple,copy_complex_real_double
  end interface vcopy
  interface c_vscal
     module procedure c_scal_simple,c_scal_double
  end interface c_vscal
  interface syrk
     module procedure syrk_simple,syrk_double
  end interface syrk
  interface herk
     module procedure herk_simple,herk_double
  end interface herk
  interface trmm
     module procedure trmm_simple,trmm_double
  end interface trmm
  interface c_trmm
     module procedure c_trmm_simple,c_trmm_double
  end interface c_trmm
  interface axpy
     module procedure axpy_simple,axpy_double,axpy_simple_to_double
  end interface axpy
  interface swap
     module procedure swap_double
  end interface swap
  interface c_axpy
     module procedure c_axpy_simple,c_axpy_double
  end interface c_axpy

contains

  subroutine linalg_initialize_timing_categories()
    use time_profiling, only: f_timing_category,f_timing_category_group
    implicit none
    character(len=*), parameter :: flibgrp='Flib LowLevel' !<this will be moved
    character(len=*), parameter :: linalgrp='BLAS-LAPACK' !<this will be moved

    call f_timing_category('Vector copy',flibgrp,&
         'Memory copy of arrays (excluded allocations)',&
         TCAT_COPY_ARRAYS)
    call f_timing_category_group(linalgrp,&
         'Basic Linear Algebra and Matrix Manupulation Subprograms (external lib)')
    call f_timing_category('Blas (d-s-c-z)GeMM',linalgrp,&
         'Blas General Matrix-Matrix multiplications of any float type',&
         TCAT_BLAS_GEMM)
    call f_timing_category('Lapack (dsy-ssy-che-zhe)eev',linalgrp,&
         'Lapack Eigenvalue Problem',&
         TCAT_LAPACK_EV)


  end subroutine linalg_initialize_timing_categories

  !> Interfaces for LAPACK routines
  !! @warning
  !!   In these interfaces the input arrays are declared as scalars,
  !!   so the passage of the arguments by addresses is compulsory when calling
  !!   these routines
  !> Cholesky factorization of a positive definite matrix
  subroutine potrf_simple(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call spotrf(uplo,n,a,lda,info)
  end subroutine potrf_simple

  subroutine potrf_double(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call dpotrf(uplo,n,a,lda,info)
  end subroutine potrf_double

  subroutine c_potrf_simple(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call cpotrf(uplo,n,a,lda,info)
  end subroutine c_potrf_simple

  subroutine c_potrf_double(uplo,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call zpotrf(uplo,n,a,lda,info)
  end subroutine c_potrf_double

  !TRiangular matrix Inverse
  subroutine trtri_simple(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call strtri(uplo,diag,n,a,lda,info)
  end subroutine trtri_simple

  subroutine trtri_double(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call dtrtri(uplo,diag,n,a,lda,info)
  end subroutine trtri_double

  subroutine c_trtri_simple(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a
    !call to LAPACK routine
    call ctrtri(uplo,diag,n,a,lda,info)
  end subroutine c_trtri_simple

  subroutine c_trtri_double(uplo,diag,n,a,lda,info)
    implicit none
    character(len=1), intent(in) :: uplo,diag
    integer, intent(in) :: lda,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a
    !call to LAPACK routine
    call ztrtri(uplo,diag,n,a,lda,info)
  end subroutine c_trtri_double

  subroutine syev_simple(jobz,uplo,n,a,lda,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,work
    real(kind=4), intent(out) :: w
    call f_timer_interrupt(TCAT_LAPACK_EV)
    !call to LAPACK routine
    call ssyev(jobz,uplo,n,a,lda,w,work,lwork,info)
    call f_timer_resume()
  end subroutine syev_simple

  subroutine syev_double(jobz,uplo,n,a,lda,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,work
    real(kind=8), intent(out) :: w
    call f_timer_interrupt(TCAT_LAPACK_EV)
    !call to LAPACK routine
    call dsyev(jobz,uplo,n,a,lda,w,work,lwork,info)
    call f_timer_resume()
  end subroutine syev_double

  subroutine heev_simple(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,work,rwork
    real(kind=4), intent(out) :: w
    call f_timer_interrupt(TCAT_LAPACK_EV)
    !call to LAPACK routine
    call cheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
    call f_timer_resume()
  end subroutine heev_simple

  subroutine heev_double(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,work,rwork
    real(kind=8), intent(out) :: w
    call f_timer_interrupt(TCAT_LAPACK_EV)
    !call to LAPACK routine
    call zheev(jobz,uplo,n,a,lda,w,work,lwork,rwork,info)
    call f_timer_resume()
  end subroutine heev_double

  subroutine sygv_simple(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: itype,lda,ldb,lwork,n
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,b,work
    real(kind=4), intent(out) :: w
    !call to LAPACK routine
    call ssygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
  end subroutine sygv_simple

  subroutine sygv_double(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: itype,lda,ldb,lwork,n
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,b,work
    real(kind=8), intent(out) :: w
    !call to LAPACK routine
    call dsygv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,info)
  end subroutine sygv_double

  subroutine hegv_simple(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n,itype,ldb
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,work,rwork,b
    real(kind=4), intent(out) :: w
    !call to LAPACK routine
    call chegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
  end subroutine hegv_simple

  subroutine gesv_double(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call dgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine gesv_double

  subroutine gesv_simple(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call sgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine gesv_simple

  subroutine hegv_double(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
    implicit none
    character(len=1), intent(in) :: jobz,uplo
    integer, intent(in) :: lda,lwork,n,itype,ldb
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,work,rwork,b
    real(kind=8), intent(out) :: w
    !call to LAPACK routine
    call zhegv(itype,jobz,uplo,n,a,lda,b,ldb,w,work,lwork,rwork,info)
  end subroutine hegv_double

  subroutine c_gesv_double(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=8), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call zgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine c_gesv_double

  subroutine c_gesv_simple(n,nrhs,a,lda,ipiv,b,ldb,info)
    implicit none
    integer, intent(in) :: n,lda,nrhs,ldb
    integer, intent(out) :: info
    real(kind=4), intent(inout) :: a,b
    integer, intent(out) :: ipiv
    !call to LAPACK routine
    call cgesv(n,nrhs,a,lda,ipiv,b,ldb,info)
  end subroutine c_gesv_simple


  !> Interfaces for BLAS routines
  !! @warning
  !!         In these interfaces the input arrays are declared as scalars,
  !!         so the passage of the arguments by addresses is compulsory when calling
  !!         these routines

  !SCALe a vector by a constant
  subroutine scal_simple(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(inout) :: dx
    !call to BLAS routine
    call SSCAL(n,da,dx,incx)
  end subroutine scal_simple

  subroutine scal_double(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(inout) :: dx
    !call to BLAS routine
    call DSCAL(n,da,dx,incx)
  end subroutine scal_double

  subroutine c_scal_simple(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(out) :: dx
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CSCAL(n,da,dx,incx)
    else
       !call to BLAS routine
       call CSCAL(n,da,dx,incx)
    end if
  end subroutine c_scal_simple

  subroutine c_scal_double(n,da,dx,incx)
    implicit none
    integer, intent(in) :: incx,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(out) :: dx
    !call to BLAS routine
    call ZSCAL(n,da,dx,incx)
  end subroutine c_scal_double

  !copy the vector
  subroutine copy_complex_real_simple(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    complex(kind=4), intent(in) :: dx
    real(kind=4), intent(out) :: dy
    !call to BLAS routine
    call f_timer_interrupt(TCAT_COPY_ARRAYS) 
    call SCOPY(n,dx,incx,dy,incy)
    call f_timer_resume() 
  end subroutine copy_complex_real_simple

  subroutine copy_complex_real_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    complex(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: dy
    !call to BLAS routine
    call f_timer_interrupt(TCAT_COPY_ARRAYS) 
    call DCOPY(n,dx,incx,dy,incy)
    call f_timer_resume() 
  end subroutine copy_complex_real_double

  subroutine copy_integer(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    integer, intent(in) :: dx
    integer, intent(out) :: dy
    !custom blas routine
    call f_timer_interrupt(TCAT_COPY_ARRAYS) 
    call icopy(n,dx,incx,dy,incy)
    call f_timer_resume() 
  end subroutine copy_integer

  subroutine copy_simple(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(out) :: dy
    !call to BLAS routine
    call f_timer_interrupt(TCAT_COPY_ARRAYS) 
    call SCOPY(n,dx,incx,dy,incy)
    call f_timer_resume() 
  end subroutine copy_simple

  subroutine copy_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(out) :: dy
    !call to BLAS routine
    call f_timer_interrupt(TCAT_COPY_ARRAYS) 
    call DCOPY(n,dx,incx,dy,incy)
    call f_timer_resume() 
  end subroutine copy_double

  subroutine copy_double_to_simple(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: dx
    real(kind=4), intent(out) :: dy
    !call to custom routine
    call f_timer_interrupt(TCAT_COPY_ARRAYS) 
    call dscopy(n,dx,incx,dy,incy)
    call f_timer_resume() 
  end subroutine copy_double_to_simple

  subroutine trmm_simple(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    real(kind=4), intent(in) :: alpha
    real(kind=4), intent(in) :: a
    real(kind=4), intent(inout) :: b
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_STRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    else
       !call to BLAS routine
       call STRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    end if
  end subroutine trmm_simple

  subroutine trmm_double(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    real(kind=8), intent(in) :: alpha
    real(kind=8), intent(in) :: a
    real(kind=8), intent(inout) :: b
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_DTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    else
       !call to BLAS routine
       call DTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    end if
  end subroutine trmm_double

  subroutine c_trmm_simple(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    complex(kind=4), intent(in) :: alpha
    real(kind=4), intent(in) :: a
    real(kind=4), intent(inout) :: b
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    else
       !call to BLAS routine
       call CTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    end if
  end subroutine c_trmm_simple

  subroutine c_trmm_double(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
    implicit none
    character(len=1), intent(in) :: side,uplo,transa,diag
    integer, intent(in) :: lda,ldb,m,n
    complex(kind=8), intent(in) :: alpha
    real(kind=8), intent(in) :: a
    real(kind=8), intent(inout) :: b
    !call to BLAS routine
    call ZTRMM(side,uplo,transa,diag,m,n,alpha,a,lda,b,ldb)
  end subroutine c_trmm_double

  subroutine axpy_simple(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(inout) :: dy
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_SAXPY(n,da,dx,incx,dy,incy)
    else
       !call to BLAS routine
       call SAXPY(n,da,dx,incx,dy,incy)
    end if
  end subroutine axpy_simple

  subroutine axpy_double(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(inout) :: dy
    !call to BLAS routine
    call DAXPY(n,da,dx,incx,dy,incy)
  end subroutine axpy_double

  subroutine axpy_simple_to_double(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: da
    real(kind=4), intent(in) :: dx
    real(kind=8), intent(inout) :: dy
    !call to custom routine, for mixed precision sum
    call dasxpdy(n,da,dx,incx,dy,incy)
  end subroutine axpy_simple_to_double

  subroutine c_axpy_simple(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=4), intent(in) :: da
    real(kind=4), intent(in) :: dx
    real(kind=4), intent(inout) :: dy
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CAXPY(n,da,dx,incx,dy,incy)
    else
       !call to BLAS routine
       call CAXPY(n,da,dx,incx,dy,incy)
    end if
  end subroutine c_axpy_simple

  subroutine c_axpy_double(n,da,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: incx,incy,n
    real(kind=8), intent(in) :: da
    real(kind=8), intent(in) :: dx
    real(kind=8), intent(inout) :: dy
    !call to BLAS routine
    call ZAXPY(n,da,dx,incx,dy,incy)
  end subroutine c_axpy_double

  subroutine swap_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    real(kind=8), intent(inout) :: dx
    real(kind=8), intent(inout) :: dy
    !call to BLAS routine
    call DSWAP(n,dx,incx,dy,incy)
  end subroutine swap_double


  !euclidean dot product
  function dot_simple(n,sx,incx,sy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    real(kind=4) :: sx,sy
    real(kind=4) :: dot_simple
    !local variables
    real(kind=4) :: cublas_sdot,sdot
    if (GPUblas) then
       !call to CUBLAS function
       dot_simple=cublas_sdot(n,sx,incx,sy,incy)
    else
       !call to BLAS function
       dot_simple=sdot(n,sx,incx,sy,incy)
    end if
  end function dot_simple

  !euclidean dot product
  function dotc_simple(n,sx,incx,sy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    complex(kind=4), intent(inout) :: sx,sy
    complex(kind=4) :: dotc_simple
    !local variables
    complex(kind=4) :: cdotc
    !call to BLAS function
    dotc_simple=cdotc(n,sx,incx,sy,incy)
  end function dotc_simple

  function dot_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    real(kind=8) :: dx,dy
    real(kind=8) :: dot_double
    !local variables
    real(kind=8) :: cublas_ddot,ddot
    if (GPUblas) then
       !call to CUBLAS function
       dot_double=cublas_ddot(n,dx,incx,dy,incy)
    else
       !call to BLAS function
       dot_double=ddot(n,dx,incx,dy,incy)
    end if
  end function dot_double

  function dotc_double(n,dx,incx,dy,incy)
    implicit none
    integer, intent(in) :: n,incx,incy
    complex(kind=8), intent(inout) :: dx,dy
    complex(kind=8) :: dotc_double
    !local variables
    complex(kind=8) :: zdotc
    !call to BLAS function
    dotc_double=zdotc(n,dx,incx,dy,incy)
  end function dotc_double

  !euclidean NoRM of a vector
  function nrm2_simple(n,x,incx)
    implicit none
    integer, intent(in) :: n,incx
    real(kind=4) :: x
    real(kind=4) :: nrm2_simple
    !local variables
    real(kind=4) :: cublas_snrm2,snrm2
    if (GPUblas .and. n>10000) then
       !call to CUBLAS function
       nrm2_simple=cublas_snrm2(n,x,incx)
    else
       !call to BLAS function
       nrm2_simple=snrm2(n,x,incx)
    end if
  end function nrm2_simple

  function nrm2_double(n,x,incx)
    implicit none
    integer, intent(in) :: n,incx
    real(kind=8) :: x
    real(kind=8) :: nrm2_double
    !local variables
    real(kind=8) :: cublas_dnrm2,dnrm2
    if (GPUblas .and. n>10000) then
       !call to CUBLAS function
       nrm2_double=cublas_dnrm2(n,x,incx)
    else
       !call to BLAS routine
       nrm2_double=dnrm2(n,x,incx)
    end if
  end function nrm2_double

  !GEneral Matrix-Matrix multiplication routines
  subroutine gemm_simple(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    real(kind=4), intent(in) :: alpha,beta
    real(kind=4) :: a
    real(kind=4) :: b
    real(kind=4) :: c
    call f_timer_interrupt(TCAT_BLAS_GEMM) 
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_SGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    else
       !call to BLAS routine
       call SGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    end if
    call f_timer_resume()
  end subroutine gemm_simple

  subroutine gemm_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8) :: a
    real(kind=8) :: b
    real(kind=8) :: c
    call f_timer_interrupt(TCAT_BLAS_GEMM)
    !call to BLAS routine
    if (GPUblas) then
       call cublas_DGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    else
       call DGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    end if
    call f_timer_resume()
  end subroutine gemm_double

  subroutine gemmsy_double_wrap(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8) :: a
    real(kind=8) :: b
    real(kind=8) :: c
    !call to BLAS routine
    call gemmsy_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
  end subroutine gemmsy_double_wrap

  subroutine c_gemm_simple(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    complex(kind=4), intent(in) :: alpha,beta
    real(kind=4) :: a
    real(kind=4) :: b
    real(kind=4) :: c
    call f_timer_interrupt(TCAT_BLAS_GEMM)
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    else
       !call to BLAS routine
       call CGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    end if
    call f_timer_resume()
  end subroutine c_gemm_simple

  subroutine c_gemm_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: transa,transb
    integer, intent(in) :: k,lda,ldb,ldc,m,n
    complex(kind=8), intent(in) :: alpha,beta
    real(kind=8) :: a
    real(kind=8) :: b
    real(kind=8) :: c
    !call to BLAS routine
    call f_timer_interrupt(TCAT_BLAS_GEMM)
    call ZGEMM(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,c,ldc)
    call f_timer_resume()
  end subroutine c_gemm_double

  !SYmmetric Rank K operation
  subroutine syrk_simple(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=4), intent(in) :: alpha,beta
    real(kind=4), intent(in) :: a
    real(kind=4), intent(out) :: c 
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_SSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    else
       !call to BLAS routine
       call SSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    end if
  end subroutine syrk_simple

  subroutine syrk_double(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8), intent(in) :: a
    real(kind=8), intent(out) :: c 
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_DSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    else
       !call to BLAS routine
       call DSYRK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    end if
  end subroutine syrk_double

  !HErmitian Rank K operation
  subroutine herk_simple(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=4), intent(in) :: alpha,beta
    real(kind=4), intent(in) :: a
    real(kind=4), intent(out) :: c 
    if (GPUblas) then
       !call to CUBLAS routine
       call cublas_CHERK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    else
       !call to BLAS routine
       call CHERK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    end if
  end subroutine herk_simple

  subroutine herk_double(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
    implicit none
    character(len=1), intent(in) :: trans,uplo
    integer, intent(in) :: k,lda,ldc,n
    real(kind=8), intent(in) :: alpha,beta
    real(kind=8), intent(in) :: a
    real(kind=8), intent(out) :: c 
    !call to BLAS routine
    call ZHERK(uplo,trans,n,k,alpha,a,lda,beta,c,ldc)
  end subroutine herk_double

  !> Determinant of a 3x3 matrix
  pure function det_3x3(a) result(det)
    implicit none
    real(f_double), dimension(3,3), intent(in) :: a
    real(f_double) :: det

    det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) &
         + a(1,2)*(a(3,1)*a(2,3) - a(2,1)*a(3,3))  &
         + a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
  end function det_3x3


end module wrapper_linalg
