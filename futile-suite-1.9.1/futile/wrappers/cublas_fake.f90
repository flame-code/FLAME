!> @file 
!! Fake routines for cublas
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
  subroutine cublas_SSCAL()
    implicit none
    stop 'FAKE SSCAL'
  END SUBROUTINE cublas_SSCAL

  subroutine cublas_CSCAL()
    implicit none
    stop 'FAKE CSCAL'
  END SUBROUTINE cublas_CSCAL

  subroutine cublas_STRMM()
    implicit none
    stop 'FAKE STRMM'
  END SUBROUTINE cublas_STRMM

  subroutine cublas_DTRMM()
    implicit none
    stop 'FAKE DTRMM'
  END SUBROUTINE cublas_DTRMM

  subroutine cublas_CTRMM()
    implicit none
    stop 'FAKE CTRMM'
  END SUBROUTINE cublas_CTRMM

  subroutine cublas_SAXPY()
    implicit none
    stop 'FAKE SAXPY'
  END SUBROUTINE cublas_SAXPY

  subroutine cublas_CAXPY()
    implicit none
    stop 'FAKE CAXPY'
  END SUBROUTINE cublas_CAXPY

  subroutine cublas_SGEMM()
    implicit none
    stop 'FAKE SGEMM'
  END SUBROUTINE cublas_SGEMM

  subroutine cublas_CGEMM()
    implicit none
    stop 'FAKE CGEMM'
  END SUBROUTINE cublas_CGEMM

  subroutine cublas_SSYRK()
    implicit none
    stop 'FAKE SSYRK'
  END SUBROUTINE cublas_SSYRK

  subroutine cublas_CHERK()
    implicit none
    stop 'FAKE CHERK'
  END SUBROUTINE cublas_CHERK

  subroutine cublas_SNRM2()
    implicit none
    stop 'FAKE SNRM2'
  END SUBROUTINE cublas_SNRM2

  subroutine cublas_DNRM2()
    implicit none
    stop 'FAKE DNRM2'
  END SUBROUTINE cublas_DNRM2

  subroutine cublas_SDOT()
    implicit none
    stop 'FAKE SDOT'
  END SUBROUTINE cublas_SDOT

  subroutine cublas_DDOT()
    implicit none
    stop 'FAKE DDOT'
  END SUBROUTINE cublas_DDOT

  subroutine cublas_DGEMM()
    implicit none
    stop 'FAKE DGEMM'
  END SUBROUTINE cublas_DGEMM

  subroutine cublas_DSYRK()
    implicit none
    stop 'FAKE DSYRK'
  END SUBROUTINE cublas_DSYRK

subroutine poisson_cublas_daxpy()
   implicit none
   stop 'poisson_cublas_daxpy'
 END SUBROUTINE poisson_cublas_daxpy
