!> @file
!!  Fake routine when gemmsy_double is not activated
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Fake routine when gemmsy_double is not activated
subroutine gemmsy_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,y,ldy)
  use wrapper_linalg
  implicit none
  character(len=1), intent (in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldy
  real(kind=8), intent(in) :: alpha,beta
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: b
  real(kind=8), intent(inout) :: y
  call gemm_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,y,ldy)

END SUBROUTINE gemmsy_double
