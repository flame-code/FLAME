!> @file
!!  Symmetric Multiplication Matrix-Matrix
!! @author 
!!    From The work of Éric Bainville (2010)
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 


!>  Symmetric matrix matrix multiplication (@f$ <\Psi|H|\Psi>@f$)
subroutine gemmsy_double(transa,transb,m,n,k,alpha,a,lda,b,ldb,beta,y,ldy)
  implicit none
  character(len=1), intent (in) :: transa, transb
  integer, intent(in) :: m, n, k, lda, ldb, ldy
  real(kind=8), intent(in) :: alpha,beta
  real(kind=8), intent(in) :: a
  real(kind=8), intent(in) :: b
  real(kind=8), intent(inout) :: y
  character(len=1) :: tra
  
  if ( m /= n ) then
    write(*,*)'ERROR (dgemmsy): the m and n dimensions differ: ',m,n
    stop 
  end if
  if ( transa == 't' .or. transa == 'T' ) tra = 'n'
  if ( transa == 'n' .or. transa == 'N' ) tra = 't'
  call dgemmsy(tra,transb,k,m,alpha,a,lda,b,ldb,beta,y,ldy,0)

END SUBROUTINE gemmsy_double
