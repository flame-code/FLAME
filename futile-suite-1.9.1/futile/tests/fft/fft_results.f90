!> @file
!!   Test the results for the FFT API 
!! @copyright
!!   Copyright (C) Stefan Goedecker, CEA Grenoble, 2002, Basel University, 2009
!!   Copyright (C) 2009-2013 BigDFT group
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS
program fft_results
  use futile
  !use numerics
  use f_random
  implicit none
  character(len=*), parameter :: inputs=&
       "- {name: ndims, shortname: n, default: 36,"//&
       "  help_string: Array of dimension 3 for the points of the simulation box,"//&
       "  help_dict: {Allowed values: list of integers}}"//f_cr//&
       "- {name: sizes, shortname: s, default: 10.0,"//&
       "  help_string: Array of dimension 3 for the size of the simulation box,"//&
       "  help_dict: {Allowed values: list of floating point numbers}}"//f_cr//&
       "- {name: pref, shortname: p, default: 2,"//&
       "  help_string: Array of dimension 3 for the reference impulses,"//&
       "  help_dict: {Allowed values: list of integers. You may also write 'auto' for random number repetition}}"//f_cr//&
       "- {name: nrep, shortname: r, default: 1,"//&
       "  help_string: Number of repetitions of the bench,"//&
       "  help_dict: {Allowed values: integer}}"//f_cr//&
       "- {name: quiet, shortname: q, default: No,"//&
       "  help_string: Quiet output, show only progress bar when applicable,"//&
       "  help_dict: {Allowed values: Boolean}}"


 
  logical :: quiet
  integer :: inzee,nrep,irep
  character(len=4) :: auto
  integer, dimension(3) :: n,pref
  real(f_double), dimension(3) :: L,h,sigma
  type(dictionary), pointer :: options
  real(f_double), dimension(:,:,:,:,:), allocatable :: zinout
  type(f_progress_bar) :: bar

  call f_lib_initialize()

  call yaml_argparse(options,inputs)
  
  n=options//'ndims'
  L=options//'sizes'
  auto=' ' 
  call f_err_open_try()
  pref=options//'pref'
  !perform a separate treatment for the reference impulses
  if (f_err_check()) then
     call f_err_close_try()
     auto=options//'pref'
     if (.not. (auto .eqv. 'auto')) call f_err_throw('Invalid reference impulses, input:'//auto)
  end if
  call f_err_close_try()

  nrep=options//'nrep'
  quiet=options//'quiet'

  call dict_free(options)
  zinout=f_malloc0([2,n(1),n(2),n(3),2],id='zinout')

  h=L/n
  sigma=L/10.0_f_double

  if (len_trim(auto) >0) then
     if (quiet) bar=f_progress_bar_new(nstep=nrep)
     do irep=1,nrep
        call f_random_number(harvest=pref,ranges=n/2)
        call verify_single_harmonic(pref,n,zinout,inzee,.true.,.not. quiet)
        if (quiet) call dump_progress_bar(bar,step=irep)
     end do
  else

!!$     call verify_1d_buildingblocks(pref(1),n(1),n(2)*n(3),zinout,inzee,.true.,.true.)
!!$     call verify_1d_buildingblocks(pref(1),n(1),n(2)*n(3),zinout,inzee,.false.,.false.)

     call verify_single_harmonic(pref,n,zinout,inzee,.false.,.not. quiet)
  end if

  call take_timings(nrep,n,zinout,inzee)

  call f_free(zinout)
  call f_lib_finalize()

end program fft_results

subroutine take_timings(nrep,n,zinout,inzee) 
  use futile
  implicit none
  integer, intent(in) :: nrep
  integer, dimension(3) :: n
  integer, intent(inout) :: inzee
  real(f_double), dimension(2,n(1),n(2),n(3),2), intent(inout) :: zinout
  !local variables
  integer :: irep
  integer(f_long) :: t0,t1
  !then take timing for the number of repetitions

  if (nrep > 100) call f_zero(zinout)
  t0=f_time()
  do irep=1,nrep
     !test the 3d FFT
     !out-of-place case of traditional output
     call FFT_3d(n(1),n(2),n(3),n(1),n(2),n(3),zinout,1,inzee)
  end do
  t1=f_time()
  call yaml_map('No. of repetitions',nrep)
  call yaml_map('Time per repetition per element (ns)',&
       1.e-3_f_double*(t1-t0)/nrep/product(n),fmt='(f12.7)')

end subroutine take_timings


subroutine verify_1d_buildingblocks(pref,n,ndat,zinout,inzee,transpose,realfft)
  use futile
  use wrapper_linalg
  implicit none
  logical, intent(in) :: transpose,realfft
  integer, intent(in) :: pref,n,ndat
  integer, intent(inout) :: inzee
  real(f_double), dimension(2,n*ndat,2), intent(out) :: zinout
  !local variables
  integer :: idat,i1,ind
  real(f_double) :: twopi

  twopi=6.283185307179586d0
  inzee=.if. realfft .then. 2 .else. 1
  do i1=1,n
     do idat=1,ndat
        ind=idat+(i1-1)*ndat
        zinout(2,ind,inzee)=0.0_f_double
        zinout(1,ind,inzee)=cos((i1-1)*twopi/n*pref)
     end do
  end do
  do i1=1,n
     print *,'input',zinout(1,(i1-1)*ndat+1,inzee),zinout(1,(i1-1)*ndat+ndat,inzee)
  end do

  if (realfft) then
     call vcopy(n*ndat,zinout(1,1,2),2,zinout(1,1,1),1)
  end if

  print *,'testcopy',zinout(:,:,2)
  print *,'test2',zinout(:,:,1)

  call FFT_1d(n,ndat/2,zinout,1,inzee,transpose,.false.)

  do i1=1,n/2
     print *,'output',zinout(1,i1,inzee),zinout(1,i1+(ndat-1)*(n/2),inzee)
  end do
  print *,'testout',zinout(:,:,inzee)

  call validate(pref,ndat,n,zinout(1,1,inzee),transpose)
    
end subroutine verify_1d_buildingblocks

subroutine validate(pref,ndat,n,zinout,transpose)
  use futile
!  use numerics
  implicit none
  logical, intent(in) :: transpose 
  integer, intent(in) :: ndat,n,pref
  real(f_double), dimension(2,n*ndat) :: zinout
  !local variables
  integer :: i,idat,ind1,ind2
  real(f_double) :: sum_tot,test1,test2,onehalf

  onehalf = 0.5d0
  !verification of the result. Only two values should be different from zero
  do idat=1,ndat
     if (transpose) then
        ind1=j(pref,n)+1+(idat-1)*n
        ind2=j(-pref,n)+1+(idat-1)*n
     else
        ind1=idat+j(pref,n)*ndat
        ind2=idat+j(-pref,n)*ndat
     end if
     test1=zinout(1,ind1)
     test2=zinout(1,ind2)
     call f_assert(test1-onehalf*real(n,f_double),id='test1')
     call f_assert(test2-onehalf*real(n,f_double),id='test2')

     zinout(1,ind1)=zinout(1,ind1)-test1
     zinout(1,ind2)=zinout(1,ind2)-test2

     call f_zero(sum_tot)
     do i=1,n
        if (transpose) then
           ind1=i+(idat-1)*n
        else
           ind1=idat+(i-1)*ndat
        end if
        sum_tot=sum_tot+abs(zinout(1,ind1))+abs(zinout(2,ind1))
     end do
     call f_assert(sum_tot,id='total_sum')
  end do

contains

  !>real space coordinate, from 0,...,n-1
  pure function j(p,n)
    implicit none
    integer, intent(in) :: p,n
    integer :: j
    j=p-((p+n)/n-1)*n
  end function j
end subroutine validate

subroutine verify_single_harmonic(pref,n,zinout,inzee,assert,dump)
  use futile
  implicit none
  logical, intent(in) :: assert !<if true only assert the validity of the test
  logical, intent(in) :: dump
  integer, dimension(3), intent(in) :: pref,n
  integer, intent(inout) :: inzee
  real(f_double), dimension(2,n(1),n(2),n(3),2), intent(inout) :: zinout
  !local variables
  real(f_double), parameter :: tolerance=1.e-12_f_double
  integer :: i1,i2,i3,i
  real(f_double) :: refres,renorm,maxdiff,p2
  real(f_double), dimension(3) :: r,pval
  real(f_double), dimension(3,3) :: gmunu
  real(f_double) :: twopi

  twopi=6.283185307179586d0
  inzee=1
  !fill the inout array
  do i3=1,n(3)
     r(3)=(i3-1)!*h(3)
     do i2=1,n(2)
        r(2)=(i2-1)!*h(2)
        do i1=1,n(1)
           r(1)=(i1-1)!*h(1)
           !r2=sum(r**2/sigma)
           !f=safe_exp(-onehalf*r2)
           zinout(1,i1,i2,i3,inzee)=product(cos(r*twopi/n*pref))
           zinout(2,i1,i2,i3,inzee)=0.0_f_double
        end do
     end do
  end do

  if (f_debug) print *,'testin',zinout(1,:,1,1,inzee)

  !test the 3d FFT
  !out-of-place case of traditional output
  call FFT_3d(n(1),n(2),n(3),n(1),n(2),n(3),zinout,1,inzee)
  
  if (f_debug) print *,'testout',maxloc(zinout)

  if (dump) call yaml_map('Reference impulse (p_ref)',pref,advance=.if. assert .then. 'no' .else. 'yes')
  refres=sum(abs(zinout(:,:,:,:,inzee)))-product(n)
  if (assert) then
     call f_assert(abs(refres/product(n)) < tolerance,id=yaml_toa(refres))
     if (dump) call yaml_comment('...OK')
  else
     if (dump) call yaml_map('Value at p_ref',&
          zinout(:,j(pref(1),n(1))+1,j(pref(2),n(2))+1,j(pref(3),n(3))+1,inzee)/product(n)*8.0_f_double)
     if (dump) call yaml_map('Symmetric value',&
          zinout(:,j(-pref(1),n(1))+1,j(-pref(2),n(2))+1,j(-pref(3),n(3))+1,inzee)/product(n)*8.0_f_double)
     if (dump) call yaml_map('Other values',refres)
  end if
     !then inspect results of the modification of coordinates
  do i1=1,n(1)
     call f_assert(j(p(i1,n(1)),n(1))+1 == i1,id=trim(yaml_toa(i1,fmt='(i6)')))
     !print *,'p',i1,p(i1,n(1)),n(1)+2-i1,p(n(1)+2-i1,n(1))
     call f_assert(-p(i1,n(1)) == p(n(1)+2-i1,n(1)) .or. i1==n(1)/2+1 ,id=trim(yaml_toa(n(1)+2-i1,fmt='(i6)')))
     !print *,'j',i1,j(-p(i1,n(1)),n(1))+1,n(1)+2-i1
     call f_assert( j(-p(i1,n(1)),n(1)) == n(1)+1-i1 .or. i1==1 ,id=trim(yaml_toa(n(1)+2-i1,fmt='(i6)')))
  end do

  call f_zero(gmunu)
  forall(i=1:3) gmunu(i,i)=1.0_f_double

  !calculate the derivatives of the input function
  do i3=1,n(3)
     pval(3)=real(p(i3,n(3)),f_double)/n(3)
     do i2=1,n(2)
        pval(2)=real(p(i2,n(2)),f_double)/n(2)
        do i1=1,n(1)
           pval(1)=real(p(i1,n(1)),f_double)/n(1)
           p2=gmunu(1,1)*pval(1)**2+&
                gmunu(2,2)*pval(2)**2+&
                gmunu(3,3)*pval(3)**2+&
                2.0_f_double*gmunu(1,2)*pval(1)*pval(2)+&
                2.0_f_double*gmunu(1,3)*pval(1)*pval(3)+&
                2.0_f_double*gmunu(2,3)*pval(2)*pval(3)
           zinout(1,i1,i2,i3,inzee)=&
                zinout(1,i1,i2,i3,inzee)*p2
           zinout(2,i1,i2,i3,inzee)=&
                zinout(2,i1,i2,i3,inzee)*p2
        end do
     end do
  end do
  !now coming back
  !out-of-place case of traditional output
  call FFT_3d(n(1),n(2),n(3),n(1),n(2),n(3),zinout,-1,inzee)

  renorm=1.0_f_double/real(product(n),f_double)

  if (f_debug) print *,'test',zinout(1,:,1,1,inzee)*renorm

  !fill the inout array
  call f_zero(maxdiff)
  pval=real(pref,f_double)
  do i3=1,n(3)
     r(3)=(i3-1)!*h(3)
     do i2=1,n(2)
        r(2)=(i2-1)!*h(2)
        do i1=1,n(1)
           r(1)=(i1-1)!*h(1)
           maxdiff=max(maxdiff,&
                abs(zinout(1,i1,i2,i3,inzee)*renorm-&
                sum((pval/n)**2)*&
                product(cos(r*twopi/n*pref))))
           if (f_debug .and. i2==1 .and. i3==1) then
              print *,'test0',zinout(1,i1,i2,i3,inzee)*renorm,&
                   sum((abs(pval)/n)**2)*&
                   product(cos(r*twopi/n*pref))
           end if
           maxdiff=max(maxdiff,abs(zinout(2,i1,i2,i3,inzee)))
        end do
     end do
  end do
  
  if (dump) call yaml_map('Maximum difference',maxdiff)
  if (assert) call f_assert(maxdiff,id='Maximum difference',&
       tol=1.e-10_f_double*maxval(n))

  contains

    !>impulse coordinate, from 0,...,n/2+1,-n/2+1,...-1
    pure function p(i,n)
      implicit none
      integer, intent(in) :: i,n
      integer :: p
      p=i-(i/(n/2+2))*n-1
    end function p

    !>real space coordinate, from 0,...,n-1
    pure function j(p,n)
      implicit none
      integer, intent(in) :: p,n
      integer :: j
      j=p-((p+n)/n-1)*n
    end function j


end subroutine verify_single_harmonic
