!> @file
!!  Define routines for performance information
!! @author
!!    Copyright (C) 2010-2013 BigDFT group (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
module f_perfs
  use dictionaries
  use f_precisions
  implicit none

  private

  character(len=*), parameter, public :: F_PERF_GFLOPS='Flops '
  character(len=*), parameter, public :: F_PERF_READB ='ReadB '
  character(len=*), parameter, public :: F_PERF_WRITEB='WriteB'

  character(len=len(F_PERF_GFLOPS)), dimension(3), parameter :: keys=[&
       F_PERF_GFLOPS,F_PERF_READB,F_PERF_WRITEB]

  type, public :: f_perf
     type(dictionary), pointer :: model => null() !performance model
  end type f_perf

  interface f_perf_set_model
     module procedure f_perf_set_model,f_perf_set_model_i
  end interface f_perf_set_model

  public :: f_perf_set_model,f_perf_aggregate,f_perf_free

  contains

    pure function f_perf_null()
      type(f_perf) :: f_perf_null
    end function f_perf_null

    subroutine ensure_initialized(perf)
      implicit none
      type(f_perf), intent(inout) :: perf
      !local variables
      integer :: ikey

      if (associated(perf%model)) return
      
      call dict_init(perf%model)
      do ikey=1,size(keys)
         call set(perf%model//keys(ikey),0)
      end do
      
    end subroutine ensure_initialized
   
    subroutine f_perf_set_model_i(perf,quantity,val)
      implicit none
      type(f_perf), intent(inout) :: perf
      character(len=*), intent(in) :: quantity
      integer(f_integer), intent(in) :: val
      call f_perf_set_model(perf,quantity,int(val,f_long))
    end subroutine f_perf_set_model_i

    subroutine f_perf_set_model(perf,quantity,val)
      !use yaml_output
      implicit none
      type(f_perf), intent(inout) :: perf
      character(len=*), intent(in) :: quantity
      integer(f_long), intent(in) :: val

      call ensure_initialized(perf)
      !call yaml_map('test1',trim(quantity))
      !call yaml_map('test',perf%model)
      if (quantity .notin. perf%model) then
         call f_err_throw('Quantity '//trim(quantity)//' not allowed in the performance model')
         return
      end if
      !ival_tmp=perf%model//quantity
      call set(perf%model//quantity,val)
      
    end subroutine f_perf_set_model
 
    subroutine f_perf_aggregate(dest,src)
      !use yaml_output
      implicit none
      type(dictionary), pointer :: dest
      type(f_perf), intent(in) :: src
      !local variables
      integer :: ikey
      integer(f_long) :: ival_tmp,ival

      !aggregate performances in the dictionary. Create them if absent
      do ikey=1,size(keys)
         ival_tmp = 0
         ival_tmp = dest .get. keys(ikey)
         !call yaml_map('test',src%model)
         ival=src%model//keys(ikey)
         call set(dest//keys(ikey),ival+ival_tmp)           
      end do
    end subroutine f_perf_aggregate

    subroutine f_perf_free(perf)
      type(f_perf), intent(inout) :: perf
      call dict_free(perf%model)
      perf=f_perf_null()
    end subroutine f_perf_free

end module f_perfs
