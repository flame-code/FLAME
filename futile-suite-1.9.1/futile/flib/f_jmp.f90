!> @file
!! Manage the setjmp / longjmp control
!! like operations on external files and basic operations in memory
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_jmp
  use f_precisions
  implicit none
  
  !points to the data of the jmpbuf environment
  type, public :: f_jmpbuf
     integer :: signal=0 !<present jmp signal
     integer :: destroy_signal=0 !<signal of the last jmp return  
     character(len=128) :: id=' '!<identification of the jmp buffer
     character, dimension(:), pointer :: jmp_buf=>null()
     integer(f_long) :: t0=0 !<jmpbuf creation time
     integer(f_address) :: callback =0
  end type f_jmpbuf

  contains

    subroutine nullify_jmpbuf(jb)
      implicit none
      type(f_jmpbuf), intent(out) :: jb
      jb%signal=0
      jb%destroy_signal=0
      jb%id=' '
      nullify(jb%jmp_buf)
      jb%t0=int(0,f_long)
      jb%callback=int(0,f_long)
    end subroutine nullify_jmpbuf

    !> set the setjmp call and retrieve the data associated to the
    !!jmp_buf environment. Save them in the jmp_buf for later use
    subroutine f_jmpbuf_set(jmpbuf)
      use dynamic_memory
      implicit none
      type(f_jmpbuf), intent(inout) :: jmpbuf
      !local variables
      integer :: bufsize

      if (.not. associated(jmpbuf%jmp_buf)) then
         call getjmpbufsize(bufsize)
         !print *,'bufsize,set',bufsize
         jmpbuf%jmp_buf=f_malloc_str_ptr(len(jmpbuf%jmp_buf),bufsize,id='jmp_buf')
      end if
      !print *,'signal',jmpbuf%signal,jmpbuf%callback,jmpbuf%destroy_signal
      call setandcpyjmpbuf(jmpbuf%signal,jmpbuf%callback,jmpbuf,&
           jmpbuf%destroy_signal)

      !print *,'this part has to be written at the last'
      !print *,'sgn222',jmpbuf%signal!,jmpbuf%jmp_buf

  end subroutine f_jmpbuf_set

  subroutine f_longjmp(jmpbuf,signal)
    implicit none
    type(f_jmpbuf), intent(inout) :: jmpbuf
    integer, intent(in), optional :: signal
    !local variables
    integer :: sgn
    
    if (jmpbuf%signal == jmpbuf%destroy_signal) then
       call f_jmpbuf_free(jmpbuf)
    else
       sgn=jmpbuf%signal+1
       if (present(signal)) sgn=signal
       !print *,'signalbeforelj',sgn
       call longjmpwrap(jmpbuf%jmp_buf,sgn)
    end if
      
  end subroutine f_longjmp

    subroutine f_jmpbuf_free(jmpbuf)
      use dynamic_memory
      implicit none
      type(f_jmpbuf), intent(inout) :: jmpbuf
      call f_free_str_ptr(len(jmpbuf%jmp_buf),jmpbuf%jmp_buf)
      call nullify_jmpbuf(jmpbuf)
    end subroutine f_jmpbuf_free
    
end module f_jmp
