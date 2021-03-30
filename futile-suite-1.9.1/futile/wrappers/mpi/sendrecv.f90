!> @file
!! Wrapper for sendreceive operations
!! @author
!!    Copyright (C) 2012-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module f_sendrecv
  use time_profiling, only: TIMING_UNINITIALIZED
  use dictionaries, only: f_err_throw
  use f_precisions
  use fmpi_types
  use yaml_strings
  implicit none

  private

  integer, public, save :: TCAT_SEND         = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_RECV         = TIMING_UNINITIALIZED
  integer, public, save :: TCAT_WAIT         = TIMING_UNINITIALIZED

  interface fmpi_send
     module procedure mpisend_d0, mpisend_gpu,mpisend_i1
  end interface fmpi_send

  interface fmpi_recv
     module procedure mpirecv_d0,mpirecv_gpu,mpirecv_i1
  end interface fmpi_recv

  public :: fmpi_send,fmpi_recv,fmpi_waitall,fmpi_wait

  contains

    subroutine fmpi_waitall(ncount, array_of_requests,array_of_statuses,simulate)
      implicit none
      ! Local variables
      integer, intent(in) :: ncount
      integer, dimension(ncount),intent(in) :: array_of_requests
      logical, intent(in), optional :: simulate
      integer, dimension(FMPI_STATUS_SIZE,ncount), intent(out), optional :: array_of_statuses
      ! Local variables
      logical :: sim
      integer :: ierr,tcat

      !no wait if no requests
      if (ncount==0) return

      sim=.false.
      if (present(simulate)) sim=simulate
      if (sim) return

      tcat=TCAT_WAIT
      ! Synchronize the communication
      call f_timer_interrupt(tcat)
      if (present(array_of_statuses)) then
         call MPI_WAITALL(ncount, array_of_requests,array_of_statuses, ierr)
      else
         call MPI_WAITALL(ncount, array_of_requests, FMPI_STATUSES_IGNORE, ierr)
      end if
      call f_timer_resume()
      if (ierr/=0) then
         call f_err_throw('An error in calling to MPI_WAITALL occured',&
              err_id=ERR_MPI_WRAPPERS)
         return
      end if

    end subroutine fmpi_waitall

    subroutine fmpi_wait(request)
      implicit none
      ! Local variables
      integer,intent(in) :: request
      ! Local variables
      integer :: ierr,tcat

      if (request /= FMPI_REQUEST_NULL) then
         tcat=TCAT_WAIT
         ! Synchronize the communication
         call f_timer_interrupt(tcat)
         call MPI_WAIT(request, FMPI_STATUSES_IGNORE, ierr)
         call f_timer_resume()
         if (ierr/=0) then
            call f_err_throw('An error in calling to MPI_WAIT occured',&
                 err_id=ERR_MPI_WRAPPERS)
         end if
      end if
    end subroutine fmpi_wait

    subroutine mpisend_i1(buf,dest,tag,comm,request,simulate,verbose)
      use yaml_output
      implicit none
      integer(f_integer), dimension(:),  intent(in) :: buf
      integer, intent(in) :: dest
      integer, intent(in), optional :: tag
      integer, intent(in), optional :: comm
      integer, intent(out), optional :: request !<toggle the isend operatio
      logical, intent(in), optional :: simulate,verbose
      !local variables
      logical :: verb,sim
      integer :: mpi_comm,ierr,tag_,tcat,count

      mpi_comm=fmpi_comm(comm)
      if (present(tag)) then
         tag_=tag
      else
         tag_=mpirank(mpi_comm)
      end if

      verb=.false.
      if (present(verbose)) verb=verbose .and. dest /=mpirank_null()

      count=size(buf)

      if (verb) then
         call yaml_mapping_open('MPI_(I)SEND')
         call yaml_map('Elements',count)
         call yaml_map('Source',mpirank(mpi_comm))
         call yaml_map('Dest',dest)
         call yaml_map('Tag',tag_)
         call yaml_mapping_close()
      end if

      sim=.false.
      if (present(simulate)) sim=simulate
      if (sim) return

      tcat=TCAT_SEND
      ! Synchronize the communication
      call f_timer_interrupt(tcat)

      if (present(request)) then
         call MPI_ISEND(buf,count,mpitype(buf),dest,tag,mpi_comm,request,ierr)
      else
         call MPI_SEND(buf,count,mpitype(buf),dest,tag,mpi_comm,ierr)
      end if
      call f_timer_resume()
      if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)SEND occured',&
           err_id=ERR_MPI_WRAPPERS)

    end subroutine mpisend_i1

    subroutine mpisend_d0(buf,count,dest,tag,comm,request,simulate,verbose)
      use yaml_output
      implicit none
      real(f_double) :: buf !fake intent(in)
      integer, intent(in) :: count
      integer, intent(in) :: dest
      integer, intent(in), optional :: tag
      integer, intent(in), optional :: comm
      integer, intent(out), optional :: request !<toggle the isend operation
      logical, intent(in), optional :: simulate,verbose
      !local variables
      logical :: verb,sim
      integer :: mpi_comm,ierr,tag_,tcat

      mpi_comm=fmpi_comm(comm)
      if (present(tag)) then
         tag_=tag
      else
         tag_=mpirank(mpi_comm)
      end if

      verb=.false.
      if (present(verbose)) verb=verbose .and. dest /=mpirank_null()

      if (verb) then
         call yaml_mapping_open('MPI_(I)SEND')
         call yaml_map('Elements',count)
         call yaml_map('Source',mpirank(mpi_comm))
         call yaml_map('Dest',dest)
         call yaml_map('Tag',tag_)
         call yaml_mapping_close()
      end if

      sim=.false.
      if (present(simulate)) sim=simulate
      if (sim) return

      tcat=TCAT_SEND
      ! Synchronize the communication
      call f_timer_interrupt(tcat)

      if (present(request)) then
         call MPI_ISEND(buf,count,mpitype(buf),dest,tag,mpi_comm,request,ierr)
      else
         call MPI_SEND(buf,count,mpitype(buf),dest,tag,mpi_comm,ierr)
      end if
      call f_timer_resume()
      if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)SEND occured',&
           err_id=ERR_MPI_WRAPPERS)

    end subroutine mpisend_d0

    subroutine mpisend_gpu(buf,count,dest,tag,comm,request,simulate,verbose,type,offset)
      use yaml_output
      use iso_c_binding
      use f_precisions, only: f_address
      implicit none
      type(c_ptr) :: buf !fake intent(in)
      integer, intent(in) :: count
      integer, intent(in) :: dest
      integer, intent(in), optional :: tag
      integer, intent(in), optional :: comm
      integer, intent(out), optional :: request !<toggle the isend operation
      logical, intent(in), optional :: simulate,verbose
      integer, intent(in) :: type
      integer, intent(in), optional :: offset
      !local variables
      real(f_double),pointer :: a 
      logical :: verb,sim
      integer :: mpi_comm,ierr,tag_,tmpsize,tcat
      integer(f_address) tmpint
      type(c_ptr) :: tmpaddr

      mpi_comm=fmpi_comm(comm)
      if (present(tag)) then
         tag_=tag
      else
         tag_=mpirank(mpi_comm)
      end if

      verb=.false.
      if (present(verbose)) verb=verbose .and. dest /=mpirank_null()

      if (verb) then
         call yaml_mapping_open('MPI_(I)SEND')
         call yaml_map('Elements',count)
         call yaml_map('Source',mpirank(mpi_comm))
         call yaml_map('Dest',dest)
         call yaml_map('Tag',tag_)
         call yaml_mapping_close()
      end if

      sim=.false.
      if (present(simulate)) sim=simulate
      if (sim) return

      tcat=TCAT_SEND
      ! Synchronize the communication
      call f_timer_interrupt(tcat)

      !if(present(offset) .and. offset/=0)then
      !this is not guaranteed to work by fortran spec
      if (present(offset)) then
         if (offset /=0) then
            tmpint = TRANSFER(buf, tmpint)
            call mpi_type_size(type, tmpsize, ierr)
            tmpint = tmpint + offset*tmpsize
            tmpaddr= TRANSFER(tmpint, tmpaddr)
            call c_f_pointer(tmpaddr, a)
         else
            call c_f_pointer(buf, a)
         end if
      else
         call c_f_pointer(buf, a)
      end if
      if (present(request)) then
         call MPI_ISEND(a,count,type,dest,tag,mpi_comm,request,ierr)
      else
         call MPI_SEND(a,count,type,dest,tag,mpi_comm,ierr)
      end if
      call f_timer_resume()

      if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)SEND occured',&
           err_id=ERR_MPI_WRAPPERS)

    end subroutine mpisend_gpu

    subroutine mpirecv_i1(buf,source,tag,comm,status,request,simulate,verbose)
      use yaml_output
      implicit none
      integer(f_integer), dimension(:), intent(inout) :: buf !to be checked only after request is completed if needed
      integer, intent(in), optional :: source
      integer, intent(in), optional :: tag
      integer, intent(in), optional :: comm
      integer, intent(out), optional :: request !<toggle the isend operation
      integer, dimension(FMPI_STATUS_SIZE), intent(out), optional :: status !<for the blocking operation
      logical, intent(in), optional :: simulate,verbose
      !local variables
      logical :: verb,sim
      integer :: mpi_comm,ierr,mpi_source,mpi_tag,mpi_type,tcat,count

      mpi_comm=fmpi_comm(comm)
      mpi_source=FMPI_ANY_SOURCE
      mpi_tag=FMPI_ANY_TAG
      if (present(source)) then
         mpi_source=source
         mpi_tag=source
      end if
      if (present(tag)) mpi_tag=tag
      count=size(buf)
      verb=.false.
      if (present(verbose)) verb=verbose .and. mpi_source /= mpirank_null()
      if (verb) call yaml_comment('Receiving'//count//'elements from'//source//'in'//mpirank(mpi_comm))

      if (verb) then
         call yaml_mapping_open('MPI_(I)RECV')
         call yaml_map('Elements',count)
         call yaml_map('Source',source)
         call yaml_map('Dest',mpirank(mpi_comm))
         call yaml_map('Tag',mpi_tag)
         call yaml_mapping_close()
      end if

      sim=.false.
      if (present(simulate)) sim=simulate
      if (sim) return

      tcat=TCAT_RECV
      ! Synchronize the communication
      call f_timer_interrupt(tcat)

      mpi_type=mpitype(buf)
      if (present(request)) then
         call MPI_IRECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,request,ierr)
      else
         if (present(status)) then
            call MPI_RECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,status,ierr)
         else
            call MPI_RECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,FMPI_STATUS_IGNORE,ierr)
         end if
      end if
      call f_timer_resume()
      if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)RECV occured',&
           err_id=ERR_MPI_WRAPPERS)

    end subroutine mpirecv_i1


    subroutine mpirecv_d0(buf,count,source,tag,comm,status,request,simulate,verbose)
      use yaml_output
      implicit none
      real(f_double), intent(inout) :: buf !fake intent(out)
      integer, intent(in) :: count
      integer, intent(in), optional :: source
      integer, intent(in), optional :: tag
      integer, intent(in), optional :: comm
      integer, intent(out), optional :: request !<toggle the isend operation
      integer, dimension(FMPI_STATUS_SIZE), intent(out), optional :: status !<for the blocking operation
      logical, intent(in), optional :: simulate,verbose
      !local variables
      logical :: verb,sim
      integer :: mpi_comm,ierr,mpi_source,mpi_tag,mpi_type,tcat

      mpi_comm=fmpi_comm(comm)
      mpi_source=FMPI_ANY_SOURCE
      mpi_tag=FMPI_ANY_TAG
      if (present(source)) then
         mpi_source=source
         mpi_tag=source
      end if
      if (present(tag)) mpi_tag=tag
      verb=.false.
      if (present(verbose)) verb=verbose .and. mpi_source /= mpirank_null()
      if (verb) call yaml_comment('Receiving'//count//'elements from'//source//'in'//mpirank(mpi_comm))

      if (verb) then
         call yaml_mapping_open('MPI_(I)RECV')
         call yaml_map('Elements',count)
         call yaml_map('Source',source)
         call yaml_map('Dest',mpirank(mpi_comm))
         call yaml_map('Tag',mpi_tag)
         call yaml_mapping_close()
      end if

      sim=.false.
      if (present(simulate)) sim=simulate
      if (sim) return

      tcat=TCAT_RECV
      ! Synchronize the communication
      call f_timer_interrupt(tcat)

      mpi_type=mpitype(buf)
      if (present(request)) then
         call MPI_IRECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,request,ierr)
      else
         if (present(status)) then
            call MPI_RECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,status,ierr)
         else
            call MPI_RECV(buf,count,mpi_type,mpi_source,mpi_tag,mpi_comm,FMPI_STATUS_IGNORE,ierr)
         end if
      end if
      call f_timer_resume()
      if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)RECV occured',&
           err_id=ERR_MPI_WRAPPERS)

    end subroutine mpirecv_d0

    subroutine mpirecv_gpu(buf,count,source,tag,comm,status,request,simulate,verbose,type)
      use yaml_output
      use iso_c_binding
      implicit none
      type(c_ptr) :: buf !fake intent(in)
      real(f_double),pointer:: a !fake intent(out)
      integer, intent(in) :: count
      integer, intent(in), optional :: source
      integer, intent(in), optional :: tag
      integer, intent(in), optional :: comm
      integer, intent(out), optional :: request !<toggle the isend operation
      integer, dimension(FMPI_STATUS_SIZE), intent(out), optional :: status !<for the blocking operation
      logical, intent(in), optional :: simulate,verbose
      integer, intent(in) :: type
      !local variables
      logical :: verb,sim
      integer :: mpi_comm,ierr,mpi_source,mpi_tag,mpi_type,tcat

      mpi_comm=fmpi_comm(comm)
      mpi_source=FMPI_ANY_SOURCE
      mpi_tag=FMPI_ANY_TAG
      if (present(source)) then
         mpi_source=source
         mpi_tag=source
      end if
      if (present(tag)) mpi_tag=tag
      verb=.false.
      if (present(verbose)) verb=verbose .and. mpi_source /= mpirank_null()
      if (verb) call yaml_comment('Receiving'//count//'elements from'//source//'in'//mpirank(mpi_comm))

      if (verb) then
         call yaml_mapping_open('MPI_(I)RECV')
         call yaml_map('Elements',count)
         call yaml_map('Source',source)
         call yaml_map('Dest',mpirank(mpi_comm))
         call yaml_map('Tag',mpi_tag)
         call yaml_mapping_close()
      end if

      sim=.false.
      if (present(simulate)) sim=simulate
      if (sim) return

      tcat=TCAT_RECV
      ! Synchronize the communication
      call f_timer_interrupt(tcat)

      mpi_type=type
      call c_f_pointer(buf, a)
      if (present(request)) then
         call MPI_IRECV(a,count,mpi_type,mpi_source,mpi_tag,mpi_comm,request,ierr)
      else
         if (present(status)) then
            call MPI_RECV(a,count,mpi_type,mpi_source,mpi_tag,mpi_comm,status,ierr)
         else
            call MPI_RECV(a,count,mpi_type,mpi_source,mpi_tag,mpi_comm,FMPI_STATUS_IGNORE,ierr)
         end if
      end if
      call f_timer_resume()
      if (ierr/=0) call f_err_throw('An error in calling to MPI_(I)RECV occured',&
           err_id=ERR_MPI_WRAPPERS)

    end subroutine mpirecv_gpu

end module f_sendrecv
