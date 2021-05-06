!> @file
!!  Module for the test OP2P (overlap point to point)
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module for the test OP2P (overlap point to point)
module mpi_layer
  use module_base
  implicit none
  private

  public :: send_mpi_profile,receive_mpi_profile,wait_mpi_profile

contains

  !> fake receiving of the arrays
  subroutine receive_mpi_profile(istep,isendproc,irecvproc,ncount,itag,irequest,recvbuf)
    implicit none
    integer, intent(in) :: istep,isendproc,irecvproc,ncount,itag,irequest
    real(kind=8), intent(in) :: recvbuf
    !local variables
    integer :: ierr

    !here we can add something to trap the IRECV call
    !print '(3(a,i4),i4)','NON_BLOCKING RECV, from',isendproc,' to',irecvproc,', step, elems:',istep,ncount

    call MPI_IRECV(recvbuf,ncount,MPI_DOUBLE_PRECISION,isendproc,&
         itag,MPI_COMM_WORLD,irequest,ierr)

    !output error signal
    if (ierr /=0) then
       write(*,*)'ERROR in IRECV, iproc, istep',irecvproc,istep
    end if

  end subroutine receive_mpi_profile

  !> fake sending of the arrays
  subroutine send_mpi_profile(istep,isendproc,irecvproc,ncount,itag,irequest,sendbuf)
    implicit none
    integer, intent(in) :: istep,isendproc,irecvproc,ncount,itag,irequest
    real(kind=8), intent(in) :: sendbuf
    !local variables
    integer :: ierr

    !here we can add something to trap the ISEND call
    !print '(3(a,i4),i4)','NON_BLOCKING SEND, from',isendproc,' to',irecvproc,', step, elems:',istep,ncount

    call MPI_ISEND(sendbuf,ncount,MPI_DOUBLE_PRECISION,irecvproc,&
         itag,MPI_COMM_WORLD,irequest,ierr)

    !output error signal
    if (ierr /=0) then
       write(*,*)'ERROR in ISEND, iproc, istep',irecvproc,istep
    end if
  end subroutine send_mpi_profile

  subroutine wait_mpi_profile(iproc,istep,nreq,requests)
    implicit none
    integer, intent(in) :: iproc,istep,nreq
    integer, dimension(nreq), intent(in) :: requests
    !local variables
    logical :: error_found
    integer :: jproc,isr
    integer, dimension(MPI_STATUS_SIZE,4) :: mpistat
    !local variables
    integer :: ierr

    !verify that the messages have been passed
    call MPI_WAITALL(nreq,requests,mpistat,ierr)
    if (ierr /=0)  then
       write(*,*),'ERROR WAITALL, iproc,step,ierr:',iproc,istep,ierr,mpistat,MPI_STATUSES_IGNORE
    end if
  end subroutine wait_mpi_profile

end module mpi_layer
