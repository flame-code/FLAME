!> @file
!! Include fortran file for allgather operations
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  integer, intent(in), optional :: sendcount
  integer, intent(in), optional :: recvcount
  integer, dimension(:), intent(in), optional :: recvcounts
  integer, dimension(:), intent(in), target, optional :: displs
  integer, intent(in), optional :: comm
  !local variables
  logical :: allv,ldispl,in_place
  integer :: ierr,mpi_comm,me,nprc,sendsize,recvsize,ntot,jprc
  integer :: recvcnt,sendcnt,dspl,ourselves
  integer, dimension(:), pointer :: displs_
  external :: MPI_ALLGATHERV,MPI_ALLGATHER

  !check of the arguments
  sendsize=kind(sendbuf)
  recvsize=sendsize
  if (present(recvbuf)) recvsize=kind(recvbuf)

  if (present(comm)) then
     mpi_comm=comm
  else
     mpi_comm=MPI_COMM_WORLD !or bigdft_mpi%mpi_comm?
  end if

  me=mpirank(mpi_comm)
  ourselves=mpisize(mpi_comm)

  !identify recvcnt
  !decide also if it is a allgatherv or not
  allv=present(recvcounts)
  if (allv) then
     if (size(recvcounts) /= ourselves) &
          call f_err_throw('Size of recvcounts is not of length group size',&
          err_name='ERR_MPI_WRAPPERS')
     recvcnt=recvcounts(me+1)
     allv=any(recvcounts/=recvcnt)
     if (.not. allv .and. present(displs)) then
        nprc=size(displs)
        if (nprc /= size(recvcounts)) &
             call f_err_throw('Error in mpiallgather, sizes of displs and recvcounts do not coincide ('//&
             trim(yaml_toa(nprc))//' /= '//&
             trim(yaml_toa(size(recvcounts)))//')',&
             err_name='ERR_MPI_WRAPPERS')
        dspl=0
        do jprc=1,nprc
           if (dspl /= displs(jprc)) then
              allv=.true.
              exit
           end if
           dspl=dspl+recvcnt
        end do
     end if
  else if (present(recvcount)) then
     recvcnt=recvcount
  else if (present(sendcount)) then
     recvcnt=sendcount
  else
     call f_err_throw('mpiallgather: Unable to determine recvcount. '//&
          'one of recvcounts, recvcount or sendcount should be present',&
          err_name='ERR_MPI_WRAPPERS')
  end if

  ldispl=allv .and. .not. present(displs)

  !determine sendcnt
  if (present(sendcount)) then
     sendcnt=sendcount
  else if (present(recvcount)) then
     sendcnt=recvcount
  else
     sendcnt=recvcounts(me+1)
  end if

  !determine if in_place or not
  in_place=have_mpi2 .and. .not. present(recvbuf)

  !in the allv case, determine the displacements in the case they are not given
  if (ldispl) then
     nprc=size(recvcounts)
     displs_=f_malloc_ptr(nprc,id='displs')
     dspl=0
     do jprc=1,nprc
        displs_(jprc)=dspl
        dspl=dspl+recvcounts(jprc)
     end do
  else if (present(displs)) then
     displs_ => displs
  else
     nullify(displs_)
  end if

  !check the validity of the buffers
  if (int(sendcnt,kind=8)*sendsize /= int(recvcnt,kind=8)*recvsize) &
       call f_err_throw(&
       'allgatherv error, send and receive sizes are not consistent ('&
       //trim(yaml_toa([sendcnt,sendsize,recvcnt,recvsize]))//')',&
       err_name='ERR_MPI_WRAPPERS')

  !then launch the communication
  if (in_place) then
     if (allv) then
        call f_timer_interrupt(TCAT_ALLGATHERV)
        !case with MPI_IN_PLACE
        call MPI_ALLGATHERV(MPI_IN_PLACE,sendcnt,mpitype(sendbuf),&
             sendbuf,recvcounts,displs_,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
     else
        call f_timer_interrupt(TCAT_ALLGATHER)
        !case with MPI_IN_PLACE
        call MPI_ALLGATHER(MPI_IN_PLACE,sendcnt,mpitype(sendbuf),&
             sendbuf,recvcnt,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
     end if
  else if (present(recvbuf)) then
     if (allv) then
        call f_timer_interrupt(TCAT_ALLGATHERV)
        call MPI_ALLGATHERV(sendbuf,sendcnt,mpitype(sendbuf),&
             recvbuf,recvcounts,displs_,mpitype(recvbuf),mpi_comm,ierr)
        call f_timer_resume()
     else
        call f_timer_interrupt(TCAT_ALLGATHER)
        call MPI_ALLGATHER(sendbuf,sendcnt,mpitype(sendbuf),&
             recvbuf,recvcnt,mpitype(recvbuf),mpi_comm,ierr)
        call f_timer_resume()
     end if
  else
     !here the receive buffer is not given, therefore 
     !we are in the case of a mpi1 implementation
     !not supporting MPI_IN_PLACE
     !therefore a temporary buffer has to be created at the place of 
     !the sendbuffer
     !its size must be bigger as the sendbuf
     if (present(recvcounts)) then
        ntot=sum(recvcounts)
     else
        ntot=recvcnt*ourselves
     end if
     copybuf = f_malloc(ntot,id='copybuf')
     call f_memcpy(n=ntot,src=sendbuf,dest=copybuf(1))
     if (allv) then
        call f_timer_interrupt(TCAT_ALLGATHERV)
        call MPI_ALLGATHERV(copybuf(1+displs_(me+1)),sendcnt,mpitype(sendbuf),&
             sendbuf,recvcounts,displs_,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
     else
        call f_timer_interrupt(TCAT_ALLGATHER)
        call MPI_ALLGATHER(copybuf(1+displs_(me)),sendcnt,mpitype(sendbuf),&
             sendbuf,recvcnt,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
     end if
     call f_free(copybuf)
  end if

  if (ldispl) call f_free_ptr(displs_)
  !otherwise nullify the pointer
  nullify(displs_)

  !retrieve the error if it is still possible
  if (ierr /=0) &
       call f_err_throw('Error in the exectution of gather operation',&
       err_name='ERR_MPI_WRAPPERS')

!!$
!!$#ifdef HAVE_MPI2
!!$    call f_timer_interrupt(TCAT_ALLGATHERV)
!!$    !case with MPI_IN_PLACE
!!$    call MPI_ALLGATHERV(MPI_IN_PLACE,counts(me),mpitype(buffer),&
!!$         buffer,counts,displs,mpitype(buffer),mpi_comm,ierr)
!!$    call f_timer_resume()
!!$#else
!!$    !local variables
!!$    real(kind=8), dimension(:), allocatable :: copybuf
!!$
!!$    !Here we have a performance penalty by copying all buffer, instead of
!!$    !just the send part, but I don't see how to get buffer(displs(me))
!!$    copybuf = f_malloc(sum(counts),id='copybuf')
!!$
!!$    call dcopy(sum(counts),buffer,1,copybuf,1) 
!!$    ierr=0 !put just for MPIfake compatibility
!!$    call f_timer_interrupt(TCAT_ALLGATHERV)
!!$    call MPI_ALLGATHERV(copybuf(1+displs(me+1)),counts(me),mpitype(buffer),&
!!$         buffer,counts,displs,mpitype(buffer),mpi_comm,ierr)
!!$    call f_timer_resume()
!!$    call f_free(copybuf)
!!$#endif
!!$
!!$    if (ierr /=0) call f_err_throw('Error in mpi_allgatherv',&
!!$         err_id=ERR_MPI_WRAPPERS)
!!$    !stop 'MPIALLGATHERV_DBL'
