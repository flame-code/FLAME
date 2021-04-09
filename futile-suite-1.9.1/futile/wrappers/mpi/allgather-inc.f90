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
  type(fmpi_win), intent(out), optional :: win
  type(f_enumerator), intent(in), optional :: algorithm
  !local variables
  logical :: in_place
  integer :: recvcnt,sendcnt,algo,jproc,iproc,nproc
  integer :: recvsize,sendsize,ierr,mpi_comm
  integer(f_long) :: buf_size,ntot,offset
  type(fmpi_counts) :: ctr
  type(fmpi_win) :: window_


  !check of the arguments
  sendsize=kind(sendbuf)
  recvsize=sendsize
  if (present(recvbuf)) recvsize=kind(recvbuf)
  mpi_comm=fmpi_comm(comm)
  iproc=mpirank(mpi_comm)
  nproc=mpisize(mpi_comm)

  if (present(recvbuf)) then
     buf_size=f_size(sendbuf)
  else
     buf_size=f_size(sendbuf)/nproc !assuming that the sendbuf is not malformed
  end if

  ctr=get_control_buffers(buf_size,recvcount,recvcounts,displs,comm)
  call get_srcounts(iproc,buf_size,sendcnt,recvcnt,sendcount,recvcount,recvcounts)

  !check the validity of the buffers
  if (int(sendcnt,f_long)*sendsize /= int(recvcnt,f_long)*recvsize) &
       call f_err_throw(&
       'allgatherv error, send and receive sizes are not consistent ('&
       //trim(yaml_toa([sendcnt,sendsize,recvcnt,recvsize]))//')',&
       err_name='ERR_MPI_WRAPPERS')

  !determine if in_place or not
  in_place=have_mpi2 .and. .not. present(recvbuf)

  !determine the algorithm
  algo=AUTOMATIC_ALGO
  if (present(algorithm)) algo=toi(algorithm)
  if (algo==AUTOMATIC_ALGO) then
     if (ctr%count==-1) then
        algo=VARIABLE_ALGO
     else
        algo=NOT_VARIABLE_ALGO
     end if
  end if
  if (present(win)) algo=VARIABLE_ONE_SIDED_GET_ALGO

  select case(algo)
  case(NOT_VARIABLE_ALGO)
     if (in_place) then
        call f_timer_interrupt(TCAT_ALLGATHER)
        !case with MPI_IN_PLACE
        call MPI_ALLGATHER(FMPI_IN_PLACE,sendcnt,mpitype(sendbuf),&
             sendbuf,recvcnt,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
     else if (present(recvbuf)) then
        call f_timer_interrupt(TCAT_ALLGATHER)
        call MPI_ALLGATHER(sendbuf,sendcnt,mpitype(sendbuf),&
             recvbuf,recvcnt,mpitype(recvbuf),mpi_comm,ierr)
        call f_timer_resume()
     else
        !here the receive buffer is not given, therefore 
        !we are in the case of a mpi1 implementation
        !not supporting MPI_IN_PLACE
        !therefore a temporary buffer has to be created at the place of 
        !the sendbuffer
        !its size must be bigger as the sendbuf
        call get_ntot_and_offset(iproc,nproc,ctr,ntot,offset)
        copybuf = f_malloc(ntot,id='copybuf')
        call f_memcpy(n=int(ntot,f_integer),src=sendbuf,dest=copybuf(1))
        call f_timer_interrupt(TCAT_ALLGATHER)
        call MPI_ALLGATHER(copybuf(1+offset),sendcnt,mpitype(sendbuf),&
             sendbuf,recvcnt,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
        call f_free(copybuf)
     end if
  case(VARIABLE_ALGO)
     if (in_place) then
        call f_timer_interrupt(TCAT_ALLGATHERV)
        !case with MPI_IN_PLACE
        call MPI_ALLGATHERV(FMPI_IN_PLACE,sendcnt,mpitype(sendbuf),&
             sendbuf,ctr%counts,ctr%displs,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
     else if (present(recvbuf)) then
        call f_timer_interrupt(TCAT_ALLGATHERV)
        call MPI_ALLGATHERV(sendbuf,sendcnt,mpitype(sendbuf),&
             recvbuf,ctr%counts,ctr%displs,mpitype(recvbuf),mpi_comm,ierr)
        call f_timer_resume()
     else
        !here the receive buffer is not given, therefore 
        !we are in the case of a mpi1 implementation
        !not supporting MPI_IN_PLACE
        !therefore a temporary buffer has to be created at the place of 
        !the sendbuffer
        !its size must be as large as the sendbuf
        call get_ntot_and_offset(iproc,nproc,ctr,ntot,offset)
        copybuf = f_malloc(ntot,id='copybuf')
        call f_memcpy(n=int(ntot,f_integer),src=sendbuf,dest=copybuf(1))
        call f_timer_interrupt(TCAT_ALLGATHERV)
        call MPI_ALLGATHERV(copybuf(1+offset),sendcnt,mpitype(sendbuf),&
             sendbuf,ctr%counts,ctr%displs,mpitype(sendbuf),mpi_comm,ierr)
        call f_timer_resume()
        call f_free(copybuf)
     end if
  case(VARIABLE_ONE_SIDED_GET_ALGO)
     !here the arguments have to be passed and are used natively
     call fmpi_win_create(window_,sendbuf,sendcnt,comm=comm)
     call fmpi_win_fence(window_,FMPI_WIN_OPEN)
     if (ctr%count == -1) then
        do jproc=0,nproc-1
           if (ctr%counts(jproc) > 0) &
                call fmpi_get(recvbuf,jproc,window_,ctr%counts(jproc),&
                origin_displ=ctr%displs(jproc))
        end do
     else if (ctr%count > 0) then
        do jproc=0,nproc-1
           call fmpi_get(recvbuf,jproc,window_,ctr%count,&
                origin_displ=jproc*ctr%count)
        end do
     end if
     if (.not. present(win)) then
        call fmpi_win_shut(window_)
     else
        win=window_
     end if
     !if we are here the mpi_get has been lauched correctly
     ierr = FMPI_SUCCESS
  end select

  call free_fmpi_counts(ctr)

!!$  logical :: allv,ldispl,in_place
!!$  integer :: ierr,mpi_comm,me,nprc,sendsize,recvsize,ntot,jprc
!!$  integer :: recvcnt,sendcnt,dspl,ourselves,jproc,algo
!!$  integer, dimension(:), pointer :: displs_
!!$  type(fmpi_win) :: window_
!!$  external :: MPI_ALLGATHERV,MPI_ALLGATHER
!!$
!!$  !check of the arguments
!!$  sendsize=kind(sendbuf)
!!$  recvsize=sendsize
!!$  if (present(recvbuf)) recvsize=kind(recvbuf)
!!$
!!$  mpi_comm=fmpi_comm(comm)
!!$  me=mpirank(mpi_comm)
!!$  ourselves=mpisize(mpi_comm)
!!$
!!$  !identify recvcnt
!!$  !decide also if it is a allgatherv or not
!!$  allv=present(recvcounts)
!!$  if (allv) then
!!$     if (size(recvcounts) /= ourselves) &
!!$          call f_err_throw('Size of recvcounts is not of length group size',&
!!$          err_name='ERR_MPI_WRAPPERS')
!!$     recvcnt=recvcounts(me+1)
!!$     allv=any(recvcounts/=recvcnt)
!!$     if (.not. allv .and. present(displs)) then
!!$        nprc=size(displs)
!!$        if (nprc /= size(recvcounts)) &
!!$             call f_err_throw('Error in mpiallgather, sizes of displs and recvcounts do not coincide ('//&
!!$             trim(yaml_toa(nprc))//' /= '//&
!!$             trim(yaml_toa(size(recvcounts)))//')',&
!!$             err_name='ERR_MPI_WRAPPERS')
!!$        dspl=0
!!$        do jprc=1,nprc
!!$           if (dspl /= displs(jprc)) then
!!$              allv=.true.
!!$              exit
!!$           end if
!!$           dspl=dspl+recvcnt
!!$        end do
!!$     end if
!!$  else if (present(recvcount)) then
!!$     recvcnt=recvcount
!!$  else if (present(sendcount)) then
!!$     recvcnt=sendcount
!!$  else
!!$     call f_err_throw('mpiallgather: Unable to determine recvcount. '//&
!!$          'one of recvcounts, recvcount or sendcount should be present',&
!!$          err_name='ERR_MPI_WRAPPERS')
!!$  end if
!!$
!!$  ldispl=allv .and. .not. present(displs)
!!$
!!$  !determine sendcnt
!!$  if (present(sendcount)) then
!!$     sendcnt=sendcount
!!$  else if (present(recvcount)) then
!!$     sendcnt=recvcount
!!$  else
!!$     sendcnt=recvcounts(me+1)
!!$  end if
!!$
!!$  !determine if in_place or not
!!$  in_place=have_mpi2 .and. .not. present(recvbuf)
!!$
!!$  !in the allv case, determine the displacements in the case they are not given
!!$  if (ldispl) then
!!$     nprc=size(recvcounts)
!!$     displs_=f_malloc_ptr(nprc,id='displs')
!!$     dspl=0
!!$     do jprc=1,nprc
!!$        displs_(jprc)=dspl
!!$        dspl=dspl+recvcounts(jprc)
!!$     end do
!!$  else if (present(displs)) then
!!$     displs_ => displs
!!$  else
!!$     nullify(displs_)
!!$  end if
!!$
!!$  !check the validity of the buffers
!!$  if (int(sendcnt,f_long)*sendsize /= int(recvcnt,f_long)*recvsize) &
!!$       call f_err_throw(&
!!$       'allgatherv error, send and receive sizes are not consistent ('&
!!$       //trim(yaml_toa([sendcnt,sendsize,recvcnt,recvsize]))//')',&
!!$       err_name='ERR_MPI_WRAPPERS')
!!$
!!$  !then launch the communication
!!$  if (in_place) then
!!$     if (allv) then
!!$        call f_timer_interrupt(TCAT_ALLGATHERV)
!!$        !case with MPI_IN_PLACE
!!$        call MPI_ALLGATHERV(FMPI_IN_PLACE,sendcnt,mpitype(sendbuf),&
!!$             sendbuf,recvcounts,displs_,mpitype(sendbuf),mpi_comm,ierr)
!!$        call f_timer_resume()
!!$     else
!!$        call f_timer_interrupt(TCAT_ALLGATHER)
!!$        !case with MPI_IN_PLACE
!!$        call MPI_ALLGATHER(FMPI_IN_PLACE,sendcnt,mpitype(sendbuf),&
!!$             sendbuf,recvcnt,mpitype(sendbuf),mpi_comm,ierr)
!!$        call f_timer_resume()
!!$     end if
!!$  else if (present(recvbuf)) then
!!$     if (allv) then
!!$        call f_timer_interrupt(TCAT_ALLGATHERV)
!!$        call MPI_ALLGATHERV(sendbuf,sendcnt,mpitype(sendbuf),&
!!$             recvbuf,recvcounts,displs_,mpitype(recvbuf),mpi_comm,ierr)
!!$        call f_timer_resume()
!!$     else
!!$        call f_timer_interrupt(TCAT_ALLGATHER)
!!$        call MPI_ALLGATHER(sendbuf,sendcnt,mpitype(sendbuf),&
!!$             recvbuf,recvcnt,mpitype(recvbuf),mpi_comm,ierr)
!!$        call f_timer_resume()
!!$     end if
!!$  else
!!$     !here the receive buffer is not given, therefore 
!!$     !we are in the case of a mpi1 implementation
!!$     !not supporting MPI_IN_PLACE
!!$     !therefore a temporary buffer has to be created at the place of 
!!$     !the sendbuffer
!!$     !its size must be bigger as the sendbuf
!!$     if (present(recvcounts)) then
!!$        ntot=sum(recvcounts)
!!$     else
!!$        ntot=recvcnt*ourselves
!!$     end if
!!$     copybuf = f_malloc(ntot,id='copybuf')
!!$     call f_memcpy(n=ntot,src=sendbuf,dest=copybuf(1))
!!$     if (allv) then
!!$        call f_timer_interrupt(TCAT_ALLGATHERV)
!!$        call MPI_ALLGATHERV(copybuf(1+displs_(me+1)),sendcnt,mpitype(sendbuf),&
!!$             sendbuf,recvcounts,displs_,mpitype(sendbuf),mpi_comm,ierr)
!!$        call f_timer_resume()
!!$     else
!!$        call f_timer_interrupt(TCAT_ALLGATHER)
!!$        call MPI_ALLGATHER(copybuf(1+displs_(me)),sendcnt,mpitype(sendbuf),&
!!$             sendbuf,recvcnt,mpitype(sendbuf),mpi_comm,ierr)
!!$        call f_timer_resume()
!!$     end if
!!$     call f_free(copybuf)
!!$  end if
!!$
!!$  if (ldispl) call f_free_ptr(displs_)
!!$  !otherwise nullify the pointer
!!$  nullify(displs_)

  !retrieve the error if it is still possible
  if (ierr /=FMPI_SUCCESS) &
       call f_err_throw('Error in the exectution of gather operation',&
       err_name='ERR_MPI_WRAPPERS')
