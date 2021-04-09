!> @file
!!  Test of the overlap point to point (mpi)
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to test the mpi behaviour
program mpi_check
  use wrapper_MPI
  use yaml_output
  use yaml_strings
  use f_utils
  use f_enums
  use test_mpi_wrappers
  use dynamic_memory
  implicit none

  logical :: failed
  integer :: ntot,ierr,i,iproc,nproc,nspin,nother,ibuf
  integer, dimension(:,:), allocatable  :: isizes
  real(kind=8), dimension(:), allocatable :: copybuffer,buffer
  logical(f_byte), dimension(:,:,:), allocatable :: lbuf

  call f_lib_initialize()
  call mpiinit()
  iproc=mpirank()
  nproc=mpisize()

  ntot=107
  nspin=2
  nother=10

  allocate(buffer(ntot),copybuffer(ntot),isizes(nspin,nother))

  lbuf = f_malloc([ntot,1,1],id='ldata')

  do ibuf=1,ntot
    lbuf(ibuf,1,1) = modulo(ibuf,nproc)==iproc
  end do
  buffer=1.0d0
  copybuffer=0.d0
  isizes=1
  !last component indicates the size
  isizes(nspin,nother)=ntot

  !reference (supposed)
  call MPI_ALLREDUCE(buffer,copybuffer,ntot,&
       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)

!!$  call MPI_ALLREDUCE(MPI_IN_PLACE,buffer,ntot,&
!!$       MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)
  call fmpi_allreduce(buffer(1),isizes(nspin,nother),FMPI_SUM)

  call fmpi_allreduce(lbuf,op=FMPI_BOR)

  call f_assert(all(lbuf),'Reduction for the byte array did not went through')

  failed=.false.
  do i=1,ntot
     if (abs(copybuffer(i)-buffer(i))>1.d-12 .or. abs(copybuffer(i)-real(nproc,kind=8))> 1.d-12) then
        write(*,'(a,i7,2(1pe25.17),i6,i6)')'Failed',i,buffer(i),copybuffer(i),iproc,nproc
        failed=.true.
     end if
  end do
  call f_assert(.not. failed,'At least one proc ('+iproc+') did not reduce correctly')

!!$  call f_open_file_mpi(unit,'buffer_file.dat',comm=mpicomm())
!!$  call MPI_FILE_SEEK(unit, offset+jd*kd*ld*nq*4, MPI_SEEK_CUR, ierr)
!!$  CALL MPI_FILE_SET_VIEW
!!$
!!$  call MPI_FILE_CLOSE(unit,ierr)
!!$

  deallocate(buffer,copybuffer,isizes)

  call f_free(lbuf)

  call test_mpi_alltoallv(iproc, nproc,fmpi_comm(),ntot,1)

  !wait all processes before finalisation
  call fmpi_barrier()

  !verify the mpiallgather
  call test_mpi_allgather(ntot,fmpi_comm())

  if (iproc==0) call yaml_map('MPI Test succeeded',.true.)
  call MPI_FINALIZE(ierr)

  call f_lib_finalize_noreport()

contains

subroutine f_open_file_mpi(unit,file,status,position,action,binary,comm)
  use yaml_strings, only: f_strcpy
  use dictionaries, only: f_err_throw
  implicit none
  !> integer of the unit. On entry, it indicates the
  !! suggested unit number. On exit, it indicates the free unit
  !! which has been used for the file opening
  integer, intent(inout) :: unit
  !> filename
  character(len=*), intent(in) :: file
  !> status
  character(len=*), intent(in), optional :: status
  !> position
  character(len=*), intent(in), optional :: position
  !> action
  character(len=*), intent(in), optional :: action
  !> if true, the file will be opened in the unformatted i/o
  !! if false or absent, the file will be opened for formatted i/o
  logical, intent(in), optional :: binary
  !> if present, the file is opened at the level of a mpi_communicator
  integer, intent(in), optional :: comm
  !local variables
  logical :: onsite
  integer :: ierr,amode,iact,ist,iacc,ipos

  if (.not. present(comm)) then
     onsite=.true.
  else
     onsite=comm == MPI_COMM_SELF
  end if

  !still not clear what to do in the MPI_COMM_SELF case
  if (onsite) then
     call f_open_file(unit,file,status,position,action,binary)
  else
     !we should check if the filename is common to all the processes
     iact=MPI_MODE_RDWR
     if (present(action)) then
        select case(trim(action))
        case('read','READ')
           iact=MPI_MODE_RDONLY
        case('write','WRITE')
           iact=MPI_MODE_WRONLY
        case('readwrite','READWRITE')
           iact=MPI_MODE_RDWR
        case default
           call f_err_throw('unrecognized action='+action)
        end select
     end if
     ist=MPI_MODE_CREATE
     if (present(status)) then
        select case(trim(status))
        case('old','OLD','unknown','UNKNOWN','replace','REPLACE')
           ist=MPI_MODE_CREATE
        case('new','NEW')
           ist=MPI_MODE_EXCL
        case('scratch','SCRATCH')
           ist=MPI_MODE_DELETE_ON_CLOSE
        case default
           call f_err_throw('unrecognized status='+action)
        end select
     end if
     iacc=MPI_MODE_SEQUENTIAL !direct not yet in the api
     ipos=0
     if (present(position)) then
        select case(trim(position))
        case('asis','ASIS','rewind','REWIND')
        case('append','APPEND')
           ipos=MPI_MODE_APPEND
        end select
     end if
     amode=iact+ist+iacc+ipos
     call mpi_file_open(comm,file,amode,MPI_INFO_NULL,unit,ierr)
     if (ierr /= MPI_SUCCESS) then
        call f_err_throw('Error in opening mpi file='//&
             trim(file)//', iostat='//trim(yaml_toa(ierr)),&
             err_name='INPUT_OUTPUT_ERROR')
     end if
  end if
end subroutine f_open_file_mpi

end program mpi_check
