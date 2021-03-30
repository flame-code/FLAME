!> @file
!!    Fake functions for MPI in the case of serial version
!! @author
!!    Copyright (C) 2007-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


subroutine  MPI_INIT(ierr)
  implicit none
  integer, intent(out) :: ierr
  ierr=0
END SUBROUTINE MPI_INIT

subroutine MPI_INITIALIZED(init,ierr)
  implicit none
  logical, intent(out) :: init
  integer, intent(out) :: ierr
  init=.false.
  ierr=0
END SUBROUTINE  MPI_INITIALIZED

subroutine  MPI_COMM_CREATE(MPI_COMM_WORLD,MPI_GROUP,MPI_COMM,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: MPI_GROUP,MPI_COMM,ierr
  MPI_GROUP=1
  MPI_COMM=1
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_CREATE

!>Used to have the mpimaxtag
!! if COMM_KEYVAL=1, gives success and 1 for attribute_wal
subroutine  MPI_COMM_GET_ATTR(COMM,COMM_KEYVAL,ATTRIBUTE_VAL,FLAG,ierr)
  implicit none
  integer, intent(in) :: COMM, COMM_KEYVAL
  logical, intent(out) :: FLAG
  !integer(kind=MPI_ADDRESS_KIND), intent(out) :: ATTRIBUTE_VAL
  integer(kind=8), intent(out) :: ATTRIBUTE_VAL
  integer, intent(out) :: ierr
  ATTRIBUTE_VAL=1
  ierr=COMM*COMM_KEYVAL*0
  FLAG=(COMM_KEYVAL==1)
END SUBROUTINE MPI_COMM_GET_ATTR

subroutine  MPI_COMM_GROUP(MPI_COMM_WORLD,MPI_GROUP,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: MPI_GROUP,ierr
  MPI_GROUP=1
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_GROUP

subroutine  MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: iproc,ierr
  iproc=0
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_RANK

subroutine  MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: nproc,ierr
  nproc=1
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_COMM_SIZE

subroutine  MPI_GROUP_INCL(GROUP,N,NRANKS,NEWGROUP,ierr)
  implicit none
  integer, intent(in) :: GROUP,N
  integer, intent(in) :: NRANKS(N)
  integer, intent(out) :: NEWGROUP,ierr
  NEWGROUP=size(NRANKS)
  ierr=GROUP*0
END SUBROUTINE MPI_GROUP_INCL

subroutine  MPI_GROUP_FREE(GROUP,ierr)
  implicit none
  integer, intent(in) :: GROUP
  integer, intent(out) :: ierr
  ierr=GROUP*0
END SUBROUTINE MPI_GROUP_FREE

subroutine mpi_test(request,flag,MPI_Status)
  implicit none
  integer, intent(in) :: request
  integer, intent(out) :: flag
  integer, intent(out) :: MPI_Status
  flag = 1 + 0*request
  MPI_Status = 1
end subroutine mpi_test

subroutine mpi_wait(request,MPI_Status)
  implicit none
  integer, intent(in) :: request
  integer, intent(out) :: MPI_Status
  MPI_Status = 1 + 0*request
end subroutine mpi_wait

subroutine mpi_file_close()
stop 'mpi_file_close'
end

subroutine mpi_file_open()
stop 'mpi_file_open'
end
subroutine mpi_file_read()
stop 'mpi_file_read'
end
subroutine mpi_file_write()
stop 'mpi_file_write'
end



subroutine mpi_file_set_view()
stop 'mpi_file_set_view'
end


!here we have routines which do not transform the argument for nproc==1
!these routines can be safely called also in the serial version
subroutine  MPI_FINALIZE(ierr)
  implicit none
  integer, intent(out) :: ierr
  ierr=0
END SUBROUTINE MPI_FINALIZE

subroutine MPI_BCAST(buffer,n,mpitype,root,comm,ierr)
  implicit none
  integer :: buffer,n,mpitype,root,comm,ierr
  ierr=0
END SUBROUTINE MPI_BCAST

subroutine  MPI_BARRIER(MPI_COMM_WORLD,ierr)
  implicit none
  integer, intent(in) :: MPI_COMM_WORLD
  integer, intent(out) :: ierr
  ierr=MPI_COMM_WORLD*0
END SUBROUTINE MPI_BARRIER


!> These routines in serial version should not be called.
!! A stop is added when necessary, otherwise for copying routines, the corresponding copy
!! is implemented whenever possible
subroutine MPI_REDUCE()
  implicit none
  stop 'MPIFAKE: REDUCE'
END SUBROUTINE MPI_REDUCE

subroutine  MPI_ALLREDUCE(i1,i2,i3,i4,op,mpi_comm,ierr)
  implicit none
  integer :: i1,i2,i3,i4
  integer, intent(in) :: op, mpi_comm
  integer, intent(out) :: ierr
  !Initialize ierr in case when MPI_ALLREDUCE is called
  ierr = 0
  !stop 'MPIFAKE: ALLREDUCE' eliminated due to ABINIT module
END SUBROUTINE MPI_ALLREDUCE

subroutine  MPI_IALLREDUCE(i1,i2,i3,i4,op,mpi_comm,ierr)
  implicit none
  integer :: i1,i2,i3,i4
  integer, intent(in) :: op, mpi_comm
  integer, intent(out) :: ierr
  !Initialize ierr in case when MPI_ALLREDUCE is called
  ierr = 0
  stop 'MPIFAKE: IALLREDUCE'
END SUBROUTINE MPI_IALLREDUCE

subroutine  MPI_ALLGatherV()
  implicit none
  stop 'MPIFAKE: ALLGATHERV'
END SUBROUTINE  MPI_ALLGatherV

subroutine  MPI_ALLGATHER()
  implicit none
  stop 'MPIFAKE: ALLGATHER'
END SUBROUTINE  MPI_ALLGATHER

subroutine  MPI_GatherV()
  implicit none
  stop 'MPIFAKE: GATHERV'
END SUBROUTINE  MPI_GatherV

subroutine  MPI_Gather()
  implicit none
  stop 'MPIFAKE: GATHER'
END SUBROUTINE  MPI_Gather


subroutine  MPI_ALLTOALL()
  implicit none
  stop 'MPIFAKE: ALLTOALL'
END SUBROUTINE  MPI_ALLTOALL

subroutine  MPI_ALLTOALLV()
  implicit none
  stop 'MPIFAKE: ALLTOALLV'
END SUBROUTINE  MPI_ALLTOALLV

subroutine  MPI_IALLTOALLV()
  implicit none
  stop 'MPIFAKE: IALLTOALLV'
END SUBROUTINE  MPI_IALLTOALLV

subroutine  MPI_REDUCE_SCATTER()
  implicit none
  stop 'MPIFAKE: REDUCE_SCATTER'
END SUBROUTINE  MPI_REDUCE_SCATTER

subroutine  MPI_ABORT()
  implicit none
  stop ' MPIFAKE: MPI_ABORT'
END SUBROUTINE  MPI_ABORT

subroutine  MPI_IRECV()
  implicit none
  stop 'MPIFAKE: IRECV'
END SUBROUTINE  MPI_IRECV

subroutine  MPI_RECV()
  implicit none
  stop 'MPIFAKE: RECV'
END SUBROUTINE  MPI_RECV

subroutine  MPI_ISEND()
  implicit none
  stop 'MPIFAKE: ISEND'
END SUBROUTINE  MPI_ISEND

subroutine  MPI_SEND()
  implicit none
  stop 'MPIFAKE: SEND'
END SUBROUTINE  MPI_SEND

subroutine  MPI_WAITALL()
  implicit none
  stop 'MPIFAKE: WAITALL'
END SUBROUTINE  MPI_WAITALL

subroutine MPI_GET_PROCESSOR_NAME(nodename_local,namelen,ierr)
  implicit none
  integer, intent(out) :: namelen,ierr
  character(len=*), intent(inout) :: nodename_local
  ierr=0
  namelen=9
  nodename_local(1:9)='localhost'
  nodename_local(10:len(nodename_local))=repeat(' ',max(len(nodename_local)-10,0))
END SUBROUTINE  MPI_GET_PROCESSOR_NAME

subroutine  mpi_error_string()
  implicit none
  stop 'MPIFAKE: mpi_error_string'
END SUBROUTINE  MPI_ERROR_STRING

subroutine  MPI_SCATTER()
  implicit none
  stop 'MPIFAKE: SCATTER'
END SUBROUTINE  MPI_SCATTER

subroutine  MPI_SCATTERV()
  implicit none
  stop 'MPIFAKE: SCATTERV'
END SUBROUTINE  MPI_SCATTERV

subroutine mpi_attr_get ()
  implicit none
  stop 'MPIFAKE: mpi_attr_get'
END SUBROUTINE  MPI_ATTR_GET

subroutine mpi_type_size()
  implicit none
  stop 'MPIFAKE: mpi_type_size'
END SUBROUTINE  MPI_TYPE_SIZE

subroutine mpi_comm_free ()
  implicit none
  stop 'MPIFAKE: mpi_comm_free'
END SUBROUTINE  MPI_COMM_FREE

subroutine mpi_waitany ()
  implicit none
  return !stop 'MPIFAKE: mpi_waitany'
END SUBROUTINE  MPI_WAITANY

subroutine mpi_irsend()
  implicit none
  stop 'MPIFAKE: mpi_irsend'
END SUBROUTINE  MPI_IRSEND

subroutine mpi_rsend()
  implicit none
  stop 'MPIFAKE: mpi_rsend'
END SUBROUTINE  MPI_RSEND

subroutine mpi_win_free()
  implicit none
  stop 'MPIFAKE: mpi_win_free'
END SUBROUTINE  MPI_WIN_FREE

subroutine mpi_win_fence()
  implicit none
  stop 'MPIFAKE: mpi_win_fence'
END SUBROUTINE  MPI_WIN_FENCE

subroutine mpi_win_create()
  implicit none
  stop 'MPIFAKE: mpi_win_create'
END SUBROUTINE  MPI_WIN_CREATE

subroutine mpi_win_lock()
  implicit none
  stop 'MPIFAKE: mpi_win_lock'
END SUBROUTINE  MPI_WIN_LOCK

subroutine mpi_win_unlock()
  implicit none
  stop 'MPIFAKE: mpi_win_unlock'
END SUBROUTINE  MPI_WIN_UNLOCK

subroutine mpi_win_wait()
  implicit none
  stop 'MPIFAKE: mpi_win_wait'
END SUBROUTINE  MPI_WIN_WAIT

subroutine mpi_win_complete()
  implicit none
  stop 'MPIFAKE: mpi_win_complete'
END SUBROUTINE  MPI_WIN_COMPLETE

subroutine mpi_win_post()
  implicit none
  stop 'MPIFAKE: mpi_win_post'
END SUBROUTINE  MPI_WIN_POST

subroutine mpi_win_start()
  implicit none
  stop 'MPIFAKE: mpi_win_start'
END SUBROUTINE  MPI_WIN_START

subroutine mpi_get()
  implicit none
  stop 'MPIFAKE: mpi_get'
END SUBROUTINE  MPI_GET

subroutine mpi_accumulate()
  implicit none
  stop 'MPIFAKE: mpi_accumulate'
END SUBROUTINE  MPI_ACCUMULATE

subroutine mpi_get_address()
  implicit none
  stop 'MPIFAKE: mpi_get_address'
END SUBROUTINE  MPI_GET_ADDRESS

subroutine mpi_type_create_struct()
  implicit none
  stop 'MPIFAKE: mpi_type_create_structure'
END SUBROUTINE  MPI_TYPE_CREATE_STRUCT

subroutine mpi_type_vector()
  implicit none
  stop 'MPIFAKE: mpi_type_vector'
END SUBROUTINE  MPI_TYPE_VECTOR

subroutine mpi_type_create_hvector()
  implicit none
  stop 'MPIFAKE: mpi_type_create_hvector'
END SUBROUTINE  MPI_TYPE_CREATE_HVECTOR

subroutine mpi_type_commit()
  implicit none
  stop 'MPIFAKE: mpi_type_commit'
END SUBROUTINE  MPI_TYPE_COMMIT

subroutine mpi_type_contiguous()
  implicit none
  stop 'MPIFAKE: mpi_type_contiguous'
END SUBROUTINE  MPI_TYPE_CONTIGUOUS

subroutine mpi_type_free()
  implicit none
  stop 'MPIFAKE: mpi_type_free'
END SUBROUTINE  MPI_TYPE_FREE

subroutine mpi_testall()
  implicit none
  stop 'MPIFAKE: mpi_testall'
END SUBROUTINE  MPI_TESTALL

subroutine mpi_info_create()
  implicit none
  stop 'MPIFAKE: mpi_info_create'
END SUBROUTINE  MPI_INFO_CREATE

subroutine mpi_info_set()
  implicit none
  stop 'MPIFAKE: mpi_info_set'
END SUBROUTINE  MPI_INFO_SET

subroutine mpi_info_free()
  implicit none
  stop 'MPIFAKE: mpi_info_free'
END SUBROUTINE  MPI_INFO_FREE

subroutine mpi_type_get_extent
  implicit none
  stop 'MPIFAKE: mpi_type_get_extent'
END SUBROUTINE mpi_type_get_extent

real(kind=8) function mpi_wtime()
  implicit none
  integer(kind=8) :: itns
  call nanosec(itns)
  mpi_wtime=real(itns,kind=8)*1.d-9
end function mpi_wtime

real(kind=8) function mpi_wtick()
  implicit none
  mpi_wtick=1.d-9
end function mpi_wtick

subroutine mpi_errhandler_get(comm,errhandler,ierr)
  implicit none
  integer, intent(in) :: comm
  integer, intent(out) :: errhandler,ierr
  errhandler = comm
  ierr = 0
end subroutine mpi_errhandler_get

subroutine mpi_errhandler_set(comm,errhandler,ierr)
  implicit none
  integer, intent(in) :: comm,errhandler
  integer, intent(out) :: ierr
  ierr = comm*errhandler
end subroutine mpi_errhandler_set

subroutine mpi_request_free(request,ierr)
  implicit none
  integer, intent(inout) :: request
  integer, intent(out) :: ierr
  ierr = request
end subroutine mpi_request_free

subroutine mpi_comm_dup(comm, newcomm, ierr)
  implicit none
  integer,intent(in) :: comm
  integer,intent(out) :: newcomm
  integer,intent(out) :: ierr
  newcomm = comm
  ierr = 0
end subroutine mpi_comm_dup

subroutine mpi_iprobe()
  implicit none
end subroutine mpi_iprobe

subroutine mpi_group_translate_ranks()
  implicit none
end subroutine mpi_group_translate_ranks

subroutine mpi_comm_split()
  implicit none
end subroutine mpi_comm_split

subroutine mpi_group_rank()
  implicit none
end subroutine mpi_group_rank

subroutine mpi_pack_size()
  implicit none
end subroutine mpi_pack_size

subroutine mpi_pack()
  implicit none
end subroutine mpi_pack

subroutine mpi_unpack()
  implicit none
end subroutine mpi_unpack

subroutine mpi_error_class()
  implicit none
end subroutine mpi_error_class

subroutine mpi_put()
  implicit none
end subroutine mpi_put
