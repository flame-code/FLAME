!> @file
!! Test the MPI wrappers
!! @author
!!    Copyright (C) 2017-2017 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


module test_mpi_wrappers
  use f_precisions
  use dynamic_memory
  use wrapper_mpi
  use yaml_output
  use f_alltoall
  use fmpi_types
  implicit none

  private


  public :: test_mpi_alltoallv,test_mpi_allgather

  contains

    subroutine test_mpi_alltoallv(iproc, nproc, comm, maxsize_local, ntest)
      implicit none
      ! Calling arguments
      integer,intent(in) :: iproc, nproc, comm, maxsize_local, ntest

      ! Local variables
      integer :: itest, jproc, iel, irep
      integer(f_long) :: nstride, nsize_local, nGB, nMB, nkB, nB, i
      real(f_double) :: fac, t1, t2
      integer,dimension(:),allocatable :: sendcounts, recvcounts, senddispls, recvdispls
      real(f_double),dimension(:),allocatable :: sendbuf, recvbuf
      integer,parameter :: nrep = 5
      real(f_double),dimension(nrep,2) :: times
      character(len=512) :: memstring
      logical,dimension(2) :: correct

      nstride = maxsize_local/(ntest*nproc)

      fac = 1.0_f_long/real(ceiling(log(real(nproc+1,kind=f_double))),f_double)

      sendcounts = f_malloc(0.to.nproc-1,id='sendcounts')
      recvcounts = f_malloc(0.to.nproc-1,id='recvdispls')
      senddispls = f_malloc(0.to.nproc-1,id='senddispls')
      recvdispls = f_malloc(0.to.nproc-1,id='recvdispls')

      if (iproc==0) then
          call yaml_sequence_open('Performing mpi_alltoallv benchmarks')
      end if
      do itest=1,ntest
          nsize_local = itest*nstride*nproc
          sendbuf = f_malloc(nsize_local,id='sendbuf')
          recvbuf = f_malloc(nsize_local,id='recvbuf')
          do jproc=0,nproc-1
              sendcounts(jproc) = nstride*itest
              recvcounts(jproc) = nstride*itest
              senddispls(jproc) = jproc*nstride*itest
              recvdispls(jproc) = jproc*nstride*itest
          end do
          iel = 0
          do jproc=0,nproc-1
              do i=1,itest*nstride
                  iel = iel + 1
                  sendbuf(iel) = real(iproc,kind=f_double) + fac*real(jproc,kind=f_double)
              end do
          end do

          if (iproc==0) then
              call yaml_sequence(advance='no')
               call yaml_map('Local number of elements',nsize_local)
               nB = nsize_local*kind(sendbuf)
               nGB = nB/1024/1024/1024
               nB = nB-nGB*1024*1024*1024
               nMB = nB/1024/1024
               nB = nB-nMB*1024*1024
               nkB=nB/1024
               nB = nB-nkB*1024
               memstring = trim(yaml_toa(nGB))//'GB'//&
                          &trim(yaml_toa(nMB))//'MB'//&
                          &trim(yaml_toa(nkB))//'kB'//&
                          &trim(yaml_toa(nB))//'B'
               call yaml_map('Local size send/receive buffer',trim(memstring))
               call yaml_sequence_open('Timings')
          end if

          do irep=1,nrep
              t1 = mpi_wtime()
              !call mpialltoallv(sendbuf, sendcounts, senddispls, &
              !     recvbuf, recvcounts, recvdispls, comm, algorithm='mpi_get')
              call fmpi_alltoall(sendbuf,recvbuf=recvbuf,&
                   sendcounts=sendcounts,sdispls=senddispls,&
                   recvcounts=recvcounts,rdispls=recvdispls,comm=comm,algorithm=ONESIDED_ENUM)
              t2 = mpi_wtime()
              times(irep,1) = t2 - t1
              correct(1) = check_result(iproc, nproc, comm, itest, nstride, nsize_local, fac, recvbuf)
          end do
          if (iproc==0) then
              call yaml_sequence(advance='no')
              call yaml_map('Algorithm','mpi_get')
              call yaml_map('Result correct', correct(1))
              call yaml_map('All measurements',times(:,1),fmt='(es10.3)')
              call yaml_map('Min/Max/Median',[minval(times(:,1)),maxval(times(:,1)),median(nrep,times(:,1))],fmt='(es10.3)')
          end if

          do irep=1,nrep
              t1 = mpi_wtime()
              !call mpialltoallv(sendbuf, sendcounts, senddispls, &
              !     recvbuf, recvcounts, recvdispls, comm, algorithm='native')
              call fmpi_alltoall(sendbuf,recvbuf=recvbuf,&
                   sendcounts=sendcounts,sdispls=senddispls,&
                   recvcounts=recvcounts,rdispls=recvdispls,comm=comm,algorithm=VARIABLE_ENUM)

              t2 = mpi_wtime()
              times(irep,2) = t2 - t1
              correct(2) = check_result(iproc, nproc, comm, itest, nstride, nsize_local, fac, recvbuf)
          end do
          if (iproc==0) then
              call yaml_sequence(advance='no')
              call yaml_map('Algorithm','native')
              call yaml_map('Result correct', correct(2))
              call yaml_map('All measurements',times(:,2),fmt='(es10.3)')
              call yaml_map('Min/Max/Median',[minval(times(:,2)),maxval(times(:,2)),median(nrep,times(:,2))],fmt='(es10.3)')
          end if

          if (iproc==0) then
               call yaml_sequence_close()
          end if

          call f_free(sendbuf)
          call f_free(recvbuf)

      end do

      call f_free(sendcounts)
      call f_free(recvcounts)
      call f_free(senddispls)
      call f_free(recvdispls)

      if (iproc==0) then
          call yaml_sequence_close()
      end if



    end subroutine test_mpi_alltoallv

    subroutine test_mpi_allgather(nvctrp,comm)
      use dictionaries, only: f_err_throw
      use f_utils, only: f_pause
      implicit none
      integer, intent(in) :: comm,nvctrp
      !local variables
      integer :: nproc,iproc,jproc
      integer, dimension(:), allocatable :: counts,displs
      real(f_double), dimension(:), allocatable :: matp,mat

      nproc=mpisize(comm)
      iproc=mpirank(comm)
      matp=f_malloc(nvctrp,id='matp')
      mat=f_malloc(nvctrp*nproc,id='mat')

      !initialize with test data
      matp=iproc

      counts=f_malloc(0.to.nproc-1,id='counts')
      displs=f_malloc0(0.to.nproc-1,id='displs')
      counts=nvctrp
      do jproc=0,nproc-2
         displs(jproc+1)=displs(jproc)+counts(jproc)
      end do

      call allgather_with_algo()
      !call allgather_with_algo(VARIABLE_ENUM) this cannot be forced
      call allgather_with_algo(NOT_VARIABLE_ENUM)
      call allgather_with_algo(ONESIDED_ENUM)

      call f_free(mat,matp)
      call f_free(counts,displs)

    contains

      subroutine allgather_with_algo(algo)
        use f_enums
        implicit none
        type(f_enumerator), intent(in), optional :: algo
        call fmpi_allgather(sendbuf=matp,sendcount=nvctrp,recvbuf=mat,&
             recvcount=nvctrp,comm=comm,algorithm=algo)

        !verify test data
        if (.not. verify_mat(nproc,nvctrp,mat)) then
           call f_err_throw('Mpi all gather not working')
        end if
      end subroutine allgather_with_algo

    end subroutine test_mpi_allgather



    function verify_mat(nproc,nvctrp,mat) result(ok)
      implicit none
      integer, intent(in) :: nproc,nvctrp
      real(f_double), dimension(nvctrp,nproc), intent(in) :: mat
      logical :: ok
      !local variables
      integer :: jproc,i

      do jproc=1,nproc
         do i=1,nvctrp
            if (nint(mat(i,jproc)) /= jproc-1) then
               call yaml_map('Error',[i,jproc-1,nint(mat(i,jproc))])
               ok=.false.
               return
            end if
         end do
      end do
      ok=.true.
      
    end function verify_mat

    function check_result(iproc, nproc, comm, itest, nstride, nsize_local, fac, recvbuf)
      implicit none

      ! Calling arguments
      logical :: check_result
      integer,intent(in) :: iproc, nproc, comm, itest
      integer(f_long),intent(in) :: nstride, nsize_local
      real(f_double),intent(in) :: fac
      real(f_double),dimension(nsize_local),intent(in) :: recvbuf

      ! Local variables
      integer :: iel, jproc
      integer(f_long) :: i
      real(f_double) :: val
      
      check_result = .true.
      iel = 0
      do jproc=0,nproc-1
          do i=1,itest*nstride
              val = fac*real(iproc,kind=f_double) + real(jproc,kind=f_double)
              iel = iel + 1
              if (abs(val-recvbuf(iel))>1.e-10_f_double) then
                  check_result = .false.
              end if
          end do
      end do
      !call mpiallred(check_result, 1, mpi_land, comm=comm)
      call fmpi_allreduce(check_result, 1, FMPI_LAND, comm=comm)
    end function check_result


    subroutine insertion_sort(n, arr)
      implicit none

      ! Calling variables
      integer,intent(in) :: n
      real(f_double),dimension(n),intent(inout) :: arr

      ! Local variables
      integer :: i, j
      real(f_double) :: temp

      do i=2,n
          j = i-1
          temp = arr(i)
          !do while (j>=1 .and. arr(j)>temp)
          do
             if (j<1) exit
             if (arr(j)<=temp) exit
             arr(j+1) = arr(j)
             j = j-1
          end do
          arr(j+1) = temp
       end do

    end subroutine insertion_sort



    function median(n, a)
      implicit none
      integer,intent(in) :: n
      real(f_double),dimension(n),intent(in) :: a
      real(f_double) :: median

      ! Local variables
      real(f_double),dimension(:),allocatable :: arr

      arr = f_malloc(n,id='arr')
      call f_memcpy(src=a, dest=arr)
      ! Sort the data (insertion sort)
      call insertion_sort(n, arr)

       ! Take the median
       if (mod(n,2)==1) then
           ! Odd n
           median = arr(n-n/2)
       else
           ! Even n
           median = 0.5_f_double*(arr(n/2)+arr(n/2+1))
       end if

      call f_free(arr)

    end function median


end module test_mpi_wrappers
