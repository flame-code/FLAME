!> @file
!!  Test of the overlap point to point
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Program to test the OP2P (overlap point to point)
program example_op2p
  use BigDFT_API
  use overlap_point_to_point
  use mpi_layer
  implicit none
  !define the array of the partitioning of the orbital
  character(len=*), parameter :: subname='example_op2p'
  integer :: ndim,nproc,istep,norb,iproc,iorb,ierr,jproc,jorb,isorb,nelempr,ngroups,igroup,norb_res,ielem
  integer :: nelem,ioff,iadd,norbgr,iadd_result,i_stat,i_all,icount
  real(kind=4) :: tt
  real(kind=8) :: maxerr
  type(OP2P_descriptors) :: OP2P
  integer, dimension(:), allocatable :: norb_par
  integer, dimension(:,:), allocatable :: orbs_attributes
  real(kind=8), dimension(:), allocatable :: allarr,psi,results,expected_results
  real(kind=8), dimension(:), pointer :: allptr
  !try with a fake interface
  interface
     subroutine fake_operation(istep,iproc,igroup,remote_result,&
          idata_glob,jdata_glob,idata_loc,jdata_loc,ndatai,ndataj,&
          nvctri,nvctrj,nvctri_results,nvctrj_results,&
          objects_data_i,objects_data_j,&
          results_data_i,results_data_j)
       implicit none
       logical, intent(in) :: remote_result
       integer, intent(in) :: istep,iproc,igroup,idata_glob,jdata_glob,idata_loc,jdata_loc,ndatai,ndataj
       integer, intent(in) :: nvctri,nvctrj,nvctri_results,nvctrj_results
       real(kind=8), dimension(nvctri), intent(in) :: objects_data_i
       real(kind=8), dimension(nvctrj), intent(in) :: objects_data_j 
       real(kind=8), dimension(nvctri_results), intent(inout) :: results_data_i
       real(kind=8), dimension(nvctrj_results), intent(inout) :: results_data_j 
     end subroutine fake_operation
  end interface

  ! Start MPI in parallel version
  call MPI_INIT(ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,iproc,ierr)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ierr)

  !decide the number of groups to run with (maximum 10 groups)
  call random_integer(10,ngroups)
  ngroups=2 !for the moment only one group
  !decide the total number of elements for each group (maximum 1000 elements)       
  call random_integer(1000,norb)
  !the number of orbitals should be bigger than the number of groups
  norb=16!max(norb,ngroups)
  !allocate the corresponding arrays
  allocate(orbs_attributes(norb,3+ndebug),stat=i_stat)
  call memocc(i_stat,orbs_attributes,'orbs_attributes',subname)
  allocate(psi(norb+ndebug),stat=i_stat)
  call memocc(i_stat,psi,'psi',subname)
  allocate(results(norb+ndebug),stat=i_stat)
  call memocc(i_stat,results,'results',subname)
  allocate(expected_results(norb+ndebug),stat=i_stat)
  call memocc(i_stat,expected_results,'expected_results',subname)

  !decide artificially the orbital repartition and groups
  !objects_attributes(:,1) <= group to which the object belongs
  !objects_attributes(:,2) <= size in number of elements of the object
  !objects_attributes(:,3) <= size in number of elements of the results of the operations associated to the object

  !decide how many elements belong to each group
!!$  norb_res=norb
!!$  do igroup=ngroups,2,-1
!!$     call random_integer(10,nelem)
!!$     !at least one element per group
!!$     nelem=max(norb/ngroups+nelem,1)
!!$     do ielem=1,nelem
!!$        orbs_attributes(norb_res,1)=igroup
!!$        norb_res=norb_res-1
!!$     end do
!!$  end do
!!$
!!$  !first group for all the other elements
!!$  do iorb=1,norb_res
!!$     orbs_attributes(iorb,1)=1
!!$  end do

  orbs_attributes(1:10,1)=1
  orbs_attributes(11:16,1)=2
  !number of components for each of the elements (and for the results)
  do iorb=1,norb
     orbs_attributes(iorb,2)=1
     orbs_attributes(iorb,3)=1
  end do
  !initialize objects
  do iorb=1,norb
     call random_number(tt) !the seed is equal for each mpi process
     psi(iorb)=real(tt,kind=8)
     results(iorb)=0.0d0
     expected_results(iorb)=0.0d0
  end do

!!$  !calculate expected results (could be done with the same fake_operation )
!!$  do iorb=1,norb
!!$     expected_results(iorb)=0.0d0
!!$     do jorb=1,norb
!!$        expected_results(iorb)=expected_results(iorb)+psi(iorb)*psi(jorb)
!!$     end do
!!$  end do
  icount=0
  !call the reference operation to calculate the expected result
  ioff=0
  iadd=1
  iadd_result=1
  do igroup=1,ngroups
     norbgr=0
     do iorb=1,norb
        if (orbs_attributes(iorb,1) == igroup) then
           norbgr=norbgr+1
        end if
     end do
     !print *,'reference',ioff,iadd,norbgr,igroup,iadd_result
     call fake_operation(0,0,igroup,.false.,ioff,ioff,ioff,ioff,norbgr,norbgr,norbgr,norbgr,norbgr,norbgr,&
          psi(iadd),psi(iadd),&
          expected_results(iadd_result),expected_results(iadd_result))
     !count the number of elements for each object
     do iorb=1,norb
        if (orbs_attributes(iorb,1) == igroup) then
           ioff=ioff+1
           iadd=iadd+orbs_attributes(iorb,2)
           iadd_result=iadd_result+orbs_attributes(iorb,3)
        end if
     end do
  end do

  print *,'iproc,icount',iproc,icount

  icount=0
  !allocation of the data should be done accordingly to the number of objects
  allocate(norb_par(0:nproc-1+ndebug),stat=i_stat)
  call memocc(i_stat,norb_par,'norb_par',subname)

  isorb=0
  norb_par(0:3)=3
  norb_par(0:7)=2
  do jproc=0,nproc-2
     if (jproc < iproc) then
        isorb=isorb+norb_par(jproc)
     end if
  end do

!!$  norb_res=norb
!!$  do jproc=0,nproc-2
!!$     call random_integer(5,nelempr)
!!$     nelempr=min(max(norb/nproc+nelempr,0),norb_res)
!!$     if (norb_res==0) nelempr=0
!!$     norb_par(jproc)=nelempr
!!$     if (jproc < iproc) then
!!$        isorb=isorb+norb_par(jproc)
!!$     end if
!!$     norb_res=norb_res-nelempr
!!$  end do
!!$  !last process
!!$  norb_par(nproc-1)=norb_res


!!$norb_par(0)=102
!!$if (iproc==0) isorb=0
!!$if(iproc ==1) isorb=102
!!$norb_par(1)=98
  !print *,'a',norb_par,orbs_attributes

  !here any processor will initialise the global communications arrays needed for executing the op2p
  call initialize_OP2P_descriptors(.true.,iproc,nproc,norb,orbs_attributes,norb_par,OP2P)

  !this simulates the concurrent call of the op2p routine by all the processors
  !print *,'starting',iproc
  !adjust isorb followng the maximum
  if (norb_par(iproc)==0) isorb=norb-1
  call OP2P_communication(iproc,nproc,OP2P,psi(isorb+1),results(isorb+1),fake_operation)!,send_mpi_profile,receive_mpi_profile,&
  !wait_mpi_profile)
  call MPI_BARRIER(MPI_COMM_WORLD,ierr)
  !print *,'iproc,icount2',iproc,icount
  call free_OP2P_descriptors(OP2P,subname)
  !print *,'barrier',ierr
  !if (iproc==0) flush(unit=6)
  maxerr=0.0d0

  do iorb=1,norb_par(iproc)
     maxerr=max(maxerr,abs(results(iorb+isorb)-expected_results(iorb+isorb)))
     if (abs(results(iorb+isorb)-expected_results(iorb+isorb)) > 0.0d0) then
        !print *,'Exception:',iproc,iorb+isorb,psi(iorb+isorb),results(iorb+isorb),expected_results(iorb+isorb),&
        !     results(iorb+isorb)-expected_results(iorb+isorb)
        write(*,'(1x,a,i5,i5,1pe19.11)')'Exception:',iproc,iorb+isorb,results(iorb+isorb)-expected_results(iorb+isorb)
     end if
  end do
  !if (maxerr > epsilon(1.e0)) then
     !write(*,*)'error,iproc',iproc,maxerr
  !   stop 'ERROR'
  !end if
  
  if (nproc > 1) call mpiallred(maxerr,1,comm=MPI_SUM)

  if (iproc==0) print *,'FINAL result,maxerr',iproc,maxerr
  
  i_all=-product(shape(norb_par))*kind(norb_par)
  deallocate(norb_par,stat=i_stat)
  call memocc(i_stat,i_all,'norb_par',subname)
  i_all=-product(shape(psi))*kind(psi)
  deallocate(psi,stat=i_stat)
  call memocc(i_stat,i_all,'psi',subname)
  i_all=-product(shape(results))*kind(results)
  deallocate(results,stat=i_stat)
  call memocc(i_stat,i_all,'results',subname)
  i_all=-product(shape(orbs_attributes))*kind(orbs_attributes)
  deallocate(orbs_attributes,stat=i_stat)
  call memocc(i_stat,i_all,'orbs_attributes',subname)
  i_all=-product(shape(expected_results))*kind(expected_results)
  deallocate(expected_results,stat=i_stat)
  call memocc(i_stat,i_all,'expected_results',subname)

  !finalize memory counting
  call memocc(0,0,'count','stop')

  call MPI_FINALIZE(ierr)

end program example_op2p

!apply a fake operation as if all the data were already sent or received
subroutine fake_operation(istep,iproc,igroup,remote_result,&
     idata_glob,jdata_glob,idata_loc,jdata_loc,ndatai,ndataj,&
     nvctri,nvctrj,nvctri_results,nvctrj_results,&
     objects_data_i,objects_data_j,&
     results_data_i,results_data_j)
  implicit none
  logical, intent(in) :: remote_result
  integer, intent(in) :: istep,iproc,igroup,idata_glob,jdata_glob,idata_loc,jdata_loc,ndatai,ndataj
  integer, intent(in) :: nvctri,nvctrj,nvctri_results,nvctrj_results
  real(kind=8), dimension(nvctri), intent(in) :: objects_data_i
  real(kind=8), dimension(nvctrj), intent(in) :: objects_data_j 
  real(kind=8), dimension(nvctri_results), intent(inout) :: results_data_i
  real(kind=8), dimension(nvctrj_results), intent(inout) :: results_data_j 
  !local variables
  integer :: iorb,jorb,iorb_glob,jorb_glob,iorb_loc,jorb_loc
  
  print '(a,6i3)','fake_operation,istep,iproc,ndatai,ndataj,idata_glob,jdata_glob',&
      istep,iproc,ndatai,ndataj,idata_glob,jdata_glob
  !fill the results data with the sum of the products of the objects positions
  !put to zero the sending element
  if (remote_result) then
     do jorb=1,ndataj
        results_data_j(jorb)=0.0d0
     end do
  end if
  do iorb=1,ndatai
     iorb_loc=iorb+idata_loc
     iorb_glob=iorb+idata_glob
     do jorb=1,ndataj
        jorb_loc=jorb+jdata_loc
        jorb_glob=jorb+jdata_glob
        !in the final passage we should take the data from the objects_data arrays
        !print *,'objects_data',objects_data_i(iorb),objects_data_j(jorb)
        results_data_i(iorb)=&
             results_data_i(iorb)+objects_data_i(iorb)*objects_data_j(jorb)
        !icount=icount+1
        !print '(a,7(i4))','prova valori',iproc,istep,iorb_glob,jorb_glob,&
        !     int(objects_data_i(iorb_loc)),int(objects_data_j(jorb_loc)),int(results_data_i(iorb_loc))
        if (remote_result) then
           results_data_j(jorb)=&
                results_data_j(jorb)+objects_data_j(jorb)*objects_data_i(iorb)
           !icount=icount+1
        end if
     end do
  end do
  !print *,' end' 
end subroutine fake_operation


subroutine random_integer(nrange,irand)
  implicit none
  integer, intent(in) :: nrange
  integer, intent(out) :: irand
  !local variables
  real(kind=4) :: tt
  
  call random_number(tt)

  irand=nint(tt*real(nrange,kind=4))

end subroutine random_integer


subroutine fake_extra_operation(icount)
  implicit none
  integer, intent(inout) :: icount
  
  icount=icount+1
end subroutine fake_extra_operation
