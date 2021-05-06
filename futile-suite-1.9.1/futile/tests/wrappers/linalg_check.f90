!> @file
!!  Test of (some) linear algebra operations
!! @author
!!    Copyright (C) 2016-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
program linalg_check
  use futile
  use wrapper_linalg
  !use wrapper_mpi
  use f_random
  implicit none
  logical :: symm
  integer :: norb,nvctr,ncplx,nthreads
  integer :: iorb,jorb,icplx,ivctr!,idum
  !integer, parameter :: f_double = selected_real_kind(15, 307)
  real(f_double) :: diff
  real(f_double), dimension(2) :: dotp
  real(f_double), dimension(:,:,:), allocatable :: psi,hpsi
  real(f_double), dimension(:,:,:), allocatable :: ref_mat,mat 
  type(dictionary), pointer :: options
  !$ integer :: omp_get_max_threads

  call f_lib_initialize()
  !call mpiinit() !shall we use mpi here? do not think so...
  !call f_malloc_set_status(iproc=mpirank())
  !read command line
  call linalg_check_command_line_options(options)

  call yaml_new_document()
  
  !retrieve command line arguments
  norb=options//'norb'
  nvctr=options//'nvctr'
  ncplx=options//'ncplx'
  symm=options//'symm'
  call yaml_map('Options',options)
  call dict_free(options)

  !allocate arrays
  psi=f_malloc([ncplx,nvctr,norb],id='psi')
  hpsi=f_malloc([ncplx,nvctr,norb],id='hpsi')

  ref_mat=f_malloc([ncplx,norb,norb],id='ref_mat')
  mat=f_malloc([ncplx,norb,norb],id='mat')

!!$  idum=11
!!$  it=f_range([ncplx,nvctr,norb])
!!$  do while(iterating(it,on=psi,stride=4))
!!$     it%z(1)=real(builtin_rand(idum),f_double)
!!$     it%z(2)=real(builtin_rand(idum),f_double)
!!$     it%z(3)=real(builtin_rand(idum),f_double)
!!$     it%z(4)=real(builtin_rand(idum),f_double)
!!$     psi(it(1),it(2),it(3))=real(builtin_rand(idum),f_double)
!!$  end do

  !now fill the arrays with random numbers (use omp parallelisation)
  !$omp parallel do private(ivctr,iorb,icplx)
  do iorb=1,norb
     do ivctr=1,nvctr
        do icplx=1,ncplx
           psi(icplx,ivctr,iorb)=1.0_f_double !real(builtin_rand(idum),f_double)
           hpsi(icplx,ivctr,iorb)=1.0_f_double !real(builtin_rand(idum),f_double)
        end do
     end do
  end do
  !$omp end parallel do

  !then calculate the expected matrix
  do iorb=1,norb
     do jorb=1,norb
        dotp=0.0_f_double
        do ivctr=1,nvctr
           do icplx=1,ncplx
              dotp(1)=dotp(1)+&
                   psi(icplx,ivctr,iorb)*hpsi(icplx,ivctr,jorb)
              if (ncplx==2) then
                 dotp(2)=dotp(2)+(-1)**(icplx)*&
                      psi(3-icplx,ivctr,iorb)*hpsi(icplx,ivctr,jorb)
              end if
           end do
        end do
        ref_mat(:,iorb,jorb)=dotp(1:ncplx)
     end do
  end do

  !and compare with the result coming from dgemm
  call subspace_matrix(symm,psi,hpsi,ncplx,nvctr,norb,mat)

  !then compare the matrices
  call f_diff(int(size(mat),f_long),mat,ref_mat,diff)

  nthreads=1
  !$ nthreads=omp_get_max_threads()

  call yaml_map('Number of omp threads',nthreads)
  call yaml_map('The maximum difference is',diff)

  if (diff > epsilon(1.0_f_double) .and. norb <=5) then
     call yaml_map('Matrix (Real part)',mat(1,:,:))
     call yaml_map('Reference (Real part)',ref_mat(1,:,:))
  end if

  call f_free(psi)
  call f_free(hpsi)
  call f_free(mat)
  call f_free(ref_mat)

  !call mpifinalize()
  call f_lib_finalize()

contains

  !> Identify the options from command line
  !! and write the result in options dict
  subroutine linalg_check_command_line_options(options)
    implicit none
    !> dictionary of the options of the run
    !! on entry, it contains the options for initializing
    !! on exit, it contains in the key "BigDFT", a list of the
    !! dictionaries of each of the run that the local instance of BigDFT
    !! code has to execute
    type(dictionary), pointer :: options
    !local variables
    type(yaml_cl_parse) :: parser !< command line parser

    !define command-line options
    parser=yaml_cl_parse_null()
    !between these lines, for another executable using BigDFT as a blackbox,
    !other command line options can be specified
    !then the bigdft options can be specified
    call yaml_cl_parse_option(parser,&
         '{name: norb,'//&
         'shortname: o,'//&
         'default: 10,'//&
         'help_string: Number of orbitals for each set of vector,'//&
         'help_dict: {Allowed values: list of integers}}')

    call yaml_cl_parse_option(parser,&
         '{name: nvctr,'//&
         'shortname: v,'//&
         'default: 1000,'//&
         'help_string: Number of components for each vector,'//&
         'help_dict: {Allowed values: list of integers}}')

    call yaml_cl_parse_option(parser,&
         '{name: ncplx,'//&
         'shortname: c,'//&
         'default: 1,'//&
         'help_string: Number of complex components,'//&
         'help_dict: {Allowed values: 1, 2}}')

    call yaml_cl_parse_option(parser,&
         '{name: symm,'//&
         'shortname: s,'//&
         'default: No,'//&
         'help_string: Suppose symmetric result,'//&
         'help_dict: {Allowed values: Boolean}}')


    !parse command line, and retrieve arguments
    call yaml_cl_parse_cmd_line(parser,args=options)
    !free command line parser information
    call yaml_cl_parse_free(parser)

  end subroutine linalg_check_command_line_options

end program linalg_check

subroutine subspace_matrix(symm,psi,hpsi,ncplx,nvctrp,norb,lambda)
  use futile, only: f_double
  use wrapper_linalg, only: gemm, c_gemm
  implicit none
  logical, intent(in) :: symm !<symmetrize the result
  integer, intent(in) :: ncplx,nvctrp,norb
  real(f_double), dimension(nvctrp*ncplx,norb), intent(in) :: psi,hpsi
  real(f_double), dimension(ncplx,norb,norb), intent(out) :: lambda

  if (nvctrp==0) return
  if(ncplx==1) then
     if (symm) then
        !call gemmsy('T','N',norb,norb,nvctrp,1.0_f_double,psi(1,1),&
        call gemm('T','N',norb,norb,nvctrp,1.0_f_double,psi(1,1),&
             nvctrp,hpsi(1,1),nvctrp,0.0_f_double,&
             lambda(1,1,1),norb)
     else
        call gemm('T','N',norb,norb,nvctrp,1.0_f_double,psi(1,1),&
             nvctrp,hpsi(1,1),nvctrp,0.0_f_double,&
             lambda(1,1,1),norb)
     end if
  else
     !this part should be recheck in the case of nspinor == 2
     call c_gemm('C','N',norb,norb,nvctrp,(1.0_f_double,0.0_f_double),psi(1,1),&
          nvctrp,hpsi(1,1),nvctrp,(0.0_f_double,0.0_f_double),&
          lambda(1,1,1),norb)
  end if

end subroutine subspace_matrix
