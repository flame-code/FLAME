!> @file
!!  Test of the overlap point to point, modern version
!! @author
!!    Copyright (C) 2015-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


program OP2P_check
  use futile
  use overlap_point_to_point
  use wrapper_MPI
  implicit none
  logical :: symmetric,nearest_neighbor
  integer :: iproc,jproc,nproc,norbp,ngroup,igroup,ndim,norb,iobj,jobj,kobj
  integer, dimension(:), allocatable :: nobj,nobj_p
  integer, dimension(:,:), allocatable :: nobj_par
  type(dictionary), pointer :: options

  call f_lib_initialize()

  call mpiinit()
  iproc = mpirank()

  call f_malloc_set_status(iproc=iproc)

  !read command line
  call OP2P_check_command_line_options(options)

  if (iproc==0) then
    call yaml_new_document()
    call yaml_dict_dump(options)
  end if

  nproc=mpisize()
  ngroup=dict_len(options//'objects')
  if (iproc==0 .and. ngroup <= 0) call f_err_throw('Error, number of groups must be more than one')

  nobj=f_malloc(ngroup,id='nobj')
  nobj=options//'objects'
  ndim=options//'ndim'
  symmetric=options//'symmetric'
  nearest_neighbor=options//'nn-pattern'
  call dict_free(options)

  !construct the number of objects per processor
  norb=sum(nobj)
  norbp=norb/nproc

  nobj_p=f_malloc(0.to.nproc-1,id='nobj_p')
  nobj_p=norbp
  !the first processes have more orbitals
  jproc=norb-norbp*nproc-1
  call f_increment(nobj_p(:jproc))

  if (iproc==0 .and. sum(nobj_p) /= norb) &
       call f_err_throw('Error in orbital repartition; norb is'+norb+' and nobj_p is'+yaml_toa(nobj_p))

  !construct the OP2P scheme and test it
  nobj_par = f_malloc0((/ 0.to.nproc-1, 1.to.ngroup /),id='nobj_par')
  iobj=0
  igroup=1
  do jproc=0,nproc-1
     kobj=0
     do jobj=1,nobj_p(jproc)
        if (iobj == nobj(igroup)) then
           nobj_par(jproc,igroup)=kobj
           iobj=0
           kobj=0
           call f_increment(igroup)
        end if
        call f_increment(iobj)
        call f_increment(kobj)
     end do
     nobj_par(jproc,igroup)=kobj
  end do

  if (iproc==0 .and. any(sum(nobj_par,dim=1) /= nobj)) &
        call f_err_throw('Error in orbital repartition'+yaml_toa(mpirank())+';'+yaml_toa(sum(nobj_par,dim=1)))

  if (iproc==0) then
     call yaml_map('Orbital repartition per group',nobj)
     call yaml_map('Orbital repartition per mpi',nobj_p)
     call yaml_map('Groups per proc',nobj_par)
     call yaml_map('Starting simulation for the operator, symmetricity activated',symmetric)
  end if

  call f_free(nobj)
  call f_free(nobj_p)


  call OP2P_unitary_test(mpiworld(),mpirank(),nproc,ngroup,ndim,nobj_par,symmetric,nearest_neighbor)

  call f_free(nobj_par)
  call mpifinalize()
  call f_lib_finalize()

  contains

    !> Identify the options from command line
    !! and write the result in options dict
    subroutine OP2P_check_command_line_options(options)
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
      call OP2P_check_options(parser)
      !parse command line, and retrieve arguments
      call yaml_cl_parse_cmd_line(parser,args=options)
      !free command line parser information
      call yaml_cl_parse_free(parser)

    end subroutine OP2P_check_command_line_options

end program OP2P_check


!> Check the options
subroutine OP2P_check_options(parser)
  use futile
  implicit none
  type(yaml_cl_parse), intent(inout) :: parser

  call yaml_cl_parse_option(parser,'ndim','10',&
       'Size of the object','n',&
       dict_new('Usage' .is. &
       'Sizes of the unitary object of the check',&
       'Allowed values' .is. &
       'Yaml list of integers. If a scalar integer is given, all the dimensions will have this size.'))

  call yaml_cl_parse_option(parser,'objects','[10, 10]',&
       'Objects per group','o',&
       dict_new('Usage' .is. &
       'Set the number of objects per group',&
       'Allowed values' .is. &
       'Yaml list of integers. The orbitals are then distributed among the processors.'))

  call yaml_cl_parse_option(parser,'symmetric','Yes',&
       'Symmetricity','s',&
       dict_new('Usage' .is. &
       'Boolean, set the symmetricity of the operation.'))

  call yaml_cl_parse_option(parser,'nn-pattern','No',&
       'Nearest-Neigbor communication','c',&
       dict_new('Usage' .is. &
       'Boolean, adjust the communication pattern of the operation.'))

end subroutine OP2P_Check_options
