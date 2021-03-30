!> @file
!!  High level routines which needs more medium-level modules of the f_lib.
!!  They should be external, in the sense that no interface should be needed to call them.

!! @author Luigi Genovese
!!    Copyright (C) 2012-2015 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

!> Dump all errors which were handled
subroutine f_dump_all_errors(unit)
  use dictionaries, only: f_get_error_dict,f_get_past_error,&
       max_field_length,f_get_no_of_errors
  use yaml_output, only: yaml_dict_dump,yaml_map,yaml_get_default_stream
  implicit none
  !> if positive, write also the errors in the corresponding stream
  integer, intent(in) :: unit
  !local variables
  integer :: ierr,iunit_def,err_past
  character(len=max_field_length) :: add_msg

  !retrieve current unit
  call yaml_get_default_stream(iunit_def)

  if (unit > 0) iunit_def=unit

  do ierr=0,f_get_no_of_errors()-1
     err_past=f_get_past_error(ierr,add_msg)
     call yaml_dict_dump(f_get_error_dict(err_past),unit=iunit_def)
     if (trim(add_msg)/= 'UNKNOWN') &
          call yaml_map('Additional Info',add_msg,unit=iunit_def)
  end do
end subroutine f_dump_all_errors


!> Initialize all arrays and dictionaries for handling the errors
subroutine initialize_flib_errors()
  use dictionaries, only: dictionaries_errors
  use yaml_output, only: yaml_output_errors
  use f_utils, only: f_utils_errors
  use yaml_parse, only: yaml_parse_errors
  use dynamic_memory, only: dynamic_memory_errors
  use f_refcnts, only: refcnts_errors
  use time_profiling, only: timing_errors
  use f_input_file, only: input_file_errors
  implicit none

  call dictionaries_errors()
  call f_utils_errors()
  call yaml_output_errors()
  !Initialize the error to parse yaml documents
  call yaml_parse_errors()
  call refcnts_errors()
  call dynamic_memory_errors()
  call timing_errors()
  call input_file_errors()

end subroutine initialize_flib_errors


subroutine initialize_flib_timing_categories()
  use time_profiling, only: f_timing_category,f_timing_category_group
  use dynamic_memory, only: TCAT_ARRAY_ALLOCATIONS,&
        TCAT_ROUTINE_PROFILING
  use f_utils, only: TCAT_INIT_TO_ZERO
  implicit none
  character(len=*), parameter :: flibgrp='Flib LowLevel'
  !group of f_lib operations, separate category
  call f_timing_category_group(flibgrp,&
       'Low Level operations of flib module collection')
  !extract from BigDFT categories the ones which are associated to flib
  call f_timing_category('Init to Zero',flibgrp,'Memset of storage space',&
       TCAT_INIT_TO_ZERO)
  call f_timing_category('Array allocations',flibgrp,&
       'Heap storage allocation and associated profiling',&
       TCAT_ARRAY_ALLOCATIONS)
  call f_timing_category('Routine Profiling',flibgrp,&
       'Profiling performances for debugging',&
       TCAT_ROUTINE_PROFILING)


end subroutine initialize_flib_timing_categories


!> Routine which initializes f_lib global pointers, to be called before any action
!! is taken
subroutine f_lib_initialize()
  use f_environment
  use dictionaries, only: f_err_initialize
  use dynamic_memory, only: f_malloc_initialize
  use time_profiling, only: f_timing_initialize
  use f_bibliography, only: f_bib_initialize
  implicit none

  call f_environment_acquire()
  !general initialization, for lowest level f_lib calling
  call f_err_initialize()
  call initialize_flib_errors()
  !initializtion of memory allocation profiling
  call f_malloc_initialize()
  !initialization of timing module
  call f_timing_initialize()
  !initialization of internal timing categories of f_lib
  call initialize_flib_timing_categories()
  !initialization of bibliography database
  call f_bib_initialize()

end subroutine f_lib_initialize


!> Calls f_err_severe from outside the module
subroutine f_lib_err_severe_external(message)
  use dictionaries, only: f_err_severe
  implicit none
  character(len=*), intent(in) :: message
  write(0,*)trim(message)
  call f_err_severe()
end subroutine f_lib_err_severe_external


!> Routine which finalize f_lib and dump the finalization information
subroutine f_lib_finalize()
!  use dictionaries_base, only: dictionary_check_leak
  use dictionaries, only: f_err_finalize,dict_get_num,dictionary_check_leak
  use dynamic_memory, only: f_malloc_finalize
  use yaml_output, only: yaml_close_all_streams,yaml_map,yaml_comment!,yaml_walltime_toa
  use yaml_parse, only: yaml_parse_errors_finalize
  use time_profiling, only: f_timing_finalize,f_clock
  use yaml_strings, only: operator(//)
  use f_utils, only: f_humantime
  use module_f_objects, only: f_object_finalize
  use f_bibliography, only: f_bib_finalize
  implicit none
  !local variables
  integer :: ndict,ndict_max,iproc,nlibs,nlibs_max
  call f_object_finalize()
  call f_malloc_finalize(process_id=iproc)
  !print maximal value of dictionary usage
  if (iproc == 0) then
     call yaml_map('Walltime since initialization',f_humantime(f_clock()))
     call dict_get_num(ndict,ndict_max,nlibs,nlibs_max)
     call yaml_map('Max No. of dictionaries used',ndict_max, advance='no')
     call yaml_comment('('//ndict//' still in use)')
     !general finalization, the f_lib should come back to uninitialized status
     call yaml_map('Number of dictionary folders allocated',nlibs_max)
  end if
  call yaml_close_all_streams()
  call f_bib_finalize()
  call yaml_parse_errors_finalize()
  call f_err_finalize()
  call f_timing_finalize()
  !debug, once again
  call dictionary_check_leak()

end subroutine f_lib_finalize


!> Finalize f_lib but do not print out report information
subroutine f_lib_finalize_noreport()
  use dynamic_memory, only: f_malloc_finalize
  use dictionaries, only: f_err_finalize
  use yaml_output, only: yaml_close_all_streams
  use yaml_parse, only: yaml_parse_errors_finalize
  use time_profiling, only: f_timing_finalize
  use f_bibliography, only: f_bib_finalize
  implicit none
  call f_malloc_finalize(dump=.false.)
  call yaml_close_all_streams()
  call f_bib_finalize()
  call yaml_parse_errors_finalize()
  call f_err_finalize()
  call f_timing_finalize()
end subroutine f_lib_finalize_noreport

!> Some constants about Futile.
module f_lib_package
include 'configure.inc'
end module f_lib_package

!> complete module of futile, to be used for high-level programs
module futile
  use f_environment
  use f_precisions
  use yaml_strings
  use dictionaries
  use yaml_output
  use f_utils
  use yaml_parse
  use dynamic_memory
  use f_refcnts
  use time_profiling
  use f_input_file
  use f_enums
  use f_bibliography
  use f_ternary
  use f_random
end module futile
