!> @file
!! Include fortran file for quick return in case
!! of omp activated
!! This file is useful in case of non-thread safe routine
!! that will become deactivated inside a parallel region
!! file included in several modules
!! this file should be included as !$ include to be considered only for OMP pragmas
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  logical :: in_omp
  logical, external :: omp_get_nested!,omp_in_parallel
  integer, external :: omp_get_thread_num
  in_omp=omp_get_thread_num() /=0 
  if (.not. in_omp) in_omp=omp_get_nested() !in_omp=omp_in_parallel()


  !!disable everything if we are into a OMP section
  !!timing routines are not thread-safe
  if(in_omp) return
  
