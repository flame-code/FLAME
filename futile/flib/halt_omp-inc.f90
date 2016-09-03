!> @file
!! Include fortran file for halting program in case
!! of omp activated
!! This file is useful in case of non-thread safe routine
!! that will raise a severe error if called inside a parallel region
!! this file should be included as !$ include to be considered only for OMP pragmas
!! @author
!!    Copyright (C) 2012-2013 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

  logical :: in_omp
  logical, external :: omp_in_parallel,omp_get_nested
  in_omp=omp_in_parallel() .or. omp_get_nested()

  !!disable everything if we are into a OMP section
  !!timing routines are not thread-safe
  if(in_omp) then
     print *,'#SEVERE: critical f_lib routines cannot be called from within a parallel region! Exiting...'
     call f_err_severe_internal()
  end if
