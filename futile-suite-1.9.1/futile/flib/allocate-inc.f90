!> @file
!! Include fortran file for allocation template
!! 
!! @author
!!    Copyright (C) 2012-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
include 'allocate-base-inc.f90'
if (m%srcdata_add > int(0,kind=8)) &
     call c_memcopy(array,m%srcdata_add,f_sizeof(array))
include 'allocate-end-inc.f90'
