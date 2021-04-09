!> @file
!! Include fortran file for mpi_get

!! @author
!!    Copyright (C) 2017-2017 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
  include 'get-decl-inc.f90'
  origin_ptr=>f_subptr(origin_addr,size=count,from=from)
  include 'get-end-inc.f90'
