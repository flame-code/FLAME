!> @file
!! Include fortran file for deallocation template, multiple freeing of homogeneous and inhomogeneous kinds
!! @author
!!    Copyright (C) 2012-2014 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
 
  call f_free(arrayA)
  call f_free(arrayB)
  if (present(arrayC)) then
     call f_free(arrayC)
  end if
  if (present(arrayD)) then
     call f_free(arrayD)
  end if
  if (present(arrayE)) then
     call f_free(arrayE)
  end if
  if (present(arrayF)) then
     call f_free(arrayF)
  end if
  if (present(arrayG)) then
     call f_free(arrayG)
  end if
  if (present(arrayH)) then
     call f_free(arrayH)
  end if
