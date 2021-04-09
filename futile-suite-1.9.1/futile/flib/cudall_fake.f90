!> @file
!!    Routines to bind fake argument for cuda low-level operations
!! @author
!!    Copyright (C) 2016-2016 BigDFT group  (LG)
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
subroutine cudafree()
   implicit none
   stop 'free'
 END SUBROUTINE cudafree
subroutine cudamalloc()
   implicit none
   stop 'allocation'
 END SUBROUTINE cudamalloc
subroutine cudamemset()
   implicit none
   stop 'memset'
 END SUBROUTINE cudamemset
subroutine get_gpu_data()
   implicit none
   stop 'get'
END SUBROUTINE get_gpu_data
subroutine synchronize()
   implicit none
   stop 'synchronize'
END SUBROUTINE synchronize
subroutine copy_gpu_data()
   implicit none
   stop 'copy_gpu_data'
END SUBROUTINE copy_gpu_data
subroutine reset_gpu_data()
  implicit none
  stop 'reset'
END SUBROUTINE reset_gpu_data
subroutine cuda_get_mem_info()
   implicit none
   stop 'cuda_get_mem_info'
 END SUBROUTINE cuda_get_mem_info
subroutine cudagetdevicecount()
   implicit none
   stop 'cudagetdevicecount'
 END SUBROUTINE cudagetdevicecount
subroutine cudasetdevice()
   implicit none
   stop 'cudasetdevice'
 END SUBROUTINE cudasetdevice
subroutine cudaresetdevice()
   implicit none
   stop 'cudaresetdevice'
 END SUBROUTINE cudaresetdevice
