!> @file
!! Fake routines in order to compile without OpenCL
!! @author 
!!    Copyright (C) 2010-2011 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 

subroutine release_acceleration_OCL()
    implicit none
    stop 'FAKE release_acceleration_OCL'
END SUBROUTINE release_acceleration_OCL

subroutine init_acceleration_OCL()
    implicit none
    stop 'FAKE init_acceleration_OCL'
END SUBROUTINE init_acceleration_OCL

subroutine allocate_data_OCL()
    implicit none
    stop 'FAKE allocate_data_OCL'
END SUBROUTINE allocate_data_OCL


subroutine free_gpu_OCL()
    implicit none
    stop 'FAKE free_gpu_OCL'
END SUBROUTINE free_gpu_OCL


subroutine local_hamiltonian_OCL()
    implicit none
    stop 'FAKE local_hamiltonian_OCL'
END SUBROUTINE local_hamiltonian_OCL

subroutine finish_hamiltonian_OCL()
    implicit none
    stop 'FAKE finish_hamiltonian_OCL'
END SUBROUTINE finish_hamiltonian_OCL

subroutine preconditionall_OCL()
    implicit none
    stop 'FAKE preconditionall_OCL'
END SUBROUTINE preconditionall_OCL

subroutine local_partial_density_OCL()
    implicit none
    stop 'FAKE local_partial_density_OCL'
END SUBROUTINE local_partial_density_OCL

subroutine daub_to_isf_OCL()
    implicit none
    stop 'FAKE daub_to_isf_OCL'
END SUBROUTINE daub_to_isf_OCL


