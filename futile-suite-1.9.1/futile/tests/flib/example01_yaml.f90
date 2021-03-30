!------------------------------------------------------------------------------
! Futile:
!------------------------------------------------------------------------------
!
! MODULE: Module Name
!
!> @author
!> Luigi Genovese and Damien Caliste
!
! DESCRIPTION: 
!> Small program to show how yaml_map works.
!
! REVISION HISTORY:
! 17/11/2016 - Initial Version
!------------------------------------------------------------------------------
program example01_yaml
    use yaml_output
    implicit none
    call f_lib_initialize()
    call yaml_map("Logical value",.true.)
    call yaml_map("Integer value",123)
    call yaml_map("Real value",4.5)
    call f_lib_finalize()
end program example01_yaml
