!*****************************************************************************************
module mod_symfunc
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    private
    !public:: 
    type, public:: typ_symfunc
        integer:: ng=-1
        integer:: nat=-1
        real(8), allocatable:: y(:,:)
        real(8), allocatable:: y0d_bond(:,:)
        real(8), allocatable:: y0d(:,:,:)
        real(8), allocatable:: y0dr(:,:,:)
        type(typ_linked_lists):: linked_lists
    end type typ_symfunc
    type, public:: typ_symfunc_arr
        integer:: nconf=-1
        type(typ_symfunc), allocatable:: symfunc(:)
    end type typ_symfunc_arr
!contains
!*****************************************************************************************
end module mod_symfunc
!*****************************************************************************************
