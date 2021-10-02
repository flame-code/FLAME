!*****************************************************************************************
module mod_symfunc_data
    use mod_linked_lists, only: typ_linked_lists
    implicit none
    public
    type typ_symfunc_data
        integer:: ng=-1
        integer:: nat=-1
        real(8), allocatable:: y(:,:)
        real(8), allocatable:: y0d_bond(:,:)
        real(8), allocatable:: y0d(:,:,:)
        real(8), allocatable:: y0dr(:,:,:)
        type(typ_linked_lists):: linked_lists
    end type typ_symfunc_data
!contains
!*****************************************************************************************
end module mod_symfunc_data
!*****************************************************************************************
