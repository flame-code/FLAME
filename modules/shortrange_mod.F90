!*****************************************************************************************
module mod_shortrange
    implicit none
    type typ_tosifumi
        integer:: ntypinter, ntypat
        integer,allocatable :: interaction(:,:)
        real(8),allocatable :: aaa(:)
        real(8),allocatable :: bbb(:)
        real(8),allocatable :: ccc(:)
        real(8),allocatable :: ddd(:)
        real(8),allocatable :: eee(:)
        real(8),allocatable :: fff(:)
        !character(5):: 
        !real(8), allocatable:: 
    end type typ_tosifumi
    type typ_shortrange
        integer:: ntypinter
        integer:: ntypat
        integer,allocatable:: interaction(:,:)
        real(8):: alpha
        real(8),allocatable:: qq(:)
        type(typ_tosifumi):: tosifumi
    end type typ_shortrange
end module mod_shortrange
!*****************************************************************************************
