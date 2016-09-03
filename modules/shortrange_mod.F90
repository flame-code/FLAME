!*****************************************************************************************
module mod_shortrange
    implicit none
    type typ_tosifumi
        integer:: ntypinter, ntypat, interaction(10,10)
        real(8):: aaa(10)
        real(8):: bbb(10)
        real(8):: ccc(10)
        real(8):: ddd(10)
        real(8):: eee(10)
        real(8):: fff(10)
        !character(5):: 
        !real(8), allocatable:: 
    end type typ_tosifumi
    type typ_shortrange
        integer:: ntypinter
        integer:: ntypat
        integer:: interaction(10,10)
        real(8):: alpha
        real(8):: qq(10)
        type(typ_tosifumi):: tosifumi
    end type typ_shortrange
end module mod_shortrange
!*****************************************************************************************
