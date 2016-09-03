!*****************************************************************************************
module mod_spline
    implicit none
    type typ_spline
        !variables for spline of function "erf over r" in realspace
        integer:: nsp !number of spline nodes.
        real(8):: hsp !size of intervals for splines.
        real(8), allocatable:: fsp(:,:,:)
        real(8), allocatable:: fdsp(:,:,:) 
        logical:: do_coulomb=.false.
        logical:: do_tosifumi=.false.
    end type typ_spline
end module mod_spline
!*****************************************************************************************
