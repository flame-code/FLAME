!*****************************************************************************************
module mod_poisson_neargrid
    implicit none 
    type typ_poisson
        integer:: ngp(3) !number of grid points in x,y,z direction.
        real(8):: h(3) !grid spacing in x,y,z direction.
        real(8), allocatable:: rho(:,:,:) !charge density array.
        integer, allocatable:: irho(:,:,:)
        integer, allocatable:: krho(:,:,:)
        real(8), allocatable:: volu(:,:) 
        integer, allocatable:: path(:,:)
        integer:: vdim, pdim, vnum, nvol, pnum
    end type typ_poisson 
    type(typ_poisson):: poisson
end module mod_poisson_neargrid
!*****************************************************************************************
module mod_poisson_ongrid
    implicit none 
    type typ_poisson
        integer:: ngpx !number of grid points in x direction.
        integer:: ngpy !number of grid points in y direction.
        integer:: ngpz !number of grid points in z direction.
        real(8):: hx !grid spacing in x direction.
        real(8):: hy !grid spacing in y direction.
        real(8):: hz !grid spacing in z direction.
        real(8), allocatable:: rho(:,:,:) !charge density array.
        integer, allocatable:: irho(:,:,:)
    end type typ_poisson 
    type(typ_poisson):: poisson
end module mod_poisson_ongrid
!*****************************************************************************************
module mod_poisson_weight
    implicit none 
    type, public :: weight
    real(8),allocatable:: w(:)
    end type
    type typ_poisson
        type(weight), allocatable:: weight(:,:,:)
        integer:: ngp(3) !number of grid points in x,y,z direction.
        real(8):: h(3) !grid spacing in x,y,z direction.
        real(8), allocatable:: rho(:,:,:) !charge density array.
        integer, allocatable:: irho(:,:,:)
        integer, allocatable:: krho(:,:,:)
        real(8), allocatable:: volu(:,:) 
        integer, allocatable:: path(:,:)
        integer:: vdim, pdim, vnum, nvol, pnum
        real(8)::step_size
    end type typ_poisson 
    type(typ_poisson):: poisson
end module mod_poisson_weight
!*****************************************************************************************
