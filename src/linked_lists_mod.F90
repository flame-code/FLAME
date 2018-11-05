!*****************************************************************************************
module mod_linked_lists
    use mod_atoms, only: typ_atoms
    implicit none
    type typ_pia
        real(8):: r
        real(8):: fc
        real(8):: fcd
        real(8):: dr(3)
    end type typ_pia
    type typ_pia_arr
        type(typ_pia), allocatable:: pia(:)
    end type typ_pia_arr
    type, extends(typ_atoms):: typ_linked_lists
        !limits of subcells to be considered in short range calculation.
        integer, allocatable:: limnbx(:,:,:) 
        integer, allocatable:: limnby(:,:)
        integer, allocatable:: head(:,:,:) !heads of cells
        integer, allocatable:: list(:) !lists of particles in cells
        integer, allocatable:: prime(:,:,:)
        integer, allocatable:: last(:,:,:)
        integer, allocatable:: perm(:)
        integer, allocatable:: maincell(:)
        integer:: mx, my, mz !length of head array in x,y,z directions.
        !mlimnb is maximum number of subcells in one direction to be considered in 
        !short range calculation and it is related to the rcut.
        integer:: mlimnb
        integer:: mlimnb1
        integer:: mlimnb2
        integer:: mlimnb3
        integer:: next 
        real(8):: avgnndis !average nearest neighbor distance.
        real(8):: scl !subcell length.
        real(8):: rcut !cutoff radius in real space.
        logical:: triplex=.false.
        integer :: maxbound_rad
        integer :: maxbound_ang
        integer, allocatable:: bound_ang(:,:)
        integer, allocatable:: bound_rad(:,:)
        integer, allocatable:: prime_bound(:)
    end type typ_linked_lists
end module mod_linked_lists
!*****************************************************************************************
