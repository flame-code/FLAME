!*****************************************************************************************
module mod_electrostatics
    use mod_linked_lists, only: typ_linked_lists, typ_pia_arr
    use mod_spline, only: typ_spline
#if defined(HAVE_BPS)
    use Poisson_Solver, only: coulomb_operator
#endif
    implicit none
    type:: typ_poisson
        integer:: nat
        integer:: lda=0 !leading dimension of array pot
        integer:: ngpx !number of grid points in x direction.
        integer:: ngpy !number of grid points in y direction.
        integer:: ngpz !number of grid points in z direction.
        real(8):: hx !grid spacing in x direction.
        real(8):: hy !grid spacing in y direction.
        real(8):: hz !grid spacing in z direction.
        real(8):: hgrid(3,3) !grid spacing in all directions.
        real(8), allocatable:: rho(:,:,:) !charge density array.
        real(8), allocatable:: pot(:,:,:) !potential array.
        real(8), allocatable:: pots(:,:,:) !surface potential array.
        real(8), allocatable:: dpot(:,:,:) !differential of surface potential array.
        real(8), allocatable:: qgrad(:)
        real(8), allocatable:: qgrad_real(:)
        real(8), allocatable:: gw_ewald(:)
        real(8), allocatable:: q(:)
        real(8), allocatable:: gw(:)
        real(8), allocatable:: rcart(:,:)
        real(8), allocatable:: rgrad(:,:)
        integer:: npl, npu
        real(8):: dpm(3)
        real(8):: qpm(3,3)
        real(8):: beta !This is exactly the beta in the P3D paper.
        real(8):: cv(3,3)
        logical:: point_particle= .false.
        !logical:: cal_poisson= .false.
        !logical:: cal_qgrad= .false.
        !logical:: cal_force= .false.
        !logical:: set_grid= .true.
        logical:: reset_rho= .true.
        logical:: initialized= .false.
        character(256):: task_finit=""
        character(20):: bc='unknown'
        integer(8), allocatable:: plan_f(:) !Plans of forward fftw with size ngpz
        integer(8), allocatable:: plan_b(:) !Plans of inverse fftw with size ngpz
        integer(8), allocatable:: plan_fs(:) !Plans of forward fftw with size 2
        !bounds to assign charge density to grid points within a shpere.
        !integer, allocatable:: mboundg(:,:,:)
        !integer:: nbgpx !number of bound grid points in x direction.
        !integer:: nbgpy !number of bound grid points in y direction.
        !integer:: nbgpz !number of bound grid points in z direction.
        !integer:: nagpx !number of additional grid points in x direction.
        !integer:: nagpy !number of additional grid points in y direction.
        !integer:: nagpz !number of additional grid points in z direction.
        real(8):: epotfixed !fixed part of electrostatic energt in ewald.
        !real(8):: rcut !cutoff radius in real space.
        real(8):: alpha =-1 !splitting ewald parameter.
        real(8):: rgcut
        real(8):: ecut
        logical:: gw_identical= .false. !if True, all gaussian width assumed to be identical
        real(8):: efield !external electric field
        real(8):: vu, vl !voltage on upper and lower plane.
        real(8):: cell(3) !cell size in x,y,z direction.
        !type(typ_poisson_p3d):: poisson_p3d
        type(typ_linked_lists):: linked_lists
        type(typ_pia_arr):: pia_arr
        type(typ_spline):: spline
#if defined(HAVE_BPS)
        type(coulomb_operator):: pkernel
#endif
    end type typ_poisson
    !type, extends(typ_ewald):: typ_ewald_p3d
end module mod_electrostatics
!*****************************************************************************************
