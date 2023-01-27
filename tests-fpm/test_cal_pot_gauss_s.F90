!*****************************************************************************************
subroutine test_cal_pot_gauss_s
    use iso_fortran_env, only: error_unit, output_unit
    use mod_fit_bf_cent2, only: cal_pot_gauss_s
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_rat, atom_allocate, atom_deallocate_old
    use mod_electrostatics, only: typ_poisson
    use mod_colors, only: green_passed, red_failed
    implicit none
    !integer:: nx, ny, nz
    real(8):: gw, q, xyz(3), ener, gausswidth(1), errmax
    real(8):: errmax_vgrad
    !logical::
    type(typ_parini):: parini
    type(typ_poisson):: poisson
    type(typ_atoms):: atoms
    real(8), allocatable:: pot(:,:,:)
    real(8), allocatable:: vgrad(:,:,:)
    real(8), allocatable:: pot_t(:,:,:)
    real(8), allocatable:: vgrad_t(:,:,:)
    integer:: ix, iy, iz
    real(8):: x, y, z
    parini%screening_factor=0.3d0
    parini%cal_scn=.true.
    parini%psolver='bigdft'
    poisson%bc='free'
    xyz(1)=12.5d0
    xyz(2)=12.5d0
    xyz(3)=12.5d0
    poisson%task_finit="alloc_rho"
    poisson%ngpx=50
    poisson%ngpy=50
    poisson%ngpz=50
    poisson%hgrid(1:3,1:3)=0.d0
    poisson%hgrid(1,1)=0.5d0
    poisson%hgrid(2,2)=0.5d0
    poisson%hgrid(3,3)=0.5d0
    poisson%xyz111(1:3)=0.d0
    gw=2.d0
    gausswidth(1)=gw
    q=1.2d0
    allocate(pot(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(pot_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    call cal_pot_gauss_s(parini,poisson,pot,.true.,xyz,gw,q,vgrad)
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    call atom_allocate(atoms,1,0,0)
    atoms%ratp(1:3,1)=xyz(1:3)
    call update_rat(atoms)
    atoms%boundcond='free'
    call init_hartree(parini,atoms,poisson,gausswidth)
    call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,1),q,gausswidth(1), &
        5.d0*gausswidth(1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    call get_hartree(parini,poisson,atoms,gausswidth,ener)
    errmax=maxval(abs(poisson%pot-pot))
    if(errmax<1.d-7) then
        write(output_unit,'(2a)') green_passed,' in test_cal_pot_gauss_s: errmax'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_cal_pot_gauss_s: errmax=  '//achar(27),errmax
        call exit(1)
    end if
    gw=gw+1.d-5
    call cal_pot_gauss_s(parini,poisson,pot_t,.true.,xyz,gw,q,vgrad_t)
    vgrad_t=(pot_t-pot)/1.d-5
    !do iz=poisson%ngpz/2,poisson%ngpz/2
    !do iy=poisson%ngpy/2,poisson%ngpy/2
    !do ix=1,poisson%ngpx
    !    x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
    !    y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
    !    z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
    !    write(83,'(3f8.3,5es19.10)') x,y,z,poisson%pot(ix,iy,iz),pot(ix,iy,iz),poisson%rho(ix,iy,iz),vgrad_t(ix,iy,iz),vgrad(ix,iy,iz)
    !enddo
    !enddo
    !enddo
    errmax_vgrad=maxval(abs(vgrad-vgrad_t))
    if (errmax_vgrad<1.d-6) then
        write(output_unit,'(2a)') green_passed,' in test_cal_pot_gauss_s: errmax_vgrad'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_cal_pot_gauss_s: errmax_vgrad=',errmax_vgrad
        call exit(1)
    end if
    call fini_hartree(parini,atoms,poisson)
    call atom_deallocate_old(atoms)
    deallocate(pot)
    deallocate(vgrad)
end subroutine test_cal_pot_gauss_s
!*****************************************************************************************
