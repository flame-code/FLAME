!*****************************************************************************************
subroutine test_cent2_analytic
    use iso_fortran_env, only: error_unit, output_unit
    use mod_fit_bf_cent2, only: cal_pot_gauss_s, cal_pot_gauss_p, cal_pot_r2gauss_s
    use mod_cent2, only: cent2_analytic
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_rat, atom_allocate, atom_deallocate_old
    use mod_electrostatics, only: typ_poisson
    use mod_colors, only: green_passed, red_failed
    implicit none
    integer:: nat
    real(8):: gw, qr0(2), pat(3,2), qr2(2), ener, gausswidth(2), gwr0(2), gwp1(2), gwr2(2), errmax
    real(8):: errmax_vgrad, voxel
    !logical::
    type(typ_parini):: parini
    type(typ_poisson):: poisson
    type(typ_atoms):: atoms
    real(8), allocatable:: vgrad(:,:,:)
    integer:: ix, iy, iz
    real(8):: x, y, z
    nat=2
    call atom_allocate(atoms,nat,0,0)
    parini%screening_factor=0.3d0
    parini%cal_scn=.true.
    parini%psolver='bigdft'
    poisson%bc='free'
    atoms%ratp(1,1)=12.5d0 ; atoms%ratp(2,1)=12.5d0 ; atoms%ratp(3,1)=12.5d0
    atoms%ratp(1,2)=17.5d0 ; atoms%ratp(2,2)=12.5d0 ; atoms%ratp(3,2)=12.5d0
    call update_rat(atoms)
    poisson%task_finit="alloc_rho"
    poisson%ngpx=60
    poisson%ngpy=50
    poisson%ngpz=50
    poisson%hgrid(1:3,1:3)=0.d0
    poisson%hgrid(1,1)=0.5d0
    poisson%hgrid(2,2)=0.5d0
    poisson%hgrid(3,3)=0.5d0
    poisson%xyz111(1:3)=0.d0
    gwr0(1)=1.5d0
    gwr0(2)=1.7d0
    gwp1(1)=1.6d0
    gwp1(2)=1.4d0
    gwr2(1)=1.2d0
    gwr2(2)=1.3d0
    qr0(1)=-1.2d0
    qr0(2)= 1.4d0
    pat(1,1)=0.9d0 ; pat(2,1)=0.8d0 ; pat(3,1)=1.3d0
    pat(1,2)=0.4d0 ; pat(2,2)=0.6d0 ; pat(3,2)=0.5d0
    qr2(1)=-1.1d0
    qr2(2)= 0.9d0
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    atoms%boundcond='free'
    call init_hartree(parini,atoms,poisson,gausswidth)
    allocate(vgrad(poisson%ngpx,poisson%ngpy,poisson%ngpz))

    poisson%pot=0.d0
    call cal_pot_gauss_s(parini,poisson,poisson%pot,.false.,atoms%ratp(1,1),gwr0(1),qr0(1),vgrad)
    call cal_pot_gauss_s(parini,poisson,poisson%pot,.false.,atoms%ratp(1,2),gwr0(2),qr0(2),vgrad)
    call cal_pot_gauss_p(parini,poisson,poisson%pot,.false.,atoms%ratp(1,1),gwp1(1),pat(1,1),vgrad)
    call cal_pot_gauss_p(parini,poisson,poisson%pot,.false.,atoms%ratp(1,2),gwp1(2),pat(1,2),vgrad)
    call cal_pot_r2gauss_s(parini,poisson,poisson%pot,.false.,atoms%ratp(1,1),gwr2(1),qr2(1),vgrad)
    call cal_pot_r2gauss_s(parini,poisson,poisson%pot,.false.,atoms%ratp(1,2),gwr2(2),qr2(2),vgrad)

    poisson%rho=0.d0
    call put_gto_sym_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,qr0,gwr0, &
        5.d0*maxval(gwr0),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    call put_gto_p_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,pat,gwp1, &
        5.d0*maxval(gwp1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    call put_r2gto_sym_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,qr2,gwr2, &
        5.d0*maxval(gwr2),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)

    ener=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        ener=ener+poisson%rho(ix,iy,iz)*poisson%pot(ix,iy,iz)
    enddo
    enddo
    enddo
    voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    ener=ener*0.5d0*voxel
    call cent2_analytic(parini,atoms,gwr0,gwp1,gwr2,qr0,pat,qr2)
    write(*,'(a,2f20.10)') 'epot,ener= ',atoms%epot,ener
    errmax=abs(atoms%epot-ener)
    if(errmax<1.d-11) then
        write(output_unit,'(2a)') green_passed,' in test_cent2_analytic: errmax'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_cent2_analytic: errmax=',errmax
        call exit(1)
    end if
    call fini_hartree(parini,atoms,poisson)
    call atom_deallocate_old(atoms)
    deallocate(vgrad)
end subroutine test_cent2_analytic
!*****************************************************************************************
