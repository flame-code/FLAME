!*****************************************************************************************
subroutine test_optimize_basis_functions
    use iso_fortran_env, only: error_unit, output_unit
    use mod_fit_bf_cent2, only: cal_pot_gauss_s, cal_pot_gauss_p, cal_pot_r2gauss_s
    use mod_fit_bf_cent2, only: typ_fitpar, optimize_basis_functions
    use mod_cent2, only: cent2_analytic
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, update_rat, atom_allocate_old, atom_deallocate_old
    use mod_electrostatics, only: typ_poisson
    use mod_colors, only: green_passed, red_failed
    use mod_flm_futile
    implicit none
    integer:: nat, iat, itypat
    real(8):: gw, ener, errmax
    real(8):: errmax_vgrad, voxel, qtot
    !logical::
    type(typ_parini):: parini
    type(typ_poisson), allocatable:: poisson_ref(:)
    real(8), allocatable:: gausswidth(:)
    real(8), allocatable:: qr0(:), qr2(:), pat(:,:), gwr0(:), gwr2(:), gwp1(:)
    type(typ_atoms_arr):: atoms_arr
    type(typ_fitpar):: fitpar
    !real(8), allocatable:: vgrad(:,:,:)
    integer:: ix, iy, iz, nfiles, iconf
    real(8):: x, y, z
    nfiles=1
    atoms_arr%nconf=nfiles
    allocate(atoms_arr%atoms(atoms_arr%nconf))
    allocate(poisson_ref(nfiles))
    call atom_allocate_old(atoms_arr%atoms(1),4,0,0)
    associate(atoms=>atoms_arr%atoms(1),poisson=>poisson_ref(1))
    atoms%boundcond='free'
    atoms%ratp(1,1)= 9.d0 ; atoms%ratp(2,1)=9.d0 ; atoms%ratp(3,1)=9.d0
    atoms%ratp(1,2)=17.d0 ; atoms%ratp(2,2)=9.d0 ; atoms%ratp(3,2)=9.d0
    atoms%ratp(1,3)=13.d0 ; atoms%ratp(2,3)=9.d0 ; atoms%ratp(3,3)=9.d0
    atoms%ratp(1,4)=21.d0 ; atoms%ratp(2,4)=9.d0 ; atoms%ratp(3,4)=9.d0
    atoms%sat(1)='Li'
    atoms%sat(2)='Li'
    atoms%sat(3)='S'
    atoms%sat(4)='S'
    call update_rat(atoms)
    allocate(gausswidth(atoms%nat))
    parini%types_main='Li S'
    call set_atomc_types_info(parini)
    do iconf=1,nfiles
        do iat=1,atoms_arr%atoms(iconf)%nat
            do itypat=1,parini%ntypat
                if(trim(atoms_arr%atoms(iconf)%sat(iat))==trim(parini%stypat(itypat))) then
                    atoms_arr%atoms(iconf)%itypat(iat)=parini%ltypat(itypat)
                    exit
                endif
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    parini%iverbose=0
    parini%screening_factor=0.25d0
    parini%cal_scn=.true.
    parini%psolver='bigdft'
    poisson%bc='free'
    poisson%task_finit="alloc_rho"
    poisson%ngpx=280
    poisson%ngpy=80
    poisson%ngpz=80
    poisson%hgrid(1:3,1:3)=0.d0
    poisson%hgrid(1,1)=0.25d0
    poisson%hgrid(2,2)=0.25d0
    poisson%hgrid(3,3)=0.25d0
    poisson%xyz111(1:3)=0.d0
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    poisson%rho=f_malloc([1.to.poisson%ngpx,1.to.poisson%ngpy,1.to.poisson%ngpz], &
        id='poisson%rho')
    allocate(qr0(atoms%nat),source=0.d0)
    allocate(qr2(atoms%nat),source=0.d0)
    allocate(pat(3,atoms%nat),source=0.d0)
    atoms%qat(1)=2.8d0
    atoms%qat(2)=2.9d0
    atoms%qat(3)=7.0d0
    atoms%qat(4)=6.1d0
    qr0(1)=0.6d0*atoms%qat(1)  ; qr2(1)=0.4d0*atoms%qat(1)  ; pat(1,1)= 0.1d0
    qr0(2)=0.6d0*atoms%qat(2)  ; qr2(2)=0.4d0*atoms%qat(2)  ; pat(1,2)= 0.2d0
    qr0(3)=0.4d0*atoms%qat(3)  ; qr2(3)=0.6d0*atoms%qat(3)  ; pat(1,3)=-0.2d0
    qr0(4)=0.4d0*atoms%qat(4)  ; qr2(4)=0.6d0*atoms%qat(4)  ; pat(1,4)= 0.3d0
    allocate(gwr0(atoms%nat),source=1.d0)
    allocate(gwr2(atoms%nat),source=1.d0)
    allocate(gwp1(atoms%nat),source=1.d0)
    do iconf=1,nfiles
        do iat=1,atoms_arr%atoms(iconf)%nat
            if(trim(atoms_arr%atoms(iconf)%sat(iat))=='Li') then
                gwr0(iat)=1.1d0
                gwr2(iat)=1.1d0
            elseif(trim(atoms_arr%atoms(iconf)%sat(iat))=='S') then
                gwr0(iat)=0.9d0
                gwr2(iat)=0.9d0
            endif
        enddo
    enddo
    poisson%rho=0.d0
    !do iat=1,atoms%nat
    call put_gto_sym_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,qr0,gwr0, &
        5.d0*maxval(gwr0),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    !call put_gto_p_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,pat,gwp1, &
    !    5.d0*maxval(gwp1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    call put_r2gto_sym_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,qr2,gwr2, &
        5.d0*maxval(gwr2),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    !enddo
    !qtot=0.d0
    !do iz=1,poisson%ngpz
    !do iy=1,poisson%ngpy
    !do ix=1,poisson%ngpx
    !    qtot=qtot+poisson%rho(ix,iy,iz)
    !enddo
    !enddo
    !enddo
    !voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    !qtot=qtot*voxel
    !write(*,'(a,es24.15)') 'total electronic charge=',qtot
    call fitpar%init_fit_bf(parini%ntypat)
    fitpar%gwz_s=1.d0
    fitpar%bz_s=1.d0
    fitpar%relaxcore=.true.
    fitpar%applycore=.true.
    fitpar%qcore_type(1)=1.d0
    fitpar%qcore_type(2)=2.d0
    fitpar%gwc_s1(1)=1.0d0
    fitpar%gwc_s1(2)=1.0d0
    fitpar%bc_s1(1)=1.d0
    fitpar%bc_s1(2)=1.d0
    fitpar%gwv_s1(1)=1.059d0
    fitpar%gwv_s1(2)=1.029d0
    fitpar%bv_s1(1)=-1.213d0
    fitpar%bv_s1(2)=0.55d0
    fitpar%gwv_p1=1.3d0
    fitpar%gwv_p2=1.0d0

    fitpar%gwv_p1(1)=fitpar%gwv_s1(1)
    fitpar%gwv_p1(2)=fitpar%gwv_s1(2)

    parini%paropt_geopt%nit=1
    parini%paropt_geopt%lprint=.true.
    parini%paropt_geopt%alphax=1.d-2
    parini%paropt_geopt%fmaxtol=1.d-5
    parini%paropt_geopt%condnum=100.d0
    call optimize_basis_functions(parini,fitpar,nfiles,atoms_arr,poisson_ref,4.0d0)
    if(parini%iverbose>0) then
        write(*,'(es24.15)') fitpar%gwc_s1(1)
        write(*,'(es24.15)') fitpar%gwc_s1(2)
        write(*,'(es24.15)') fitpar%bc_s1(1)
        write(*,'(es24.15)') fitpar%bc_s1(2)
        write(*,'(es24.15)') fitpar%gwv_s1(1)
        write(*,'(es24.15)') fitpar%gwv_s1(2)
        write(*,'(es24.15)') fitpar%bv_s1(1)
        write(*,'(es24.15)') fitpar%bv_s1(2)
    endif
    errmax=0.d0
    errmax=max(errmax,abs(fitpar%gwc_s1(1)-( 1.001691124592907d+00)))
    errmax=max(errmax,abs(fitpar%gwc_s1(2)-( 1.001782021797540d+00)))
    errmax=max(errmax,abs(fitpar%bc_s1(1) -( 9.994362918023644d-01)))
    errmax=max(errmax,abs(fitpar%bc_s1(2) -( 9.994059927341535d-01)))
    errmax=max(errmax,abs(fitpar%gwv_s1(1)-( 1.055574096465681d+00)))
    errmax=max(errmax,abs(fitpar%gwv_s1(2)-( 1.025453475513902d+00)))
    errmax=max(errmax,abs(fitpar%bv_s1(1) -(-1.212473511885337d+00)))
    errmax=max(errmax,abs(fitpar%bv_s1(2) -( 5.509507082974783d-01)))
    call fitpar%fini_fit_bf()
    deallocate(gausswidth)
    deallocate(gwr0,gwp1,gwr2,qr0,qr2,pat)
    end associate
    do iconf=1,nfiles
        call fini_hartree(parini,atoms_arr%atoms(iconf),poisson_ref(iconf))
    enddo
    deallocate(poisson_ref)
    do iconf=1,atoms_arr%nconf
        call atom_deallocate_old(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
    if(errmax<1.d-11) then
        write(output_unit,'(2a)') green_passed,' in test_optimize_basis_functions: errmax'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_optimize_basis_functions: errmax=',errmax
        call exit(1)
    end if
end subroutine test_optimize_basis_functions
!*****************************************************************************************
