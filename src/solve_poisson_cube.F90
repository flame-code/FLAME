!*****************************************************************************************
subroutine solve_poisson(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, get_rat_iat
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    if(trim(parini%psolver)=='p3d') then
        call solve_poisson_cube_p3d(parini)
    elseif(trim(parini%psolver)=='bigdft') then
        call solve_poisson_cube_bigdft(parini)
    else
        write(*,*) 'ERROR: Use this subtask when you have the following:'
        write(*,*) '                            (1) bigdft with free BC'
        write(*,*) '                            (2) P3D with slab BC'
        stop
    endif
end subroutine solve_poisson
!*****************************************************************************************
subroutine solve_poisson_cube_p3d(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson
    type(typ_poisson):: poisson_ion
    type(typ_atoms):: atoms
    integer:: istat, igpx, igpy, igpz, iat
    real(8):: epot, rgcut_a, t1, t2, t3, t4, pi
    real(8),allocatable::  gausswidth(:)
    integer:: nbgpx, nbgpy, nbgpz
    pi=4.d0*atan(1.d0)
    call cube_read('rho.cube',atoms,poisson)
    atoms%boundcond='slab'
    allocate(gausswidth(atoms%nat))
    poisson_ion%alpha=0.4d0 !atoms%rcov(iat)
    gausswidth=0.4d0 !atoms%rcov(iat)
    poisson%task_finit=""
    call init_hartree(parini,atoms,poisson,gausswidth)
    poisson%cell(1)=poisson%hgrid(1,1)*poisson%ngpx
    poisson%cell(2)=poisson%hgrid(2,2)*poisson%ngpy
    poisson%cell(3)=poisson%hgrid(3,3)*poisson%ngpz
    !-------------------------------------------------------
    poisson_ion%ngpx=poisson%ngpx
    poisson_ion%ngpy=poisson%ngpy
    poisson_ion%ngpz=poisson%ngpz
    poisson_ion%hgrid(1,1)=poisson%hgrid(1,1)
    poisson_ion%hgrid(2,2)=poisson%hgrid(2,2)
    poisson_ion%hgrid(3,3)=poisson%hgrid(3,3)
    if(.not. parini%gaussian_width>0.d0) then
        stop 'ERROR: gaussian_width must be set.'
    endif
    rgcut_a=parini%gaussian_width !3.d0
    nbgpx=int(rgcut_a/poisson_ion%hgrid(1,1))+2
    nbgpy=int(rgcut_a/poisson_ion%hgrid(2,2))+2
    nbgpz=int(rgcut_a/poisson_ion%hgrid(3,3))+2
    poisson_ion%ngpz=poisson_ion%ngpz+2*nbgpz
    write(*,'(a,3i5)') 'ngpx, ngpy, ngpz  ',poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz
    write(*,'(a,3i5)') 'nbgpx,nbgpy,nbgpz ',nbgpx,nbgpy,nbgpz
    poisson_ion%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson_ion,gausswidth)
    poisson_ion%reset_rho=.true.
    poisson_ion%nat=atoms%nat
    poisson_ion%cv=atoms%cellvec
    poisson_ion%bc=atoms%boundcond
    poisson_ion%q(1:poisson_ion%nat)=atoms%zat(1:atoms%nat)
    poisson_ion%gw(1:poisson_ion%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson_ion%rcart)
    call put_charge_density(parini,poisson_ion)
    t1=0.d0
    t2=0.d0
    t3=0.d0
    t4=0.d0
    do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
            do igpx=1,poisson%ngpx
                t1=t1+poisson%rho(igpx,igpy,igpz)
                t2=t2+poisson_ion%rho(igpx,igpy,igpz+nbgpz)
                poisson%rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)-poisson_ion%rho(igpx,igpy,igpz+nbgpz)
                t4=t4+(igpz-1)*poisson_ion%hgrid(3,3)*poisson%rho(igpx,igpy,igpz)
                t3=t3+poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    t1=t1*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    t2=t2*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    t3=t3*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    t4=t4*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    write(*,*) 't1=',t1
    write(*,*) 't2=',t2
    write(*,*) 't3=',t3
    write(*,*) 't4=',t4
    call fini_hartree(parini,atoms,poisson_ion)
    call cube_write('total_rho.cube',atoms,poisson,'rho')
    !-------------------------------------------------------
    call update_ratp(atoms)
    t1=0.d0
    do iat=1,atoms%nat
        t1=t1+atoms%zat(iat)*atoms%ratp(3,iat)
    enddo
    !t1=t1*(2*pi)/(cell(1)*cell(2))
    write(*,*) 't1=',t1
    if(parini%ewald) then
        write(*,*) 'ERROR: ewald=True is wrong when reading from cube file.'
        stop
    endif
    if(trim(parini%psolver)=='kwald') then
        write(*,*) 'ERROR: psolver=kwald is wrong for grid base charge density.'
        stop
    endif
    call get_hartree(parini,poisson,atoms,gausswidth,epot)
    call cube_write('pot_p3d.cube',atoms,poisson,'pot')
    call fini_hartree(parini,atoms,poisson)
end subroutine solve_poisson_cube_p3d
!*****************************************************************************************
subroutine solve_poisson_cube_bigdft(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp!, get_rat_iat
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson
    type(typ_poisson):: poisson_scn
    type(typ_poisson):: poisson_ion
    type(typ_atoms):: atoms
    integer:: istat, igpx, igpy, igpz, iat
    real(8):: epot, rgcut_a, qtot, pi!, qtot_e, qtot_i
    real(8):: epot_scn, ehartree_scn_excl, tt1, tt2
!    real(8):: xyz(3), dxyz(3)
    real(8),allocatable::  gausswidth(:)
    integer:: nbgpx, nbgpy, nbgpz, ix, iy, iz
    !real(8) :: xyz(3)
    !integer:: ny,nz
    !real(8):: max_ion, max_ele
    pi=4.d0*atan(1.d0)
    if(.not. parini%gaussian_width>0.d0) then
        stop 'ERROR: gaussian_width must be set.'
    endif
    if(parini%ewald) then
        write(*,*) 'ERROR: ewald=True is wrong when reading from cube file.'
        stop
    endif
    if(trim(parini%psolver)=='kwald') then
        write(*,*) 'ERROR: psolver=kwald is wrong for grid base charge density.'
        stop
    endif
    if(parini%cal_scn .and. parini%screening_factor==0.d0) then
        stop 'ERROR: cal_scn is TRUE and screening_factor is 0, MEANINGLESS!'
    endif
    !-------------------------------------------------------
    call cube_read('rho.cube',atoms,poisson)
    atoms%boundcond='free'
    allocate(gausswidth(atoms%nat))
    gausswidth=parini%gaussian_width
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    poisson%task_finit=""
    call init_hartree(parini,atoms,poisson,gausswidth)
    poisson%cell(1)=poisson%hgrid(1,1)*poisson%ngpx
    poisson%cell(2)=poisson%hgrid(2,2)*poisson%ngpy
    poisson%cell(3)=poisson%hgrid(3,3)*poisson%ngpz
    !-------------------------------------------------------
    poisson_ion%alpha=parini%gaussian_width
    poisson_ion%rgcut=parini%rgcut_ewald*poisson_ion%alpha
    poisson_ion%ngpx=poisson%ngpx
    poisson_ion%ngpy=poisson%ngpy
    poisson_ion%ngpz=poisson%ngpz
    poisson_ion%hgrid(1,1)=poisson%hgrid(1,1)
    poisson_ion%hgrid(2,2)=poisson%hgrid(2,2)
    poisson_ion%hgrid(3,3)=poisson%hgrid(3,3)
    rgcut_a=parini%rgcut_ewald*maxval(gausswidth) !parini%gaussian_width !3.d0
    nbgpx=int(rgcut_a/poisson_ion%hgrid(1,1))+2
    nbgpy=int(rgcut_a/poisson_ion%hgrid(2,2))+2
    nbgpz=int(rgcut_a/poisson_ion%hgrid(3,3))+2
    poisson_ion%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson_ion,gausswidth)
    poisson_ion%reset_rho=.true.
    poisson_ion%nat=atoms%nat
    poisson_ion%cv=atoms%cellvec
    poisson_ion%bc=atoms%boundcond
    poisson_ion%q(1:poisson_ion%nat)=atoms%zat(1:atoms%nat)
    poisson_ion%gw(1:poisson_ion%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson_ion%rcart)
    call put_charge_density(parini,poisson_ion)
    !-------------------------------------------------------
    !call get_rat_iat(atoms,1,xyz)
    !ny=int(xyz(2)/poisson%hgrid(2,2)) 
    !nz=int(xyz(3)/poisson%hgrid(3,3)) 
    !max_ion=maxval(poisson_ion%rho(:,:,:))
    !max_ele=maxval(poisson%rho(:,:,:))
    !do igpz=1,poisson%ngpz
    !    do igpy=1,poisson%ngpy
    !        do igpx=1,poisson%ngpx
    !            if (poisson_ion%rho(igpx,igpy,igpz)==max_ion) write(*,*) 'ind max_ion', igpx,igpy,igpz
    !            if (poisson%rho(igpx,igpy,igpz)==max_ele) then
    !                write(*,*) 'ind max_ele', igpx,igpy,igpz
    !                ny = igpy
    !                nz = igpz
    !            end if
    !        end do
    !    end do
    !end do
    !do igpx=1,poisson%ngpx
    !    write(1370,'(3es14.6)') igpx*poisson%hgrid(1,1),poisson%rho(igpx,ny,nz),poisson_ion%rho(igpx,ny,nz)
    !end do
    !-------------------------------------------------------
    qtot=0.d0
    !qtot_e=0.d0
    !qtot_i=0.d0
    do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
            do igpx=1,poisson%ngpx
    !            qtot_e=qtot_e+poisson%rho(igpx,igpy,igpz)
    !            qtot_i=qtot_i+poisson_ion%rho(igpx,igpy,igpz)
                poisson%rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)-poisson_ion%rho(igpx,igpy,igpz)
                qtot=qtot+poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    !do igpx=1,poisson%ngpx
    !    write(1371,'(2es14.6)') igpx*poisson%hgrid(1,1),poisson%rho(igpx,ny,nz)
    !end do
    qtot=qtot*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    write(*,*) 'qtot= ',qtot
    !write(*,*) 'qtot_e= ',qtot_e
    !write(*,*) 'qtot_i= ',qtot_i
    call fini_hartree(parini,atoms,poisson_ion)
    !-------------------------------------------------------
    call update_ratp(atoms)
    call get_hartree(parini,poisson,atoms,gausswidth,epot)
!    write(*,'(a,es24.15,es14.5)') 'ehartree_scn_excl ',epot,poisson%screening_factor
!    !-------------------------------------------------------
!    do iat=1,atoms%nat
!        do iz=0,6
!        do iy=0,6
!        do ix=0,6
!        dxyz(1)=ix*0.3d0
!        dxyz(2)=iy*0.3d0
!        dxyz(3)=iz*0.3d0
!        if(ix>3) dxyz(1)=(ix-7)*0.3d0
!        if(iy>3) dxyz(2)=(iy-7)*0.3d0
!        if(iz>3) dxyz(3)=(iz-7)*0.3d0
!        poisson_ion%q(1:poisson_ion%nat)=0.d0
!        poisson_ion%q(iat)=1.d0
!        poisson_ion%rho=0.d0
!        xyz(1:3)=poisson_ion%rcart(1:3,iat)
!        poisson_ion%rcart(1,iat)=xyz(1)+dxyz(1)
!        poisson_ion%rcart(2,iat)=xyz(2)+dxyz(2)
!        poisson_ion%rcart(3,iat)=xyz(3)+dxyz(3)
!        call put_charge_density(parini,poisson_ion)
!        poisson_ion%rcart(1:3,iat)=xyz(1:3)
!        epot_trial=0.d0
!        do igpz=1,poisson%ngpz
!        do igpy=1,poisson%ngpy
!        do igpx=1,poisson%ngpx
!            epot_trial=epot_trial+poisson_ion%rho(igpx,igpy,igpz)*poisson%pot(igpx,igpy,igpz)
!        enddo
!        enddo
!        enddo
!        epot_trial=epot_trial*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
!        !if(ix==0 .and. iy==0 .and. iz==0) then
!        !    epot_trial0=epot_trial
!        !    !write(*,'(a,i5,es24.15,es14.5)') 'iat,epot_trial ',iat,epot_trial,poisson%screening_factor
!        !else
!        !    write(*,'(a,i5,4es24.15,es14.5)') 'iat,epot_trial ',iat,dxyz(1),dxyz(2),dxyz(3),epot_trial-epot_trial0,poisson%screening_factor
!            write(*,'(a,i5,4es24.15,es14.5)') 'iat,epot_trial ',iat,dxyz(1),dxyz(2),dxyz(3),epot_trial,poisson%screening_factor
!        !endif
!        enddo
!        enddo
!        enddo
!    enddo
    !-------------------------------------------------------
    write(*,'(a,es24.15,es14.5)') 'ehartree_scn_excl ',epot,poisson%screening_factor
    !call cube_write('total_rho.cube',atoms,poisson,'rho')
    !call cube_write('total_pot.cube',atoms,poisson,'pot')
    call fini_hartree(parini,atoms,poisson)
end subroutine solve_poisson_cube_bigdft
!*****************************************************************************************
