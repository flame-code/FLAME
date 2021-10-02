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
    integer:: istat, igpx, igpy, igpz, iat, ntrial, itrial, nsegx, nsegy, nsegz
    integer:: jgpx, jgpy, jgpz
    real(8):: epot, rgcut_a, qtot, pi!, qtot_e, qtot_i
    real(8):: epot_scn, ehartree_scn_excl, tt1, tt2
    real(8):: xyz(3), dxyz(3), epot_trial, gwt
    real(8):: dx, dy, dz, r2, coeff, rloc, c1, c2
    real(8):: xmin, ymin, zmin, xmax, ymax, zmax
    real(8):: q_one(1), gw_one(1)
    real(8), allocatable::  gausswidth(:)
    real(8), allocatable::  rat_trial(:,:)
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
    poisson_ion%xyz111=poisson%xyz111
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
    !call put_charge_density(parini,poisson_ion)
    nbgpx=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(1,1))+2
    nbgpy=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(2,2))+2
    nbgpz=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(3,3))+2
    poisson_ion%rho=0.d0
    do iat=1,atoms%nat
    !if(trim(atoms%sat(iat))=='Mg') gwt=0.6d0
    !if(trim(atoms%sat(iat))=='O' ) gwt=0.3d0
    gwt=gausswidth(iat)
    coeff=atoms%zat(iat)/(gwt**3*pi**1.5d0)
    !coeff=atoms%zat(iat)*2.d0/(3.d0*gwt**5*pi**1.5d0)
    !coeff=atoms%zat(iat)*4.d0/(15.d0*gwt**7*pi**1.5d0)
    if(trim(atoms%sat(iat))=='Mg') then
        !0.65406138674  2 -5.223929095  0.913704167481045 rloc nloc c1 .. cnloc
        rloc=0.65406138674d0
        c1=-5.223929095d0
        c2=0.913704167481045d0
    endif
    if(trim(atoms%sat(iat))=='O') then
        !0.3454999999    2 -11.7435870154  1.90653967947 rloc nloc c1 .. cnloc
        rloc=0.3454999999d0
        c1=-11.7435870154d0
        c2=1.90653967947d0
    endif
    !c1=-c1
    !c2=-c2
    jgpx=int(poisson_ion%rcart(1,iat)/poisson%hgrid(1,1))
    jgpy=int(poisson_ion%rcart(2,iat)/poisson%hgrid(2,2))
    jgpz=int(poisson_ion%rcart(3,iat)/poisson%hgrid(3,3))
    do igpz=jgpz-nbgpz,jgpz+nbgpz
    do igpy=jgpy-nbgpy,jgpy+nbgpy
    do igpx=jgpx-nbgpx,jgpx+nbgpx
        dx=(igpx-1)*poisson%hgrid(1,1)-poisson_ion%rcart(1,iat)
        dy=(igpy-1)*poisson%hgrid(2,2)-poisson_ion%rcart(2,iat)
        dz=(igpz-1)*poisson%hgrid(3,3)-poisson_ion%rcart(3,iat)
        r2=dx**2+dy**2+dz**2
        if(r2<10.d0**2*gwt**2) then
        !if(r2<10.d0**2*rloc**2) then
            poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*exp(-r2/gwt**2)
            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*r2*exp(-r2/gwt**2)
            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*r2**2*exp(-r2/gwt**2)
            !tt1=exp(-r2/(2.d0*rloc**2))
            !tt2=-3.d0*rloc**4*(c1-2.d0*c2)+rloc**2*(c1-7.d0*c2)*r2+c2*r2**2+rloc**3*sqrt(2.d0/pi)*atoms%zat(iat)
            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+tt1*tt2/(4.d0*rloc**6*pi)
        endif
    enddo
    enddo
    enddo
    enddo
    tt1=0.d0
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        tt1=tt1+poisson_ion%rho(igpx,igpy,igpz)
        !write(33,'(3i5,f8.5)') igpx,igpy,igpz,poisson_ion%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    tt1=tt1*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    write(*,*) 'TT1 ',tt1
    !stop 'WWWWWWWWWWWWWWW'
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
                poisson%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)-poisson%rho(igpx,igpy,igpz)
                !poisson%rho(igpx,igpy,igpz)=-poisson%rho(igpx,igpy,igpz)
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
    !-------------------------------------------------------
    call update_ratp(atoms)
    call get_hartree(parini,poisson,atoms,gausswidth,epot)
    write(*,'(a,es24.15,es14.5)') 'ehartree_scn_excl ',epot,poisson%screening_factor
    atoms%fat=0.d0
    call force_gto_sym_ortho(parini,poisson_ion%bc,atoms%nat,poisson_ion%rcart, &
        poisson_ion%q,gausswidth,6.d0,poisson_ion%xyz111, &
        poisson_ion%ngpx,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
        poisson_ion%hgrid,poisson%pot,atoms%fat)
    do iat=1,atoms%nat
        write(*,'(a,i4,3es19.10)') 'FAT ',iat,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
    enddo
    !-------------------------------------------------------
    poisson_ion%gw(1:poisson_ion%nat)=1.d0
    xmin= huge(1.d0)
    ymin= huge(1.d0)
    zmin= huge(1.d0)
    xmax=-huge(1.d0)
    ymax=-huge(1.d0)
    zmax=-huge(1.d0)
    do iat=1,atoms%nat
        if(poisson_ion%rcart(1,iat)<xmin) xmin=poisson_ion%rcart(1,iat)
        if(poisson_ion%rcart(2,iat)<ymin) ymin=poisson_ion%rcart(2,iat)
        if(poisson_ion%rcart(3,iat)<zmin) zmin=poisson_ion%rcart(3,iat)
        if(poisson_ion%rcart(1,iat)>xmax) xmax=poisson_ion%rcart(1,iat)
        if(poisson_ion%rcart(2,iat)>ymax) ymax=poisson_ion%rcart(2,iat)
        if(poisson_ion%rcart(3,iat)>zmax) zmax=poisson_ion%rcart(3,iat)
    enddo
    dxyz(1)=1.5d0
    dxyz(2)=1.5d0
    dxyz(3)=1.5d0
    nsegx=int((poisson%hgrid(1,1)*poisson%ngpx-16.d0)/dxyz(1))+1
    nsegy=int((poisson%hgrid(2,2)*poisson%ngpy-16.d0)/dxyz(2))+1
    nsegz=int((poisson%hgrid(3,3)*poisson%ngpz-16.d0)/dxyz(3))+1
    dxyz(1)=(poisson%hgrid(1,1)*poisson%ngpx-16.d0)/real(nsegx,kind=8)
    dxyz(2)=(poisson%hgrid(2,2)*poisson%ngpy-16.d0)/real(nsegy,kind=8)
    dxyz(3)=(poisson%hgrid(3,3)*poisson%ngpz-16.d0)/real(nsegz,kind=8)
    !nseg=10
    !ntrial=(nseg+1)**3
    ntrial=(nsegx+1)*(nsegy+1)*(nsegz+1)
    write(*,'(a,3i3,i6)') 'nsegx,nsegy,nsegz,ntrial ',nsegx,nsegy,nsegz,ntrial
    !ntrial=atoms%nat
    allocate(rat_trial(3,ntrial))
    itrial=0
    do iz=0,nsegz
    do iy=0,nsegy
    do ix=0,nsegx
        itrial=itrial+1
        rat_trial(1,itrial)=8.d0+dxyz(1)*ix
        rat_trial(2,itrial)=8.d0+dxyz(2)*iy
        rat_trial(3,itrial)=8.d0+dxyz(3)*iz
    enddo
    enddo
    enddo
    !itrial=0
    !do iat=1,atoms%nat
    !    itrial=itrial+1
    !    rat_trial(1,itrial)=poisson_ion%rcart(1,iat)
    !    rat_trial(2,itrial)=poisson_ion%rcart(2,iat)
    !    rat_trial(3,itrial)=poisson_ion%rcart(3,iat)
    !enddo
    do itrial=1,ntrial
        xyz(1)=rat_trial(1,itrial)-poisson_ion%rcart(1,1)
        xyz(2)=rat_trial(2,itrial)-poisson_ion%rcart(2,1)
        xyz(3)=rat_trial(3,itrial)-poisson_ion%rcart(3,1)
        !    put_gto_sym_ortho(parini,bc,reset,nat,rxyz,qat,gw,rgcut,xyz111,ngx,ngy,ngz,hgrid,rho)
        q_one(1)=1.d0
        gw_one(1)=1.d0
        call put_gto_sym_ortho(parini,poisson_ion%bc,.true.,1,rat_trial(1,itrial),q_one,gw_one, &
            6.d0,poisson_ion%xyz111,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz,poisson_ion%hgrid,poisson_ion%rho)
        epot_trial=0.d0
        do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
        do igpx=1,poisson%ngpx
            epot_trial=epot_trial+poisson_ion%rho(igpx,igpy,igpz)*poisson%pot(igpx,igpy,igpz)
        enddo
        enddo
        enddo
        epot_trial=epot_trial*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
        !if(ix==0 .and. iy==0 .and. iz==0) then
        !    epot_trial0=epot_trial
        !    !write(*,'(a,i5,es24.15,es14.5)') 'iat,epot_trial ',iat,epot_trial,poisson%screening_factor
        !else
        !    write(*,'(a,i5,4es24.15,es14.5)') 'iat,epot_trial ',iat,dxyz(1),dxyz(2),dxyz(3),epot_trial-epot_trial0,poisson%screening_factor
            write(*,'(a,i5,4es24.15,es14.5)') 'iat,epot_trial ',1,xyz(1),xyz(2),xyz(3),epot_trial,poisson%screening_factor
            write(71,'(a,i3,4(a2,es24.15),a,es14.5,a)') '  - [',1,', ',xyz(1),', ',xyz(2),', ',xyz(3),', ',epot_trial,', ',poisson%screening_factor,']'
        !endif
    enddo
    write(*,'(a,6f8.1)') 'MINMAX ',xmin,ymin,zmin,xmax,ymax,zmax
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(1,1)*poisson%ngpx
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(2,2)*poisson%ngpy
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(3,3)*poisson%ngpz
    !-------------------------------------------------------
    !call cube_write('total_rho.cube',atoms,poisson,'rho')
    !call cube_write('total_pot.cube',atoms,poisson,'pot')
    call fini_hartree(parini,atoms,poisson)
    call fini_hartree(parini,atoms,poisson_ion)
end subroutine solve_poisson_cube_bigdft
!*****************************************************************************************
