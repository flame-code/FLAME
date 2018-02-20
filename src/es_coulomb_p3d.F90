!*****************************************************************************************
!CAUTION: Why is this subroutine called in ANN cent1 when system is bulk
!and kwald is the psolver.
subroutine construct_poisson(parini,atoms,poisson)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson):: poisson_rough
    type(typ_poisson), intent(inout):: poisson
    !local variables
    include 'fftw3.f'
    real(8):: pi, shortrange_at_rcut
    call f_routine(id='construct_poisson')
    pi=4.d0*atan(1.d0)
    poisson_rough%hx=parini%hx_ewald
    poisson_rough%hy=parini%hy_ewald
    poisson_rough%hz=parini%hz_ewald
    if (parini%ewald .and. parini%alpha_ewald>0.d0) then
        poisson%alpha=parini%alpha_ewald
    else if (poisson%alpha< 0.d0 .and. parini%alpha_ewald<= 0.d0) then
            write(*,*) "ERROR : alpha is undefined"
            stop
    endif
    poisson%linked_lists%rcut=parini%rcut_ewald
    poisson_rough%rgcut=parini%rgcut_ewald*poisson%alpha
    poisson%spline%nsp=parini%nsp_ewald
    poisson%cell(1)=atoms%cellvec(1,1)
    poisson%cell(2)=atoms%cellvec(2,2)
    poisson%cell(3)=atoms%cellvec(3,3)
    poisson%vu=parini%vu_ewald
    poisson%vl=parini%vl_ewald
    call calparam(parini,atoms,poisson_rough,poisson)
    associate(ngpx=>poisson%ngpx)
    associate(ngpy=>poisson%ngpy)
    associate(ngpz=>poisson%ngpz)
    associate(nbgpy=>poisson%nbgpy)
    associate(nbgpz=>poisson%nbgpz)
    if(trim(atoms%boundcond)=='bulk') then
        poisson%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
            id='poisson%rho')
        poisson%pot=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
            id='poisson%pot')
    elseif(trim(atoms%boundcond)=='slab') then
        poisson%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
            id='poisson%rho')
        poisson%pot=f_malloc([1.to.ngpx+2,1.to.ngpy,1.to.poisson%ngpztot], &
            id='poisson%pot')
    else
        write(*,*) 'ERROR: other BCs are not yet considered.'
    endif
    poisson%mboundg=f_malloc([1.to.2,-nbgpy.to.nbgpy,-nbgpz.to.nbgpz],id='poisson%mboundg')
    if(trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)=='bigdft') then
            call construct_ewald_bps(parini,atoms,poisson)
        endif
    elseif(trim(atoms%boundcond)=='slab') then
        call ps2dp1df_construction(poisson)
    endif
    call determine_glimitsphere(poisson)
    poisson%epotfixed=dot_product(atoms%qat,atoms%qat)/(sqrt(2.d0*pi)*poisson%alpha)
    shortrange_at_rcut=erfc(poisson%linked_lists%rcut/(sqrt(2.d0)*poisson%alpha))/(poisson%linked_lists%rcut)
    if(parini%iverbose>=2) then
        write(*,*) 'real part in rcut',shortrange_at_rcut
    endif
    end associate
    end associate
    end associate
    end associate
    end associate
    call f_release_routine()
end subroutine construct_poisson
!*****************************************************************************************
subroutine destruct_poisson(parini,atoms,poisson)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    call f_routine(id='destruct_poisson')
    if(trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)=='bigdft') then
            call destruct_ewald_bps(poisson)
        endif
    elseif(trim(atoms%boundcond)=='slab') then
        call ps2dp1df_destruction(poisson)
    endif
    call f_free(poisson%rho)
    call f_free(poisson%pot)
    if(trim(parini%bias_type)=='p3dbias') then
     !   deallocate(poisson%pots)
    endif
    call f_free(poisson%mboundg)
    call f_release_routine()
end subroutine destruct_poisson
!*****************************************************************************************
subroutine calculate_forces_energy(parini,poisson,atoms)
    use mod_interface
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    !local variables
    real(8):: epotlong, epotplane !, epotshort
    real(8):: time(10)
    !The following dummy varables definition to be deleted later.
    !real(8):: totrho
    real(8):: beta, pi, charge,c,charge0,E,tmp
    integer:: igpx, igpy, igpz, iat
    integer:: ix, iy, iz, jx, jy, jz, kx, ky, kz
    integer:: npl, npu, nlayer, ngpx, ngpy, ngpz
    real(8),allocatable :: pots_layer(:,:,:,:)
    real(8):: vl, vu, A, d, rl, ru, dipole_correction, dipole
    real(8),allocatable :: gausswidth(:)  
    call f_routine(id='calculate_forces_energy')
    ngpz=poisson%ngpz
    ngpy=poisson%ngpy
    ngpx=poisson%ngpx
    poisson%point_particle= .true.

    pi=4.d0*atan(1.d0)
    beta=0.d0
    do iat=1,atoms%nat
        beta=beta+atoms%qat(iat)*atoms%rat(3,iat)
    enddo
    beta=beta*2.d0*pi*poisson%ngpx*poisson%ngpy/(poisson%cell(1)*poisson%cell(2))
    gausswidth=f_malloc([1.to.atoms%nat],id='gausswidth')
    gausswidth(:)=poisson%alpha

    !write(*,*) 'total momentum z component',beta
    !write(*,*) 'total momentum z component',0.13074051987178871d5/beta
    call cpu_time(time(1))
    call putgaussgrid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%qat,gausswidth,poisson)
    call cpu_time(time(2))
    !-----------------------------------------------------------------------
    !totrho=0.d0
    !do iz=1,ngpz;do iy=1,ngpy;do ix=1,poisson%ngpx
    !    totrho=totrho+rho(ix,iy,iz)
    !enddo;enddo;enddo
    !write(*,*) 'totrho',totrho
    !-----------------------------------------------------------------------
    !pot=0.d0
    call solve_poisson_slab_p3d(parini,poisson,poisson%cell,poisson%hx,poisson%hy,poisson%hz,epotlong,beta)
    call cpu_time(time(3))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call longerange_forces(parini,atoms,poisson,gausswidth)
    call cpu_time(time(4))
    !call shortenergy(atoms,poisson%linked_lists,poisson%spline,poisson%alpha,poisson%cell,epotshort)
    call cpu_time(time(5))
    epotplane=0.d0
    if(trim(parini%bias_type)=='p3dbias') then
        call bias_potener_forces(parini,poisson,atoms,epotplane) 
    end if

    if(trim(parini%bias_type)=='fixed_efield' .or. trim(parini%bias_type)=='fixed_potdiff') then
        call bias_field_potener_forces(parini,poisson,atoms,epotplane) 
    endif
    call cpu_time(time(6))
    !atoms%epot=epotlong+epotshort-poisson%epotfixed+epotplane
    atoms%epot=epotlong+-poisson%epotfixed+epotplane
    write(*,*) '-----------------------------------------------------------'
    write(*,'(a50,e32.15)') 'epotfixed',poisson%epotfixed
    write(*,'(a50,e32.15)') 'epotlong',epotlong
    write(*,'(a50,e32.15)') 'epotplane',epotplane
    write(*,'(a50,e32.15)') 'epottotal',atoms%epot
    !write(*,*) '-----------------------------------------------------------'
!    write(*,'(a50,f32.15)') 'Time for putgaussgrid ',time(2)-time(1)
!    write(*,'(a50,f32.15)') 'Time for long range ', time(3)-time(2)
!    write(*,'(a50,f32.15)') 'Time for long range forces',time(4)-time(3)
!    write(*,'(a50,f32.15)') 'Time for short range',time(5)-time(4)
!    write(*,'(a50,f32.15)') 'Time for plane ',time(6)-time(5)
!    write(*,'(a50,f32.15)') 'Time for Total without plane',time(5)-time(1)
!    write(*,'(a50,f32.15)') 'Time for Total',time(6)-time(1)
    call f_release_routine()
end subroutine calculate_forces_energy
!*****************************************************************************************
subroutine calparam(parini,atoms,poisson_rough,poisson)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(in):: poisson_rough
    type(typ_poisson), intent(inout):: poisson
    !local variables
    real(8):: tt1, tt2
    integer:: ngptot
    associate(ngpx=>poisson%ngpx)
    associate(ngpy=>poisson%ngpy)
    associate(ngpz=>poisson%ngpz)
    if(trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)=='bigdft') then
            call set_ngp_bps(parini,atoms,poisson_rough,poisson)
            !write(*,*) poisson%ngpx,poisson%ngpy, &
            !    poisson%ngpz
            !stop 'AFTER CALL TO set_ngp_bps'
        elseif(trim(parini%psolver_ann)=='kwald') then
            return
        endif
    elseif(trim(atoms%boundcond)=='slab') then
        ngpx=int(poisson%cell(1)/poisson_rough%hx)+1
        ngpy=int(poisson%cell(2)/poisson_rough%hx)+1
        ngpz=int(poisson%cell(3)/poisson_rough%hz)+1
        if(mod(ngpx,2)/=0) ngpx=ngpx+1
        if(mod(ngpy,2)/=0) ngpy=ngpy+1
        poisson%hx=poisson%cell(1)/real(ngpx,8)
        poisson%hy=poisson%cell(2)/real(ngpy,8)
        poisson%hz=poisson%cell(3)/real(ngpz-1,8)
    elseif(trim(atoms%boundcond)=='wire') then
        stop 'ERROR: wire BCs is not complete yet.'
    elseif(trim(atoms%boundcond)=='free') then
        stop 'ERROR: free BCs is not complete yet.'
    else
        write(*,'(2a)') 'ERROR: unknown BC in calparam ',trim(atoms%boundcond)
        stop
    endif
    poisson%nbgpx=int(poisson_rough%rgcut/poisson%hx)+2
    poisson%nbgpy=int(poisson_rough%rgcut/poisson%hy)+2
    poisson%nbgpz=int(poisson_rough%rgcut/poisson%hz)+2
    if(trim(atoms%boundcond)=='bulk') then
        poisson%nagpx=poisson%nbgpx+1
        poisson%nagpy=poisson%nbgpy+1
        poisson%nagpz=poisson%nbgpz+1
    elseif(trim(atoms%boundcond)=='slab') then
        poisson%nagpx=poisson%nbgpx+1
        poisson%nagpy=poisson%nbgpy+1
        poisson%nagpz=0
        ngpz=ngpz+2*poisson%nbgpz
    elseif(trim(atoms%boundcond)=='wire') then
        stop 'ERROR: wire BCs is not complete yet.'
    elseif(trim(atoms%boundcond)=='free') then
        stop 'ERROR: free BCs is not complete yet.'
    else
        write(*,'(2a)') 'ERROR: unknown BC in calparam ',trim(atoms%boundcond)
        stop
    endif
    tt1=real((ngpx+2*poisson%nbgpx)*(ngpy+2*poisson%nbgpy),8)
    tt2=real((ngpx+2)*(ngpy),8)
    poisson%ngpztot=ngpz*(int(tt1/tt2)+2)
    ngptot=ngpx*ngpy*ngpz
    write(*,'(a50,4i)') 'ngpx,ngpy,ngpz,ngptot',ngpx,ngpy,ngpz,ngptot
    write(*,'(a50,3i)') 'nbgpx,nbgpy,nbgpz',poisson%nbgpx,poisson%nbgpy,poisson%nbgpz
    write(*,'(a50,3i)') 'nagpx,nagpy,nagpz',poisson%nagpx,poisson%nagpy,poisson%nagpz
    write(*,'(a50,3f14.7)') 'hgx,hgy,hgz',poisson%hx,poisson%hy,poisson%hz
    write(*,'(a50,i)') 'ngpztot',poisson%ngpztot
    end associate
    end associate
    end associate
end subroutine calparam
!*****************************************************************************************
!This subroutine determines the limits of grids in a sphere.
subroutine determine_glimitsphere(poisson)
    use mod_interface
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: ix, iy, iz
    real(8):: rgcut, rgcutsq
    rgcut=max(poisson%hx*poisson%nbgpx,poisson%hy*poisson%nbgpy,poisson%hz*poisson%nbgpz)
    rgcutsq=rgcut**2
    do iz=-poisson%nbgpz,poisson%nbgpz
        do iy=-poisson%nbgpy,poisson%nbgpy
            poisson%mboundg(1,iy,iz)=1
            poisson%mboundg(2,iy,iz)=0
        enddo
    enddo
    do iz=0,poisson%nbgpz
    do iy=-poisson%nbgpy,poisson%nbgpy
    do ix=0,poisson%nbgpx
        if(ix**2*poisson%hx**2+iy**2*poisson%hy**2+iz**2*poisson%hz**2<=rgcutsq) then
            poisson%mboundg(1,iy,iz)=-ix
            poisson%mboundg(2,iy,iz)=ix
        endif
    enddo
    enddo
    enddo
    do iz=-poisson%nbgpz,-1
        do iy=-poisson%nbgpy,poisson%nbgpy
            poisson%mboundg(1,iy,iz)=poisson%mboundg(1,iy,-iz)
            poisson%mboundg(2,iy,iz)=poisson%mboundg(2,iy,-iz)
        enddo
    enddo
end subroutine determine_glimitsphere
!*****************************************************************************************
subroutine putgaussgrid(parini,bc,reset,nat,rxyz,qat,gausswidth,poisson)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gausswidth(nat)
    type(typ_poisson), intent(inout):: poisson
    !work array that is bigger than rho array, big enough to include of 
    !grid points that are outside of box.
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8), allocatable:: wx(:), wy(:), wz(:)
    real(8):: rhoz, rhoyz, pi
    real(8):: hgxinv, hgyinv, hgzinv
    real(8):: width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz
    real(8):: xat, yat, zat, facqiat, fac, width
    integer:: iat, iw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz
    integer:: ngpx, ngpy, ngpz, iii
    real(8), allocatable:: wa(:,:,:)
    call f_routine(id='putgaussgrid')
    associate(nagpx=>poisson%nagpx,nagpy=>poisson%nagpy,nagpz=>poisson%nagpz)
    ngpx=poisson%ngpx
    ngpy=poisson%ngpy
    ngpz=poisson%ngpz
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>poisson%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-poisson%nbgpx.to.poisson%nbgpx],id='wx')
    wy=f_malloc([-poisson%nbgpy.to.poisson%nbgpy],id='wy')
    wz=f_malloc([-poisson%nbgpz.to.poisson%nbgpz],id='wz')
    pi=4.d0*atan(1.d0)
    hgxinv=1.d0/poisson%hx
    hgyinv=1.d0/poisson%hy
    hgzinv=1.d0/poisson%hz
!    write(*,*) hgxinv,hgyinv,hgzinv
    !assaigning the gaussian charge density on grid points, saving in the 
    !work array which is bigger than rho array.
    !write(*,*) '********************************** ',associated(wa)
    if(trim(bc)=='bulk') then
        iii=0
    elseif(trim(bc)=='slab') then
        iii=1
    endif
    wa=0.d0
    if (parini%ewald) then
        width= poisson%alpha
        width_inv=1.d0/width
        fac=1.d0/(width*sqrt(pi))**3
        width_inv_hgx=width_inv*poisson%hx
        width_inv_hgy=width_inv*poisson%hy
        width_inv_hgz=width_inv*poisson%hz
    endif
    do iat=1,nat
        !shift the gaussian centers
        iatox=nint(rxyz(1,iat)*hgxinv)+1
        iatoy=nint(rxyz(2,iat)*hgyinv)+1
        iatoz=nint(rxyz(3,iat)*hgzinv)+1+poisson%nbgpz*iii
        xat=rxyz(1,iat)-(iatox-1)*poisson%hx
        yat=rxyz(2,iat)-(iatoy-1)*poisson%hy
        zat=rxyz(3,iat)-(iatoz-1-poisson%nbgpz*iii)*poisson%hz
        !construct the one-dimensional gaussians

        if (.not.parini%ewald) then
            width=gausswidth(iat)
            width_inv=1.d0/width
            fac=1.d0/(width*sqrt(pi))**3
            width_inv_hgx=width_inv*poisson%hx
            width_inv_hgy=width_inv*poisson%hy
            width_inv_hgz=width_inv*poisson%hz
        endif

        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        do iw=-poisson%nbgpx,poisson%nbgpx
            wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
        enddo
        do iw=-poisson%nbgpy,poisson%nbgpy
            wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
        enddo
        do iw=-poisson%nbgpz,poisson%nbgpz
            wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
        enddo
        facqiat=fac*qat(iat)
        do iz=-poisson%nbgpz,poisson%nbgpz
            rhoz=facqiat*wz(iz)
            jz=iatoz+iz
            do iy=-poisson%nbgpy,poisson%nbgpy
                rhoyz=rhoz*wy(iy)
                jy=iatoy+iy
                do ix=poisson%mboundg(1,iy,iz),poisson%mboundg(2,iy,iz)
                    jx=iatox+ix
                    !write(*,'(5i5)') iat,iatox,ix,iatoy,iy
                    wa(jx,jy,jz)=wa(jx,jy,jz)+rhoyz*wx(ix)
                enddo
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    do iz=1-nagpz,ngpz+nagpz
        do iy=1-nagpy,0
            do ix=1-nagpx,0
                wa(ix+ngpx,iy+ngpy,iz)=wa(ix+ngpx,iy+ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=1,ngpx
                wa(ix,iy+ngpy,iz)=wa(ix,iy+ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix-ngpx,iy+ngpy,iz)=wa(ix-ngpx,iy+ngpy,iz)+wa(ix,iy,iz)
            enddo
        enddo
        do iy=1,ngpy
            do ix=1-nagpx,0
                wa(ix+ngpx,iy,iz)=wa(ix+ngpx,iy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix-ngpx,iy,iz)=wa(ix-ngpx,iy,iz)+wa(ix,iy,iz)
            enddo
        enddo
        do iy=ngpy+1,ngpy+nagpy
            do ix=1-nagpx,0
                wa(ix+ngpx,iy-ngpy,iz)=wa(ix+ngpx,iy-ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=1,ngpx
                wa(ix,iy-ngpy,iz)=wa(ix,iy-ngpy,iz)+wa(ix,iy,iz)
            enddo
            do ix=ngpx+1,ngpx+nagpx
                wa(ix-ngpx,iy-ngpy,iz)=wa(ix-ngpx,iy-ngpy,iz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    do iz=1-nagpz,0
        do iy=1,ngpy
            do ix=1,ngpx
                wa(ix,iy,iz+ngpz)=wa(ix,iy,iz+ngpz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    do iz=ngpz+1,ngpz+nagpz
        do iy=1,ngpy
            do ix=1,ngpx
                wa(ix,iy,iz-ngpz)=wa(ix,iy,iz-ngpz)+wa(ix,iy,iz)
            enddo
        enddo
    enddo
    if(reset) then
        do iz=1,ngpz
            do iy=1,ngpy
                do ix=1,ngpx
                    poisson%rho(ix,iy,iz)=wa(ix,iy,iz)
                enddo
            enddo
        enddo
    else
        do iz=1,ngpz
            do iy=1,ngpy
                do ix=1,ngpx
                    poisson%rho(ix,iy,iz)=poisson%rho(ix,iy,iz)+wa(ix,iy,iz)
                enddo
            enddo
        enddo
    endif
    !do iz=1,ngpz
    !    do iy=1,ngpy
    !        do ix=1,ngpx
    !            write(61,'(3i4,es20.10)') ix,iy,iz,poisson%rho(ix,iy,iz)
    !        enddo
    !    enddo
    !enddo
    !stop
    call f_free(wx)
    call f_free(wy)
    call f_free(wz)
    call f_free(wa)
    end associate
    call f_release_routine()
end subroutine putgaussgrid
!*****************************************************************************************
subroutine longerange_forces(parini,atoms,poisson,gausswidth)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: gausswidth(atoms%nat)
    !work array that is bigger than rho array, big enough to include of 
    !grid points that are outside of box.
    !local variables
    real(8), allocatable:: wx(:), wy(:), wz(:) !values of one dimensional Gaussian functions
    real(8), allocatable:: vx(:), vy(:), vz(:) !derivatives of wx,wy,wz arrays
    real(8):: rhoz, rhoyz, pi
    real(8):: hgxinv, hgyinv, hgzinv, hgxhgyhgz
    real(8):: width, width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz
    real(8):: xat, yat, zat, facqiat, fac
    integer:: iat, ivw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz, iyt, izt, iii
    real(8):: fx, fy, fz, dx, dy, dz, derrhoz, derrhoyz, derrhozy
    real(8), allocatable:: wa(:,:,:)
    call f_routine(id='longerange_forces')
    associate(nagpx=>poisson%nagpx,nagpy=>poisson%nagpy,nagpz=>poisson%nagpz)
    associate(ngpx=>poisson%ngpx)
    associate(ngpy=>poisson%ngpy)
    associate(ngpz=>poisson%ngpz)
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>poisson%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-poisson%nbgpx.to.poisson%nbgpx],id='wx')
    wy=f_malloc([-poisson%nbgpy.to.poisson%nbgpy],id='wy')
    wz=f_malloc([-poisson%nbgpz.to.poisson%nbgpz],id='wz')
    vx=f_malloc([-poisson%nbgpx.to.poisson%nbgpx],id='vx')
    vy=f_malloc([-poisson%nbgpy.to.poisson%nbgpy],id='vy')
    vz=f_malloc([-poisson%nbgpz.to.poisson%nbgpz],id='vz')
    pi=4.d0*atan(1.d0)
    hgxhgyhgz=poisson%hx*poisson%hy*poisson%hz
    hgxinv=1.d0/poisson%hx
    hgyinv=1.d0/poisson%hy
    hgzinv=1.d0/poisson%hz
    !---------------------------------------------------------------------------
    do iz=1-nagpz,ngpz+nagpz
        izt=iz+(sign(ngpz,-iz)+sign(ngpz,ngpz-iz))/2
        do iy=1-poisson%nagpy,ngpy+poisson%nagpy
            iyt=iy+(sign(ngpy,-iy)+sign(ngpy,ngpy-iy))/2
            do ix=1-poisson%nagpx,0
                wa(ix,iy,iz)=poisson%rho(ix+ngpx,iyt,izt)
            enddo
            do ix=1,ngpx
                wa(ix,iy,iz)=poisson%rho(ix,iyt,izt)
            enddo
            do ix=ngpx+1,ngpx+poisson%nagpx
                wa(ix,iy,iz)=poisson%rho(ix-ngpx,iyt,izt)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    if(trim(atoms%boundcond)=='bulk') then
        iii=0
    elseif(trim(atoms%boundcond)=='slab') then
        iii=1
    endif
    !initialize the density 
    do iat=1,atoms%nat  
        !shift the Gaussian centers
        iatox=nint(atoms%rat(1,iat)*hgxinv)+1
        iatoy=nint(atoms%rat(2,iat)*hgyinv)+1
        iatoz=nint(atoms%rat(3,iat)*hgzinv)+1+poisson%nbgpz*iii
        xat=atoms%rat(1,iat)-(iatox-1)*poisson%hx
        yat=atoms%rat(2,iat)-(iatoy-1)*poisson%hy
        zat=atoms%rat(3,iat)-(iatoz-1-poisson%nbgpz*iii)*poisson%hz
        width=gausswidth(iat)
        width_inv=1.d0/width
        fac=2.d0/(width*(width*sqrt(pi))**3)
        width_inv_hgx=width_inv*poisson%hx
        width_inv_hgy=width_inv*poisson%hy
        width_inv_hgz=width_inv*poisson%hz
        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        !construct the one-dimensional gaussians
        do ivw=-poisson%nbgpx,poisson%nbgpx
            dx=width_inv_hgx*ivw-width_inv_xat
            wx(ivw)=exp(-dx*dx)
            vx(ivw)=-wx(ivw)*dx
        enddo
        do ivw=-poisson%nbgpy,poisson%nbgpy
            dy=width_inv_hgy*ivw-width_inv_yat
            wy(ivw)=exp(-dy*dy)
            vy(ivw)=-wy(ivw)*dy
        enddo
        do ivw=-poisson%nbgpz,poisson%nbgpz
            dz=width_inv_hgz*ivw-width_inv_zat
            wz(ivw)=exp(-dz*dz)
            vz(ivw)=-wz(ivw)*dz
        enddo
        fx=0.d0;fy=0.d0;fz=0.d0
        facqiat=fac*atoms%qat(iat)
        do iz=-poisson%nbgpz,poisson%nbgpz
            rhoz=facqiat*wz(iz)
            derrhoz=facqiat*vz(iz)
            jz=iatoz+iz
            do iy=-poisson%nbgpy,poisson%nbgpy
                rhoyz=rhoz*wy(iy)
                derrhozy=rhoz*vy(iy)
                derrhoyz=derrhoz*wy(iy)
                jy=iatoy+iy
                do ix=poisson%mboundg(1,iy,iz),poisson%mboundg(2,iy,iz)
                    jx=iatox+ix
                    fx=fx+rhoyz*vx(ix)*wa(jx,jy,jz)
                    fy=fy+derrhozy*wx(ix)*wa(jx,jy,jz)
                    fz=fz+derrhoyz*wx(ix)*wa(jx,jy,jz)
                enddo
            enddo
        enddo
        atoms%fat(1,iat)=atoms%fat(1,iat)+fx*hgxhgyhgz
        atoms%fat(2,iat)=atoms%fat(2,iat)+fy*hgxhgyhgz
        atoms%fat(3,iat)=atoms%fat(3,iat)+fz*hgxhgyhgz
    enddo
    if((.not. poisson%point_particle) .and. trim(parini%bias_type)=='fixed_efield') then
        do iat=1,atoms%nat
            atoms%fat(3,iat)=atoms%fat(3,iat)-parini%efield*0.5d0*atoms%qat(iat)
        enddo
    endif
    call f_free(wa)
    call f_free(wx)
    call f_free(wy)
    call f_free(wz)
    call f_free(vx)
    call f_free(vy)
    call f_free(vz)
    end associate
    end associate
    end associate
    end associate
    call f_release_routine()
end subroutine longerange_forces
!*****************************************************************************************
