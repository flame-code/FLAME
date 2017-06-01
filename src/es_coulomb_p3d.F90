!*****************************************************************************************
!CAUTION: Why is this subroutine called in ANN eem1 when system is bulk
!and kwald is the psolver.
subroutine construct_ewald_p3d(parini,atoms,ewald_p3d)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d):: ewald_p3d_rough
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    include 'fftw3.f'
    real(8):: pi, shortrange_at_rcut
    call f_routine(id='construct_ewald_p3d')
    pi=4.d0*atan(1.d0)
    ewald_p3d_rough%hgx=parini%hx_ewald
    ewald_p3d_rough%hgy=parini%hy_ewald
    ewald_p3d_rough%hgz=parini%hz_ewald
    if (parini%ewald .and. parini%alpha_ewald>0.d0) then
        ewald_p3d%alpha=parini%alpha_ewald
    else if (ewald_p3d%alpha< 0.d0 .and. parini%alpha_ewald<= 0.d0) then
            write(*,*) "ERROR : alpha is undefined"
            stop
    endif
    ewald_p3d%linked_lists%rcut=parini%rcut_ewald
    ewald_p3d_rough%rgcut=parini%rgcut_ewald*ewald_p3d%alpha
    ewald_p3d%spline%nsp=parini%nsp_ewald
    ewald_p3d%cell(1)=atoms%cellvec(1,1)
    ewald_p3d%cell(2)=atoms%cellvec(2,2)
    ewald_p3d%cell(3)=atoms%cellvec(3,3)
    ewald_p3d%vu=parini%vu_ewald
    ewald_p3d%vl=parini%vl_ewald
    call calparam(parini,atoms,ewald_p3d_rough,ewald_p3d)
    associate(ngpx=>ewald_p3d%poisson_p3d%ngpx)
    associate(ngpy=>ewald_p3d%poisson_p3d%ngpy)
    associate(ngpz=>ewald_p3d%poisson_p3d%ngpz)
    associate(nbgpy=>ewald_p3d%nbgpy)
    associate(nbgpz=>ewald_p3d%nbgpz)
    if(trim(atoms%boundcond)=='bulk') then
        ewald_p3d%poisson_p3d%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
            id='ewald_p3d%poisson_p3d%rho')
        ewald_p3d%poisson_p3d%pot=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
            id='ewald_p3d%poisson_p3d%pot')
    elseif(trim(atoms%boundcond)=='slab') then
        ewald_p3d%poisson_p3d%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
            id='ewald_p3d%poisson_p3d%rho')
        ewald_p3d%poisson_p3d%pot=f_malloc([1.to.ngpx+2,1.to.ngpy,1.to.ewald_p3d%ngpztot], &
            id='ewald_p3d%poisson_p3d%pot')
    else
        write(*,*) 'ERROR: other BCs are not yet considered.'
    endif
    ewald_p3d%mboundg=f_malloc([1.to.2,-nbgpy.to.nbgpy,-nbgpz.to.nbgpz],id='ewald_p3d%mboundg')
    if(trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)=='bigdft') then
            call construct_ewald_bps(parini,atoms,ewald_p3d)
        endif
    elseif(trim(atoms%boundcond)=='slab') then
        call ps2dp1df_construction(ewald_p3d%poisson_p3d)
    endif
    call determine_glimitsphere(ewald_p3d)
    ewald_p3d%epotfixed=dot_product(atoms%qat,atoms%qat)/(sqrt(2.d0*pi)*ewald_p3d%alpha)
    shortrange_at_rcut=erfc(ewald_p3d%linked_lists%rcut/(sqrt(2.d0)*ewald_p3d%alpha))/(ewald_p3d%linked_lists%rcut)
    if(parini%iverbose>=2) then
        write(*,*) 'real part in rcut',shortrange_at_rcut
    endif
    end associate
    end associate
    end associate
    end associate
    end associate
    call f_release_routine()
end subroutine construct_ewald_p3d
!*****************************************************************************************
subroutine destruct_ewald_p3d(parini,atoms,ewald_p3d)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_potential, only: bias  
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    call f_routine(id='destruct_ewald_p3d')
    if(trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)=='bigdft') then
            call destruct_ewald_bps(ewald_p3d)
        endif
    elseif(trim(atoms%boundcond)=='slab') then
        call ps2dp1df_destruction(ewald_p3d%poisson_p3d)
    endif
    call f_free(ewald_p3d%poisson_p3d%rho)
    call f_free(ewald_p3d%poisson_p3d%pot)
    if(trim(bias)=='yes') then
        call f_free(ewald_p3d%poisson_p3d%pots)
    endif
    call f_free(ewald_p3d%mboundg)
    call f_release_routine()
end subroutine destruct_ewald_p3d
!*****************************************************************************************
subroutine calculate_forces_energy(parini,ewald_p3d,atoms)
    use mod_interface
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    use mod_potential, only: bias 
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    !local variables
    real(8):: epotlong, epotplane !, epotshort
    real(8):: time(10)
    !The following dummy varables definition to be deleted later.
    !real(8):: totrho
    real(8):: beta, pi, charge,c,charge0,E,tmp
    integer:: igpx, igpy, igpz, iat, iz
    integer:: npl, npu, nlayer, ngpx, ngpy, ngpz
    real(8),allocatable :: pots_layer(:,:,:,:)
    real(8):: vl, vu, A, d, rl, ru, dipole_correction, dipole
    real(8),allocatable :: gausswidth(:)  
    call f_routine(id='calculate_forces_energy')
    ngpz=ewald_p3d%poisson_p3d%ngpz
    ngpy=ewald_p3d%poisson_p3d%ngpy
    ngpx=ewald_p3d%poisson_p3d%ngpx

    pi=4.d0*atan(1.d0)
    beta=0.d0
    do iat=1,atoms%nat
        beta=beta+atoms%qat(iat)*atoms%rat(3,iat)
    enddo
    beta=beta*2.d0*pi*ewald_p3d%poisson_p3d%ngpx*ewald_p3d%poisson_p3d%ngpy/(ewald_p3d%cell(1)*ewald_p3d%cell(2))
    gausswidth=f_malloc([1.to.atoms%nat],id='gausswidth')
    gausswidth(:)=ewald_p3d%alpha

    !write(*,*) 'total momentum z component',beta
    !write(*,*) 'total momentum z component',0.13074051987178871d5/beta
    call cpu_time(time(1))
    call putgaussgrid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%qat,gausswidth,ewald_p3d)
    call cpu_time(time(2))
    !-----------------------------------------------------------------------
    !totrho=0.d0
    !do iz=1,ngpz;do iy=1,ngpy;do ix=1,ewald_p3d%poisson_p3d%ngpx
    !    totrho=totrho+rho(ix,iy,iz)
    !enddo;enddo;enddo
    !write(*,*) 'totrho',totrho
    !-----------------------------------------------------------------------
    !pot=0.d0
    call calculate_potener_pot(ewald_p3d%poisson_p3d,ewald_p3d%cell,ewald_p3d%hgx,ewald_p3d%hgy,ewald_p3d%hgz,epotlong,beta)
    call cpu_time(time(3))
    do igpz=1,ewald_p3d%poisson_p3d%ngpz
    do igpy=1,ewald_p3d%poisson_p3d%ngpy
    do igpx=1,ewald_p3d%poisson_p3d%ngpx
        ewald_p3d%poisson_p3d%rho(igpx,igpy,igpz)=ewald_p3d%poisson_p3d%pot(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call longerange_forces(atoms,ewald_p3d,gausswidth)
    call cpu_time(time(4))
    !call shortenergy(atoms,ewald_p3d%linked_lists,ewald_p3d%spline,ewald_p3d%alpha,ewald_p3d%cell,epotshort)
    call cpu_time(time(5))
    epotplane=0.d0
    if (trim(bias)=='yes') then
        vl=parini%vl_ewald
        vu=parini%vu_ewald+parini%vu_ac_ewald*sin(parini%frequency_ewald*parini%time_dynamics)
        d= ewald_p3d%cell(3)
        write(*,*)'--------------------------------------------------------'
        write(*,*)"distance between planes",d
        A= ewald_p3d%poisson_p3d%ngpx*ewald_p3d%poisson_p3d%ngpy*ewald_p3d%hgx*ewald_p3d%hgy
        write(*,*)"A= ", A
        c= A/(4.d0*pi*d)
        write(*,*) "C=A/(4pid) = ", A/(4.d0*pi*d)
        dipole = beta*(ewald_p3d%hgx*ewald_p3d%hgy)
        charge0= -dipole/(2*pi*d)
        charge = -dipole/(2*pi*d)+c*(vu-vl)
        write(88,*)vu , beta, charge
        write(89,*)vu , beta, charge0
        write(*,*)'dipole = ', dipole/(2*pi)
        dipole_correction = 3/(4*pi)*dipole**2/(ewald_p3d%cell(3)*ewald_p3d%cell(2)*ewald_p3d%cell(1))
        write(*,*)'dipole correction ', dipole_correction
        write(*,*)'pot correction ', charge0*(vu-vl)
        !dipole_correction =0.d0
        dipole_correction =dipole_correction +0.5*charge0*(vu-vl)+0.5*c*(vu-vl)**2
        write(*,*)'charge on upper  plate  ', charge
        write(*,*)'--------------------------------------------------------'
        ewald_p3d%poisson_p3d%npu=ewald_p3d%poisson_p3d%ngpz-ewald_p3d%nbgpz
        ewald_p3d%poisson_p3d%npl=1+ewald_p3d%nbgpz  
        write(*,*) "min rat_z " ,minval(atoms%rat(3,:)),"max rat_z ",maxval(atoms%rat(3,:))
        npl=ewald_p3d%poisson_p3d%npl
        npu=ewald_p3d%poisson_p3d%npu

        ewald_p3d%poisson_p3d%pots=f_malloc([1.to.ewald_p3d%poisson_p3d%ngpx+2,1.to.ewald_p3d%poisson_p3d%ngpy,npl.to.npu],id='ewald_p3d%poisson_p3d%pots')
        write(*,*)"npu,npl",ewald_p3d%poisson_p3d%npu,ewald_p3d%poisson_p3d%npl
        nlayer=1
        if (parini%cal_charge) then 
            nlayer=5
            pots_layer=f_malloc([1.to.ewald_p3d%poisson_p3d%ngpx,1.to.ewald_p3d%poisson_p3d%ngpy,1.to.2,1.to.nlayer-1],id='pots_layer')
        endif
        ewald_p3d%poisson_p3d%pots=0.d0
        call erfc_surface_zero(parini,atoms,ewald_p3d,nlayer)
        if (parini%cal_charge) then
            pots_layer(1:ngpx,1:ngpy,1,:)=ewald_p3d%poisson_p3d%pots(1:ngpx,1:ngpy,npl+1:npl+nlayer-1)
            pots_layer(1:ngpx,1:ngpy,2,:)=ewald_p3d%poisson_p3d%pots(1:ngpx,1:ngpy,npu-1:npu-(nlayer-1):-1)
            ewald_p3d%poisson_p3d%pots(:,:,npl+1:npl+nlayer-1)  =0.d0
            ewald_p3d%poisson_p3d%pots(:,:,npu-(nlayer-1):npu-1)=0.d0
        endif


        !    force repulsion on planes
!        do iat=1,atoms%nat
!            rl=abs(atoms%rat(3,iat))
!            ru=abs( ewald_p3d%cell(3)-atoms%rat(3,iat))
!            atoms%fat(3,iat)=160.d0*exp(-1.6d0*rl)-160.d0*exp(-1.6d0*ru)
!        enddo
!
        call sollaplaceq(ewald_p3d%poisson_p3d,ewald_p3d%hgz,ewald_p3d%cell,vl,vu)
        call calculate_force_ener_plane(atoms,ewald_p3d,epotplane)

        if (parini%cal_charge) then 
            call surface_charge(ewald_p3d,pots_layer,vl,vu)
            call f_free(pots_layer)
        endif
        epotplane = epotplane+dipole_correction
        call f_free(ewald_p3d%poisson_p3d%pots)
    end if
    if (trim(parini%bias_field)=='yes') then
        vl=parini%vl_ewald
        vu=parini%vu_ewald
        d = ewald_p3d%cell(3)
        E =- (vu-vl)/d
        write(*,*)'--------------------------------------------------------'
        write(*,*)"distance between planes",d
        A= ewald_p3d%poisson_p3d%ngpx*ewald_p3d%poisson_p3d%ngpy*ewald_p3d%hgx*ewald_p3d%hgy
        c= A/(4.d0*pi*d)
        write(*,*) "C=A/(4pid) = ", A/(4.d0*pi*d)
        dipole = beta*(ewald_p3d%hgx*ewald_p3d%hgy)
        charge0= -dipole/(2*pi*d)
        charge = -dipole/(2*pi*d)+c*(vu-vl)
        write(88,*)vu , beta, charge
        write(89,*)vu , beta, charge0
        write(*,*)'dipole = ', dipole/(2*pi)
        write(*,*)'charge on upper  plate  ', charge
        write(*,*)'--------------------------------------------------------'
        ewald_p3d%poisson_p3d%npu=ewald_p3d%poisson_p3d%ngpz-ewald_p3d%nbgpz
        ewald_p3d%poisson_p3d%npl=1+ewald_p3d%nbgpz  
        write(*,*) "min rat_z " ,minval(atoms%rat(3,:)),"max rat_z ",maxval(atoms%rat(3,:))
        npl=ewald_p3d%poisson_p3d%npl
        npu=ewald_p3d%poisson_p3d%npu
        epotplane = 0.d0
        do iat=1,atoms%nat
            epotplane = epotplane - E * atoms%qat(iat)*atoms%rat(3,iat)
            atoms%fat(3,iat)=atoms%fat(3,iat)+ E * atoms%qat(iat)
        enddo
    endif
    call cpu_time(time(6))
    !atoms%epot=epotlong+epotshort-ewald_p3d%epotfixed+epotplane
    atoms%epot=epotlong+-ewald_p3d%epotfixed+epotplane
    write(*,*) '-----------------------------------------------------------'
    write(*,'(a50,e32.15)') 'epotfixed',ewald_p3d%epotfixed
    write(*,'(a50,e32.15)') 'epotlong',epotlong
    write(*,'(a50,e32.15)') 'epotplane',epotplane
    write(*,'(a50,e32.15)') 'epottotal',atoms%epot
    !write(*,*) '-----------------------------------------------------------'
    write(*,'(a50,f32.15)') 'Time for putgaussgrid ',time(2)-time(1)
    write(*,'(a50,f32.15)') 'Time for long range ', time(3)-time(2)
    write(*,'(a50,f32.15)') 'Time for long range forces',time(4)-time(3)
    write(*,'(a50,f32.15)') 'Time for short range',time(5)-time(4)
    write(*,'(a50,f32.15)') 'Time for plane ',time(6)-time(5)
    write(*,'(a50,f32.15)') 'Time for Total without plane',time(5)-time(1)
    write(*,'(a50,f32.15)') 'Time for Total',time(6)-time(1)
    call f_release_routine()
end subroutine calculate_forces_energy
!*****************************************************************************************
subroutine calparam(parini,atoms,ewald_p3d_rough,ewald_p3d)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ewald_p3d), intent(in):: ewald_p3d_rough
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    real(8):: tt1, tt2
    integer:: ngptot
    associate(ngpx=>ewald_p3d%poisson_p3d%ngpx)
    associate(ngpy=>ewald_p3d%poisson_p3d%ngpy)
    associate(ngpz=>ewald_p3d%poisson_p3d%ngpz)
    if(trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)=='bigdft') then
            call set_ngp_bps(atoms,ewald_p3d_rough,ewald_p3d)
            !write(*,*) ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy, &
            !    ewald_p3d%poisson_p3d%ngpz
            !stop 'AFTER CALL TO set_ngp_bps'
        elseif(trim(parini%psolver_ann)=='kwald') then
            return
        endif
    elseif(trim(atoms%boundcond)=='slab') then
        ngpx=int(ewald_p3d%cell(1)/ewald_p3d_rough%hgx)+1
        ngpy=int(ewald_p3d%cell(2)/ewald_p3d_rough%hgx)+1
        ngpz=int(ewald_p3d%cell(3)/ewald_p3d_rough%hgz)+1
        if(mod(ngpx,2)/=0) ngpx=ngpx+1
        if(mod(ngpy,2)/=0) ngpy=ngpy+1
        ewald_p3d%hgx=ewald_p3d%cell(1)/real(ngpx,8)
        ewald_p3d%hgy=ewald_p3d%cell(2)/real(ngpy,8)
        ewald_p3d%hgz=ewald_p3d%cell(3)/real(ngpz-1,8)
    elseif(trim(atoms%boundcond)=='wire') then
        stop 'ERROR: wire BCs is not complete yet.'
    elseif(trim(atoms%boundcond)=='free') then
        stop 'ERROR: free BCs is not complete yet.'
    else
        write(*,'(2a)') 'ERROR: unknown BC in calparam ',trim(atoms%boundcond)
        stop
    endif
    ewald_p3d%nbgpx=int(ewald_p3d_rough%rgcut/ewald_p3d%hgx)+2
    ewald_p3d%nbgpy=int(ewald_p3d_rough%rgcut/ewald_p3d%hgy)+2
    ewald_p3d%nbgpz=int(ewald_p3d_rough%rgcut/ewald_p3d%hgz)+2
    if(trim(atoms%boundcond)=='bulk') then
        ewald_p3d%nagpx=ewald_p3d%nbgpx+1
        ewald_p3d%nagpy=ewald_p3d%nbgpy+1
        ewald_p3d%nagpz=ewald_p3d%nbgpz+1
    elseif(trim(atoms%boundcond)=='slab') then
        ewald_p3d%nagpx=ewald_p3d%nbgpx+1
        ewald_p3d%nagpy=ewald_p3d%nbgpy+1
        ewald_p3d%nagpz=0
        ngpz=ngpz+2*ewald_p3d%nbgpz
    elseif(trim(atoms%boundcond)=='wire') then
        stop 'ERROR: wire BCs is not complete yet.'
    elseif(trim(atoms%boundcond)=='free') then
        stop 'ERROR: free BCs is not complete yet.'
    else
        write(*,'(2a)') 'ERROR: unknown BC in calparam ',trim(atoms%boundcond)
        stop
    endif
    tt1=real((ngpx+2*ewald_p3d%nbgpx)*(ngpy+2*ewald_p3d%nbgpy),8)
    tt2=real((ngpx+2)*(ngpy),8)
    ewald_p3d%ngpztot=ngpz*(int(tt1/tt2)+2)
    ngptot=ngpx*ngpy*ngpz
    write(*,'(a50,4i)') 'ngpx,ngpy,ngpz,ngptot',ngpx,ngpy,ngpz,ngptot
    write(*,'(a50,3i)') 'nbgpx,nbgpy,nbgpz',ewald_p3d%nbgpx,ewald_p3d%nbgpy,ewald_p3d%nbgpz
    write(*,'(a50,3i)') 'nagpx,nagpy,nagpz',ewald_p3d%nagpx,ewald_p3d%nagpy,ewald_p3d%nagpz
    write(*,'(a50,3f14.7)') 'hgx,hgy,hgz',ewald_p3d%hgx,ewald_p3d%hgy,ewald_p3d%hgz
    write(*,'(a50,i)') 'ngpztot',ewald_p3d%ngpztot
    end associate
    end associate
    end associate
end subroutine calparam
!*****************************************************************************************
!This subroutine determines the limits of grids in a sphere.
subroutine determine_glimitsphere(ewald_p3d)
    use mod_interface
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    integer:: ix, iy, iz
    real(8):: rgcut, rgcutsq
    rgcut=max(ewald_p3d%hgx*ewald_p3d%nbgpx,ewald_p3d%hgy*ewald_p3d%nbgpy,ewald_p3d%hgz*ewald_p3d%nbgpz)
    rgcutsq=rgcut**2
    do iz=-ewald_p3d%nbgpz,ewald_p3d%nbgpz
        do iy=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
            ewald_p3d%mboundg(1,iy,iz)=1
            ewald_p3d%mboundg(2,iy,iz)=0
        enddo
    enddo
    do iz=0,ewald_p3d%nbgpz
    do iy=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
    do ix=0,ewald_p3d%nbgpx
        if(ix**2*ewald_p3d%hgx**2+iy**2*ewald_p3d%hgy**2+iz**2*ewald_p3d%hgz**2<=rgcutsq) then
            ewald_p3d%mboundg(1,iy,iz)=-ix
            ewald_p3d%mboundg(2,iy,iz)=ix
        endif
    enddo
    enddo
    enddo
    do iz=-ewald_p3d%nbgpz,-1
        do iy=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
            ewald_p3d%mboundg(1,iy,iz)=ewald_p3d%mboundg(1,iy,-iz)
            ewald_p3d%mboundg(2,iy,iz)=ewald_p3d%mboundg(2,iy,-iz)
        enddo
    enddo
end subroutine determine_glimitsphere
!*****************************************************************************************
subroutine putgaussgrid(parini,bc,reset,nat,rxyz,qat,gausswidth,ewald_p3d)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
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
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
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
    associate(nagpx=>ewald_p3d%nagpx,nagpy=>ewald_p3d%nagpy,nagpz=>ewald_p3d%nagpz)
    ngpx=ewald_p3d%poisson_p3d%ngpx
    ngpy=ewald_p3d%poisson_p3d%ngpy
    ngpz=ewald_p3d%poisson_p3d%ngpz
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>ewald_p3d%poisson_p3d%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-ewald_p3d%nbgpx.to.ewald_p3d%nbgpx],id='wx')
    wy=f_malloc([-ewald_p3d%nbgpy.to.ewald_p3d%nbgpy],id='wy')
    wz=f_malloc([-ewald_p3d%nbgpz.to.ewald_p3d%nbgpz],id='wz')
    pi=4.d0*atan(1.d0)
    hgxinv=1.d0/ewald_p3d%hgx
    hgyinv=1.d0/ewald_p3d%hgy
    hgzinv=1.d0/ewald_p3d%hgz
    write(*,*) hgxinv,hgyinv,hgzinv
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
        width= ewald_p3d%alpha
        width_inv=1.d0/width
        fac=1.d0/(width*sqrt(pi))**3
        width_inv_hgx=width_inv*ewald_p3d%hgx
        width_inv_hgy=width_inv*ewald_p3d%hgy
        width_inv_hgz=width_inv*ewald_p3d%hgz
    endif
    do iat=1,nat
        !shift the gaussian centers
        iatox=nint(rxyz(1,iat)*hgxinv)+1
        iatoy=nint(rxyz(2,iat)*hgyinv)+1
        iatoz=nint(rxyz(3,iat)*hgzinv)+1+ewald_p3d%nbgpz*iii
        xat=rxyz(1,iat)-(iatox-1)*ewald_p3d%hgx
        yat=rxyz(2,iat)-(iatoy-1)*ewald_p3d%hgy
        zat=rxyz(3,iat)-(iatoz-1-ewald_p3d%nbgpz*iii)*ewald_p3d%hgz
        !construct the one-dimensional gaussians

        if (.not.parini%ewald) then
            width=gausswidth(iat)
            width_inv=1.d0/width
            fac=1.d0/(width*sqrt(pi))**3
            width_inv_hgx=width_inv*ewald_p3d%hgx
            width_inv_hgy=width_inv*ewald_p3d%hgy
            width_inv_hgz=width_inv*ewald_p3d%hgz
        endif

        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        do iw=-ewald_p3d%nbgpx,ewald_p3d%nbgpx
            wx(iw)=exp(-(width_inv_hgx*iw-width_inv_xat)**2)
        enddo
        do iw=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
            wy(iw)=exp(-(width_inv_hgy*iw-width_inv_yat)**2)
        enddo
        do iw=-ewald_p3d%nbgpz,ewald_p3d%nbgpz
            wz(iw)=exp(-(width_inv_hgz*iw-width_inv_zat)**2)
        enddo
        facqiat=fac*qat(iat)
        do iz=-ewald_p3d%nbgpz,ewald_p3d%nbgpz
            rhoz=facqiat*wz(iz)
            jz=iatoz+iz
            do iy=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
                rhoyz=rhoz*wy(iy)
                jy=iatoy+iy
                do ix=ewald_p3d%mboundg(1,iy,iz),ewald_p3d%mboundg(2,iy,iz)
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
                    ewald_p3d%poisson_p3d%rho(ix,iy,iz)=wa(ix,iy,iz)
                enddo
            enddo
        enddo
    else
        do iz=1,ngpz
            do iy=1,ngpy
                do ix=1,ngpx
                    ewald_p3d%poisson_p3d%rho(ix,iy,iz)=ewald_p3d%poisson_p3d%rho(ix,iy,iz)+wa(ix,iy,iz)
                enddo
            enddo
        enddo
    endif
    !do iz=1,ngpz
    !    do iy=1,ngpy
    !        do ix=1,ngpx
    !            write(61,'(3i4,es20.10)') ix,iy,iz,ewald_p3d%poisson_p3d%rho(ix,iy,iz)
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
subroutine longerange_forces(atoms,ewald_p3d,gausswidth)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
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
    associate(nagpx=>ewald_p3d%nagpx,nagpy=>ewald_p3d%nagpy,nagpz=>ewald_p3d%nagpz)
    associate(ngpx=>ewald_p3d%poisson_p3d%ngpx)
    associate(ngpy=>ewald_p3d%poisson_p3d%ngpy)
    associate(ngpz=>ewald_p3d%poisson_p3d%ngpz)
    !wa(1-nagpx:ngpx+nagpx,1-nagpy:ngpy+nagpy,1-nagpz:ngpz+nagpz)=>ewald_p3d%poisson_p3d%pot
    wa=f_malloc([1-nagpx.to.ngpx+nagpx,1-nagpy.to.ngpy+nagpy,1-nagpz.to.ngpz+nagpz],id='wa')
    wx=f_malloc([-ewald_p3d%nbgpx.to.ewald_p3d%nbgpx],id='wx')
    wy=f_malloc([-ewald_p3d%nbgpy.to.ewald_p3d%nbgpy],id='wy')
    wz=f_malloc([-ewald_p3d%nbgpz.to.ewald_p3d%nbgpz],id='wz')
    vx=f_malloc([-ewald_p3d%nbgpx.to.ewald_p3d%nbgpx],id='vx')
    vy=f_malloc([-ewald_p3d%nbgpy.to.ewald_p3d%nbgpy],id='vy')
    vz=f_malloc([-ewald_p3d%nbgpz.to.ewald_p3d%nbgpz],id='vz')
    pi=4.d0*atan(1.d0)
    hgxhgyhgz=ewald_p3d%hgx*ewald_p3d%hgy*ewald_p3d%hgz
    hgxinv=1.d0/ewald_p3d%hgx
    hgyinv=1.d0/ewald_p3d%hgy
    hgzinv=1.d0/ewald_p3d%hgz
    !---------------------------------------------------------------------------
    do iz=1-nagpz,ngpz+nagpz
        izt=iz+(sign(ngpz,-iz)+sign(ngpz,ngpz-iz))/2
        do iy=1-ewald_p3d%nagpy,ngpy+ewald_p3d%nagpy
            iyt=iy+(sign(ngpy,-iy)+sign(ngpy,ngpy-iy))/2
            do ix=1-ewald_p3d%nagpx,0
                wa(ix,iy,iz)=ewald_p3d%poisson_p3d%rho(ix+ngpx,iyt,izt)
            enddo
            do ix=1,ngpx
                wa(ix,iy,iz)=ewald_p3d%poisson_p3d%rho(ix,iyt,izt)
            enddo
            do ix=ngpx+1,ngpx+ewald_p3d%nagpx
                wa(ix,iy,iz)=ewald_p3d%poisson_p3d%rho(ix-ngpx,iyt,izt)
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
        iatoz=nint(atoms%rat(3,iat)*hgzinv)+1+ewald_p3d%nbgpz*iii
        xat=atoms%rat(1,iat)-(iatox-1)*ewald_p3d%hgx
        yat=atoms%rat(2,iat)-(iatoy-1)*ewald_p3d%hgy
        zat=atoms%rat(3,iat)-(iatoz-1-ewald_p3d%nbgpz*iii)*ewald_p3d%hgz
        width=gausswidth(iat)
        width_inv=1.d0/width
        fac=2.d0/(width*(width*sqrt(pi))**3)
        width_inv_hgx=width_inv*ewald_p3d%hgx
        width_inv_hgy=width_inv*ewald_p3d%hgy
        width_inv_hgz=width_inv*ewald_p3d%hgz
        width_inv_xat=width_inv*xat
        width_inv_yat=width_inv*yat
        width_inv_zat=width_inv*zat
        !construct the one-dimensional gaussians
        do ivw=-ewald_p3d%nbgpx,ewald_p3d%nbgpx
            dx=width_inv_hgx*ivw-width_inv_xat
            wx(ivw)=exp(-dx*dx)
            vx(ivw)=-wx(ivw)*dx
        enddo
        do ivw=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
            dy=width_inv_hgy*ivw-width_inv_yat
            wy(ivw)=exp(-dy*dy)
            vy(ivw)=-wy(ivw)*dy
        enddo
        do ivw=-ewald_p3d%nbgpz,ewald_p3d%nbgpz
            dz=width_inv_hgz*ivw-width_inv_zat
            wz(ivw)=exp(-dz*dz)
            vz(ivw)=-wz(ivw)*dz
        enddo
        fx=0.d0;fy=0.d0;fz=0.d0
        facqiat=fac*atoms%qat(iat)
        do iz=-ewald_p3d%nbgpz,ewald_p3d%nbgpz
            rhoz=facqiat*wz(iz)
            derrhoz=facqiat*vz(iz)
            jz=iatoz+iz
            do iy=-ewald_p3d%nbgpy,ewald_p3d%nbgpy
                rhoyz=rhoz*wy(iy)
                derrhozy=rhoz*vy(iy)
                derrhoyz=derrhoz*wy(iy)
                jy=iatoy+iy
                do ix=ewald_p3d%mboundg(1,iy,iz),ewald_p3d%mboundg(2,iy,iz)
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
