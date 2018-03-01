!*****************************************************************************************
subroutine init_hartree(parini,atoms,poisson)
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
    call f_routine(id='init_hartree')
    if(trim(parini%psolver_ann)=='p3d') then
        call init_hartree_p3d(parini,atoms,poisson)
    elseif(trim(parini%psolver_ann)=='bigdft') then
        call init_hartree_bps(parini,atoms,poisson)
    elseif(trim(parini%psolver_ann)=='kwald') then
        !It seems that nothing needs to be done for kwald method in the
        !init_hartree routine.
    endif
    if(allocated(poisson%qgrad)) deallocate(poisson%qgrad)
    allocate(poisson%qgrad(1:atoms%nat))
    call f_release_routine()
end subroutine init_hartree
!*****************************************************************************************
subroutine fini_hartree(parini,atoms,poisson)
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
    call f_routine(id='fini_hartree')
    if(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk') then
        if(trim(atoms%boundcond)=='bulk') then
            if(trim(parini%psolver_ann)=='bigdft') then
                call destruct_ewald_bps(poisson)
            endif
        elseif(trim(atoms%boundcond)=='slab') then
            call fini_psolver_p3d_slab(poisson)
        endif
        if(trim(parini%psolver_ann)/='kwald') then
            call f_free(poisson%rho)
            call f_free(poisson%pot)
        endif
        if(trim(parini%bias_type)=='p3dbias') then
         !   deallocate(poisson%pots)
        endif
    endif
    call f_release_routine()
end subroutine fini_hartree
!*****************************************************************************************
subroutine init_hartree_bps(parini,atoms,poisson)
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
    include 'fftw3.f'
    type(typ_poisson):: poisson_rough
    real(8):: pi, shortrange_at_rcut
    real(8):: tt1, tt2
    integer:: ngptot
    call f_routine(id='init_hartree_bps')
    associate(ngpx=>poisson%ngpx)
    associate(ngpy=>poisson%ngpy)
    associate(ngpz=>poisson%ngpz)
    if(trim(atoms%boundcond)/='bulk') then
        write(*,*) 'ERROR: init_hartree_bps currently assumes BC=bulk'
        stop
    endif
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
    poisson%rgcut=poisson_rough%rgcut
    poisson%spline%nsp=parini%nsp_ewald
    poisson%cell(1)=atoms%cellvec(1,1)
    poisson%cell(2)=atoms%cellvec(2,2)
    poisson%cell(3)=atoms%cellvec(3,3)
    poisson%vu=parini%vu_ewald
    poisson%vl=parini%vl_ewald
    !---------------------------------------------------------------------------
    !call calparam(parini,atoms,poisson_rough,poisson)
    call set_ngp_bps(parini,atoms,poisson_rough,poisson)
    ngptot=ngpx*ngpy*ngpz
    write(*,'(a50,4i)') 'ngpx,ngpy,ngpz,ngptot',ngpx,ngpy,ngpz,ngptot
    !write(*,'(a50,3i)') 'nbgpx,nbgpy,nbgpz',poisson%nbgpx,poisson%nbgpy,poisson%nbgpz
    !write(*,'(a50,3i)') 'nagpx,nagpy,nagpz',poisson%nagpx,poisson%nagpy,poisson%nagpz
    write(*,'(a50,3f14.7)') 'hgx,hgy,hgz',poisson%hx,poisson%hy,poisson%hz
    !---------------------------------------------------------------------------
    poisson%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
        id='poisson%rho')
    poisson%pot=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
        id='poisson%pot')
    call construct_ewald_bps(parini,atoms,poisson)
    poisson%epotfixed=dot_product(atoms%qat,atoms%qat)/(sqrt(2.d0*pi)*poisson%alpha)
    shortrange_at_rcut=erfc(poisson%linked_lists%rcut/(sqrt(2.d0)*poisson%alpha))/(poisson%linked_lists%rcut)
    if(parini%iverbose>=2) then
        write(*,*) 'real part in rcut',shortrange_at_rcut
    endif
    end associate
    end associate
    end associate
    call f_release_routine()
end subroutine init_hartree_bps
!*****************************************************************************************
subroutine init_hartree_p3d(parini,atoms,poisson)
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
    integer:: ngptot
    integer:: nbgpx, nbgpy, nbgpz
    call f_routine(id='init_hartree_p3d')
    associate(ngpx=>poisson%ngpx)
    associate(ngpy=>poisson%ngpy)
    associate(ngpz=>poisson%ngpz)
    if(trim(atoms%boundcond)/='slab') then
        write(*,*) 'ERROR: init_hartree_p3d currently assumes BC=slab'
        stop
    endif
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
    poisson%rgcut=poisson_rough%rgcut
    poisson%spline%nsp=parini%nsp_ewald
    poisson%cell(1)=atoms%cellvec(1,1)
    poisson%cell(2)=atoms%cellvec(2,2)
    poisson%cell(3)=atoms%cellvec(3,3)
    poisson%vu=parini%vu_ewald
    poisson%vl=parini%vl_ewald
    !---------------------------------------------------------------------------
    !call calparam(parini,atoms,poisson_rough,poisson)
    ngpx=int(poisson%cell(1)/poisson_rough%hx)+1
    ngpy=int(poisson%cell(2)/poisson_rough%hx)+1
    ngpz=int(poisson%cell(3)/poisson_rough%hz)+1
    if(mod(ngpx,2)/=0) ngpx=ngpx+1
    if(mod(ngpy,2)/=0) ngpy=ngpy+1
    poisson%hx=poisson%cell(1)/real(ngpx,8)
    poisson%hy=poisson%cell(2)/real(ngpy,8)
    poisson%hz=poisson%cell(3)/real(ngpz-1,8)
    nbgpx=int(poisson_rough%rgcut/poisson%hx)+2
    nbgpy=int(poisson_rough%rgcut/poisson%hy)+2
    nbgpz=int(poisson_rough%rgcut/poisson%hz)+2
    ngpz=ngpz+2*nbgpz
    ngptot=ngpx*ngpy*ngpz
    write(*,'(a50,4i)') 'ngpx,ngpy,ngpz,ngptot',ngpx,ngpy,ngpz,ngptot
    write(*,'(a50,3i)') 'nbgpx,nbgpy,nbgpz',nbgpx,nbgpy,nbgpz
    !write(*,'(a50,3i)') 'nagpx,nagpy,nagpz',poisson%nagpx,poisson%nagpy,poisson%nagpz
    write(*,'(a50,3f14.7)') 'hgx,hgy,hgz',poisson%hx,poisson%hy,poisson%hz
    !---------------------------------------------------------------------------
    poisson%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
        id='poisson%rho')
    poisson%pot=f_malloc([1.to.ngpx+2,1.to.ngpy,1.to.ngpz], &
        id='poisson%pot')
    call init_psolver_p3d_slab(poisson)
    poisson%epotfixed=dot_product(atoms%qat,atoms%qat)/(sqrt(2.d0*pi)*poisson%alpha)
    shortrange_at_rcut=erfc(poisson%linked_lists%rcut/(sqrt(2.d0)*poisson%alpha))/(poisson%linked_lists%rcut)
    if(parini%iverbose>=2) then
        write(*,*) 'real part in rcut',shortrange_at_rcut
    endif
    end associate
    end associate
    end associate
    call f_release_routine()
end subroutine init_hartree_p3d
!*****************************************************************************************
subroutine get_hartree_simple(parini,poisson,atoms,gausswidth,ehartree)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree
    !local variables
    real(8):: dpm, pi !, gtot, epotreal
    integer:: iat
    !real(8), allocatable:: gwsq(:), ratred(:,:), gg(:) 
    !real(8), allocatable::  ewaldwidth(:)
    !real(8):: stress(3,3), kmax, c, vol, talpha
    pi=4.d0*atan(1.d0)
! !   if (parini%ewald .and. parini%alpha_ewald<0.d0) then
!        call getvol_alborz(atoms%cellvec,vol)
!        c=2
!        write(*,*)"alpha optimize", 1.d0/(sqrt(pi)*(atoms%nat/vol**2)**(1.d0/6.d0))
!         
! !   end if


    if(trim(parini%psolver_ann)/='kwald') then
        dpm=0.d0
        do iat=1,atoms%nat
            dpm=dpm+atoms%qat(iat)*atoms%rat(3,iat)
        enddo
        dpm=dpm*2.d0*pi*poisson%ngpx*poisson%ngpy/(poisson%cell(1)*poisson%cell(2))
        !do iat=1,atoms%nat
        !    write(33,'(2i4,3es14.5)') iter,iat,atoms%qat(iat),atoms%rat(3,iat),dpm
        !enddo
    endif
    !-----------------------------------------------------------------
    if(poisson%cal_rho) then
        select case(trim(parini%psolver_ann))
            case('kwald')
                !do nothing
            case('bigdft','p3d')
                call putgaussgrid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat, &
                    atoms%qat,gausswidth,poisson)
            case default
                write(*,*) 'ERROR: unknown method for hartree calculation.'
                stop
        end select
    endif
    !-----------------------------------------------------------------
    !Even if cal_poisson is false, psolver_bulk_fourier must be called
    !once more because fat is set to zero after dU/dq=0 in CENT
    if(poisson%cal_poisson .or. trim(parini%psolver_ann)=='kwald') then
        select case(trim(parini%psolver_ann))
            case('kwald')
                call psolver_bulk_fourier(parini,poisson,atoms,gausswidth, &
                    ehartree,poisson%qgrad)
            case('bigdft')
                call psolver_bps(poisson,atoms,ehartree)
            case('p3d')
                    call psolver_p3d(parini,poisson,atoms,ehartree,dpm)
            case default
                write(*,*) 'ERROR: unknown method for hartree calculation.'
                stop
        end select
    endif
    !-----------------------------------------------------------------
    if(poisson%cal_qgrad) then
        select case(trim(parini%psolver_ann))
            case('kwald')
                !do nothing
            case('bigdft','p3d')
                call get_g_from_pot(parini,atoms,poisson,gausswidth,poisson%qgrad)
                call apply_external_field(parini,atoms,poisson,ehartree,poisson%qgrad)
            case default
                write(*,*) 'ERROR: unknown method for hartree calculation.'
                stop
        end select
    endif
    !-----------------------------------------------------------------
    if(poisson%cal_force) then
        select case(trim(parini%psolver_ann))
            case('kwald')
                !do nothing
            case('bigdft','p3d')
                call longerange_forces(parini,atoms,poisson,gausswidth)
            case default
                write(*,*) 'ERROR: unknown method for hartree calculation.'
                stop
        end select
    endif
    !-----------------------------------------------------------------
end subroutine get_hartree_simple
!*****************************************************************************************
subroutine get_hartree(parini,poisson,atoms,gausswidth,ehartree)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree
    !local variables
    real(8):: epotreal !, pi
    !integer:: iat
    real(8), allocatable:: gg(:)
    real(8), allocatable::  ewaldwidth(:)
    real(8):: stress(3,3), kmax !, talpha
    call f_routine(id='get_hartree')
    call f_timing(TCAT_PSOLVER,'ON')

    epotreal=0.d0
    ewaldwidth=f_malloc([1.to.atoms%nat],id='ewaldwidth')
    if(parini%ewald) then
        gg=f_malloc([1.to.atoms%nat],id='gg')
       ewaldwidth(:)=poisson%alpha
    end if

    poisson%qgrad(1:atoms%nat)=0.d0
    if(parini%ewald) then
        !kmax=2.d0/poisson%alpha*sqrt(-log(1.d-3))
        kmax=2.d0/poisson%alpha*sqrt(-log(1.d3*parini%tolerance_ewald))
        poisson%ecut=kmax**2/2.d0
        poisson%gw_identical=.true.
        write(*,*)"ecut", poisson%ecut, "alpha",poisson%alpha
    else
        poisson%ecut=parini%ecut_ewald
    !!!!!  else
    !!!!!      talpha=minval(gausswidth)
    !!!!!      kmax=2.d0/talpha*sqrt(-log(1.d3*parini%tolerance_ewald))
    !!!!!      !kmax=2.d0/talpha*sqrt(-log(1.d-3))
    !!!!!      ecut=kmax**2/2.d0*(4.d0*pi**2)
    !!!!!      write(*,*)"ecut", ecut, "alpha",poisson%alpha
    endif

    if(.not. parini%ewald) then
        ewaldwidth=gausswidth
    endif
    
    call get_hartree_simple(parini,poisson,atoms,ewaldwidth,ehartree)
    if(parini%ewald) then
        call real_part(parini,atoms,gausswidth,poisson%alpha,epotreal,gg,stress)
        atoms%stress=atoms%stress+stress
        ehartree=ehartree+epotreal
        poisson%qgrad=poisson%qgrad+gg
    end if

    if(parini%ewald) then
        call f_free(gg)
        call f_free(ewaldwidth)
    endif
    call f_timing(TCAT_PSOLVER,'OF')
    call f_release_routine()
end subroutine get_hartree
!*****************************************************************************************
subroutine apply_external_field(parini,atoms,poisson,ehartree,g)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    !use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: ehartree, g(atoms%nat)
    !local variables
    !real(8):: dpm, pi, gtot, ecut, epotreal, alphasq
    integer:: iat, igpx, igpy, igpz
    !real(8), allocatable:: gwsq(:), ratred(:,:), gg(:) 
    !real(8), allocatable::  ewaldwidth(:)
    real(8):: dipole !,ext_pot
    if(trim(parini%psolver_ann)/='p3d') then
        write(*,*) 'ERROR: currently external electric field is supposed to be'            
        write(*,*) '       used together with the P3D method.'
        stop
    endif
    if((.not. poisson%point_particle) .and. trim(parini%bias_type)=='fixed_efield') then
        do igpz=1,poisson%ngpz
            !igpz=1 is not necessarily z=0, now in this way external potential is not
            !zero at z=0 but since shift the potential, such that external potential
            !to be zero at z=0, has no effect on energy, we keep it this way.
            !ext_pot=parini%efield*igpz*poisson%hz
            !do igpy=1,poisson%ngpy
            !    do igpx=1,poisson%ngpx
            !        poisson%pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)+ext_pot
            !    enddo
            !enddo
        enddo
        dipole=0.d0
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%rat(3,iat)
        enddo
        ehartree=ehartree+parini%efield*0.5d0*dipole
        do iat=1,atoms%nat
            !atoms%fat(3,iat)=atoms%fat(3,iat)-parini%efield*0.5d0*atoms%qat(iat)
            g(iat)=g(iat)+parini%efield*0.5d0*atoms%rat(3,iat)
        enddo
    elseif((.not. poisson%point_particle) .and. trim(parini%bias_type)=='fixed_potdiff') then
        !efield=0.d0 !to be corrected by Samare
        !do igpz=1,poisson%ngpz
        !    !igpz=1 is not necessarily z=0, now in this way external potential is not
        !    !zero at z=0 but since shift the potential, such that external potential
        !    !to be zero at z=0, has no effect on energy, we keep it this way.
        !    ext_pot=efield*igpz*hz
        !    do igpy=1,poisson%ngpy
        !        do igpx=1,poisson%ngpx
        !            poisson%pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)+ext_pot
        !        enddo
        !    enddo
        !enddo
    endif
end subroutine apply_external_field
!*****************************************************************************************
subroutine real_part(parini,atoms,gausswidth,alpha,epotreal,gg,stress)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_linked_lists, only: typ_linked_lists
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_linked_lists):: linked_lists
    real(8):: dx, dy, dz, r, rsq, xiat, yiat, ziat
    real(8):: t, tt1, tt2, tt3, ttt
    real(8):: rcutsq, fx, fy, fz, pi, qiat, qiatjat
    integer:: ip, jp, jpt, jl, il, iat_maincell, jat_maincell
    integer:: ix, iy, iz, jy, jz, iat, jat, maincell, maincell_iat
    real(8), intent(in):: gausswidth(atoms%nat),alpha
    real(8):: gg(atoms%nat),rr
    real(8):: cell(3) , vol, stress(3,3)
    real(8)::epotreal, alphatwoinv, rbetainv, alphasq, betainv
    pi = 4.d0 *atan(1.d0)
    call getvol_alborz(atoms%cellvec,vol)
    rr=linked_lists%rcut
    linked_lists%rcut =sqrt(2.d0)*max(maxval(gausswidth(:)),alpha)*sqrt(-log(parini%tolerance_ewald))
    call linkedlists_init(parini,atoms,cell,linked_lists)
    rcutsq=linked_lists%rcut**2
    alphatwoinv =1.d0/(sqrt(2.d0)*alpha)
    alphasq=(alphatwoinv)**2
    epotreal=0.d0
    gg=0.d0
    tt2=2.d0/sqrt(pi)
    stress = 0.d0
    include 'act1_cell_linkedlist.inc'
    do  iat=ip,il
        qiat=linked_lists%qat(iat)
        xiat=linked_lists%rat(1,iat)
        yiat=linked_lists%rat(2,iat)
        ziat=linked_lists%rat(3,iat)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        maincell_iat=linked_lists%maincell(iat)
        do  jat=jp,jl
            dx=xiat-linked_lists%rat(1,jat)
            dy=yiat-linked_lists%rat(2,jat)
            dz=ziat-linked_lists%rat(3,jat)
            rsq= dx*dx+dy*dy+dz*dz
            maincell=maincell_iat+linked_lists%maincell(jat)
            if(rsq<rcutsq .and. maincell >-1) then
                iat_maincell=linked_lists%perm(iat)
                jat_maincell=linked_lists%perm(jat)
                qiatjat=linked_lists%qat(jat)*qiat
                r=sqrt(rsq)
                betainv=1.d0/sqrt(gausswidth(iat_maincell)**2+gausswidth(jat_maincell)**2)
                rbetainv=r*betainv
                tt1= 1.d0/r *(-erfc(rbetainv)+erfc(r*alphatwoinv))
                epotreal = epotreal + tt1*qiatjat
                gg(iat_maincell) = gg(iat_maincell)+tt1*linked_lists%qat(jat)
                gg(jat_maincell) = gg(jat_maincell)+tt1*qiat

                tt3=tt2*(exp(-rsq*alphasq)*alphatwoinv-exp(-rbetainv**2)*betainv)
                ttt=(tt1+tt3)/rsq*qiatjat
                fx=ttt*dx;fy=ttt*dy;fz=ttt*dz
                !-----------------------------------
                atoms%fat(1,iat_maincell)=atoms%fat(1,iat_maincell)+fx
                atoms%fat(2,iat_maincell)=atoms%fat(2,iat_maincell)+fy
                atoms%fat(3,iat_maincell)=atoms%fat(3,iat_maincell)+fz
                atoms%fat(1,jat_maincell)=atoms%fat(1,jat_maincell)-fx
                atoms%fat(2,jat_maincell)=atoms%fat(2,jat_maincell)-fy
                atoms%fat(3,jat_maincell)=atoms%fat(3,jat_maincell)-fz
                stress(1,1)=stress(1,1)+fx*dx
                stress(1,2)=stress(1,2)+fx*dy
                stress(1,3)=stress(1,3)+fx*dz
                stress(2,1)=stress(2,1)+fy*dx
                stress(2,2)=stress(2,2)+fy*dy
                stress(2,3)=stress(2,3)+fy*dz
                stress(3,1)=stress(3,1)+fz*dx
                stress(3,2)=stress(3,2)+fz*dy
                stress(3,3)=stress(3,3)+fz*dz
            endif
        enddo
    enddo
    include 'act2_cell_linkedlist.inc'
    stress=stress/vol
    tt2=0.d0
    t=sqrt(2.d0/pi)
    do iat=1,atoms%nat
        tt2=tt2+(1/gausswidth(iat)-1/alpha)*atoms%qat(iat)**2
        gg(iat)=gg(iat)+(1/gausswidth(iat)-1/alpha)*atoms%qat(iat)*t
    end do
    tt2=tt2/sqrt(2*pi)
    epotreal = epotreal+tt2
    linked_lists%rcut=rr

end subroutine real_part
!*****************************************************************************************
