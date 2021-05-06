!*****************************************************************************************
subroutine init_hartree(parini,atoms,poisson,gausswidth)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson):: poisson_rough
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: gausswidth(atoms%nat)
    !local variables
    real(8):: kmax
    call f_routine(id='init_hartree')
    if(poisson%initialized) then
        stop 'ERRRO: poisson is already initialized.'
    endif
    if(trim(parini%psolver)=='p3d') then
        call init_hartree_p3d(parini,atoms,poisson)
    elseif(trim(parini%psolver)=='bigdft') then
        call init_hartree_bps(parini,atoms,poisson)
    elseif(trim(parini%psolver)=='kwald') then
        !It seems that nothing needs to be done for kwald method in the
        !init_hartree routine.
    endif
    if(allocated(poisson%qgrad)) deallocate(poisson%qgrad)
    allocate(poisson%qgrad(1:atoms%nat))

    poisson%rgrad=f_malloc([1.to.3,1.to.atoms%nat],id='poisson%rgrad')
    poisson%rcart=f_malloc([1.to.3,1.to.atoms%nat],id='poisson%rcart')
    poisson%q=f_malloc([1.to.atoms%nat],id='poisson%q')
    poisson%gw=f_malloc([1.to.atoms%nat],id='poisson%gw')
    poisson%gw_ewald=f_malloc([1.to.atoms%nat],id='poisson%gw_ewald')
    if(parini%ewald) then
        poisson%qgrad_real=f_malloc([1.to.atoms%nat],id='poisson%qgrad_real')
        poisson%gw_ewald(:)=poisson%alpha
    end if

    if(parini%ewald .and. parini%ecut_auto) then
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
        poisson%gw_ewald=gausswidth
    endif
    poisson%initialized=.true.
    call f_release_routine()
end subroutine init_hartree
!*****************************************************************************************
subroutine fini_hartree(parini,atoms,poisson)
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
    if(.not. poisson%initialized) then
        stop 'ERRRO: poisson cannot be finalized because it is not initialized.'
    endif
    if(parini%ewald) then
        call f_free(poisson%qgrad_real)
    endif
    call f_free(poisson%rgrad)
    call f_free(poisson%rcart)
    call f_free(poisson%q)
    call f_free(poisson%gw)
    call f_free(poisson%gw_ewald)
    if(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk') then
        if(trim(atoms%boundcond)=='bulk') then
            if(trim(parini%psolver)=='bigdft') then
                call fini_psolver_bps(poisson)
            endif
        elseif(trim(atoms%boundcond)=='slab') then
            call fini_psolver_p3d(poisson)
        endif
        if(trim(parini%psolver)/='kwald') then
            call f_free(poisson%rho)
            call f_free(poisson%pot)
        endif
        if(trim(parini%bias_type)=='p3dbias') then
         !   deallocate(poisson%pots)
        endif
    endif
    poisson%initialized=.false.
    call f_release_routine()
end subroutine fini_hartree
!*****************************************************************************************
subroutine init_hartree_bps(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    include 'fftw3.f'
    type(typ_poisson):: poisson_rough
    real(8):: pi, shortrange_at_rcut
    real(8):: tt1, tt2
    integer:: nbgpx, nbgpy, nbgpz
    integer:: ngptot, ind
    call f_routine(id='init_hartree_bps')
    associate(ngpx=>poisson%ngpx)
    associate(ngpy=>poisson%ngpy)
    associate(ngpz=>poisson%ngpz)
    pi=4.d0*atan(1.d0)
    ind=index(poisson%task_finit,'set_ngp')
    if(ind>0) then
    write(*,*) 'ECUT',parini%ecut_ewald
    poisson_rough%hx=pi/sqrt(2.d0*parini%ecut_ewald)
    poisson_rough%hy=pi/sqrt(2.d0*parini%ecut_ewald)
    poisson_rough%hz=pi/sqrt(2.d0*parini%ecut_ewald)
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
        !write(*,*) 'BBBBBBBBBBBB'
        !write(*,*) ngpx,ngpy,ngpz
        !stop
    if(trim(atoms%boundcond)=='bulk') then
    poisson%hgrid(1,1)=atoms%cellvec(1,1)/ngpx ; poisson%hgrid(2,1)=atoms%cellvec(2,1)/ngpx ; poisson%hgrid(3,1)=atoms%cellvec(3,1)/ngpx
    poisson%hgrid(1,2)=atoms%cellvec(1,2)/ngpy ; poisson%hgrid(2,2)=atoms%cellvec(2,2)/ngpy ; poisson%hgrid(3,2)=atoms%cellvec(3,2)/ngpy
    poisson%hgrid(1,3)=atoms%cellvec(1,3)/ngpz ; poisson%hgrid(2,3)=atoms%cellvec(2,3)/ngpz ; poisson%hgrid(3,3)=atoms%cellvec(3,3)/ngpz
    elseif(trim(atoms%boundcond)=='free') then
    if(atoms%cellvec(2,1)/=0.d0 .or. atoms%cellvec(3,1)/=0.d0 .or. &
       atoms%cellvec(1,2)/=0.d0 .or. atoms%cellvec(3,2)/=0.d0 .or. &
       atoms%cellvec(1,3)/=0.d0 .or. atoms%cellvec(2,3)/=0.d0) then
        write(*,'(a)') 'ERROR: this routine is only for orthogonal cell:'
        write(*,'(3es14.5)') atoms%cellvec(1,1),atoms%cellvec(2,1),atoms%cellvec(3,1)
        write(*,'(3es14.5)') atoms%cellvec(1,2),atoms%cellvec(2,2),atoms%cellvec(3,2)
        write(*,'(3es14.5)') atoms%cellvec(1,3),atoms%cellvec(2,3),atoms%cellvec(3,3)
        stop
    endif
    poisson%hgrid=0.d0
    poisson%hgrid(1,1)=atoms%cellvec(1,1)/(ngpx-1)
    poisson%hgrid(2,2)=atoms%cellvec(2,2)/(ngpy-1)
    poisson%hgrid(3,3)=atoms%cellvec(3,3)/(ngpz-1)
    nbgpx=int(poisson_rough%rgcut/poisson%hgrid(1,1))+2
    nbgpy=int(poisson_rough%rgcut/poisson%hgrid(2,2))+2
    nbgpz=int(poisson_rough%rgcut/poisson%hgrid(3,3))+2
    ngpx=ngpx+2*nbgpx
    ngpy=ngpy+2*nbgpy
    ngpz=ngpz+2*nbgpz
    else
        write(*,*) 'ERROR: init_hartree_bps currently assumes BC=bulk or free'
        stop
    endif
    poisson%lda=ngpx !To be corrected for BCs other than bulk, e.g. +2 for slab+P3D
    ngptot=ngpx*ngpy*ngpz
    call yaml_mapping_open('grid info',flow=.true.)
    call yaml_map('ngpx',ngpx,fmt='(i8)')
    call yaml_map('ngpy',ngpy,fmt='(i8)')
    call yaml_map('ngpz',ngpz,fmt='(i8)')
    call yaml_map('ngptot',ngptot,fmt='(i8)')
    call yaml_map('hxx',poisson%hgrid(1,1),fmt='(f20.10)')
    call yaml_map('hyx',poisson%hgrid(2,1),fmt='(f20.10)')
    call yaml_map('hzx',poisson%hgrid(3,1),fmt='(f20.10)')
    call yaml_map('hxy',poisson%hgrid(1,2),fmt='(f20.10)')
    call yaml_map('hyy',poisson%hgrid(2,2),fmt='(f20.10)')
    call yaml_map('hzy',poisson%hgrid(3,2),fmt='(f20.10)')
    call yaml_map('hxz',poisson%hgrid(1,3),fmt='(f20.10)')
    call yaml_map('hyz',poisson%hgrid(2,3),fmt='(f20.10)')
    call yaml_map('hzz',poisson%hgrid(3,3),fmt='(f20.10)')
    !write(*,'(a50,4i)') 'ngpx,ngpy,ngpz,ngptot',ngpx,ngpy,ngpz,ngptot
    !write(*,'(a50,3i)') 'nbgpx,nbgpy,nbgpz',nbgpx,nbgpy,nbgpz
    !write(*,'(a50,3i)') 'nagpx,nagpy,nagpz',poisson%nagpx,poisson%nagpy,poisson%nagpz
    !write(*,'(a,3f20.10)') 'hxx,hyx,hzx ',poisson%hgrid(1,1),poisson%hgrid(2,1),poisson%hgrid(3,1)
    !write(*,'(a,3f20.10)') 'hxy,hyy,hzy ',poisson%hgrid(1,2),poisson%hgrid(2,2),poisson%hgrid(3,2)
    !write(*,'(a,3f20.10)') 'hxz,hyz,hzz ',poisson%hgrid(1,3),poisson%hgrid(2,3),poisson%hgrid(3,3)
    call yaml_mapping_close()
    endif
    !---------------------------------------------------------------------------
    ind=index(poisson%task_finit,'alloc_rho')
    if(ind>0) then
    poisson%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
        id='poisson%rho')
    endif
    poisson%pot=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
        id='poisson%pot')
    call init_psolver_bps(parini,atoms,poisson)
    poisson%epotfixed=dot_product(atoms%qat,atoms%qat)/(sqrt(2.d0*pi)*poisson%alpha)
    shortrange_at_rcut=erfc(poisson%linked_lists%rcut/(sqrt(2.d0)*poisson%alpha))/(poisson%linked_lists%rcut)
    if(parini%iverbose>=2) then
        call yaml_map('real part in rcut',shortrange_at_rcut)
        !write(*,*) 'real part in rcut',shortrange_at_rcut
    endif
    end associate
    end associate
    end associate
    call f_release_routine()
end subroutine init_hartree_bps
!*****************************************************************************************
subroutine init_hartree_p3d(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson):: poisson_rough
    type(typ_poisson), intent(inout):: poisson
    !local variables
    include 'fftw3.f'
    real(8):: pi, shortrange_at_rcut
    integer:: ngptot, ind
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
    ind=index(poisson%task_finit,'set_ngp')
    if(ind>0) then
    poisson_rough%hx=pi/sqrt(2.d0*parini%ecut_ewald)
    poisson_rough%hy=pi/sqrt(2.d0*parini%ecut_ewald)
    poisson_rough%hz=pi/sqrt(2.d0*parini%ecutz_ewald)
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
    poisson%hgrid=0.d0
    poisson%hgrid(1,1)=poisson%cell(1)/real(ngpx,8)
    poisson%hgrid(2,2)=poisson%cell(2)/real(ngpy,8)
    poisson%hgrid(3,3)=poisson%cell(3)/real(ngpz-1,8)
    nbgpx=int(poisson_rough%rgcut/poisson%hx)+2
    nbgpy=int(poisson_rough%rgcut/poisson%hy)+2
    nbgpz=int(poisson_rough%rgcut/poisson%hz)+2
    ngpz=ngpz+2*nbgpz
    ngptot=ngpx*ngpy*ngpz
    call yaml_mapping_open('grid info',flow=.true.)
    call yaml_map('ngpx',ngpx,fmt='(i8)')
    call yaml_map('ngpy',ngpy,fmt='(i8)')
    call yaml_map('ngpz',ngpz,fmt='(i8)')
    call yaml_map('ngptot',ngptot,fmt='(i8)')
    call yaml_map('nbgpx',nbgpx,fmt='(i8)')
    call yaml_map('nbgpy',nbgpy,fmt='(i8)')
    call yaml_map('nbgpz',nbgpz,fmt='(i8)')
    call yaml_map('hxx',poisson%hgrid(1,1),fmt='(f20.10)')
    call yaml_map('hyx',poisson%hgrid(2,1),fmt='(f20.10)')
    call yaml_map('hzx',poisson%hgrid(3,1),fmt='(f20.10)')
    call yaml_map('hxy',poisson%hgrid(1,2),fmt='(f20.10)')
    call yaml_map('hyy',poisson%hgrid(2,2),fmt='(f20.10)')
    call yaml_map('hzy',poisson%hgrid(3,2),fmt='(f20.10)')
    call yaml_map('hxz',poisson%hgrid(1,3),fmt='(f20.10)')
    call yaml_map('hyz',poisson%hgrid(2,3),fmt='(f20.10)')
    call yaml_map('hzz',poisson%hgrid(3,3),fmt='(f20.10)')
    !write(*,'(a50,4i)') 'ngpx,ngpy,ngpz,ngptot',ngpx,ngpy,ngpz,ngptot
    !write(*,'(a50,3i)') 'nbgpx,nbgpy,nbgpz',nbgpx,nbgpy,nbgpz
    !write(*,'(a50,3i)') 'nagpx,nagpy,nagpz',poisson%nagpx,poisson%nagpy,poisson%nagpz
    !write(*,'(a50,3f14.7)') 'hgx,hgy,hgz',poisson%hx,poisson%hy,poisson%hz
    !write(*,'(a,3f20.10)') 'hxx,hyx,hzx ',poisson%hgrid(1,1),poisson%hgrid(2,1),poisson%hgrid(3,1)
    !write(*,'(a,3f20.10)') 'hxy,hyy,hzy ',poisson%hgrid(1,2),poisson%hgrid(2,2),poisson%hgrid(3,2)
    !write(*,'(a,3f20.10)') 'hxz,hyz,hzz ',poisson%hgrid(1,3),poisson%hgrid(2,3),poisson%hgrid(3,3)
    call yaml_mapping_close()
    endif
    poisson%lda=ngpx+2
    !---------------------------------------------------------------------------
    ind=index(poisson%task_finit,'alloc_rho')
    if(ind>0) then
    poisson%rho=f_malloc([1.to.ngpx,1.to.ngpy,1.to.ngpz], &
        id='poisson%rho')
    endif
    poisson%pot=f_malloc([1.to.poisson%lda,1.to.ngpy,1.to.ngpz], &
        id='poisson%pot')
    call init_psolver_p3d(poisson)
    poisson%epotfixed=dot_product(atoms%qat,atoms%qat)/(sqrt(2.d0*pi)*poisson%alpha)
    shortrange_at_rcut=erfc(poisson%linked_lists%rcut/(sqrt(2.d0)*poisson%alpha))/(poisson%linked_lists%rcut)
    if(parini%iverbose>=2) then
        call yaml_map('real part in rcut',shortrange_at_rcut)
        !write(*,*) 'real part in rcut',shortrange_at_rcut
    endif
    end associate
    end associate
    end associate
    call f_release_routine()
end subroutine init_hartree_p3d
!*****************************************************************************************
subroutine put_charge_density(parini,poisson)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    !use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    !local variables
    if(trim(parini%psolver)=='kwald' .and. trim(poisson%bc)/='bulk') then
        write(*,*) 'ERROR: kwald works only with bulk BC.'
    endif
    if(trim(parini%psolver)=='p3d' .and. trim(poisson%bc)/='slab') then
        write(*,*) 'ERROR: P3D works only with slab BC.'
    endif
    !-----------------------------------------------------------------
    select case(trim(parini%psolver))
        case('kwald')
            !do nothing
        case('bigdft')
            if(parini%cell_ortho) then
                call put_gto_sym_ortho(parini,poisson%bc,poisson%reset_rho,poisson%nat, &
                    poisson%rcart,poisson%q,poisson%gw_ewald,poisson%rgcut, &
                    poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            else
                call put_gto_sym(parini,poisson%bc,poisson%reset_rho,poisson%nat, &
                    poisson%rcart,poisson%q,poisson%gw,poisson%rgcut, &
                    poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            endif
        case('p3d')
            if(parini%cell_ortho) then
                call put_gto_sym_ortho(parini,poisson%bc,poisson%reset_rho,poisson%nat, &
                    poisson%rcart,poisson%q,poisson%gw_ewald,poisson%rgcut, &
                    poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            else
                write(*,*) 'ERROR: P3D works only with orthogonal cell.'
                stop
            endif
        case default
            write(*,*) 'ERROR: unknown method for hartree calculation.'
            stop
    end select
end subroutine put_charge_density
!*****************************************************************************************
subroutine get_psolver(parini,poisson,atoms,gausswidth,ehartree)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree
    !local variables
    select case(trim(parini%psolver))
        case('kwald')
            call get_psolver_fourier(parini,poisson,atoms,gausswidth, &
                ehartree,poisson%qgrad)
        case('bigdft')
            call get_psolver_bps(poisson,atoms,ehartree)
        case('p3d')
                call get_psolver_p3d(parini,poisson,poisson%cell,poisson%hgrid(1,1),poisson%hgrid(2,2),poisson%hgrid(3,3),ehartree)
        case default
            write(*,*) 'ERROR: unknown method for hartree calculation.'
            stop
    end select
end subroutine get_psolver
!*****************************************************************************************
subroutine get_hartree_grad_rho(parini,poisson,atoms,ehartree)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    !use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
    !local variables

    select case(trim(parini%psolver))
        case('kwald')
            !do nothing
        case('bigdft')
            if(parini%cell_ortho) then
                call qgrad_gto_sym_ortho(parini,poisson%bc,poisson%nat,poisson%rcart,poisson%q,poisson%gw_ewald, &
                    poisson%rgcut,poisson%lda,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%pot,poisson%qgrad)
            else
                call rqgrad_gto_sym(parini,poisson%bc,poisson%nat,poisson%rcart,poisson%q,poisson%gw, &
                    poisson%rgcut,poisson%lda,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%pot,poisson%rgrad,poisson%qgrad)
            endif
        case('p3d')
            if(parini%cell_ortho) then
                if(trim(parini%bias_type)/='no') then
                    call apply_external_field(parini,atoms,poisson,ehartree,poisson%qgrad,"qgrad")
                endif
                call qgrad_gto_sym_ortho(parini,poisson%bc,poisson%nat,poisson%rcart,poisson%q,poisson%gw_ewald, &
                    poisson%rgcut,poisson%lda,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%pot,poisson%qgrad)
            else
                write(*,*) 'ERROR: P3D works only with orthogonal cell.'
                stop
            endif
        case default
            write(*,*) 'ERROR: unknown method for hartree calculation.'
            stop
    end select
end subroutine get_hartree_grad_rho
!*****************************************************************************************
subroutine get_hartree_force(parini,poisson,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    !use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8):: ehartree
    !local variables
    integer:: iat
    select case(trim(parini%psolver))
        case('kwald')
            !do nothing
        case('bigdft')
            if(parini%cell_ortho) then
                call force_gto_sym_ortho(parini,poisson%bc,poisson%nat,poisson%rcart, &
                    poisson%q,poisson%gw_ewald,poisson%rgcut,poisson%lda,poisson%ngpx, &
                    poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%pot,atoms%fat)
            else
                call force_gto_sym(parini,poisson%bc,poisson%nat,poisson%rcart, &
                    poisson%q,poisson%gw,poisson%rgcut,poisson%lda,poisson%ngpx, &
                    poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%pot,atoms%fat)
            endif
        case('p3d')
            if(parini%cell_ortho) then
                call apply_external_field(parini,atoms,poisson,ehartree,poisson%qgrad,"force")
                call force_gto_sym_ortho(parini,poisson%bc,poisson%nat,poisson%rcart, &
                    poisson%q,poisson%gw_ewald,poisson%rgcut,poisson%lda,poisson%ngpx, &
                    poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%pot,atoms%fat)
            else
                write(*,*) 'ERROR: P3D works only with orthogonal cell.'
                stop
            endif
        case default
            write(*,*) 'ERROR: unknown method for hartree calculation.'
            stop
    end select
end subroutine get_hartree_force
!*****************************************************************************************
subroutine get_hartree(parini,poisson,atoms,gausswidth,ehartree)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use time_profiling
    !use mod_timing , only: TCAT_PSOLVER
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gausswidth(atoms%nat)
    real(8), intent(out):: ehartree
    !local variables
    real(8):: epotreal !, pi
    integer:: ind
    real(8):: stress(3,3) !, talpha
    call f_routine(id='get_hartree')
    !call f_timing(TCAT_PSOLVER,'ON')
    
    !real(8), allocatable:: gwsq(:), ratred(:,:), gg(:) 
    !real(8):: stress(3,3), kmax, c, vol, talpha
!    pi=4.d0*atan(1.d0)
! !   if (parini%ewald .and. parini%alpha_ewald<0.d0) then
!        call getvol_alborz(atoms%cellvec,vol)
!        c=2
!        write(*,*)"alpha optimize", 1.d0/(sqrt(pi)*(atoms%nat/vol**2)**(1.d0/6.d0))
!         
! !   end if
    !-----------------------------------------------------------------
    if(trim(parini%psolver)=='kwald' .and. trim(atoms%boundcond)/='bulk') then
        write(*,*) 'ERROR: kwald works only with bulk BC.'
    endif
    if(trim(parini%psolver)=='p3d' .and. trim(atoms%boundcond)/='slab') then
        write(*,*) 'ERROR: P3D works only with slab BC.'
    endif
    !-----------------------------------------------------------------
    !Even if cal_poisson is false, get_psolver_fourier must be called
    !once more because fat is set to zero after dU/dq=0 in CENT
    call get_psolver(parini,poisson,atoms,poisson%gw_ewald,ehartree)
    !-----------------------------------------------------------------
    if(parini%ewald .and. (trim(parini%approach_ann)=='eem1' .or. trim(parini%approach_ann)=='cent1')) then
        epotreal=0.d0
        call real_part(parini,atoms,gausswidth,poisson%alpha,epotreal,poisson%qgrad_real,stress)
        atoms%stress=atoms%stress+stress
        ehartree=ehartree+epotreal
        poisson%qgrad=poisson%qgrad+poisson%qgrad_real
    end if

    !call f_timing(TCAT_PSOLVER,'OF')
    call f_release_routine()
end subroutine get_hartree
!*****************************************************************************************
subroutine apply_external_field(parini,atoms,poisson,ehartree,g,flag)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_poisson
    use yaml_output
    !use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: ehartree, g(atoms%nat)
    !local variables
    integer:: iat, igpx, igpy, igpz
    integer:: nbgpz,nlayer                                        
    real(8):: dipole !,ext_pot
    real(8):: dv, pi, beta ,vu,vl
    real(8):: c, charge, epot, ext_pot , efield, charge0
    real(8):: qv, qv_pp, Ef
    real(8),allocatable :: pot_short(:,:,:,:)
    character(5)::flag
    pi=4.d0*atan(1.d0)        
    if(trim(parini%psolver)/='p3d') then
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
        call update_ratp(atoms)
        dipole=0.d0
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%ratp(3,iat)
        enddo
        ehartree=ehartree-parini%efield*dipole
        do iat=1,atoms%nat
            atoms%fat(3,iat)=atoms%fat(3,iat)+parini%efield*atoms%qat(iat)
            g(iat)=g(iat)-parini%efield*atoms%ratp(3,iat)
        enddo
        if (flag=="force") then
            call yaml_mapping_open('dipole for fixed_efield',flow=.true.)
            call yaml_map('dipole',dipole)
            call yaml_map('K',4.d0*pi*dipole/(parini%efield*(poisson%cell(1)*poisson%cell(2)*(poisson%cell(3))))+1)
            call yaml_mapping_close()
            !write(*,*)"dipole ,K =  " , dipole , 4*pi*dipole/(parini%efield*(poisson%cell(1)*poisson%cell(2)*(poisson%cell(3))))+1
        endif
    elseif((.not. poisson%point_particle) .and. trim(parini%bias_type)=='fixed_potdiff') then
        nbgpz=int(poisson%rgcut/poisson%hz)+2
        dipole=0.d0
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%ratp(3,iat)
        enddo
        beta=-dipole*2.d0*pi/(poisson%cell(1)*poisson%cell(2)) 
        dv=(poisson%vu-poisson%vl) 
        write(*,*)'real pot = vu , vl ', poisson%vu+beta, poisson%vl-beta
        efield =- (dv+2.d0*beta)/ poisson%cell(3)

        c=poisson%cell(1)*poisson%cell(2)/(4.d0*pi*poisson%cell(3))
        charge=-dipole/poisson%cell(3)+c*dv

        ehartree = ehartree +(dv+beta)*dipole/poisson%cell(3)
        atoms%ebattery=dipole/poisson%cell(3)*dv
        !ehartree = ehartree + 0.5*c*dv**2

        do iat=1,atoms%nat
            g(iat)=g(iat)+(dv+2.d0*beta)/poisson%cell(3) * atoms%ratp(3,iat)!-atoms%rat(3,iat)/poisson%cell(3)*(dv)
            atoms%fat(3,iat)=atoms%fat(3,iat)-((dv+2.d0*beta)/poisson%cell(3))*atoms%qat(iat)!+atoms%qat(iat)/poisson%cell(3)*(dv)
        enddo
        if (flag=="force") then
            call yaml_mapping_open('dipole for fixed_potdiff',flow=.true.)
            call yaml_map('dipole',dipole)
            call yaml_map('K',-4.d0*pi*dipole/(dv*(poisson%cell(1)*poisson%cell(2)))+1)
            call yaml_mapping_close()
            !write(*,*)"dipole ,K =  " , dipole , -4*pi*dipole/(dv*(poisson%cell(1)*poisson%cell(2)))+1
        endif


    elseif((.not. poisson%point_particle) .and. trim(parini%bias_type)=='fixed_potdiff2') then
        vu=parini%vu_ewald+parini%vu_ac_ewald*sin(2*pi*parini%frequency_ewald*parini%time_dynamics)
        vl=parini%vl_ewald
        ehartree=0.d0
        nbgpz=int(poisson%rgcut/poisson%hz)+2
        poisson%npu=poisson%ngpz-nbgpz
        poisson%npl=1+nbgpz  

        allocate(poisson%pots(1:poisson%ngpx+2,1:poisson%ngpy,poisson%npl:poisson%npu))
        poisson%pots=0.d0
        !if (parini%cal_charge) then 
        !    nlayer=5
        !    allocate(pots_layer(1:poisson%ngpx,1:poisson%ngpy,1:2,1:nlayer))
        !    pots_layer = 0.d0
        !endif
        call sollaplaceq(poisson,poisson%hz,poisson%cell,vl,vu)
        !if (parini%cal_charge) then 
        !    call surface_charge(parini,poisson,pots_layer,vl,vu)
        !    deallocate(pots_layer)
        !endif

        do igpz=1,poisson%npl-1      
            do igpy=1,poisson%ngpy
                do igpx=1,poisson%ngpx
                    poisson%pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)+vl
                enddo
            enddo
        enddo
        do igpz=poisson%npl,poisson%npu       
            do igpy=1,poisson%ngpy
                do igpx=1,poisson%ngpx
                    poisson%pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)+poisson%pots(igpx,igpy,igpz)
                enddo
            enddo
        enddo
        do igpz=poisson%npu+1,poisson%ngpz      
            do igpy=1,poisson%ngpy
                do igpx=1,poisson%ngpx
                    poisson%pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)+vu
                enddo
            enddo
        enddo
        epot=0.d0
        do igpz=1,poisson%ngpz
            do igpy=1,poisson%ngpy
                do igpx=1,poisson%ngpx
                    epot=epot+poisson%pot(igpx,igpy,igpz)*poisson%rho(igpx,igpy,igpz)
                enddo
            enddo
        enddo
        epot=epot*0.5d0*poisson%hx*poisson%hy*poisson%hz
        ehartree=ehartree+epot
        dipole=0.d0
        do iat=1,atoms%nat
            dipole=dipole+atoms%qat(iat)*atoms%ratp(3,iat)
        enddo
        !write(*,*)"dipole  =  " , dipole
        beta=-dipole*2.d0*pi/(poisson%cell(1)*poisson%cell(2)) 
        dv=vu-vl
        c=poisson%cell(1)*poisson%cell(2)/(4.d0*pi*poisson%cell(3))
        charge0=-dipole/poisson%cell(3)
        charge=-dipole/poisson%cell(3)+c*dv

        ehartree=ehartree-0.5d0*charge0*dv!+0.5d0*c*dv**2!
        atoms%ebattery=-charge0*dv

        if (flag=="force") then
            call yaml_mapping_open('dipole for fixed_potdiff2',flow=.true.)
            call yaml_map('dipole',dipole)
            call yaml_map('K',-4.d0*pi*dipole/((dv)*(poisson%cell(1)*poisson%cell(2)))+1)
            call yaml_mapping_close()
            write(*,*)"dipole ,beta =  " , dipole,beta
            write(*,*)"charge =                ",charge, charge0
            write(*,*)"externalwork = ",atoms%ebattery
            !write(*,*)"dipole ,K =  " , dipole , -4*pi*dipole/((dv)*(poisson%cell(1)*poisson%cell(2)))+1
           ! allocate(pot_short(poisson%ngpx,poisson%ngpy,2,5))
           ! pot_short=0.d0
           ! poisson%pots=0.d0
           ! call surface_charge(parini,poisson,pot_short,vl,vu)
            !write(*,*)"C_0 = ",c, " K =" , charge/dv/c
           ! deallocate(pot_short)
                !open(unit=55, file="pots.txt" ,status='unknown',position='append')
                !   do igpy=1,min(10,poisson%ngpy),2
                !   do igpx=1,min(10,poisson%ngpx),2
                !   do igpz=poisson%npl,poisson%npu
                !       !write(55,*)  (iz-1-nbgpz)*poisson%hz, -poisson%pots(ix,iy,iz)+&
                !       !             (-2*beta/(poisson%ngpx*poisson%ngpy)+vu-vl)/d*(iz-1-nbgpz)*poisson%hz+beta/(poisson%ngpx*poisson%ngpy)+vl
                !       write(55,'(4es20.8)')  (igpz-1-nbgpz)*poisson%hz, poisson%pots(igpx,igpy,igpz),&
                !                    (vu-vl)/poisson%cell(3)*(igpz-1-nbgpz)*poisson%hz+vl,&
                !                    (2*beta+vu-vl)/poisson%cell(3)*(igpz-1-nbgpz)&
                !                    *poisson%hz-beta+vl

                !   enddo 
                !       write(55,*)   
                !   enddo 
                !   enddo 
                !close(55)
                open(unit=56, file="field.txt" ,status='unknown',position='append')
                !   do igpy=1,min(10,poisson%ngpy),2
                !   do igpx=1,min(10,poisson%ngpx),2
                !   do igpz=poisson%npl+1,poisson%npu-1
                !       Ef = (poisson%pots(igpx,igpy,igpz+1)-poisson%pots(igpx,igpy,igpz-1))/(2.d0*poisson%hz)
                !       write(56,'(4es20.8)')  (igpz-1-nbgpz)*poisson%hz, Ef ,&
                !                    (vu-vl)/poisson%cell(3),&
                !                    (2*beta+vu-vl)/poisson%cell(3)

                !   enddo 
                !       write(56,*)   
                !   enddo 
                !   enddo 
                       write(56,'(4es20.8)')  (vu-vl)/poisson%cell(3), (2*beta+vu-vl)/poisson%cell(3)
                close(56)
        endif

        do iat=1,atoms%nat
            g(iat)=g(iat)!-atoms%rat(3,iat)/poisson%cell(3)*(dv)
         !   g(iat)=g(iat)+(dv+2.d0*beta)/poisson%cell(3) * atoms%rat(3,iat)!-atoms%rat(3,iat)/poisson%cell(3)*(dv)
        !    atoms%fat(3,iat)=atoms%fat(3,iat)!+atoms%qat(iat)/poisson%cell(3)*(dv)
        enddo
        deallocate(poisson%pots)
    endif
end subroutine apply_external_field
!*****************************************************************************************
subroutine real_part(parini,atoms,gausswidth,alpha,epotreal,gg,stress)
    use mod_parini, only: typ_parini
    use mod_linked_lists, only: typ_linked_lists
    use mod_atoms, only: typ_atoms, update_ratp
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
    call update_ratp(linked_lists%typ_atoms)
    include 'act1_cell_linkedlist.inc'
    do  iat=ip,il
        qiat=linked_lists%qat(iat)
        xiat=linked_lists%ratp(1,iat)
        yiat=linked_lists%ratp(2,iat)
        ziat=linked_lists%ratp(3,iat)
        jpt=linked_lists%prime(ix+linked_lists%limnbx(1,jy-iy,jz-iz),jy,jz)
        jp=(iat-ip+1)*((isign(1,ip-jpt)+1)/2)+jpt
        jl=linked_lists%last(ix+linked_lists%limnbx(2,jy-iy,jz-iz),jy,jz)
        maincell_iat=linked_lists%maincell(iat)
        do  jat=jp,jl
            dx=xiat-linked_lists%ratp(1,jat)
            dy=yiat-linked_lists%ratp(2,jat)
            dz=ziat-linked_lists%ratp(3,jat)
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
