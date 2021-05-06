!*****************************************************************************************
subroutine calculate_forces_energy(parini,poisson,atoms)
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat
    use mod_parini, only: typ_parini
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    !local variables
    real(8):: epotlong, epotplane !, epotshort
    real(8):: time(10)
    !The following dummy varables definition to be deleted later.
    !real(8):: totrho
    real(8):: charge,c,charge0,E,tmp
    integer:: igpx, igpy, igpz, iat
    integer:: ix, iy, iz, jx, jy, jz, kx, ky, kz
    integer:: npl, npu, nlayer, ngpx, ngpy, ngpz
    real(8),allocatable :: pots_layer(:,:,:,:)
    real(8):: vl, vu, A, d, rl, ru, dipole_correction, dipole, epot_dielec
    real(8),allocatable :: gausswidth(:)  
    call f_routine(id='calculate_forces_energy')
    ngpz=poisson%ngpz
    ngpy=poisson%ngpy
    ngpx=poisson%ngpx
    poisson%point_particle= .true.

    gausswidth=f_malloc([1.to.atoms%nat],id='gausswidth')
    gausswidth(:)=poisson%alpha

    call cpu_time(time(1))
    if(.not. poisson%initialized) then
        stop 'ERROR: calculate_forces_energy: poisson is not initialized!'
    endif
    poisson%reset_rho=.true.
    poisson%nat=atoms%nat
    poisson%cv=atoms%cellvec
    poisson%bc=atoms%boundcond
    poisson%q(1:poisson%nat)=atoms%qat(1:atoms%nat)
    poisson%gw(1:poisson%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson%rcart)
    call put_charge_density(parini,poisson)
    call cpu_time(time(2))
    !-----------------------------------------------------------------------
    !totrho=0.d0
    !do iz=1,ngpz;do iy=1,ngpy;do ix=1,poisson%ngpx
    !    totrho=totrho+rho(ix,iy,iz)
    !enddo;enddo;enddo
    !write(*,*) 'totrho',totrho
    !-----------------------------------------------------------------------
    !pot=0.d0
    call cpu_time(time(3))
    call get_hartree(parini,poisson,atoms,gausswidth,epotlong)
    poisson%gw(1:poisson%nat)=poisson%gw_ewald(1:poisson%nat)
    call get_hartree_force(parini,poisson,atoms)
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
    epot_dielec =0.d0
    if(trim(parini%bias_type)=='dielec') then
        call dielec_potener_forces(parini,poisson,atoms,epot_dielec) 
    end if
    call cpu_time(time(6))
    !atoms%epot=epotlong+epotshort-poisson%epotfixed+epotplane
    atoms%epot=epotlong-poisson%epotfixed+epotplane + epot_dielec
    call yaml_map('epotfixed',poisson%epotfixed,fmt='(es23.15))')
    call yaml_map('epotlong',epotlong,fmt='(es23.15))')
    call yaml_map('epotplane',epotplane,fmt='(es23.15))')
    call yaml_map('epottotal',atoms%epot,fmt='(es23.15))')
    !write(*,*) '-----------------------------------------------------------'
    !write(*,'(a50,e32.15)') 'epotfixed',poisson%epotfixed
    !write(*,'(a50,e32.15)') 'epotlong',epotlong
    !write(*,'(a50,e32.15)') 'epotplane',epotplane
    !write(*,'(a50,e32.15)') 'epottotal',atoms%epot
    !write(*,*) '-----------------------------------------------------------'
!    write(*,'(a50,f32.15)') 'Time for put_gto_sym_ortho ',time(2)-time(1)
!    write(*,'(a50,f32.15)') 'Time for long range ', time(3)-time(2)
!    write(*,'(a50,f32.15)') 'Time for long range forces',time(4)-time(3)
!    write(*,'(a50,f32.15)') 'Time for short range',time(5)-time(4)
!    write(*,'(a50,f32.15)') 'Time for plane ',time(6)-time(5)
!    write(*,'(a50,f32.15)') 'Time for Total without plane',time(5)-time(1)
!    write(*,'(a50,f32.15)') 'Time for Total',time(6)-time(1)
    call f_free(gausswidth)
    call f_release_routine()
end subroutine calculate_forces_energy
!*****************************************************************************************
