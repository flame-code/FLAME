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
    call psolver_p3d_slab(parini,poisson,poisson%cell,poisson%hx,poisson%hy,poisson%hz,epotlong,beta)
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
