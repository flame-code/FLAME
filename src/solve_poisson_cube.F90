subroutine solve_poisson(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
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
    !write(*,*) poisson%ngpx,poisson%ngpy,poisson%ngpz
    allocate(poisson%pot(poisson%ngpx+2,poisson%ngpy,poisson%ngpz),stat=istat)
    if(istat/=0) stop 'ERROR: allocation of pot failed.'
    call init_psolver_p3d(poisson)
    poisson%cell(1)=poisson%hx*poisson%ngpx
    poisson%cell(2)=poisson%hy*poisson%ngpy
    poisson%cell(3)=poisson%hz*poisson%ngpz
    !-------------------------------------------------------
    poisson_ion%ngpx=poisson%ngpx
    poisson_ion%ngpy=poisson%ngpy
    poisson_ion%ngpz=poisson%ngpz
    poisson_ion%hx=poisson%hx
    poisson_ion%hy=poisson%hy
    poisson_ion%hz=poisson%hz
    if(.not. parini%gaussian_width>0.d0) then
        stop 'ERROR: gaussian_width must be set.'
    endif
    rgcut_a=parini%gaussian_width !3.d0
    nbgpx=int(rgcut_a/poisson_ion%hx)+2
    nbgpy=int(rgcut_a/poisson_ion%hy)+2
    nbgpz=int(rgcut_a/poisson_ion%hz)+2
    poisson_ion%ngpz=poisson_ion%ngpz+2*nbgpz
    allocate(poisson_ion%rho(poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz),stat=istat)
    if(istat/=0) stop 'ERROR: allocation of rho failed.'
    write(*,'(a,3i5)') 'ngpx, ngpy, ngpz  ',poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz
    write(*,'(a,3i5)') 'nbgpx,nbgpy,nbgpz ',nbgpx,nbgpy,nbgpz
    allocate(poisson_ion%pot(poisson_ion%ngpx+2,poisson_ion%ngpy,poisson_ion%ngpz),stat=istat)
    if(istat/=0) stop 'ERROR: allocation of pot failed.'
    !allocate(poisson_ion%mboundg(2,-nbgpy:nbgpy,-nbgpz:nbgpz),stat=istat)
    !if(istat/=0) stop 'Error allocating array mboundg of module ee2dp1df'
    allocate(gausswidth(atoms%nat))
    poisson_ion%alpha=0.4d0 !atoms%rcov(iat)
    gausswidth=0.4d0 !atoms%rcov(iat)
    atoms%boundcond='slab'
    call put_gto_sym_ortho(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%qat,gausswidth,poisson_ion)
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
                t4=t4+(igpz-1)*poisson_ion%hz*poisson%rho(igpx,igpy,igpz)
                t3=t3+poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    t1=t1*poisson_ion%hx*poisson_ion%hy*poisson_ion%hz
    t2=t2*poisson_ion%hx*poisson_ion%hy*poisson_ion%hz
    t3=t3*poisson_ion%hx*poisson_ion%hy*poisson_ion%hz
    t4=t4*poisson_ion%hx*poisson_ion%hy*poisson_ion%hz
    write(*,*) 't1=',t1
    write(*,*) 't2=',t2
    write(*,*) 't3=',t3
    write(*,*) 't4=',t4
    call cube_write('total_rho.cube',atoms,poisson,'rho')
    deallocate(poisson_ion%rho,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of rho failed.'
    deallocate(poisson_ion%pot,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of pot failed.'
    !deallocate(poisson_ion%mboundg,stat=istat)
    !if(istat/=0) stop 'Error deallocating array mboundg of module ee2dp1df'
    !-------------------------------------------------------
    t1=0.d0
    do iat=1,atoms%nat
        t1=t1+atoms%qat(iat)*atoms%rat(3,iat)
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
    poisson%cal_rho=.false.
    poisson%cal_poisson=.true.
    poisson%cal_qgrad=.false.
    poisson%cal_force=.false.
    call get_hartree(parini,poisson,atoms,gausswidth,epot)
    call cube_write('pot_p3d.cube',atoms,poisson,'pot')
    call fini_psolver_p3d(poisson)
    deallocate(poisson%rho,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of rho failed.'
    deallocate(poisson%pot,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of pot failed.'
end subroutine solve_poisson
!*****************************************************************************************
