subroutine solve_poisson(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson_p3d):: poisson_p3d
    type(typ_ewald_p3d):: ewald_p3d
    type(typ_atoms):: atoms
    integer:: istat, igpx, igpy, igpz, iat
    real(8):: cell(3), epot, rgcut_a, t1, t2, t3, t4, pi
    real(8),allocatable::  gausswidth(:)
    pi=4.d0*atan(1.d0)
    call cube_read('rho.cube',atoms,poisson_p3d%typ_poisson)
    !write(*,*) poisson_p3d%ngpx,poisson_p3d%ngpy,poisson_p3d%ngpz
    allocate(poisson_p3d%pot(poisson_p3d%ngpx+2,poisson_p3d%ngpy,poisson_p3d%ngpz),stat=istat)
    if(istat/=0) stop 'ERROR: allocation of pot failed.'
    call ps2dp1df_construction(poisson_p3d)
    cell(1)=poisson_p3d%hx*poisson_p3d%ngpx
    cell(2)=poisson_p3d%hy*poisson_p3d%ngpy
    cell(3)=poisson_p3d%hz*poisson_p3d%ngpz
    !-------------------------------------------------------
    ewald_p3d%poisson_p3d%ngpx=poisson_p3d%ngpx
    ewald_p3d%poisson_p3d%ngpy=poisson_p3d%ngpy
    ewald_p3d%poisson_p3d%ngpz=poisson_p3d%ngpz
    ewald_p3d%hgx=poisson_p3d%hx
    ewald_p3d%hgy=poisson_p3d%hy
    ewald_p3d%hgz=poisson_p3d%hz
    if(.not. parini%gaussian_width>0.d0) then
        stop 'ERROR: gaussian_width must be set.'
    endif
    rgcut_a=parini%gaussian_width !3.d0
    ewald_p3d%nbgpx=int(rgcut_a/ewald_p3d%hgx)+2
    ewald_p3d%nbgpy=int(rgcut_a/ewald_p3d%hgy)+2
    ewald_p3d%nbgpz=int(rgcut_a/ewald_p3d%hgz)+2
    ewald_p3d%poisson_p3d%ngpz=ewald_p3d%poisson_p3d%ngpz+2*ewald_p3d%nbgpz
    t1=real((ewald_p3d%poisson_p3d%ngpx+2*ewald_p3d%nbgpx)*(ewald_p3d%poisson_p3d%ngpy+2*ewald_p3d%nbgpy),8)
    t2=real((ewald_p3d%poisson_p3d%ngpx+2)*(ewald_p3d%poisson_p3d%ngpy),8)
    ewald_p3d%ngpztot=ewald_p3d%poisson_p3d%ngpz*(int(t1/t2)+2)
    allocate(ewald_p3d%poisson_p3d%rho(ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,ewald_p3d%poisson_p3d%ngpz),stat=istat)
    if(istat/=0) stop 'ERROR: allocation of rho failed.'
    write(*,'(a,3i5)') 'ngpx, ngpy, ngpz  ',ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,ewald_p3d%poisson_p3d%ngpz
    write(*,'(a,3i5)') 'nbgpx,nbgpy,nbgpz ',ewald_p3d%nbgpx,ewald_p3d%nbgpy,ewald_p3d%nbgpz
    write(*,'(a,1i5)') 'ngpztot           ',ewald_p3d%ngpztot
    allocate(ewald_p3d%poisson_p3d%pot(ewald_p3d%poisson_p3d%ngpx+2,ewald_p3d%poisson_p3d%ngpy,ewald_p3d%ngpztot),stat=istat)
    if(istat/=0) stop 'ERROR: allocation of pot failed.'
    allocate(ewald_p3d%mboundg(2,-ewald_p3d%nbgpy:ewald_p3d%nbgpy,-ewald_p3d%nbgpz:ewald_p3d%nbgpz),stat=istat)
    if(istat/=0) stop 'Error allocating array mboundg of module ee2dp1df'
    allocate(gausswidth(atoms%nat))
    call determine_glimitsphere(ewald_p3d)
    ewald_p3d%alpha=0.4d0 !atoms%rcov(iat)
    gausswidth=0.4d0 !atoms%rcov(iat)
    atoms%boundcond='slab'
    call putgaussgrid(parini,atoms%boundcond,.true.,atoms%nat,atoms%rat,atoms%qat,gausswidth,ewald_p3d)
    t1=0.d0
    t2=0.d0
    t3=0.d0
    t4=0.d0
    do igpz=1,poisson_p3d%ngpz
        do igpy=1,poisson_p3d%ngpy
            do igpx=1,poisson_p3d%ngpx
                t1=t1+poisson_p3d%rho(igpx,igpy,igpz)
                t2=t2+ewald_p3d%poisson_p3d%rho(igpx,igpy,igpz+ewald_p3d%nbgpz)
                poisson_p3d%rho(igpx,igpy,igpz)=poisson_p3d%rho(igpx,igpy,igpz)-ewald_p3d%poisson_p3d%rho(igpx,igpy,igpz+ewald_p3d%nbgpz)
                t4=t4+(igpz-1)*ewald_p3d%hgz*poisson_p3d%rho(igpx,igpy,igpz)
                t3=t3+poisson_p3d%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    t1=t1*ewald_p3d%hgx*ewald_p3d%hgy*ewald_p3d%hgz
    t2=t2*ewald_p3d%hgx*ewald_p3d%hgy*ewald_p3d%hgz
    t3=t3*ewald_p3d%hgx*ewald_p3d%hgy*ewald_p3d%hgz
    t4=t4*ewald_p3d%hgx*ewald_p3d%hgy*ewald_p3d%hgz
    write(*,*) 't1=',t1
    write(*,*) 't2=',t2
    write(*,*) 't3=',t3
    write(*,*) 't4=',t4
    call cube_write('total_rho.cube',atoms,poisson_p3d%typ_poisson,'rho')
    deallocate(ewald_p3d%poisson_p3d%rho,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of rho failed.'
    deallocate(ewald_p3d%poisson_p3d%pot,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of pot failed.'
    deallocate(ewald_p3d%mboundg,stat=istat)
    if(istat/=0) stop 'Error deallocating array mboundg of module ee2dp1df'
    !-------------------------------------------------------
    t1=0.d0
    do iat=1,atoms%nat
        t1=t1+atoms%qat(iat)*atoms%rat(3,iat)
    enddo
    !t1=t1*(2*pi)/(cell(1)*cell(2))
    write(*,*) 't1=',t1
    call calculate_potener_pot(poisson_p3d,cell,poisson_p3d%hx,poisson_p3d%hy,poisson_p3d%hz,epot)
    call cube_write('pot_p3d.cube',atoms,poisson_p3d%typ_poisson,'pot')
    call ps2dp1df_destruction(poisson_p3d)
    deallocate(poisson_p3d%rho,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of rho failed.'
    deallocate(poisson_p3d%pot,stat=istat)
    if(istat/=0) stop 'ERROR: deallocation of pot failed.'
end subroutine solve_poisson
!*****************************************************************************************
