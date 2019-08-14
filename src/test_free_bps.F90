!*****************************************************************************************
subroutine test_free_bps(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, update_ratp, get_rat
    use dynamic_memory
    use mod_acf, only: acf_read
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson):: poisson_grd
    type(typ_atoms):: atoms
    ! local variables
    integer:: i, j, k
    integer:: iat, ix, iy, iz 
    integer:: nbgx, nbgy, nbgz
    real(8):: ehartree_grd
    real(8):: pi
    real(8):: hx, hy, hz, dx, dy, dz , r 
    real(8):: rho_err_max, temp1,temp2, c1, d1, all_grids
    real(8),allocatable:: gausswidth(:) 
    real(8),allocatable:: pot_anl(:,:,:), pot_grd(:,:,:)
    pi=4.d0*atan(1.d0)
    associate(nat=>atoms%nat)
    associate(ngx=>poisson_grd%ngpx,ngy=>poisson_grd%ngpy,ngz=>poisson_grd%ngpz)
    call acf_read(parini,'posinp.acf',1,atoms=atoms)
    allocate(gausswidth(1:atoms%nat)) 
    atoms%qat(1:atoms%nat) = 1.d0
    poisson_grd%alpha=parini%alpha_ewald
    gausswidth(:) = poisson_grd%alpha
    poisson_grd%task_finit='set_ngp:alloc_rho'
    call init_hartree_bps(parini,atoms,poisson_grd)
    nbgx=int(poisson_grd%rgcut/poisson_grd%hgrid(1,1))+2
    nbgy=int(poisson_grd%rgcut/poisson_grd%hgrid(2,2))+2
    nbgz=int(poisson_grd%rgcut/poisson_grd%hgrid(3,3))+2
    hx=poisson_grd%hgrid(1,1)
    hy=poisson_grd%hgrid(2,2)
    hz=poisson_grd%hgrid(3,3)
    allocate(pot_grd(ngx,ngy,ngz),pot_anl(ngx,ngy,ngz))
    poisson_grd%reset_rho=.true.
    poisson_grd%nat=atoms%nat
    poisson_grd%cv=atoms%cellvec
    poisson_grd%bc=atoms%boundcond
    allocate(poisson_grd%q(1:poisson_grd%nat), poisson_grd%gw_ewald(1:poisson_grd%nat), poisson_grd%rcart(1:3,1:poisson_grd%nat))
    poisson_grd%q(1:poisson_grd%nat)=atoms%qat(1:atoms%nat)
    poisson_grd%gw_ewald(1:poisson_grd%nat)=gausswidth(1:atoms%nat)
    !poisson_grd%rcart(1:3,1:poisson_grd%nat)=atoms%rat(1:3,1:atoms%nat)
    call get_rat(atoms,poisson_grd%rcart)
    call put_charge_density(parini,poisson_grd)
    !-------------------------------------------------------
!    do iat=1,atoms%nat
!        do iz=1,ngz
!            dz=(iz-1-nbgz*1)*hz - atoms%rat(3,iat)
!            do iy=1,ngy
!                dy=(iy-1-nbgy*1)*hy - atoms%rat(2,iat)
!                do ix=1,ngx
!                    dx=(ix-1-nbgx*1)*hx - atoms%rat(1,iat)
!                    rsq=dx**2+dy**2+dz**2
!                    pot_anl(ix,iy,iz)=1.d0/(2.d0**3*pi**1.5d0)*exp(-rsq/2.d0**2) !alpha=2.d0
!                enddo
!            enddo
!        enddo
!    enddo
!    do iz=1,ngz
!        do iy=1,ngy
!            do ix=1,ngx
!                write(95,'(3i4,3es19.10)') ix,iy,iz,pot_anl(ix,iy,iz),poisson_grd%rho(ix,iy,iz),pot_anl(ix,iy,iz)-poisson_grd%rho(ix,iy,iz)
!            enddo
!        enddo
!    enddo
!    stop 'AAA '
    !-------------------------------------------------------
    call get_hartree(parini,poisson_grd,atoms,gausswidth,ehartree_grd)
    write(*,'(a,2es19.10)') "Ehartree, Ehartree_differance : ", ehartree_grd,ehartree_grd-1.d0/(gausswidth(1)*sqrt(2.d0*pi))
    pot_grd(1:ngx,1:ngy,1:ngz) = poisson_grd%pot(1:ngx,1:ngy,1:ngz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call update_ratp(atoms)
    do iat=1,atoms%nat
        do iz=1,ngz
            dz=(iz-1-nbgz*1)*hz - atoms%ratp(3,iat)
            do iy=1,ngy
                dy=(iy-1-nbgy*1)*hy - atoms%ratp(2,iat)
                do ix=1,ngx
                    dx=(ix-1-nbgx*1)*hx - atoms%ratp(1,iat)
                    r=sqrt(dx**2+dy**2+dz**2)
                    if (r <= 1.d-12) then 
                        pot_anl(ix,iy,iz)=2.d0/sqrt(pi)
                    else
                        pot_anl(ix,iy,iz)=atoms%qat(iat)*erf(r/gausswidth(iat))/r
                    endif
                enddo
            enddo
        enddo
    enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rho_err_max=0.d0
    all_grids = ngx*ngy*ngz*1.d0
    all_grids = 1.d0/all_grids
    temp2=0.d0
    do i = 1 , ngz
        do j = 1 , ngy
            do k = 1 , ngx
                c1 = pot_anl(k,j,i)
                d1 = pot_grd(k,j,i)
                temp1 = c1-d1
                temp2 = temp2+(c1-d1)**2
                if (abs(temp1) >= rho_err_max) rho_err_max=abs(temp1)
            end do
        end do
    end do
    temp2 = sqrt(temp2*all_grids)
    write(*,'(a,es19.10)') "Rho_err_max is : ",rho_err_max
    write(*,'(a,es19.10)') "Rho_err_norm is : ",temp2
    end associate
    end associate
    deallocate(poisson_grd%q, poisson_grd%gw_ewald,poisson_grd%rcart,pot_anl,pot_grd,gausswidth)
end subroutine test_free_bps
!*****************************************************************************************
