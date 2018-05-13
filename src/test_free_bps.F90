subroutine test_free_bps(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson):: poisson_grd
    !type(typ_poisson):: poisson_ana
    type(typ_atoms):: atoms
    ! local variables
    real(8), allocatable :: gausswidth(:) 
    real(8) :: ehartree_grd
    logical :: reset
    integer :: nat, i, j, k
    real(8) :: rxyz(3,1)
    real(8) :: qat(1)
    real(8) :: gw(1)
    real(8) :: rgcut
    real(8) :: hgrid(3,3)
    real(8), allocatable:: wx(:), wy(:), wz(:)
    real(8):: rhoz, rhoyz, pi
    real(8):: hgxinv, hgyinv, hgzinv
    real(8):: width_inv, width_inv_xat, width_inv_yat, width_inv_zat
    real(8):: width_inv_hgx, width_inv_hgy, width_inv_hgz
    real(8):: xat, yat, zat, facqiat, fac, width
    real(8):: hx, hy, hz, dx, dy, dz , r, rsq , r_at
    integer:: iat, iw, ix, iy, iz, iatox, iatoy, iatoz, jx, jy, jz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz
    integer:: ibcx, ibcy, ibcz
    real(8), allocatable:: pot_anl(:,:,:), pot_grd(:,:,:)
    integer, allocatable:: mboundg(:,:,:)
    real(8) :: rho_err_max, temp1,temp2, c1, d1, all_grids
    pi=4.d0*atan(1.d0)
    associate(nat=>atoms%nat)
    associate(ngx=>poisson_grd%ngpx,ngy=>poisson_grd%ngpy,ngz=>poisson_grd%ngpz)
    open(unit=1377,file='params.txt')
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
    !write(*,*) "poisson_grd reset ",poisson_grd%reset_rho 
    poisson_grd%nat=atoms%nat
    !write(*,*) "poisson_grd nat ",poisson_grd%nat
    poisson_grd%cv=atoms%cellvec
    !write(*,*) "poisson_grd CV ",poisson_grd%cv
    poisson_grd%bc=atoms%boundcond
    !write(*,*) "poisson_grd bc ",poisson_grd%bc
    allocate(poisson_grd%q(1:poisson_grd%nat), poisson_grd%gw_ewald(1:poisson_grd%nat), poisson_grd%rcart(1:3,1:poisson_grd%nat))
    poisson_grd%q(1:poisson_grd%nat)=atoms%qat(1:atoms%nat)
    !write(*,*) "poisson_grd q ",poisson_grd%q(1:poisson_grd%nat)
    poisson_grd%gw_ewald(1:poisson_grd%nat)=gausswidth(1:atoms%nat)
    !write(*,*) "poisson_grd gw ",poisson_grd%gw_ewald(1:poisson_grd%nat)
    poisson_grd%rcart(1:3,1:poisson_grd%nat)=atoms%rat(1:3,1:atoms%nat)
    write(*,*) "poisson_grd rcart ",poisson_grd%rcart(1:3,1:poisson_grd%nat)
    write(*,*) poisson_grd%bc
    write(*,*) poisson_grd%reset_rho
    write(*,*) poisson_grd%nat
    write(*,*) poisson_grd%rcart
    write(*,*) poisson_grd%q
    write(*,*) poisson_grd%gw_ewald
    write(*,*) poisson_grd%rgcut
    write(*,*) poisson_grd%ngpx,poisson_grd%ngpy,poisson_grd%ngpz
    write(*,*) poisson_grd%hgrid
    call put_charge_density(parini,poisson_grd)
    write(*,*) "put_charge_density is OK"
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
    write(*,*) "get_hartree is OK"
    write(*,*) "ehartree ", ehartree_grd,ehartree_grd-1.d0/(gausswidth(1)*sqrt(2.d0*pi))
    !write(*,*) "ngx, ngy, ngz ",size(poisson_grd%pot)
    pot_grd(1:ngx,1:ngy,1:ngz) = poisson_grd%pot(1:ngx,1:ngy,1:ngz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !pot_anl(1:ngx,1:ngy,1:ngz) = 0.d0
    !write (100,*) hx , hy , hz
    !write (100,*) ngx , ngy , ngz
    !write (100,*) (ngx/2.d0)*hx , (ngy/2.d0)*hy , (ngz/2.d0)*hz
    
    do iat=1,atoms%nat
        do iz=1,ngz
            dz=(iz-1-nbgz*1)*hz - atoms%rat(3,iat)
            do iy=1,ngy
                dy=(iy-1-nbgy*1)*hy - atoms%rat(2,iat)
                do ix=1,ngx
                    dx=(ix-1-nbgx*1)*hx - atoms%rat(1,iat)
                    r=sqrt(dx**2+dy**2+dz**2)
                    if (r <= 1.d-12) then 
                        pot_anl(ix,iy,iz)=2.d0/sqrt(pi)
                    else
                        pot_anl(ix,iy,iz)=atoms%qat(iat)*erf(r/gausswidth(iat))/r
                    endif
                    write(95,*) atoms%qat , r , erf(r/poisson_grd%gw_ewald(iat))
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
    write(*,*) "Rho_err_max is : ",rho_err_max
    write(*,*) "Rho_err_norm is : ",temp2
    write(77,*) pot_anl
    write(88,*) pot_grd
    end associate
    end associate
end subroutine test_free_bps
