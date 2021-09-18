!*****************************************************************************************
subroutine test_free_bps(parini)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, typ_atoms_arr
    use mod_atoms, only: atom_copy_old, atom_deallocate, update_ratp, get_rat
!    use dynamic_memory
    !use mod_acf, only: acf_read
    use mod_yaml_conf, only: read_yaml_conf
    implicit none
    type(typ_parini), intent(in):: parini
    ! local variables
    type(typ_poisson):: poisson_grd, poisson_scn
    type(typ_atoms):: atoms
    type(typ_atoms_arr):: atoms_arr
    integer:: iat, ix, iy, iz 
    integer:: nbgx, nbgy, nbgz
    real(8):: ehartree_grd, ehartree_scn, ehartree_ttt
    real(8):: pi
    real(8):: hx, hy, hz, dx, dy, dz , r 
    real(8):: rho_err_max, temp1,temp2, c1, d1, all_grids
    real(8),allocatable:: gausswidth(:) 
    real(8),allocatable:: pot_anl(:,:,:), pot_grd(:,:,:)
    pi=4.d0*atan(1.d0)
    associate(nat=>atoms%nat)
    associate(ngx=>poisson_grd%ngpx,ngy=>poisson_grd%ngpy,ngz=>poisson_grd%ngpz)
    associate(sf=>parini%screening_factor)
    !call acf_read(parini,'posinp.acf',1,atoms=atoms)
    call read_yaml_conf(parini,'posinp.yaml',1,atoms_arr)
    if(atoms_arr%nconf/=1) stop 'ERROR: atoms_arr%nconf/=1 in test_free_bps'
    call atom_copy_old(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(iconf)->atoms')
    call atom_deallocate(atoms_arr%atoms(1))
    deallocate(atoms_arr%atoms)
    allocate(gausswidth(1:atoms%nat)) 
    atoms%qat(1:atoms%nat) = 1.d0
    !-------------------------------------------------------
    poisson_grd%alpha=parini%alpha_ewald
    gausswidth(:) = poisson_grd%alpha
    poisson_grd%task_finit='set_ngp:alloc_rho'
    call init_hartree_bps(parini,atoms,poisson_grd)
    poisson_grd%reset_rho=.true.
    poisson_grd%nat=atoms%nat
    poisson_grd%cv=atoms%cellvec
    poisson_grd%bc=atoms%boundcond
    allocate(poisson_grd%q(1:poisson_grd%nat), poisson_grd%gw_ewald(1:poisson_grd%nat), poisson_grd%rcart(1:3,1:poisson_grd%nat))
    poisson_grd%q(1:poisson_grd%nat)=atoms%qat(1:atoms%nat)
    poisson_grd%gw_ewald(1:poisson_grd%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson_grd%rcart)
    call put_charge_density(parini,poisson_grd)
    !-------------------------------------------------------
    poisson_scn%alpha=parini%alpha_ewald
    gausswidth(:) = poisson_scn%alpha
    poisson_scn%task_finit='set_ngp:alloc_rho'
    poisson_scn%screening_factor=sf
    call init_hartree_bps(parini,atoms,poisson_scn)
    poisson_scn%reset_rho=.true.
    poisson_scn%nat=atoms%nat
    poisson_scn%cv=atoms%cellvec
    poisson_scn%bc=atoms%boundcond
    allocate(poisson_scn%q(1:poisson_scn%nat), poisson_scn%gw_ewald(1:poisson_scn%nat), poisson_scn%rcart(1:3,1:poisson_scn%nat))
    poisson_scn%q(1:poisson_scn%nat)=atoms%qat(1:atoms%nat)
    poisson_scn%gw_ewald(1:poisson_scn%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson_scn%rcart)
    call put_charge_density(parini,poisson_scn)
    !-------------------------------------------------------
    nbgx=int(poisson_grd%rgcut/poisson_grd%hgrid(1,1))+2
    nbgy=int(poisson_grd%rgcut/poisson_grd%hgrid(2,2))+2
    nbgz=int(poisson_grd%rgcut/poisson_grd%hgrid(3,3))+2
    hx=poisson_grd%hgrid(1,1)
    hy=poisson_grd%hgrid(2,2)
    hz=poisson_grd%hgrid(3,3)
    allocate(pot_grd(ngx,ngy,ngz),pot_anl(ngx,ngy,ngz))
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
    !-------------------------------------------------------
    call get_hartree(parini,poisson_scn,atoms,gausswidth,ehartree_scn)
    ehartree_ttt=0.d0
    do iz=1,ngz
    do iy=1,ngy
    do ix=1,ngx
        ehartree_ttt=ehartree_ttt+poisson_scn%rho(ix,iy,iz)*(poisson_grd%pot(ix,iy,iz)-poisson_scn%pot(ix,iy,iz))
    enddo
    enddo
    enddo
    ehartree_ttt=ehartree_ttt*(hx*hy*hz)*0.5d0
    write(*,*) 'ehartree_grd= ',ehartree_grd
    write(*,*) 'ehartree_scn= ',ehartree_scn
    write(*,*) 'ehartree_ttt= ',ehartree_ttt
    write(*,*) 'SF= ',sf

    write(*,'(a,2es19.10)') "Ehartree, Ehartree_differance : ", ehartree_grd,ehartree_grd-1.d0/(gausswidth(1)*sqrt(2.d0*pi))
    pot_grd(1:ngx,1:ngy,1:ngz) = poisson_scn%pot(1:ngx,1:ngy,1:ngz)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    call update_ratp(atoms)
    !call get_anal_potential_nonspe(atoms%nat,atoms%ratp,atoms%qat,gausswidth,ngx,ngy,ngz,nbgx,nbgy,nbgz,hx,hy,hz,pot_anl)
    call get_anal_potential_spe(atoms%nat,atoms%ratp,atoms%qat,gausswidth,ngx,ngy,ngz,nbgx,nbgy,nbgz,hx,hy,hz,sf,pot_anl)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    rho_err_max=0.d0
    all_grids = ngx*ngy*ngz*1.d0
    all_grids = 1.d0/all_grids
    temp2=0.d0
    do iz = 1 , ngz
        do iy = 1 , ngy
            do ix = 1 , ngx
                c1 = pot_anl(ix,iy,iz)
                d1 = pot_grd(ix,iy,iz)
                if(iz==ngz/2 .and. iy==ngy/2) write(21,'(3i5,3es19.10)') iz,iy,ix,c1,d1,poisson_scn%rho(ix,iy,iz)
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
    end associate
    deallocate(poisson_grd%q, poisson_grd%gw_ewald,poisson_grd%rcart,pot_anl,pot_grd,gausswidth)
end subroutine test_free_bps
!*****************************************************************************************
subroutine get_anal_potential_nonspe(nat,rat,qat,gausswidth,ngx,ngy,ngz,nbgx,nbgy,nbgz,hx,hy,hz,pot_anl)
    implicit none
    integer, intent(in):: nat, ngx, ngy, ngz , nbgx, nbgy, nbgz
    real(8), intent(in):: rat(3,nat), qat(nat), gausswidth(nat), hx, hy, hz
    real(8), intent(out):: pot_anl(ngx,ngy,ngz)
    ! local variables
    integer:: iat, ix, iy, iz
    real(8):: pi, r, dx, dy, dz
    pi=4.d0*atan(1.d0)
    pot_anl=0.d0
    do iat=1,nat
        do iz=1,ngz
            dz=(iz-1-nbgz*1)*hz - rat(3,iat)
            do iy=1,ngy
                dy=(iy-1-nbgy*1)*hy - rat(2,iat)
                do ix=1,ngx
                    dx=(ix-1-nbgx*1)*hx - rat(1,iat)
                    r=sqrt(dx**2+dy**2+dz**2)
                    if(r<=1.d-12) then 
                        pot_anl(ix,iy,iz)=pot_anl(ix,iy,iz)+qat(iat)*2.d0/sqrt(pi)/gausswidth(iat)
                    else
                        pot_anl(ix,iy,iz)=pot_anl(ix,iy,iz)+qat(iat)*erf(r/gausswidth(iat))/r
                    endif
                enddo
            enddo
        enddo
    enddo
end subroutine get_anal_potential_nonspe
!*****************************************************************************************
subroutine get_anal_potential_spe(nat,rat,qat,gausswidth,ngx,ngy,ngz,nbgx,nbgy,nbgz,hx,hy,hz,sf,pot_anl)
    implicit none
    integer, intent(in):: nat, ngx, ngy, ngz , nbgx, nbgy, nbgz
    real(8), intent(in):: rat(3,nat), qat(nat), gausswidth(nat), hx, hy, hz, sf
    real(8), intent(out):: pot_anl(ngx,ngy,ngz)
    ! local variables
    integer:: iat, ix, iy, iz
    real(8):: pi, r, dx, dy, dz, a, tt1, tt2
    pi=4.d0*atan(1.d0)
    write(*,*) 'QAT ',qat(1),qat(nat)
    write(*,*) 'RAT ',rat(1,1),rat(2,1),rat(3,1)
    write(*,*) 'RAT ',rat(1,2),rat(2,2),rat(3,2)
    pot_anl=0.d0
    write(*,*) shape(pot_anl)
    do iat=1,nat
        a=gausswidth(iat)
        do iz=1,ngz
            dz=(iz-1-nbgz*1)*hz - rat(3,iat)
            do iy=1,ngy
                dy=(iy-1-nbgy*1)*hy - rat(2,iat)
                do ix=1,ngx
                    dx=(ix-1-nbgx*1)*hx - rat(1,iat)
                    r=sqrt(dx**2+dy**2+dz**2)
                    if(r<=1.d-12) then 
                        tt1=2.d0/(a*sqrt(pi))-exp((a**2*sf**2)/4.d0)*sf*erfc((a*sf)/2.d0)
                        pot_anl(ix,iy,iz)=pot_anl(ix,iy,iz)+qat(iat)*tt1
                    else
                        tt1=exp((sf*(-4.d0*r+a**2*sf))/4.d0)
                        tt2=2.d0-erfc(r/a-(a*sf)/2.d0)-exp(2.d0*r*sf)*erfc(r/a+(a*sf)/2.d0)
                        pot_anl(ix,iy,iz)=pot_anl(ix,iy,iz)+qat(iat)*(tt1*tt2)/(2.d0*r)
                    endif
                enddo
            enddo
        enddo
    enddo
end subroutine get_anal_potential_spe
!*****************************************************************************************
