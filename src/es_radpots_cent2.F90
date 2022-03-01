!*****************************************************************************************
module mod_radpots_cent2
    implicit none
    private
    public:: typ_radpots_cent2
    public:: cal_pot_hartree

    type:: typ_radpots_cent2
        private
        integer, public:: ngp
        integer, public:: ntypat
        real(8), public:: hgp
        logical:: initialized=.false.
        real(8), allocatable, public:: linear_rho_e(:,:)
        real(8), allocatable, public:: linear_rho_n(:,:)
        real(8), allocatable, public:: linear_rho_t(:)
        real(8), allocatable, public:: linear_pot_e(:,:)
        real(8), allocatable, public:: linear_pot_n(:,:)
        real(8), allocatable, public:: gwe(:)
        real(8), allocatable, public:: gwn(:)
        contains
        procedure, public, pass(self):: init_radpots_cent2
        procedure, public, pass(self):: fini_radpots_cent2
        procedure, public, pass(self):: pot_1dto3d
        procedure, public, pass(self):: energy_1dto3d
        procedure, public, pass(self):: cal_radpots
    end type typ_radpots_cent2
contains
!*****************************************************************************************
subroutine init_radpots_cent2(self,hgp,rgcut,cellVec_max,ntypat)
    use dynamic_memory
    implicit none
    class(typ_radpots_cent2), intent(inout):: self
    real(8), intent(in):: hgp
    real(8), intent(in):: rgcut
    real(8), intent(in):: cellVec_max
    integer, intent(in):: ntypat
    !local variables
    self%hgp=hgp
    self%ntypat=ntypat
    self%ngp=ceiling(2.d0*(rgcut+cellVec_max)/hgp)
    self%linear_rho_e=f_malloc0([0.to.self%ngp,1.to.ntypat],id='linear_rho_e')
    self%linear_pot_e=f_malloc0([0.to.self%ngp,1.to.ntypat],id='linear_pot_e')
    self%linear_rho_n=f_malloc0([0.to.self%ngp,1.to.ntypat],id='linear_rho_n')
    self%linear_pot_n=f_malloc0([0.to.self%ngp,1.to.ntypat],id='linear_pot_n')
    self%linear_rho_t=f_malloc0([0.to.self%ngp],id='linear_rho_t')
    self%gwe=f_malloc0([1.to.ntypat],id='gwe')
    self%gwn=f_malloc0([1.to.ntypat],id='gwn')
    self%initialized=.true.
end subroutine init_radpots_cent2
!*****************************************************************************************
subroutine fini_radpots_cent2(self)
    use dynamic_memory
    implicit none
    class(typ_radpots_cent2), intent(inout):: self
    !local variables
    if(.not. self%initialized) return
    !if(.not. self%initialized) then
    !    stop 'ERROR: attempting finalizing radpots_cent2 while not initialized.'
    !endif
    call f_free(self%linear_rho_e)
    call f_free(self%linear_pot_e)
    call f_free(self%linear_rho_n)
    call f_free(self%linear_pot_n)
    call f_free(self%linear_rho_t)
    call f_free(self%gwe)
    call f_free(self%gwn)
    self%initialized=.false.
end subroutine fini_radpots_cent2
!*****************************************************************************************
subroutine cal_radpots(self,sf)
    implicit none
    class(typ_radpots_cent2), intent(inout):: self
    real(8), intent(in):: sf
    !local variables
    integer:: itype, ii
    real(8):: pi, tt, pref
    do itype=1,self%ntypat
        call get_scf_pot_cent2_twogauss(self%ngp,self%hgp,self%gwe(itype),sf, &
            self%linear_rho_e(0,itype),self%linear_pot_e(0,itype))
        call get_scf_pot_cent2_onegauss(self%ngp,self%hgp,self%gwn(itype),sf, &
            self%linear_rho_n(0,itype),self%linear_pot_n(0,itype))
    enddo
    pi=4.d0*atan(1.d0)
    pref=1.d0/((pi**1.5d0)*(1.d0**3))
    do ii=0,self%ngp
        tt=real(ii,kind=8)*self%hgp
        if(tt>(10.d0*1.d0)) then
            self%linear_rho_t(ii)=0.d0
        else
            self%linear_rho_t(ii)=pref*exp(-tt**2/1.d0**2)
        endif
    enddo
end subroutine cal_radpots
!*****************************************************************************************
subroutine pot_1dto3d(self,select_pot,xyz,itypat,poisson,reset,q,pot)
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_radpots_cent2), intent(in):: self
    character(*), intent(in):: select_pot
    real(8), intent(in):: xyz(3)
    integer, intent(in):: itypat
    type(typ_poisson), intent(in):: poisson
    logical, intent(in):: reset
    real(8), intent(in):: q
    real(8), intent(inout):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    select case(trim(select_pot))
        case('linear_pot_e')
            call radial_to_3d(self%ngp,self%hgp,self%linear_pot_e(0,itypat),xyz, &
                poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,reset,q,pot)
        case('linear_pot_n')
            call radial_to_3d(self%ngp,self%hgp,self%linear_pot_n(0,itypat),xyz, &
                poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,reset,q,pot)
        case default
            stop 'ERROR: unknown select_pot in pot_1dto3d'
    end select
end subroutine pot_1dto3d
!*****************************************************************************************
subroutine energy_1dto3d(self,select_pot,xyz,itypat,poisson,q,pot,ener)
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_radpots_cent2), intent(in):: self
    character(*), intent(in):: select_pot
    real(8), intent(in):: xyz(3)
    integer, intent(in):: itypat
    type(typ_poisson), intent(in):: poisson
    real(8), intent(in):: q
    real(8), intent(in):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    real(8), intent(out):: ener
    !local variables
    select case(trim(select_pot))
        case('linear_rho_e')
            call radial_to_3d_energy(self%ngp,self%hgp,self%linear_rho_e(0,itypat),xyz, &
                poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,poisson%rgcut,q,pot,ener)
        case('linear_rho_n')
            call radial_to_3d_energy(self%ngp,self%hgp,self%linear_rho_n(0,itypat),xyz, &
                poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,poisson%rgcut,q,pot,ener)
        case('linear_rho_t')
            call radial_to_3d_energy(self%ngp,self%hgp,self%linear_rho_t(0),xyz, &
                poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,poisson%rgcut,q,pot,ener)
        case default
            stop 'ERROR: unknown select_pot in energy_1dto3d'
    end select
end subroutine energy_1dto3d
!*****************************************************************************************
subroutine radial_to_3d(nrad,hrad,funrad,xyz,xyz111,ngpx,ngpy,ngpz,hgrid,reset,q,fun3d)
    implicit none
    integer, intent(in):: nrad, ngpx, ngpy, ngpz
    real(8), intent(in):: hrad, funrad(0:nrad), hgrid(3,3), xyz(3), xyz111(3), q
    logical, intent(in):: reset
    real(8), intent(out):: fun3d(ngpx,ngpy,ngpz)
    !local variables
    integer:: igpx, igpy, igpz, ii
    real(8):: dx, dy, dz, r, tt0, tt1, r0
    if(reset) fun3d(1:ngpx,1:ngpy,1:ngpz)=0.d0
    do igpz=1,ngpz
        do igpy=1,ngpy
            do igpx=1,ngpx
                dx=xyz111(1)+(igpx-1)*hgrid(1,1)-xyz(1)
                dy=xyz111(2)+(igpy-1)*hgrid(2,2)-xyz(2)
                dz=xyz111(3)+(igpz-1)*hgrid(3,3)-xyz(3)
                r=sqrt(dx**2+dy**2+dz**2)
                ii=floor(r/hrad)
                r0=real(ii,kind=8)*hrad
                tt0=funrad(ii+0)
                tt1=funrad(ii+1)
                fun3d(igpx,igpy,igpz)=fun3d(igpx,igpy,igpz)+q*((r-r0)*(tt1-tt0)/hrad+tt0)
            enddo
        enddo
    enddo
end subroutine radial_to_3d
!*****************************************************************************************
subroutine radial_to_3d_energy(nrad,hrad,funrad,xyz,xyz111,ngpx,ngpy,ngpz,hgrid,rgcut,q,pot,ener)
    implicit none
    integer, intent(in):: nrad, ngpx, ngpy, ngpz
    real(8), intent(in):: hrad, funrad(0:nrad), hgrid(3,3), xyz(3), xyz111(3), rgcut, q
    real(8), intent(in):: pot(ngpx,ngpy,ngpz)
    real(8), intent(out):: ener
    !local variables
    integer:: igpx, igpy, igpz, ii
    integer:: nbgx, nbgy, nbgz, agpx, agpy, agpz
    real(8):: dx, dy, dz, r, tt0, tt1, r0, res, dd
    nbgx=int(rgcut/hgrid(1,1))+3
    nbgy=int(rgcut/hgrid(2,2))+3
    nbgz=int(rgcut/hgrid(3,3))+3
    agpx=int(xyz(1)/hgrid(1,1))+nbgx
    agpy=int(xyz(2)/hgrid(2,2))+nbgy
    agpz=int(xyz(3)/hgrid(3,3))+nbgz
    res=0.d0
    do igpz=agpz-nbgz,agpz+nbgz
        do igpy=agpy-nbgy,agpy+nbgy
            do igpx=agpx-nbgx,agpx+nbgx
                dx=xyz111(1)+(igpx-1)*hgrid(1,1)-xyz(1)
                dy=xyz111(2)+(igpy-1)*hgrid(2,2)-xyz(2)
                dz=xyz111(3)+(igpz-1)*hgrid(3,3)-xyz(3)
                r=sqrt(dx**2+dy**2+dz**2)
                ii=floor(r/hrad)
                r0=real(ii,kind=8)*hrad
                tt0=funrad(ii+0)
                tt1=funrad(ii+1)
                dd=q*((r-r0)*(tt1-tt0)/hrad+tt0)
                res=res+dd*pot(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    ener=res*(hgrid(1,1)*hgrid(2,2)*hgrid(3,3))
end subroutine radial_to_3d_energy
!*****************************************************************************************
subroutine get_scf_pot_cent2_onegauss(ngp,hgp,gw,scf,rho,pot)
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: hgp, gw, scf
    real(8), intent(out):: rho(0:ngp)
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp
    real(8):: den_coeff
    real(8):: rrad(0:ngp), pot_scn(0:ngp),weight (0:ngp)
    real(8):: pi
    pi=4.d0*atan(1.d0)
    den_coeff=1.d0/((pi**1.5d0)*(gw**3))
    rho=0.d0
    pot=0.d0
    do igp=0,ngp
        rrad(igp)= igp*hgp
        if(rrad(igp)>(10.d0*gw)) then
            rho(igp)=0.d0
        else
            rho(igp)=den_coeff*exp(-1.d0*(rrad(igp)**2/gw**2)) !Charge density with Q=1
        endif
    enddo
    if(scf>0.d0) then
        call set_weight(ngp,rrad,weight)
        call cal_powern_screened_poisson_gaussian(ngp,rrad,weight,1.d0,gw,scf,4,pot_scn)
    endif
    call cal_pot_hartree(ngp,rrad,rho,pot)
    do igp=0,ngp
        pot(igp)=pot(igp)-pot_scn(igp)
    enddo
end subroutine get_scf_pot_cent2_onegauss
!*****************************************************************************************
subroutine get_scf_pot_cent2_twogauss(ngp,hgp,gw,scf,rho,pot)
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: hgp, gw, scf
    real(8), intent(out):: rho(0:ngp)
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp
    real(8):: den_coeff
    real(8):: rrad(0:ngp), weight(0:ngp)
    real(8):: pi, gw_t, alpha, beta !, tt
    real(8), allocatable:: pot_scn(:), pot_scn_t(:)
    pi=4.d0*atan(1.d0)
    allocate(pot_scn(0:ngp),pot_scn_t(0:ngp))
    den_coeff=1.d0/((pi**1.5d0)*(gw**3))
    rho=0.d0
    pot=0.d0
    do igp=0,ngp
        rrad(igp)= igp*hgp
        if(rrad(igp)>(10.d0*gw)) then
            rho(igp)=0.d0
        else
            rho(igp)=den_coeff*exp(-1.d0*(rrad(igp)**2/gw**2)) !Charge density with Q=1
        endif
    enddo
    gw_t=0.95d0*gw
    alpha=gw**3/(gw**3-gw_t**3)
    beta=-gw_t**3/(gw**3-gw_t**3)
    do igp=0,ngp
        rho(igp)=alpha*rho(igp)
    enddo
    den_coeff=1.d0/((pi**1.5d0)*(gw_t**3))
    do igp=0,ngp
        !rrad(igp)= igp*hgp
        if(rrad(igp)>(10.d0*gw_t)) then
            !rho(igp)=0.d0
        else
            rho(igp)=rho(igp)+beta*den_coeff*exp(-1.d0*(rrad(igp)**2/gw_t**2)) ! Charge density with Q=1
        endif
    enddo
    if(scf>0.d0) then
        call set_weight(ngp,rrad,weight)
        call cal_powern_screened_poisson_gaussian(ngp,rrad,weight,1.d0,gw,scf,4,pot_scn)
        call cal_powern_screened_poisson_gaussian(ngp,rrad,weight,1.d0,gw_t,scf,4,pot_scn_t)
    endif
    !tt=0.d0
    !do igp=0,ngp
    !    tt=tt+rrad(igp)**2*rho(igp)*weight(igp)
    !enddo
    !tt=tt*(4.d0*pi)
    !write(*,*) 'RHO ',tt
    call cal_pot_hartree(ngp,rrad,rho,pot)
    do igp=0,ngp
        pot(igp)=pot(igp)-(alpha*pot_scn(igp)+beta*pot_scn_t(igp))
    enddo
end subroutine get_scf_pot_cent2_twogauss
!*****************************************************************************************
!This routine calculates:
!4*pi/r*\int_0^r rho(r') r'^2 dr' + 4*pi*\int_r^inf rho(r') r' dr'
subroutine cal_pot_hartree(ngp,rrad,den,pot_hartree)
    !Calculates Hartree potential.
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: rrad(0:ngp), den(0:ngp)
    real(8), intent(out):: pot_hartree(0:ngp)
    !local variables
    integer:: igp
    real(8):: pi, tt
    real(8), allocatable:: wa1(:), wa2(:)
    pi=4.d0*atan(1.d0)
    allocate(wa1(0:ngp),wa2(0:ngp))
    !calculating the first integral: \int_0^r rho(r') r'^2 dr'
    wa1(0)=.0d0
    do igp=1,ngp
        !wa1(igp)=wa1(igp-1)+den(igp)*rrad(igp)*rrad(igp)*weight(igp)
        tt=den(igp-1)*rrad(igp-1)**2+den(igp)*rrad(igp)**2
        wa1(igp)=wa1(igp-1)+tt*(rrad(igp)-rrad(igp-1))*0.5d0
    enddo
    !calculating the second integral: \int_r^inf rho(r') r' dr'
    wa2(ngp)=.0d0
    do igp=ngp-1,0,-1
        tt=den(igp)*rrad(igp)+den(igp+1)*rrad(igp+1)
        wa2(igp)=wa2(igp+1)+tt*(rrad(igp+1)-rrad(igp))*0.5d0
    enddo
    !finally summing the two integrals
    !note that Hartree potential at origin is due zero not infinity,
    !this can be obtained after resolving 0/0 ambiguity
    pot_hartree(0)=(4.d0*pi)*wa2(0)
    do igp=1,ngp
        pot_hartree(igp)=4.d0*pi*(wa1(igp)/rrad(igp)+wa2(igp))
    enddo
    deallocate(wa1,wa2)
end subroutine cal_pot_hartree
!*****************************************************************************************
subroutine cal_powern_screened_poisson_gaussian(ngp,rrad,weight,q,gw,sf,npow,pot)
    !Calculates Hartree potential.
    implicit none
    integer, intent(in):: ngp, npow
    real(8), intent(in):: rrad(0:ngp), weight(0:ngp), q, gw, sf
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp, jgp, mgp
    real(8):: pi, tt0, tt1, tt2, ss1, ss2, res, slimit
    real(8), allocatable:: w1(:)
    pi=4.d0*atan(1.d0)
    slimit=-log(1.d-10) ! 45.d0 !exp(-slimit)~2.9E-20
    !write(*,*) 'slimit= ',slimit
    allocate(w1(0:ngp))
    w1(0:ngp)=0.d0
    do igp=0,ngp
        ss1=(sf*rrad(igp))**npow
        if(ss1>slimit) then
            mgp=igp
            exit
        endif
        w1(igp)=exp(-ss1)
    enddo
    !write(*,*) 'mgp,ngp= ',mgp,ngp
    tt0=q/(sqrt(pi)*gw)
    do igp=1,ngp
        res=0.d0
        do jgp=0,mgp
            ss1=(rrad(jgp)-rrad(igp))**2/gw**2
            if(ss1>slimit) then
                tt1=0.d0
            else
                tt1=exp(-ss1)
            endif
            ss2=4.d0*rrad(jgp)*rrad(igp)/gw**2
            if(ss2>slimit) then
                tt2=1.d0
            else
                tt2=(1.d0-exp(-ss2))
            endif
            res=res+w1(jgp)*tt1*tt2*weight(jgp)
        enddo
        pot(igp)=res*tt0/rrad(igp)
    enddo
    res=0.d0
    do jgp=0,mgp
        ss1=rrad(jgp)**2/gw**2
        if(ss1>slimit) then
            tt1=0.d0
        else
            tt1=exp(-ss1)
        endif
        res=res+w1(jgp)*tt1*rrad(jgp)*weight(jgp)
    enddo
    pot(0)=res*4.d0*q/(sqrt(pi)*gw**3)
    deallocate(w1)
end subroutine cal_powern_screened_poisson_gaussian
!*****************************************************************************************
subroutine set_weight(ngp,rrad,weight)
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: rrad(0:ngp)
    real(8), intent(out):: weight(0:ngp)
    !local variables
    integer:: igp
    weight(0)=(rrad(1)-rrad(0))*0.5d0
    do igp=1,ngp-1
        weight(igp)=(rrad(igp+1)-rrad(igp-1))*0.5d0
    enddo
    weight(ngp)=(rrad(ngp)-rrad(ngp-1))*0.5d0
end subroutine set_weight
!*****************************************************************************************
end module mod_radpots_cent2
!*****************************************************************************************
