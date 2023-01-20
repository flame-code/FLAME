!*****************************************************************************************
module mod_fit_bf_cent2
    use mod_trial_energy, only: typ_trial_energy
    implicit none
    private
    public:: get_basis_functions_cent2
    public:: cal_pot_gauss_s, cal_pot_gauss_p, cal_pot_r2gauss_s
    real(8):: rgcut=5.d0
    real(8), allocatable:: c_s_t(:)
    real(8), allocatable:: c_p_t(:)
    !real(8):: time_cal_pot_gauss_p=0.d0
    !real(8):: time_cal_pot_gauss_s=0.d0
    !real(8):: time_cal_pot_r2gauss_s=0.d0
    type(typ_trial_energy), allocatable:: trial_energy_all(:)
    type typ_fitpar
        private
        integer:: ntypat
        real(8):: gausswidth_ion=0.5d0
        logical, allocatable:: relaxcore(:)
        logical, allocatable:: applycore(:)
        real(8), allocatable:: gwz_s(:)
        real(8), allocatable:: gwc_s1(:)
        !real(8), allocatable:: gwc_s2(:)
        real(8), allocatable:: gwv_s1(:)
        !real(8), allocatable:: gwv_s2(:)
        real(8), allocatable:: gwv_p1(:)
        real(8), allocatable:: gwv_p2(:)
        real(8), allocatable:: grad_gwc_s1(:)
        !real(8), allocatable:: grad_gwc_s2(:)
        real(8), allocatable:: grad_gwv_s1(:)
        !real(8), allocatable:: grad_gwv_s2(:)
        real(8), allocatable:: grad_gwv_p1(:)
        real(8), allocatable:: grad_gwv_p2(:)
        real(8), allocatable:: qcore_type(:)
        real(8), allocatable:: bz_s(:)
        real(8), allocatable:: bc_s1(:)
        !real(8), allocatable:: bc_s2(:)
        real(8), allocatable:: bv_s1(:)
        !real(8), allocatable:: bv_s2(:)
        real(8), allocatable:: bv_p1(:)
        real(8), allocatable:: bv_p2(:)
        real(8), allocatable:: grad_bc_s1(:)
        !real(8), allocatable:: grad_bc_s2(:,:)
        real(8), allocatable:: grad_bv_s1(:)
        !real(8), allocatable:: grad_bv_s2(:,:)
        real(8), allocatable:: grad_bv_p1(:,:)
        real(8), allocatable:: grad_bv_p2(:,:)
        contains
        procedure, private, pass(self):: init_fit_bf
        procedure, private, pass(self):: fini_fit_bf
        !procedure, private, pass(self):: set_param
        procedure, private, pass(self):: set_bc_bcg
        procedure, private, pass(self):: report_fit_bf
    end type typ_fitpar
    type typ_bf
        private
        integer:: nbf
        real(8):: zeroder
        !real(8), allocatable:: 
        integer, allocatable:: iat_list(:)
        integer, allocatable:: ibf_list_s(:)
        integer, allocatable:: ibf_list_px(:)
        real(8), allocatable:: wa1(:,:,:)
        real(8), allocatable:: wa2(:,:,:)
        real(8), allocatable:: wa3(:,:,:)
        real(8), allocatable:: wa4(:,:,:)
        real(8), allocatable:: vgrad_1(:,:,:)
        real(8), allocatable:: vgrad_2(:,:,:)
        real(8), allocatable:: pot_1(:,:,:)
        real(8), allocatable:: pot_2(:,:,:)
        real(8), allocatable:: subsecder_grad(:,:,:,:)
        real(8), allocatable:: subfirstder_grad(:,:)
        real(8), allocatable:: subsecder(:,:,:,:)
        real(8), allocatable:: subfirstder(:,:)
        real(8), allocatable:: secder(:,:)
        real(8), allocatable:: firstder(:)
        real(8), allocatable:: csp(:)
        real(8), allocatable:: qcore(:)
        real(8), allocatable:: bv_1(:)
        real(8), allocatable:: bv_2(:)
        real(8), allocatable:: grad_bv_1(:,:)
        real(8), allocatable:: grad_bv_2(:,:)
        real(8), allocatable:: weight(:,:,:)
        character(2), allocatable:: bt(:)
        contains
        procedure, private, pass(self):: init_bf
        procedure, private, pass(self):: fini_bf
        procedure, private, pass(self):: set_bt
        procedure, private, pass(self):: set_bf_bc_bcg
        procedure, private, pass(self):: get_pot_single
        procedure, private, pass(self):: get_pot_single_core
        procedure, private, pass(self):: get_rho_single
        procedure, private, pass(self):: get_rho_single_core
        !procedure, private, pass(self):: get_linearcoeff
    end type typ_bf
contains
!*****************************************************************************************
subroutine init_bf(self,nat,rat,poisson)
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_bf), intent(inout):: self
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat)
    type(typ_poisson), intent(in):: poisson
    !local variables
    real(8):: x, y, z, rsq
    integer:: ix, iy, iz
    self%nbf=nat*2
    allocate(self%wa1(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%wa2(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%wa3(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%wa4(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%vgrad_1(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%vgrad_2(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%pot_1(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%pot_2(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%subsecder_grad(2,2,self%nbf,self%nbf))
    allocate(self%subfirstder_grad(2,self%nbf))
    allocate(self%subsecder(2,2,self%nbf,self%nbf))
    allocate(self%subfirstder(2,self%nbf))
    allocate(self%secder(self%nbf,self%nbf))
    allocate(self%firstder(self%nbf))
    allocate(self%csp(self%nbf))
    allocate(self%qcore(nat))
    allocate(self%bv_1(self%nbf))
    allocate(self%bv_2(self%nbf))
    allocate(self%grad_bv_1(2,self%nbf))
    allocate(self%grad_bv_2(2,self%nbf))
    allocate(self%bt(self%nbf))
    allocate(self%iat_list(self%nbf))
    allocate(self%ibf_list_s(nat))
    allocate(self%ibf_list_px(nat))
    allocate(self%weight(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    if(nat==1) then
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-rat(1,1)
            y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-rat(2,1)
            z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-rat(3,1)
            rsq=x**2+y**2+z**2
            self%weight(ix,iy,iz)=1.d0-exp(-(1.5d0**2*rsq)**2)
        enddo
        enddo
        enddo
    else
        self%weight=1.d0
    endif
    call self%set_bt()
end subroutine init_bf
!*****************************************************************************************
subroutine fini_bf(self)
    implicit none
    class(typ_bf), intent(inout):: self
    !integer, intent(in):: nat
    !local variables
    deallocate(self%wa1)
    deallocate(self%wa2)
    deallocate(self%wa3)
    deallocate(self%wa4)
    deallocate(self%vgrad_1)
    deallocate(self%vgrad_2)
    deallocate(self%pot_1)
    deallocate(self%pot_2)
    deallocate(self%subsecder_grad)
    deallocate(self%subfirstder_grad)
    deallocate(self%subsecder)
    deallocate(self%subfirstder)
    deallocate(self%secder)
    deallocate(self%firstder)
    deallocate(self%csp)
    deallocate(self%qcore)
    deallocate(self%bv_1)
    deallocate(self%bv_2)
    deallocate(self%grad_bv_1)
    deallocate(self%grad_bv_2)
    deallocate(self%bt)
    deallocate(self%iat_list)
    deallocate(self%ibf_list_s)
    deallocate(self%ibf_list_px)
    deallocate(self%weight)
end subroutine fini_bf
!*****************************************************************************************
subroutine set_bt(self)
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_bf), intent(inout):: self
    !local variables
    integer:: ibf, iat
    do ibf=1,self%nbf
        iat=modulo(ibf-1,self%nbf/2)+1
        self%iat_list(ibf)=iat
        if((2*(ibf-1))/self%nbf==0) then
            self%bt(ibf)='s'
            self%ibf_list_s(iat)=ibf
        elseif((2*(ibf-1))/self%nbf==1) then
            self%bt(ibf)='px'
            self%ibf_list_px(iat)=ibf
        else
            stop 'ERROR: unknon ibf in set_bt'
        endif
    enddo
end subroutine set_bt
!*****************************************************************************************
subroutine set_bf_bc_bcg(self,atoms,fitpar)
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_bf), intent(inout):: self
    type(typ_atoms), intent(in):: atoms
    type(typ_fitpar), intent(in):: fitpar
    !local variables
    integer:: ibf, iat, itypat
    do ibf=1,self%nbf
        iat=self%iat_list(ibf)
        itypat=atoms%itypat(iat)
        if(trim(self%bt(ibf))=='s') then
            self%bv_1(ibf)=fitpar%bv_s1(itypat)
            self%bv_2(ibf)=(1.d0-fitpar%bv_s1(itypat))
            self%grad_bv_1(1,ibf)=0.d0
            self%grad_bv_1(2,ibf)=0.d0
            self%grad_bv_2(1,ibf)=0.d0
            self%grad_bv_2(2,ibf)=0.d0
        elseif(trim(self%bt(ibf))=='px') then
            self%bv_1(ibf)=fitpar%bv_p1(itypat)
            self%bv_2(ibf)=fitpar%bv_p2(itypat)
            self%grad_bv_1(1,ibf)=fitpar%grad_bv_p1(1,itypat)
            self%grad_bv_1(2,ibf)=fitpar%grad_bv_p1(2,itypat)
            self%grad_bv_2(1,ibf)=fitpar%grad_bv_p2(1,itypat)
            self%grad_bv_2(2,ibf)=fitpar%grad_bv_p2(2,itypat)
        else
            stop 'ERROR: unknown bt in set_bf_bc_bcg'
        endif
    enddo
end subroutine set_bf_bc_bcg
!*****************************************************************************************
subroutine get_pot_single(self,parini,poisson,atoms,fitpar,ibf)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_bf), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    type(typ_atoms), intent(in):: atoms
    type(typ_fitpar), intent(in):: fitpar
    integer, intent(in):: ibf
    !local variables
    integer:: iat, itypat
    real(8):: q_tmp(1), p_tmp(3)
    iat=self%iat_list(ibf)
    itypat=atoms%itypat(iat)
    if(trim(self%bt(ibf))=='s') then
        q_tmp(1)=1.d0
        call cal_pot_gauss_s(parini,poisson,self%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwv_s1(itypat),q_tmp(1),self%vgrad_1)
        call cal_pot_r2gauss_s(parini,poisson,self%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwv_s1(itypat),q_tmp(1),self%vgrad_2)
    elseif(trim(self%bt(ibf))=='px') then
        p_tmp=0.d0
        p_tmp(1)=1.d0
        call cal_pot_gauss_p(parini,poisson,self%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwv_p1(itypat),p_tmp(1),self%vgrad_1)
        call cal_pot_gauss_p(parini,poisson,self%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwv_p2(itypat),p_tmp(1),self%vgrad_2)
    else
            stop 'ERROR: unknonbt in get_pot_single'
    endif
end subroutine get_pot_single
!*****************************************************************************************
subroutine get_pot_single_core(self,parini,poisson,atoms,fitpar,iat)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_bf), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    type(typ_atoms), intent(in):: atoms
    type(typ_fitpar), intent(in):: fitpar
    integer, intent(in):: iat
    !local variables
    integer:: itypat
    real(8):: q_tmp(1)
    itypat=atoms%itypat(iat)
    if(fitpar%applycore(itypat)) then
    q_tmp(1)=1.d0
    call cal_pot_gauss_s(parini,poisson,self%wa1,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),self%wa3)
    q_tmp(1)=1.d0
    call cal_pot_r2gauss_s(parini,poisson,self%wa2,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),self%wa4)
    endif
end subroutine get_pot_single_core
!*****************************************************************************************
subroutine get_rho_single(self,parini,poisson,atoms,fitpar,ibf)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_bf), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    type(typ_atoms), intent(in):: atoms
    type(typ_fitpar), intent(in):: fitpar
    integer, intent(in):: ibf
    !local variables
    integer:: iat, itypat
    real(8):: q_tmp(1), p_tmp(3)
    iat=self%iat_list(ibf)
    itypat=atoms%itypat(iat)
    if(trim(self%bt(ibf))=='s') then
        q_tmp(1)=1.d0
        !call cal_pot_gauss_s(parini,poisson,self%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwv_s1(itypat),q_tmp(1),self%vgrad_1)
    call put_gto_sym_ortho_wgrad(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp(1),fitpar%gwv_s1(itypat), &
        rgcut*fitpar%gwv_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,self%pot_1,self%vgrad_1)
        !call cal_pot_r2gauss_s(parini,poisson,self%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwv_s1(itypat),q_tmp(1),self%vgrad_2)
    call put_r2gto_sym_ortho_wgrad(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp(1),fitpar%gwv_s1(itypat), &
        rgcut*fitpar%gwv_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,self%pot_2,self%vgrad_2)
    elseif(trim(self%bt(ibf))=='px') then
        p_tmp=0.d0
        p_tmp(1)=1.d0
        !call cal_pot_gauss_p(parini,poisson,self%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwv_p1(itypat),p_tmp(1),self%vgrad_1)
    call put_gto_p_ortho_wgrad(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),p_tmp(1),fitpar%gwv_p1(itypat), &
        rgcut*fitpar%gwv_p1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,self%pot_1,self%vgrad_1)
        !call cal_pot_gauss_p(parini,poisson,self%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwv_p2(itypat),p_tmp(1),self%vgrad_2)
    call put_gto_p_ortho_wgrad(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),p_tmp(1),fitpar%gwv_p2(itypat), &
        rgcut*fitpar%gwv_p2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,self%pot_2,self%vgrad_2)
    else
            stop 'ERROR: unknonbt in get_pot_single'
    endif
end subroutine get_rho_single
!*****************************************************************************************
subroutine get_rho_single_core(self,parini,poisson,atoms,fitpar,iat)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_bf), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    type(typ_atoms), intent(in):: atoms
    type(typ_fitpar), intent(in):: fitpar
    integer, intent(in):: iat
    !local variables
    integer:: itypat
    real(8):: q_tmp(1)
    itypat=atoms%itypat(iat)
    if(fitpar%applycore(itypat)) then
    q_tmp(1)=1.d0
    !call cal_pot_gauss_s(parini,poisson,self%wa1,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),self%wa3)
    call put_gto_sym_ortho_wgrad(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp(1),fitpar%gwc_s1(itypat), &
        rgcut*fitpar%gwc_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,self%wa1,self%wa3)
    q_tmp(1)=1.d0
    !call cal_pot_r2gauss_s(parini,poisson,self%wa2,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),self%wa4)
    call put_r2gto_sym_ortho_wgrad(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp(1),fitpar%gwc_s1(itypat), &
        rgcut*fitpar%gwc_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,self%wa2,self%wa4)
    endif
end subroutine get_rho_single_core
!*****************************************************************************************
subroutine init_fit_bf(self,ntypat)
    implicit none
    class(typ_fitpar), intent(inout):: self
    integer, intent(in):: ntypat
    !local variables
    self%ntypat=ntypat
    allocate(self%relaxcore(ntypat))
    allocate(self%applycore(ntypat))
    allocate(self%gwv_s1(ntypat))
    !allocate(self%gwv_s2(ntypat))
    allocate(self%gwz_s(ntypat))
    allocate(self%gwc_s1(ntypat))
    !allocate(self%gwc_s2(ntypat))
    allocate(self%gwv_p1(ntypat))
    allocate(self%gwv_p2(ntypat))
    allocate(self%grad_gwc_s1(ntypat))
    !allocate(self%grad_gwc_s2(ntypat))
    allocate(self%grad_gwv_s1(ntypat))
    !allocate(self%grad_gwv_s2(ntypat))
    allocate(self%grad_gwv_p1(ntypat))
    allocate(self%grad_gwv_p2(ntypat))
    allocate(self%qcore_type(ntypat))
    allocate(self%bz_s(ntypat))
    allocate(self%bc_s1(ntypat))
    !allocate(self%bc_s2(ntypat))
    allocate(self%bv_s1(ntypat))
    !allocate(self%bv_s2(ntypat))
    allocate(self%bv_p1(ntypat))
    allocate(self%bv_p2(ntypat))
    allocate(self%grad_bc_s1(ntypat))
    !allocate(self%grad_bc_s2(2,ntypat))
    allocate(self%grad_bv_s1(ntypat))
    !allocate(self%grad_bv_s2(2,ntypat))
    allocate(self%grad_bv_p1(2,ntypat))
    allocate(self%grad_bv_p2(2,ntypat))
end subroutine init_fit_bf
!*****************************************************************************************
subroutine fini_fit_bf(self)
    implicit none
    class(typ_fitpar), intent(inout):: self
    !local variables
    deallocate(self%relaxcore)
    deallocate(self%applycore)
    deallocate(self%gwv_s1)
    !deallocate(self%gwv_s2)
    deallocate(self%gwz_s)
    deallocate(self%gwc_s1)
    !deallocate(self%gwc_s2)
    deallocate(self%gwv_p1)
    deallocate(self%gwv_p2)
    deallocate(self%grad_gwc_s1)
    !deallocate(self%grad_gwc_s2)
    deallocate(self%grad_gwv_s1)
    !deallocate(self%grad_gwv_s2)
    deallocate(self%grad_gwv_p1)
    deallocate(self%grad_gwv_p2)
    deallocate(self%qcore_type)
    deallocate(self%bz_s)
    deallocate(self%bc_s1)
    !deallocate(self%bc_s2)
    deallocate(self%bv_s1)
    !deallocate(self%bv_s2)
    deallocate(self%bv_p1)
    deallocate(self%bv_p2)
    deallocate(self%grad_bc_s1)
    !deallocate(self%grad_bc_s2)
    deallocate(self%grad_bv_s1)
    !deallocate(self%grad_bv_s2)
    deallocate(self%grad_bv_p1)
    deallocate(self%grad_bv_p2)
end subroutine fini_fit_bf
!*****************************************************************************************
subroutine set_bc_bcg(self)
    implicit none
    class(typ_fitpar), intent(inout):: self
    !local variables
    integer:: itypat
    do itypat=1,self%ntypat
        self%bv_p1(itypat)=get_coeff_p(self%gwv_p1(itypat),self%gwv_p2(itypat))
        self%bv_p2(itypat)=get_coeff_p(self%gwv_p2(itypat),self%gwv_p1(itypat))
        self%grad_bv_p1(1,itypat)=-get_coeff_p_grad(self%gwv_p1(itypat),self%gwv_p2(itypat))
        self%grad_bv_p1(2,itypat)= get_coeff_p_grad(self%gwv_p2(itypat),self%gwv_p1(itypat))
        self%grad_bv_p2(1,itypat)= get_coeff_p_grad(self%gwv_p1(itypat),self%gwv_p2(itypat))
        self%grad_bv_p2(2,itypat)=-get_coeff_p_grad(self%gwv_p2(itypat),self%gwv_p1(itypat))
    enddo
end subroutine set_bc_bcg
!*****************************************************************************************
subroutine report_fit_bf(self,istep,cost,cost_gw)
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_fitpar), intent(inout):: self
    integer, intent(in):: istep
    !type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: cost, cost_gw
    !local variables
    integer:: itypat, iat
    real(8):: tt1, tt2, tt3, tt4, dpm(3)
    write(*,'(a,i6,es14.5)',advance='no') 'core ',istep,cost
    do itypat=1,self%ntypat
        if(self%applycore(itypat)) then
            tt1=self%grad_gwc_s1(itypat)
            tt2=self%grad_bc_s1(itypat)
            write(*,'(2es10.1)',advance='no') tt1,tt2
        endif
    enddo
    write(*,*)
    write(*,'(a,i6,es14.5)',advance='no') 'cost ',istep,cost
    do itypat=1,self%ntypat
        tt1=self%grad_gwv_s1(itypat)
        tt2=self%grad_bv_s1(itypat)
        tt3=self%grad_gwv_p1(itypat)
        tt4=self%grad_gwv_p2(itypat)
        write(*,'(4es10.1)',advance='no') tt1,tt2,tt3,tt4
    enddo
    write(*,'(f8.3)') cost_gw
    write(*,'(a,i6)',advance='no') 'gwv ',istep
    do itypat=1,self%ntypat
        tt1=self%gwv_s1(itypat)
        tt2=self%bv_s1(itypat)
        tt3=self%gwv_p1(itypat)
        tt4=self%gwv_p2(itypat)
        write(*,'(4f8.3)',advance='no') tt1,tt2,tt3,tt4
    enddo
    write(*,*)
    write(*,'(a,i6)',advance='no') 'gwc ',istep
    do itypat=1,self%ntypat
        if(self%applycore(itypat)) then
            tt1=self%gwc_s1(itypat)
            tt2=self%bc_s1(itypat)
            write(*,'(2f8.3)',advance='no') tt1,tt2
        endif
    enddo
    write(*,*)
end subroutine report_fit_bf
!*****************************************************************************************
subroutine get_basis_functions_cent2(parini,atoms_arr_t,poisson_ref_t)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    !use mod_ann, only: typ_ann_arr
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_opt, only: typ_paropt
    use mod_trial_energy, only: get_trial_energy, trial_energy_deallocate
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), optional, target, intent(inout):: atoms_arr_t
    type(typ_poisson), optional, target, intent(inout):: poisson_ref_t(:)
    !local variables
    type(typ_poisson), pointer:: poisson_ref(:)
    !type(typ_atoms):: atoms
    type(typ_atoms_arr), pointer:: atoms_arr=>null()
    !type(typ_ann_arr):: ann_arr
    type(typ_fitpar):: fitpar, fitpar_t
    type(typ_bf):: bf
    type(typ_paropt):: paropt
    integer:: iat, iconf, istep, ibf, jbf, nwork
    integer:: ix, iy, iz, itypat, nr, ir
    !real(8):: ehartree_scn_excl
    !real(8):: dpx, dpy, dpz
    !real(8):: time1, time2, time3, time4, time5, time6, time7, time8, time9
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6, tt7, ttt, ttt_t, tt3_t, cost, cost_gw
    real(8):: ss1, ss2, ss3, ss4, ss5, fd_step
    real(8):: gausswidth(1), q_tmp(1), p_tmp(3), voxel, dpm(3), pi
    !real(8):: x, y, z, r
    real(8):: errmax, rmse, ener
    real(8), allocatable:: bz(:)
    real(8), allocatable:: gwz(:)
    real(8), allocatable:: xr(:)
    real(8), allocatable:: fr(:)
    real(8), allocatable:: work(:)
    character(256), allocatable:: fn_list(:)
    integer:: nfiles, nfiles_max, ifile
    call f_routine(id='get_basis_functions_cent2')
    if(present(atoms_arr_t)) then
        atoms_arr=>atoms_arr_t
    else
        allocate(atoms_arr)
    endif
    pi=4.d0*atan(1.d0)
    !call read_data_yaml(parini,'list_posinp_cent2.yaml',atoms_arr)
    nfiles_max=100
    allocate(fn_list(nfiles_max))
    if(present(poisson_ref_t)) then
        poisson_ref=>poisson_ref_t
        nfiles=size(poisson_ref_t)
    else
        call read_list_files_yaml('list_cubes.yaml',nfiles_max,fn_list,nfiles)
        allocate(poisson_ref(nfiles))
        atoms_arr%nconf=nfiles
        allocate(atoms_arr%atoms(nfiles))
        do ifile=1,nfiles
            call cube_read(trim(fn_list(ifile)),atoms_arr%atoms(ifile),poisson_ref(ifile))
        enddo
    endif
    do ifile=1,nfiles
        do iat=1,atoms_arr%atoms(ifile)%nat
            do itypat=1,parini%ntypat
                if(trim(atoms_arr%atoms(ifile)%sat(iat))==trim(parini%stypat(itypat))) then
                    atoms_arr%atoms(ifile)%itypat(iat)=parini%ltypat(itypat)
                    exit
                endif
            enddo
        enddo
    enddo
    !write(*,*) nfiles
    !write(*,*) fn_list(1:nfiles)
    !stop 'OOOOOOOOOOOOOO'

    call fitpar%init_fit_bf(parini%ntypat)
    call fitpar_t%init_fit_bf(parini%ntypat)

    open(unit=2134,file="gw.inp",status='old')
    do itypat=1,parini%ntypat
        read(2134,*) fitpar%gwz_s(itypat),fitpar%bz_s(itypat),fitpar%relaxcore(itypat),fitpar%applycore(itypat),tt1,tt2,ttt,tt4,ttt_t,tt6,tt7
        fitpar%qcore_type(itypat)=tt1
        fitpar%gwc_s1(itypat)=tt2
        fitpar%bc_s1(itypat)=ttt
        fitpar%gwv_s1(itypat)=tt4
        fitpar%bv_s1(itypat)=ttt_t
        fitpar%gwv_p1(itypat)=tt6
        fitpar%gwv_p2(itypat)=tt7
    enddo
    close(2134)

    allocate(trial_energy_all(nfiles))
    do ifile=1,nfiles
    !!!!!!!!!!!!!!!!
    associate(atoms=>atoms_arr%atoms(ifile))
    allocate(bz(atoms%nat),gwz(atoms%nat))
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        bz(iat)=fitpar%bz_s(itypat)
        gwz(iat)=fitpar%gwz_s(itypat)
    enddo
    call get_trial_energy(parini,atoms,poisson_ref(ifile),2*atoms%nat,bz,gwz,trial_energy_all(ifile),atoms%qtot,atoms%dpm)
    deallocate(bz,gwz)
    end associate
    enddo

    !do ifile=1,nfiles
    !    call get_poisson_ref(parini,fitpar%gausswidth_ion,poisson_ref(ifile),atoms_arr%atoms(ifile),fitpar)
    !enddo
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,4i5)') 'ngpx,ngpy,ngpz= ',poisson_ref(1)%ngpx,poisson_ref(1)%ngpy,poisson_ref(1)%ngpz
    endif
    !stop
    if(nfiles==1) then
        allocate(c_s_t(atoms_arr%atoms(1)%nat))
        allocate(c_p_t(atoms_arr%atoms(1)%nat))
    endif
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms_arr%atoms(1)%nat
        write(*,*) atoms_arr%atoms(1)%ratp(1,iat),atoms_arr%atoms(1)%ratp(2,iat),atoms_arr%atoms(1)%ratp(3,iat)
    enddo
    write(*,*) poisson_ref(1)%xyz111(1),poisson_ref(1)%xyz111(2),poisson_ref(1)%xyz111(3)
    endif
!    do ifile=1,nfiles
!    !do iz=poisson_ref(ifile)%ngpz/2,poisson_ref(ifile)%ngpz/2
!    !do iy=poisson_ref(ifile)%ngpy/2,poisson_ref(ifile)%ngpy/2
!    do iz=1+20,poisson_ref(ifile)%ngpz-20
!    do iy=1+20,poisson_ref(ifile)%ngpy-20
!    do ix=1+20,poisson_ref(ifile)%ngpx-20
!        x=poisson_ref(ifile)%xyz111(1)+(ix-1)*poisson_ref(ifile)%hgrid(1,1)
!        y=poisson_ref(ifile)%xyz111(2)+(iy-1)*poisson_ref(ifile)%hgrid(2,2)
!        z=poisson_ref(ifile)%xyz111(3)+(iz-1)*poisson_ref(ifile)%hgrid(3,3)
!        !write(70+ifile,'(3f8.3,es19.10)') x,y,z,poisson_ref(ifile)%rho(ix,iy,iz)
!        write(70,'(4f8.3,es19.10)') x,y,z,atoms_arr%atoms(ifile)%qtot,poisson_ref(ifile)%rho(ix,iy,iz)
!    enddo
!    enddo
!    enddo
!    enddo
!    stop 'PPPPPPPPPPPPPP'
    !do iat=1,atoms%nat
    !    if(trim(atoms%sat(iat))=='Li') atoms%zat(iat)=1.d0
    !    if(trim(atoms%sat(iat))=='S' ) atoms%zat(iat)=4.d0
    !enddo
    !-------------------------------------------------------

    gausswidth=fitpar%gausswidth_ion

    !q_tmp=2.d0
    !call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,1),q_tmp,gausswidth, &
    !    rgcut*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    call fitpar%set_bc_bcg()

    !do iat=1,atoms%nat
    !    if(trim(atoms%sat(iat))=='O') atoms%zat(iat)=4.d0
    !enddo

    !-----------------------------------------------------------------
    !q_tmp(1)=1.d0
    !call put_r2gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,1),q_tmp,fitpar%gwc_s2(1), &
    !    rgcut*fitpar%gwc_s2(1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !ener=0.d0
    !do iz=1,poisson%ngpz
    !do iy=1,poisson%ngpy
    !do ix=1,poisson%ngpx
    !    ener=ener+poisson%rho(ix,iy,iz)
    !enddo
    !enddo
    !enddo
    !ener=ener*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    !write(*,*) 'Total charge from put_r2gto_sym_ortho= ',ener
    !call get_hartree(parini,poisson,atoms,gausswidth,ener)
    !call cal_pot_r2gauss_s(parini,poisson,bf%wa1,.true.,atoms%ratp(1,1),fitpar%gwc_s2(1),q_tmp(1),bf%vgrad_1)
    !fitpar%gwc_s2(1)=fitpar%gwc_s2(1)+1.d-3
    !call cal_pot_r2gauss_s(parini,poisson,bf%wa2,.true.,atoms%ratp(1,1),fitpar%gwc_s2(1),q_tmp(1),bf%vgrad_2)
    !!do iz=1,poisson%ngpz
    !!do iy=1,poisson%ngpy
    !do iz=poisson%ngpz/2,poisson%ngpz/2
    !do iy=poisson%ngpy/2,poisson%ngpy/2
    !do ix=1,poisson%ngpx
    !    x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
    !    y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
    !    z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
    !    !write(61,'(3f8.3,2es19.10)') x,y,z,bf%wa1(ix,iy,iz),poisson%pot(ix,iy,iz)
    !    write(61,'(3f8.3,2es19.10)') x,y,z,bf%vgrad_1(ix,iy,iz),(bf%wa2(ix,iy,iz)-bf%wa1(ix,iy,iz))/1.d-3
    !enddo
    !enddo
    !enddo
    !stop 'UUUUUUUUUUUUUUUUU'
    !-----------------------------------------------------------------
    !gausswidth=0.5d0

    paropt%converged=.false.
    paropt%iflag=0
    paropt%lprint=parini%paropt_geopt%lprint
    paropt%dt_start=parini%paropt_geopt%dt_start
    paropt%dtmax=parini%paropt_geopt%dtmax
    paropt%alphax=parini%paropt_geopt%alphax
    paropt%dxmax=parini%paropt_geopt%dxmax
    paropt%fmaxtol=parini%paropt_geopt%fmaxtol
    paropt%nit=parini%paropt_geopt%nit
    paropt%condnum=parini%paropt_geopt%condnum
    paropt%funits=40.d0
    paropt%finc=1.2d0
    paropt%fdec=0.1d0
    paropt%ndowntol=0
    paropt%itfire=0
    nr=4*parini%ntypat !+3
    do itypat=1,parini%ntypat
        if(fitpar%relaxcore(itypat) .and. fitpar%applycore(itypat)) then
            nr=nr+1
            nr=nr+1
        endif
    enddo
    allocate(xr(nr),fr(nr))
    !nwork=3*nr !fire
    nwork=nr*nr+3*nr+3*nr*nr+3*nr !mybfgs
    allocate(work(nwork))
    ir=0
    do itypat=1,parini%ntypat
        if(fitpar%relaxcore(itypat) .and. fitpar%applycore(itypat)) then
            ir=ir+1 ; xr(ir)=fitpar%gwc_s1(itypat)
            ir=ir+1 ; xr(ir)=fitpar%bc_s1(itypat)
        endif
        ir=ir+1 ; xr(ir)=fitpar%bv_s1(itypat)
        ir=ir+1 ; xr(ir)=fitpar%gwv_s1(itypat)
        ir=ir+1 ; xr(ir)=fitpar%gwv_p1(itypat)
        ir=ir+1 ; xr(ir)=fitpar%gwv_p2(itypat)
    enddo
    !---------------------------------------------------------------------------
!    do itypat=1,parini%ntypat
!        fitpar%grad_gwc_s1(itypat)=0.d0
!        fitpar%grad_bc_s1(itypat)=0.d0
!        fitpar%grad_gwv_s1(itypat)=0.d0
!        fitpar%grad_bv_s1(itypat)=0.d0
!        fitpar%grad_gwv_p1(itypat)=0.d0
!        fitpar%grad_gwv_p2(itypat)=0.d0
!    enddo
!    call get_cost_from_pot(parini,istep,nfiles,poisson_ref,atoms_arr,fitpar,bf,cost)
!    tt1=cost
!    fd_step=1.d-6
!
!    !tt2=fitpar%grad_bv_s1(1)
!    !fitpar%bv_s1(1)=fitpar%bv_s1(1)+fd_step
!
!    !tt2=fitpar%grad_gwv_s1(2)
!    !fitpar%gwv_s1(2)=fitpar%gwv_s1(2)+fd_step
!
!    !tt2=fitpar%grad_gwv_p2(2)
!    !fitpar%gwv_p2(2)=fitpar%gwv_p2(2)+fd_step
!
!    tt2=fitpar%grad_gwv_p1(1)
!    fitpar%gwv_p1(1)=fitpar%gwv_p1(1)+fd_step
!
!    !tt2=fitpar%grad_bc_s1(1)
!    !fitpar%bc_s1(1)=fitpar%bc_s1(1)+fd_step
!
!    !tt2=fitpar%grad_gwc_s1(2)
!    !fitpar%gwc_s1(2)=fitpar%gwc_s1(2)+fd_step
!
!    do itypat=1,parini%ntypat
!        fitpar%grad_gwc_s1(itypat)=0.d0
!        fitpar%grad_bc_s1(itypat)=0.d0
!        fitpar%grad_gwv_s1(itypat)=0.d0
!        fitpar%grad_bv_s1(itypat)=0.d0
!        fitpar%grad_gwv_p1(itypat)=0.d0
!        fitpar%grad_gwv_p2(itypat)=0.d0
!    enddo
!    call get_cost_from_pot(parini,istep,nfiles,poisson_ref,atoms_arr,fitpar,bf,cost)
!    tt3=(cost-tt1)/fd_step
!    write(*,'(a,3es14.5)') 'FD= ',tt2,tt3,tt2-tt3
!    stop 'PPPPPPPPPPPPPPPPPPPPPPP'


    !---------------------------------------------------------------------------
    istep=0
    do

    do itypat=1,parini%ntypat
        fitpar%grad_gwc_s1(itypat)=0.d0
        fitpar%grad_bc_s1(itypat)=0.d0
        fitpar%grad_gwv_s1(itypat)=0.d0
        fitpar%grad_bv_s1(itypat)=0.d0
        fitpar%grad_gwv_p1(itypat)=0.d0
        fitpar%grad_gwv_p2(itypat)=0.d0
    enddo
    call get_cost_from_pot(parini,istep,nfiles,poisson_ref,atoms_arr,fitpar,bf,cost)
    !call get_cost_from_density(parini,istep,nfiles,poisson_ref,atoms_arr,fitpar,bf,cost)

    !call get_cost_gw(parini,fitpar,cost_gw)
    !cost=cost+cost_gw
!    call get_cost_gw_new(parini,fitpar,fitpar_t,cost_gw)
!    do itypat=1,parini%ntypat
!        !write(*,'(a,es14.5,i3,4es14.5)') 'SS= ',cost_gw,itypat,fitpar%grad_gwv_p1(itypat),fitpar%grad_gwv_p2(itypat),fitpar_t%grad_gwv_p1(itypat),fitpar_t%grad_gwv_p2(itypat)
!        fitpar%grad_gwc_s1(itypat)=cost_gw*fitpar%grad_gwc_s1(itypat)+cost*fitpar_t%grad_gwc_s1(itypat)
!        fitpar%grad_gwv_s1(itypat)=cost_gw*fitpar%grad_gwv_s1(itypat)+cost*fitpar_t%grad_gwv_s1(itypat)
!        !fitpar%grad_gwv_s2(itypat)=cost_gw*fitpar%grad_gwv_s2(itypat)+cost*fitpar_t%grad_gwv_s2(itypat)
!        fitpar%grad_gwv_p1(itypat)=cost_gw*fitpar%grad_gwv_p1(itypat)+cost*fitpar_t%grad_gwv_p1(itypat)
!        fitpar%grad_gwv_p2(itypat)=cost_gw*fitpar%grad_gwv_p2(itypat)+cost*fitpar_t%grad_gwv_p2(itypat)
!    enddo
!    cost=cost*cost_gw
    cost_gw=1.d0


    if(parini%mpi_env%iproc==0) then
        call fitpar%report_fit_bf(istep,cost,cost_gw)
    endif

    if(istep==parini%paropt_geopt%nit) exit
    ir=0
    do itypat=1,parini%ntypat
        if(fitpar%relaxcore(itypat) .and. fitpar%applycore(itypat)) then
            ir=ir+1 ; fr(ir)=-fitpar%grad_gwc_s1(itypat)
            ir=ir+1 ; fr(ir)=-fitpar%grad_bc_s1(itypat)
        endif
        ir=ir+1 ; fr(ir)=-fitpar%grad_bv_s1(itypat)
        ir=ir+1 ; fr(ir)=-fitpar%grad_gwv_s1(itypat)
        ir=ir+1 ; fr(ir)=-fitpar%grad_gwv_p1(itypat)
        ir=ir+1 ; fr(ir)=-fitpar%grad_gwv_p2(itypat)
    enddo
    !call fire(parini,parini%mpi_env%iproc,nr,xr,cost,fr,work,paropt)
    call mybfgs(parini%mpi_env%iproc,nr,xr,cost,fr,nwork,work,paropt)
    if(paropt%iflag<=0) exit
    ir=0
    do itypat=1,parini%ntypat
        if(fitpar%relaxcore(itypat) .and. fitpar%applycore(itypat)) then
            ir=ir+1 ; fitpar%gwc_s1(itypat)=xr(ir)
            ir=ir+1 ; fitpar%bc_s1(itypat)=xr(ir)
        endif
        ir=ir+1 ; fitpar%bv_s1(itypat)=xr(ir)
        ir=ir+1 ; fitpar%gwv_s1(itypat)=xr(ir)
        ir=ir+1 ; fitpar%gwv_p1(itypat)=xr(ir)
        ir=ir+1 ; fitpar%gwv_p2(itypat)=xr(ir)
    enddo
    call fitpar%set_bc_bcg()
    istep=istep+1
    enddo !end of loop over istep
    deallocate(xr,fr,work)
    call yaml_sequence_close()
    !-----------------------------------------------------------------
    if(parini%mpi_env%iproc==0) then
    open(unit=2134,file="gw.out",status='replace')
    do itypat=1,parini%ntypat
        tt1=fitpar%qcore_type(itypat)
        tt2=fitpar%gwc_s1(itypat)
        ttt=fitpar%bc_s1(itypat)
        tt4=fitpar%gwv_s1(itypat)
        ttt_t=fitpar%bv_s1(itypat)
        tt6=fitpar%gwv_p1(itypat)
        tt7=fitpar%gwv_p2(itypat)
        write(2134,'(2f10.6,2l2,7f10.6)') fitpar%gwz_s(itypat),fitpar%bz_s(itypat),fitpar%relaxcore(itypat),fitpar%applycore(itypat),tt1,tt2,ttt,tt4,ttt_t,tt6,tt7
    enddo
    close(2134)
    endif
    !-----------------------------------------------------------------
    if(nfiles==1) then
        call final_report(parini,fitpar,poisson_ref(1),atoms_arr%atoms(1))
    endif
    !-----------------------------------------------------------------
    do ifile=1,nfiles
    call fini_hartree(parini,atoms_arr%atoms(ifile),poisson_ref(ifile))
    enddo
    if(nfiles==1) then
        deallocate(c_s_t)
        deallocate(c_p_t)
    endif
    do iconf=1,atoms_arr%nconf
        call atom_deallocate_old(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
    call fitpar%fini_fit_bf()
    call fitpar_t%fini_fit_bf()
    do ifile=1,nfiles
        call trial_energy_deallocate(trial_energy_all(ifile))
    enddo
    deallocate(trial_energy_all)
    !write(*,'(a,f10.3)') 'time_cal_pot_gauss_p= ',time_cal_pot_gauss_p
    !write(*,'(a,f10.3)') 'time_cal_pot_gauss_s= ',time_cal_pot_gauss_s
    !write(*,'(a,f10.3)') 'time_cal_pot_r2gauss_s= ',time_cal_pot_r2gauss_s
    call f_release_routine()
end subroutine get_basis_functions_cent2
!*****************************************************************************************
subroutine final_report(parini,fitpar,poisson_ref,atoms)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_fitpar), intent(in):: fitpar
    type(typ_poisson), intent(in):: poisson_ref
    type(typ_atoms), intent(in):: atoms
    !integer, intent(in):: istep, nfiles
    !type(typ_atoms_arr), intent(in):: atoms_arr
    !type(typ_bf), intent(inout):: bf
    !real(8), intent(out):: cost
    !local variables
    type(typ_poisson):: poisson
    integer:: ix, iy, iz, iat, itypat
    real(8):: errmax, rmse, ener, x, y, z, r
    real(8):: gausswidth(1), q_tmp(1), p_tmp(3), voxel
    gausswidth=fitpar%gausswidth_ion
    voxel=poisson_ref%hgrid(1,1)*poisson_ref%hgrid(2,2)*poisson_ref%hgrid(3,3)
    poisson%ngpx=poisson_ref%ngpx
    poisson%ngpy=poisson_ref%ngpy
    poisson%ngpz=poisson_ref%ngpz
    poisson%hgrid(1:3,1:3)=0.d0
    poisson%hgrid(1,1)=poisson_ref%hgrid(1,1)
    poisson%hgrid(2,2)=poisson_ref%hgrid(2,2)
    poisson%hgrid(3,3)=poisson_ref%hgrid(3,3)
    poisson%xyz111=poisson_ref%xyz111
    poisson%cal_scn=poisson_ref%cal_scn
    poisson%screening_factor=poisson_ref%screening_factor
    poisson%bc=poisson_ref%bc
    poisson%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson,gausswidth)
    !----------------------------------------------------------
    poisson%rho=0.d0
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
!    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),atoms%zat(iat),gausswidth(1), &
!        rgcut*gausswidth(1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
!        poisson%hgrid,poisson%rho)
    if(fitpar%applycore(itypat)) then
    q_tmp(1)=fitpar%qcore_type(itypat)*fitpar%bc_s1(itypat)
    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwc_s1(itypat), &
        rgcut*fitpar%gwc_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    q_tmp(1)=fitpar%qcore_type(itypat)*(1.d0-fitpar%bc_s1(itypat))
    call put_r2gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwc_s1(itypat), &
        rgcut*fitpar%gwc_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    endif
    q_tmp(1)=c_s_t(iat)*fitpar%bv_s1(itypat)
    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwv_s1(itypat), &
        rgcut*fitpar%gwv_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    q_tmp(1)=c_s_t(iat)*(1.d0-fitpar%bv_s1(itypat))
    call put_r2gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwv_s1(itypat), &
        rgcut*fitpar%gwv_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    p_tmp=0.d0
    p_tmp(1)=c_p_t(iat)*fitpar%bv_p1(itypat)
    call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),p_tmp,fitpar%gwv_p1(itypat), &
        rgcut*fitpar%gwv_p1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    p_tmp(1)=c_p_t(iat)*fitpar%bv_p2(itypat)
    call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),p_tmp,fitpar%gwv_p2(itypat), &
        rgcut*fitpar%gwv_p2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    enddo
    !poisson%rho=poisson%rho-poisson_ref%rho
    !call cube_write('diffrho.cube',atoms,poisson,'rho')
    call get_hartree(parini,poisson,atoms,gausswidth,ener)
    if(allocated(poisson%pot)) then
    errmax=0.d0
    rmse=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        errmax=max(errmax,abs(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)))
        rmse=rmse+(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))**2
        !write(63,'(3i5,2f20.12)') ix,iy,iz,poisson_ref%pot(ix,iy,iz),poisson%pot(ix,iy,iz)
    enddo
    enddo
    enddo
    rmse=sqrt(rmse/real(poisson%ngpx*poisson%ngpy*poisson%ngpz,kind=8))
    write(*,'(a,2es14.5)') 'errmax,rmse= ',errmax,rmse
    endif
    if(allocated(poisson%pot)) then
    !do iz=1,poisson%ngpz
    do iz=poisson%ngpz/2,poisson%ngpz/2
    do iy=poisson%ngpy/2,poisson%ngpy/2
    !do ix=poisson%ngpx/2,poisson%ngpx/2
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
        write(83,'(3f8.3,4es19.10)') x,y,z,poisson_ref%pot(ix,iy,iz),poisson%pot(ix,iy,iz),poisson_ref%rho(ix,iy,iz),poisson%rho(ix,iy,iz)
    enddo
    enddo
    enddo
    endif
    if(allocated(poisson%pot)) then
    ener=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        ener=ener+poisson%rho(ix,iy,iz)*poisson%pot(ix,iy,iz)
    enddo
    enddo
    enddo
    ener=ener*0.5d0*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    !if(parini%mpi_env%iproc==0) then
    write(*,'(a,es24.10,es14.5)') 'ehartree_scn_excl ',ener,poisson%screening_factor
    !endif
    endif
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-atoms%ratp(1,1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-atoms%ratp(2,1)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-atoms%ratp(3,1)
        r=sqrt(x**2+y**2+z**2)
        if(r<4.d0) then
        write(51,'(3es14.5)') r,poisson%rho(ix,iy,iz),poisson_ref%rho(ix,iy,iz)
        endif
    enddo
    enddo
    enddo
    call fini_hartree(parini,atoms,poisson)
end subroutine final_report
!*****************************************************************************************
subroutine get_cost_from_density(parini,istep,nfiles,poisson_ref,atoms_arr,fitpar,bf,cost)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: istep, nfiles
    type(typ_poisson), intent(in):: poisson_ref(nfiles)
    type(typ_atoms_arr), intent(in):: atoms_arr
    type(typ_fitpar), intent(inout):: fitpar
    type(typ_bf), intent(inout):: bf
    real(8), intent(out):: cost
    !local variables
    integer:: ifile, ix, iy, iz, ibf, jbf, itypat, iat
    real(8):: tt1, tt2, tt3, tt4, tt5, ss1, ss2, ss3, ss4, ss5
    real(8):: dpm(3), voxel, pi, gausswidth(1), center(3), q
    type(typ_poisson):: poisson
    real(8), allocatable:: c_s(:)
    real(8), allocatable:: c_p(:)
    real(8), allocatable:: rho_ion(:,:,:)
    pi=4.d0*atan(1.d0)
    gausswidth(1)=fitpar%gausswidth_ion
    cost=0.d0
    do ifile=1,nfiles
    allocate(c_s(atoms_arr%atoms(ifile)%nat))
    allocate(c_p(atoms_arr%atoms(ifile)%nat))
    poisson%ngpx=poisson_ref(ifile)%ngpx
    poisson%ngpy=poisson_ref(ifile)%ngpy
    poisson%ngpz=poisson_ref(ifile)%ngpz
    poisson%hgrid(1:3,1:3)=0.d0
    poisson%hgrid(1,1)=poisson_ref(ifile)%hgrid(1,1)
    poisson%hgrid(2,2)=poisson_ref(ifile)%hgrid(2,2)
    poisson%hgrid(3,3)=poisson_ref(ifile)%hgrid(3,3)
    poisson%xyz111=poisson_ref(ifile)%xyz111
    poisson%cal_scn=poisson_ref(ifile)%cal_scn
    poisson%screening_factor=poisson_ref(ifile)%screening_factor
    poisson%bc=poisson_ref(ifile)%bc
    poisson%task_finit="alloc_rho"
    call init_hartree(parini,atoms_arr%atoms(ifile),poisson,gausswidth)
    allocate(rho_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    call bf%init_bf(atoms_arr%atoms(ifile)%nat,atoms_arr%atoms(ifile)%ratp,poisson)
    call bf%set_bf_bc_bcg(atoms_arr%atoms(ifile),fitpar)
    do iat=1,atoms_arr%atoms(ifile)%nat
        itypat=atoms_arr%atoms(ifile)%itypat(iat)
        bf%qcore(iat)=fitpar%qcore_type(itypat)
    enddo
    rho_ion=0.d0
!    do iat=1,atoms_arr%atoms(ifile)%nat
!    itypat=atoms_arr%atoms(ifile)%itypat(iat)
!    !call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms_arr%atoms(ifile)%ratp(1,iat),gausswidth(1),atoms_arr%atoms(ifile)%zat(iat),bf%vgrad_1)
!    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms_arr%atoms(ifile)%ratp(1,iat),atoms_arr%atoms(ifile)%zat(iat),gausswidth(1), &
!        rgcut*gausswidth(1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
!        poisson%hgrid,rho_ion)
!    enddo

    poisson%rho=rho_ion

    voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)

    do iat=1,atoms_arr%atoms(ifile)%nat
        itypat=atoms_arr%atoms(ifile)%itypat(iat)
        call bf%get_rho_single_core(parini,poisson,atoms_arr%atoms(ifile),fitpar,iat)
        tt1=bf%qcore(iat)*fitpar%bc_s1(itypat)
        tt2=bf%qcore(iat)*(1.d0-fitpar%bc_s1(itypat))
        poisson%rho=poisson%rho+tt1*bf%wa1+tt2*bf%wa2
    enddo
    call get_linearcoeff_rho(parini,fitpar,bf,poisson,poisson_ref(ifile),atoms_arr%atoms(ifile),c_s,c_p)
    do ibf=1,bf%nbf
        call bf%get_rho_single(parini,poisson,atoms_arr%atoms(ifile),fitpar,ibf)
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=bf%bv_1(ibf)*bf%pot_1(ix,iy,iz)+bf%bv_2(ibf)*bf%pot_2(ix,iy,iz)
            poisson%rho(ix,iy,iz)=poisson%rho(ix,iy,iz)+bf%csp(ibf)*tt1
        enddo
        enddo
        enddo
    enddo
    if(any(fitpar%relaxcore)) then
    do iat=1,atoms_arr%atoms(ifile)%nat
        itypat=atoms_arr%atoms(ifile)%itypat(iat)
        call bf%get_rho_single_core(parini,poisson,atoms_arr%atoms(ifile),fitpar,iat)
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        tt4=0.d0
        tt5=0.d0
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            ss1=fitpar%bc_s1(itypat)*bf%wa3(ix,iy,iz)
            ss2=(1.d0-fitpar%bc_s1(itypat))*bf%wa4(ix,iy,iz)
            ss4=bf%wa1(ix,iy,iz)-bf%wa2(ix,iy,iz)
            tt1=tt1+bf%qcore(iat)*ss1*(poisson%rho(ix,iy,iz)-poisson_ref(ifile)%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
            tt2=tt2+bf%qcore(iat)*ss2*(poisson%rho(ix,iy,iz)-poisson_ref(ifile)%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
            tt4=tt4+bf%qcore(iat)*ss4*(poisson%rho(ix,iy,iz)-poisson_ref(ifile)%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
        enddo
        enddo
        enddo
        fitpar%grad_gwc_s1(itypat)=fitpar%grad_gwc_s1(itypat)+2.d0*tt1*voxel
        fitpar%grad_gwc_s1(itypat)=fitpar%grad_gwc_s1(itypat)+2.d0*tt2*voxel
        fitpar%grad_bc_s1(itypat)=fitpar%grad_bc_s1(itypat)+2.d0*tt4*voxel
    enddo
    endif
    do ibf=1,bf%nbf
        iat=bf%iat_list(ibf)
        itypat=atoms_arr%atoms(ifile)%itypat(iat)
        call bf%get_rho_single(parini,poisson,atoms_arr%atoms(ifile),fitpar,ibf)
        tt5=0.d0
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            ss5=bf%pot_1(ix,iy,iz)-bf%pot_2(ix,iy,iz)
            tt5=tt5+c_s(iat)*ss5*(poisson%rho(ix,iy,iz)-poisson_ref(ifile)%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
        enddo
        enddo
        enddo
        fitpar%grad_bv_s1(itypat)=fitpar%grad_bv_s1(itypat)+2.d0*tt5*voxel
    enddo
    !----------------------------------------------------------
    cost=cost+bf%zeroder
    do jbf=1,bf%nbf
        do ibf=1,bf%nbf
            cost=cost+bf%csp(ibf)*bf%secder(ibf,jbf)*bf%csp(jbf)
        enddo
    enddo
    do ibf=1,bf%nbf
        cost=cost+2.d0*bf%csp(ibf)*bf%firstder(ibf)
    enddo
    do ibf=1,bf%nbf
        iat=bf%iat_list(ibf)
        itypat=atoms_arr%atoms(ifile)%itypat(iat)
        tt1=bf%csp(ibf)*(bf%grad_bv_1(1,ibf)*bf%subfirstder(1,ibf)+bf%grad_bv_2(1,ibf)*bf%subfirstder(2,ibf))
        tt2=bf%csp(ibf)*(bf%grad_bv_1(2,ibf)*bf%subfirstder(1,ibf)+bf%grad_bv_2(2,ibf)*bf%subfirstder(2,ibf))
        if(trim(bf%bt(ibf))=='s') then
            fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*tt1
            fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*tt2
        elseif(trim(bf%bt(ibf))=='px') then
            fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*tt1
            fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*tt2
        else
            stop 'ERROR: unknwon br when computing gwgrad'
        endif
        do jbf=1,bf%nbf
            tt1=bf%csp(ibf)*bf%grad_bv_1(1,ibf)*bf%subsecder(1,1,ibf,jbf)*bf%csp(jbf)*bf%bv_1(jbf)
            tt2=bf%csp(ibf)*bf%grad_bv_1(1,ibf)*bf%subsecder(1,2,ibf,jbf)*bf%csp(jbf)*bf%bv_2(jbf)
            tt3=bf%csp(ibf)*bf%grad_bv_2(1,ibf)*bf%subsecder(2,1,ibf,jbf)*bf%csp(jbf)*bf%bv_1(jbf)
            tt4=bf%csp(ibf)*bf%grad_bv_2(1,ibf)*bf%subsecder(2,2,ibf,jbf)*bf%csp(jbf)*bf%bv_2(jbf)
            if(trim(bf%bt(ibf))=='s') then
                fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            elseif(trim(bf%bt(ibf))=='px') then
                fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            else
                stop 'ERROR: unknwon br when computing gwgrad'
            endif
            tt1=bf%csp(ibf)*bf%grad_bv_1(2,ibf)*bf%subsecder(1,1,ibf,jbf)*bf%csp(jbf)*bf%bv_1(jbf)
            tt2=bf%csp(ibf)*bf%grad_bv_1(2,ibf)*bf%subsecder(1,2,ibf,jbf)*bf%csp(jbf)*bf%bv_2(jbf)
            tt3=bf%csp(ibf)*bf%grad_bv_2(2,ibf)*bf%subsecder(2,1,ibf,jbf)*bf%csp(jbf)*bf%bv_1(jbf)
            tt4=bf%csp(ibf)*bf%grad_bv_2(2,ibf)*bf%subsecder(2,2,ibf,jbf)*bf%csp(jbf)*bf%bv_2(jbf)
            if(trim(bf%bt(ibf))=='s') then
                fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            elseif(trim(bf%bt(ibf))=='px') then
                fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            else
                stop 'ERROR: unknwon br when computing gwgrad'
            endif
        enddo
        tt1=bf%csp(ibf)*(bf%bv_1(ibf)*bf%subfirstder_grad(1,ibf))
        tt2=bf%csp(ibf)*(bf%bv_2(ibf)*bf%subfirstder_grad(2,ibf))
        if(trim(bf%bt(ibf))=='s') then
            fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*tt1
            fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*tt2
        elseif(trim(bf%bt(ibf))=='px') then
            fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*tt1
            fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*tt2
        else
            stop 'ERROR: unknwon br when computing gwgrad'
        endif
        do jbf=1,bf%nbf
            tt1=bf%csp(ibf)*bf%bv_1(ibf)*bf%subsecder_grad(1,1,ibf,jbf)*bf%csp(jbf)*bf%bv_1(jbf)
            tt2=bf%csp(ibf)*bf%bv_1(ibf)*bf%subsecder_grad(1,2,ibf,jbf)*bf%csp(jbf)*bf%bv_2(jbf)
            tt3=bf%csp(ibf)*bf%bv_2(ibf)*bf%subsecder_grad(2,1,ibf,jbf)*bf%csp(jbf)*bf%bv_1(jbf)
            tt4=bf%csp(ibf)*bf%bv_2(ibf)*bf%subsecder_grad(2,2,ibf,jbf)*bf%csp(jbf)*bf%bv_2(jbf)
            if(trim(bf%bt(ibf))=='s') then
                fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*(tt1+tt2)
                fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*(tt3+tt4)
            elseif(trim(bf%bt(ibf))=='px') then
                fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*(tt1+tt2)
                fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*(tt3+tt4)
            else
                stop 'ERROR: unknwon br when computing gwgrad'
            endif
        enddo
    enddo
    !-----------------------------------------------------------------
    write(*,'(a,i6)',advance='no') 'qp ',istep
    center(1)=0.d0
    center(2)=0.d0
    center(3)=0.d0
    do iat=1,atoms_arr%atoms(ifile)%nat
        center(1)=center(1)+atoms_arr%atoms(ifile)%ratp(1,iat)
        center(2)=center(2)+atoms_arr%atoms(ifile)%ratp(2,iat)
        center(3)=center(3)+atoms_arr%atoms(ifile)%ratp(3,iat)
    enddo
    center(1)=center(1)/atoms_arr%atoms(ifile)%nat
    center(2)=center(2)/atoms_arr%atoms(ifile)%nat
    center(3)=center(3)/atoms_arr%atoms(ifile)%nat
    dpm(1)=0.d0
    dpm(2)=0.d0
    dpm(3)=0.d0
    do iat=1,atoms_arr%atoms(ifile)%nat
        itypat=atoms_arr%atoms(ifile)%itypat(iat)
        q=(atoms_arr%atoms(ifile)%zat(iat)+fitpar%qcore_type(itypat)+c_s(iat))
        dpm(1)=dpm(1)+q*(atoms_arr%atoms(ifile)%ratp(1,iat)-center(1))
        write(*,'(f8.3)',advance='no') c_s(iat)
    enddo
    dpm(1)=dpm(1)+sum(c_p(:))
    write(*,'(2f8.3)') sum(c_p(:)),dpm(1)
    write(*,'(a,i6)',advance='no') 'ps ',istep
    do iat=1,atoms_arr%atoms(ifile)%nat
        write(*,'(f8.3)',advance='no') c_p(iat)
    enddo
    write(*,*)
    write(*,'(a,i6)',advance='no') 'origin ',istep
    do itypat=1,fitpar%ntypat
        tt1=fitpar%bc_s1(itypat)/(pi**1.5d0*fitpar%gwc_s1(itypat)**3)
        write(*,'(f8.3)',advance='no') tt1
    enddo
    write(*,*)
    !-----------------------------------------------------------------

    deallocate(rho_ion)
    !if(.not.(istep==parini%paropt_geopt%nit .and. nfiles==1)) then
        deallocate(c_s)
        deallocate(c_p)
        call fini_hartree(parini,atoms_arr%atoms(ifile),poisson)
    !endif
    if(nfiles==1) then
    do iat=1,atoms_arr%atoms(1)%nat
        c_s_t(iat)=bf%csp(iat+0        )
        c_p_t(iat)=bf%csp(iat+atoms_arr%atoms(1)%nat)
    enddo
    endif
    call bf%fini_bf()
    enddo !end of loop over ifile
end subroutine get_cost_from_density
!*****************************************************************************************
subroutine get_cost_from_pot(parini,istep,nfiles,poisson_ref,atoms_arr,fitpar,bf,cost)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms_arr
    use mod_trial_energy, only: get_rmse
    use mod_cent2, only: cal_rho_pot_integral_local
    use mod_processors, only: get_proc_stake
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: istep, nfiles
    type(typ_poisson), intent(in):: poisson_ref(nfiles)
    type(typ_atoms_arr), intent(in):: atoms_arr
    type(typ_fitpar), intent(inout):: fitpar
    type(typ_bf), intent(inout):: bf
    real(8), intent(out):: cost
    !local variables
    integer:: ifile, ix, iy, iz, ibf, jbf, itypat, iat, itrial
    integer:: itrials, itriale
    integer:: icount_rate, itime0, itime1, itime2, itime3, itime4, itime5, itime6
    real(8):: time_linearcoeff, time_pot, time_grad, time_net
    real(8):: tt1, tt2, tt3, tt4, tt5, ss1, ss2, ss3, ss4, ss5
    real(8):: dpm(3), voxel, pi, gausswidth(1), center(3), q
    real(8):: rmse, xyz(3), ha2mha, q_tmp(1), p_tmp(3), pref
    type(typ_poisson):: poisson
    real(8), allocatable:: c_s(:)
    real(8), allocatable:: c_p(:)
    real(8), allocatable:: pot_ion(:,:,:)
    real(8), allocatable:: wa_gws(:,:,:,:)
    real(8), allocatable:: wa_bs(:,:,:,:)
    real(8), allocatable:: wa_gwp1(:,:,:,:)
    real(8), allocatable:: wa_gwp2(:,:,:,:)
    real(8), allocatable:: wa_gwc(:,:,:,:)
    real(8), allocatable:: wa_bc(:,:,:,:)
    call system_clock(itime0)
    pi=4.d0*atan(1.d0)
    ha2mha=1.d3
    gausswidth(1)=fitpar%gausswidth_ion
    call system_clock(count_rate=icount_rate)
    time_linearcoeff=0.d0
    time_pot=0.d0
    time_grad=0.d0
    cost=0.d0
    do ifile=1,nfiles
    associate(atoms=>atoms_arr%atoms(ifile))
    associate(trial_energy=>trial_energy_all(ifile))
    allocate(c_s(atoms%nat))
    allocate(c_p(atoms%nat))
    poisson%ngpx=poisson_ref(ifile)%ngpx
    poisson%ngpy=poisson_ref(ifile)%ngpy
    poisson%ngpz=poisson_ref(ifile)%ngpz
    poisson%hgrid(1:3,1:3)=0.d0
    poisson%hgrid(1,1)=poisson_ref(ifile)%hgrid(1,1)
    poisson%hgrid(2,2)=poisson_ref(ifile)%hgrid(2,2)
    poisson%hgrid(3,3)=poisson_ref(ifile)%hgrid(3,3)
    poisson%xyz111=poisson_ref(ifile)%xyz111
    poisson%cal_scn=poisson_ref(ifile)%cal_scn
    poisson%screening_factor=poisson_ref(ifile)%screening_factor
    poisson%bc=poisson_ref(ifile)%bc
    poisson%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson,gausswidth)
    allocate(pot_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    call bf%init_bf(atoms%nat,atoms%ratp,poisson)
    call bf%set_bf_bc_bcg(atoms,fitpar)
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        bf%qcore(iat)=fitpar%qcore_type(itypat)
    enddo
    pot_ion=0.d0
!    do iat=1,atoms%nat
!    itypat=atoms%itypat(iat)
!    call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),gausswidth(1),atoms%zat(iat),bf%vgrad_1)
!    enddo

    poisson%pot=pot_ion

    voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)

    !do iat=1,atoms%nat
    !    itypat=atoms%itypat(iat)
    !    call bf%get_pot_single_core(parini,poisson,atoms,fitpar,iat)
    !    tt1=bf%qcore(iat)*fitpar%bc_s1(itypat)
    !    tt2=bf%qcore(iat)*(1.d0-fitpar%bc_s1(itypat))
    !    poisson%pot=poisson%pot+tt1*bf%wa1+tt2*bf%wa2
    !enddo
    call system_clock(itime1)
    call get_linearcoeff(parini,fitpar,bf,poisson,poisson_ref(ifile),atoms,trial_energy,c_s,c_p)
    call system_clock(itime2)
    call get_rmse(trial_energy,bf%nbf,bf%csp,rmse)
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,es14.5)') 'RMSE= ',rmse
    endif
    cost=cost+(ha2mha*rmse)**2
    !cost=rmse**2*trial_energy%ntrial
    pref=ha2mha**2/trial_energy%ntrial
    allocate(wa_gws(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat),source=0.d0)
    allocate(wa_bs(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat),source=0.d0)
    allocate(wa_gwp1(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat),source=0.d0)
    allocate(wa_gwp2(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat),source=0.d0)
    allocate(wa_gwc(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat),source=0.d0)
    allocate(wa_bc(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat),source=0.d0)
    call system_clock(itime3)
    do iat=1,atoms%nat
        ibf=bf%ibf_list_s(iat)
        itypat=atoms%itypat(iat)
        q_tmp(1)=1.d0
        call cal_pot_gauss_s(parini,poisson,bf%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwv_s1(itypat),q_tmp(1),bf%vgrad_1)
        call cal_pot_r2gauss_s(parini,poisson,bf%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwv_s1(itypat),q_tmp(1),bf%vgrad_2)
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=bf%bv_1(ibf)*bf%vgrad_1(ix,iy,iz)+bf%bv_2(ibf)*bf%vgrad_2(ix,iy,iz)
            wa_gws(ix,iy,iz,itypat)=wa_gws(ix,iy,iz,itypat)+bf%csp(ibf)*tt1
        enddo
        enddo
        enddo
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=bf%pot_1(ix,iy,iz)-bf%pot_2(ix,iy,iz)
            wa_bs(ix,iy,iz,itypat)=wa_bs(ix,iy,iz,itypat)+bf%csp(ibf)*tt1
        enddo
        enddo
        enddo
        ibf=bf%ibf_list_px(iat)
        p_tmp=0.d0
        p_tmp(1)=1.d0
        call cal_pot_gauss_p(parini,poisson,bf%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwv_p1(itypat),p_tmp(1),bf%vgrad_1)
        call cal_pot_gauss_p(parini,poisson,bf%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwv_p2(itypat),p_tmp(1),bf%vgrad_2)
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=bf%bv_1(ibf)*bf%vgrad_1(ix,iy,iz)+bf%grad_bv_1(1,ibf)*bf%pot_1(ix,iy,iz)+bf%grad_bv_2(1,ibf)*bf%pot_2(ix,iy,iz)
            wa_gwp1(ix,iy,iz,itypat)=wa_gwp1(ix,iy,iz,itypat)+bf%csp(ibf)*tt1
        enddo
        enddo
        enddo
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=bf%bv_2(ibf)*bf%vgrad_2(ix,iy,iz)+bf%grad_bv_1(2,ibf)*bf%pot_1(ix,iy,iz)+bf%grad_bv_2(2,ibf)*bf%pot_2(ix,iy,iz)
            wa_gwp2(ix,iy,iz,itypat)=wa_gwp2(ix,iy,iz,itypat)+bf%csp(ibf)*tt1
        enddo
        enddo
        enddo
        if(any(fitpar%relaxcore)) then
        q_tmp(1)=1.d0
        call cal_pot_gauss_s(parini,poisson,bf%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),bf%vgrad_1)
        call cal_pot_r2gauss_s(parini,poisson,bf%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),bf%vgrad_2)
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=fitpar%bc_s1(itypat)*bf%vgrad_1(ix,iy,iz)+(1.d0-fitpar%bc_s1(itypat))*bf%vgrad_2(ix,iy,iz)
            wa_gwc(ix,iy,iz,itypat)=wa_gwc(ix,iy,iz,itypat)+bf%qcore(iat)*tt1
        enddo
        enddo
        enddo
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=bf%pot_1(ix,iy,iz)-bf%pot_2(ix,iy,iz)
            wa_bc(ix,iy,iz,itypat)=wa_bc(ix,iy,iz,itypat)+bf%qcore(iat)*tt1
        enddo
        enddo
        enddo
        endif
    enddo
    call system_clock(itime4)
    call get_proc_stake(parini%mpi_env,trial_energy%ntrial,itrials,itriale)
    do itrial=itrials,itriale
        xyz(1:3)=atoms%ratp(1:3,trial_energy%iat_list(itrial))+trial_energy%disp(1:3,itrial)
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,xyz,1.d0,1.d0, &
            5.d0*1.d0,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        do itypat=1,parini%ntypat
        call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,6.d0, &
            poisson%rho,wa_gws(1,1,1,itypat),tt1)
        tt2=trial_energy%E_all(itrial)-trial_energy%energy(itrial)
        fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*tt1*tt2*pref
        enddo
        do itypat=1,parini%ntypat
        call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,6.d0, &
            poisson%rho,wa_bs(1,1,1,itypat),tt1)
        tt2=trial_energy%E_all(itrial)-trial_energy%energy(itrial)
        fitpar%grad_bv_s1(itypat)=fitpar%grad_bv_s1(itypat)+2.d0*tt1*tt2*pref
        enddo
        do itypat=1,parini%ntypat
        call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,6.d0, &
            poisson%rho,wa_gwp1(1,1,1,itypat),tt1)
        tt2=trial_energy%E_all(itrial)-trial_energy%energy(itrial)
        fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*tt1*tt2*pref
        enddo
        do itypat=1,parini%ntypat
        call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,6.d0, &
            poisson%rho,wa_gwp2(1,1,1,itypat),tt1)
        tt2=trial_energy%E_all(itrial)-trial_energy%energy(itrial)
        fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*tt1*tt2*pref
        enddo
        if(any(fitpar%relaxcore)) then
        do itypat=1,parini%ntypat
        call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,6.d0, &
            poisson%rho,wa_gwc(1,1,1,itypat),tt1)
        tt2=trial_energy%E_all(itrial)-trial_energy%energy(itrial)
        fitpar%grad_gwc_s1(itypat)=fitpar%grad_gwc_s1(itypat)+2.d0*tt1*tt2*pref
        enddo
        do itypat=1,parini%ntypat
        call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,6.d0, &
            poisson%rho,wa_bc(1,1,1,itypat),tt1)
        tt2=trial_energy%E_all(itrial)-trial_energy%energy(itrial)
        fitpar%grad_bc_s1(itypat)=fitpar%grad_bc_s1(itypat)+2.d0*tt1*tt2*pref
        enddo
        endif
    enddo
    call system_clock(itime5)
    deallocate(wa_gws)
    deallocate(wa_bs)
    deallocate(wa_gwp1)
    deallocate(wa_gwp2)
    deallocate(wa_gwc)
    deallocate(wa_bc)
    !-----------------------------------------------------------------
    center(1)=0.d0
    center(2)=0.d0
    center(3)=0.d0
    do iat=1,atoms%nat
        center(1)=center(1)+atoms%ratp(1,iat)
        center(2)=center(2)+atoms%ratp(2,iat)
        center(3)=center(3)+atoms%ratp(3,iat)
    enddo
    center(1)=center(1)/atoms%nat
    center(2)=center(2)/atoms%nat
    center(3)=center(3)/atoms%nat
    dpm(1)=0.d0
    dpm(2)=0.d0
    dpm(3)=0.d0
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        q=(atoms%zat(iat)+fitpar%qcore_type(itypat)+c_s(iat))
        dpm(1)=dpm(1)+q*(atoms%ratp(1,iat)-center(1))
    enddo
    dpm(1)=dpm(1)+sum(c_p(:))
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,i6)',advance='no') 'qp ',istep
    do iat=1,atoms%nat
        write(*,'(f8.3)',advance='no') c_s(iat)
    enddo
    write(*,'(2f8.3)') sum(c_p(:)),dpm(1)
    write(*,'(a,i6)',advance='no') 'ps ',istep
    do iat=1,atoms%nat
        write(*,'(f8.3)',advance='no') c_p(iat)
    enddo
    write(*,*)
    write(*,'(a,i6)',advance='no') 'origin ',istep
    do itypat=1,fitpar%ntypat
        tt1=fitpar%bc_s1(itypat)/(pi**1.5d0*fitpar%gwc_s1(itypat)**3)
        write(*,'(f8.3)',advance='no') tt1
    enddo
    write(*,*)
    endif
    !-----------------------------------------------------------------

    deallocate(pot_ion)
    !if(.not.(istep==parini%paropt_geopt%nit .and. nfiles==1)) then
        deallocate(c_s)
        deallocate(c_p)
        call fini_hartree(parini,atoms,poisson)
    !endif
    if(nfiles==1) then
    do iat=1,atoms_arr%atoms(1)%nat
        c_s_t(iat)=bf%csp(iat+0        )
        c_p_t(iat)=bf%csp(iat+atoms_arr%atoms(1)%nat)
    enddo
    endif
    call bf%fini_bf()
    time_linearcoeff=time_linearcoeff+real(itime2-itime1,8)/real(icount_rate,8)
    time_pot=time_pot+real(itime4-itime3,8)/real(icount_rate,8)
    time_grad=time_grad+real(itime5-itime4,8)/real(icount_rate,8)
    end associate
    end associate
    enddo !end of loop over ifile
    if(parini%mpi_env%nproc>1) then
    call fmpi_allreduce(fitpar%grad_gwv_s1(1),parini%ntypat,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(fitpar%grad_bv_s1(1),parini%ntypat,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(fitpar%grad_gwv_p1(1),parini%ntypat,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(fitpar%grad_gwv_p2(1),parini%ntypat,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(fitpar%grad_gwc_s1(1),parini%ntypat,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(fitpar%grad_bc_s1(1),parini%ntypat,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    endif
    call system_clock(itime6)
    time_net=real(itime6-itime0,8)/real(icount_rate,8)
    if(parini%mpi_env%iproc==0) then
        write(*,'(a,4f8.2)') 'time_linearcoeff,time_pot,time_grad,time_net= ',time_linearcoeff,time_pot,time_grad,time_net
    endif
end subroutine get_cost_from_pot
!*****************************************************************************************
subroutine get_cost_gw_new(parini,fitpar,fitpar_t,cost_gw)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_fitpar), intent(in):: fitpar
    type(typ_fitpar), intent(inout):: fitpar_t
    real(8), intent(out):: cost_gw
    !local variables
    real(8):: dgw, dgw_min, tt, pref
    integer:: itypat
    dgw_min=0.5d0
    pref=0.5d0 !*parini%screening_factor**2
    cost_gw=1.d0
    fitpar_t%grad_gwv_s1=0.d0
    !fitpar_t%grad_gwv_s2=0.d0
    fitpar_t%grad_gwv_p1=0.d0
    fitpar_t%grad_gwv_p2=0.d0
    do itypat=1,parini%ntypat
        !dgw=fitpar%gwv_s1(itypat)-fitpar%gwv_s2(itypat)
        !if(abs(dgw)<dgw_min) then
        !    cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
        !    tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
        !    fitpar_t%grad_gwv_s1(itypat)=fitpar_t%grad_gwv_s1(itypat)-tt
        !    fitpar_t%grad_gwv_s2(itypat)=fitpar_t%grad_gwv_s2(itypat)+tt
        !endif
        !dgw=fitpar%gwv_s2(itypat)-1.1d0
        !if(abs(dgw)<dgw_min) then
        !    cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
        !    tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
        !    fitpar_t%grad_gwv_s2(itypat)=fitpar_t%grad_gwv_s2(itypat)-tt
        !endif
        dgw=fitpar%gwv_p1(itypat)-fitpar%gwv_p2(itypat)
        if(abs(dgw)<dgw_min) then
            cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
            tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
            fitpar_t%grad_gwv_p1(itypat)=fitpar_t%grad_gwv_p1(itypat)-tt
            fitpar_t%grad_gwv_p2(itypat)=fitpar_t%grad_gwv_p2(itypat)+tt
        endif
    enddo
!    dgw_min=0.0d0
!    pref=10.d0
!    do itypat=1,parini%ntypat
!        dgw=fitpar%gwv_s1(itypat)-fitpar%gwv_s2(itypat)
!        !if(abs(dgw)<dgw_min) then
!            cost_gw=cost_gw+pref*(dgw)**4
!            tt=pref*4.d0*(dgw-dgw_min)**3
!            fitpar_t%grad_gwv_s1(itypat)=fitpar_t%grad_gwv_s1(itypat)+tt
!            fitpar_t%grad_gwv_s2(itypat)=fitpar_t%grad_gwv_s2(itypat)-tt
!        !endif
!        dgw=fitpar%gwc_s1(itypat)-fitpar%gwc_s2(itypat)
!        !if(abs(dgw)<dgw_min) then
!            cost_gw=cost_gw+pref*(dgw)**4
!            tt=pref*4.d0*(dgw-dgw_min)**3
!            fitpar_t%grad_gwc_s1(itypat)=fitpar_t%grad_gwc_s1(itypat)+tt
!            fitpar_t%grad_gwc_s2(itypat)=fitpar_t%grad_gwc_s2(itypat)-tt
!        !endif
!    enddo
    dgw_min=0.5d0
    pref=0.2d0
    fitpar_t%grad_gwc_s1=0.d0
!    fitpar_t%grad_gwc_s2=0.d0
!    do itypat=1,parini%ntypat
!        if(fitpar%relaxcore(itypat) .and. fitpar%qcore_type(itypat)<0.d0 .and. fitpar%gwc_s2(itypat)>0.d0) then
!            !dgw=fitpar%gwc_s1(itypat)-fitpar%gwc_s2(itypat)
!            !if(abs(dgw)<dgw_min) then
!            !    cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
!            !    tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
!            !    fitpar_t%grad_gwc_s1(itypat)=fitpar_t%grad_gwc_s1(itypat)-tt
!            !    fitpar_t%grad_gwc_s2(itypat)=fitpar_t%grad_gwc_s2(itypat)+tt
!            !endif
!        endif
!    enddo
    !dgw_min=0.5d0
    !pref=0.1d0
    !do itypat=1,parini%ntypat
    !    if(fitpar%relaxcore(itypat) .and. fitpar%qcore_type(itypat)<0.d0 .and. fitpar%gwc_s2(itypat)>0.d0) then
    !        dgw=fitpar%gwc_s1(itypat)-fitpar%gwv_s2(itypat)
    !        if(abs(dgw)<dgw_min) then
    !            cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
    !            tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
    !            fitpar_t%grad_gwc_s1(itypat)=fitpar_t%grad_gwc_s1(itypat)-tt
    !            fitpar_t%grad_gwv_s2(itypat)=fitpar_t%grad_gwv_s2(itypat)+tt
    !        endif
    !    endif
    !enddo
end subroutine get_cost_gw_new
!*****************************************************************************************
subroutine get_linearcoeff_rho(parini,fitpar,bf,poisson,poisson_ref,atoms,c_s,c_p)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_fitpar), intent(inout):: fitpar
    type(typ_bf), intent(inout):: bf
    type(typ_poisson), intent(in):: poisson, poisson_ref
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: c_s(atoms%nat), c_p(atoms%nat)
    !local variables
    real(8):: voxel, center(3), cons(3)
    real(8):: tt1, tt2, tt3, tt4
    integer:: ix, iy, iz, iat, info, nbf, nc
    integer:: ibf, jbf, lwork
    real(8), allocatable:: secder(:,:), rhs(:)
    real(8), allocatable:: secder_t(:,:), work(:), eval(:)
    integer, allocatable:: ipiv(:)
    nbf=2*atoms%nat
    nc=2
    allocate(secder(nbf+nc,nbf+nc),source=0.d0)
    allocate(rhs(nbf+nc),source=0.d0)
    voxel=poisson_ref%hgrid(1,1)*poisson_ref%hgrid(2,2)*poisson_ref%hgrid(3,3)
    do ibf=1,bf%nbf
        call bf%get_rho_single(parini,poisson,atoms,fitpar,ibf)
        tt1=0.d0
        tt2=0.d0
        do iz=1,poisson_ref%ngpz
        do iy=1,poisson_ref%ngpy
        do ix=1,poisson_ref%ngpx
            tt1=tt1+bf%pot_1(ix,iy,iz)*(poisson%rho(ix,iy,iz)-poisson_ref%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
            tt2=tt2+bf%pot_2(ix,iy,iz)*(poisson%rho(ix,iy,iz)-poisson_ref%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
        enddo
        enddo
        enddo
        bf%subfirstder(1,ibf)=tt1*voxel
        bf%subfirstder(2,ibf)=tt2*voxel
        tt1=0.d0
        tt2=0.d0
        do iz=1,poisson_ref%ngpz
        do iy=1,poisson_ref%ngpy
        do ix=1,poisson_ref%ngpx
            tt1=tt1+bf%vgrad_1(ix,iy,iz)*(poisson%rho(ix,iy,iz)-poisson_ref%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
            tt2=tt2+bf%vgrad_2(ix,iy,iz)*(poisson%rho(ix,iy,iz)-poisson_ref%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
        enddo
        enddo
        enddo
        bf%subfirstder_grad(1,ibf)=tt1*voxel
        bf%subfirstder_grad(2,ibf)=tt2*voxel
        bf%wa1=bf%pot_1
        bf%wa2=bf%pot_2
        bf%wa3=bf%vgrad_1
        bf%wa4=bf%vgrad_2
        do jbf=ibf,bf%nbf
            call bf%get_rho_single(parini,poisson,atoms,fitpar,jbf)
            tt3=0.d0
            do iz=1,poisson_ref%ngpz
            do iy=1,poisson_ref%ngpy
            do ix=1,poisson_ref%ngpx
                tt1=bf%bv_1(ibf)*bf%wa1(ix,iy,iz)+bf%bv_2(ibf)*bf%wa2(ix,iy,iz)
                tt2=bf%bv_1(jbf)*bf%pot_1(ix,iy,iz)+bf%bv_2(jbf)*bf%pot_2(ix,iy,iz)
                tt3=tt3+tt1*tt2*bf%weight(ix,iy,iz)
            enddo
            enddo
            enddo
            secder(ibf,jbf)=tt3*voxel
            secder(jbf,ibf)=tt3*voxel
            tt1=0.d0
            tt2=0.d0
            tt3=0.d0
            tt4=0.d0
            do iz=1,poisson_ref%ngpz
            do iy=1,poisson_ref%ngpy
            do ix=1,poisson_ref%ngpx
                tt1=tt1+bf%pot_1(ix,iy,iz)*bf%wa1(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt2=tt2+bf%pot_2(ix,iy,iz)*bf%wa1(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt3=tt3+bf%pot_1(ix,iy,iz)*bf%wa2(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt4=tt4+bf%pot_2(ix,iy,iz)*bf%wa2(ix,iy,iz)*bf%weight(ix,iy,iz)
            enddo
            enddo
            enddo
            bf%subsecder(1,1,jbf,ibf)=tt1*voxel
            bf%subsecder(2,1,jbf,ibf)=tt2*voxel
            bf%subsecder(1,2,jbf,ibf)=tt3*voxel
            bf%subsecder(2,2,jbf,ibf)=tt4*voxel
            bf%subsecder(1,1,ibf,jbf)=tt1*voxel
            bf%subsecder(1,2,ibf,jbf)=tt2*voxel
            bf%subsecder(2,1,ibf,jbf)=tt3*voxel
            bf%subsecder(2,2,ibf,jbf)=tt4*voxel
            tt1=0.d0
            tt2=0.d0
            tt3=0.d0
            tt4=0.d0
            do iz=1,poisson_ref%ngpz
            do iy=1,poisson_ref%ngpy
            do ix=1,poisson_ref%ngpx
                tt1=tt1+bf%pot_1(ix,iy,iz)*bf%wa3(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt2=tt2+bf%pot_2(ix,iy,iz)*bf%wa3(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt3=tt3+bf%pot_1(ix,iy,iz)*bf%wa4(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt4=tt4+bf%pot_2(ix,iy,iz)*bf%wa4(ix,iy,iz)*bf%weight(ix,iy,iz)
            enddo
            enddo
            enddo
            bf%subsecder_grad(1,1,ibf,jbf)=tt1*voxel
            bf%subsecder_grad(1,2,ibf,jbf)=tt2*voxel
            bf%subsecder_grad(2,1,ibf,jbf)=tt3*voxel
            bf%subsecder_grad(2,2,ibf,jbf)=tt4*voxel
            tt1=0.d0
            tt2=0.d0
            tt3=0.d0
            tt4=0.d0
            do iz=1,poisson_ref%ngpz
            do iy=1,poisson_ref%ngpy
            do ix=1,poisson_ref%ngpx
                tt1=tt1+bf%vgrad_1(ix,iy,iz)*bf%wa1(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt2=tt2+bf%vgrad_2(ix,iy,iz)*bf%wa1(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt3=tt3+bf%vgrad_1(ix,iy,iz)*bf%wa2(ix,iy,iz)*bf%weight(ix,iy,iz)
                tt4=tt4+bf%vgrad_2(ix,iy,iz)*bf%wa2(ix,iy,iz)*bf%weight(ix,iy,iz)
            enddo
            enddo
            enddo
            bf%subsecder_grad(1,1,jbf,ibf)=tt1*voxel
            bf%subsecder_grad(2,1,jbf,ibf)=tt2*voxel
            bf%subsecder_grad(1,2,jbf,ibf)=tt3*voxel
            bf%subsecder_grad(2,2,jbf,ibf)=tt4*voxel
        enddo
        !-------------------------------------------------------------
        tt1=0.d0
        do iz=1,poisson_ref%ngpz
        do iy=1,poisson_ref%ngpy
        do ix=1,poisson_ref%ngpx
            tt2=bf%bv_1(ibf)*bf%wa1(ix,iy,iz)+bf%bv_2(ibf)*bf%wa2(ix,iy,iz)
            tt1=tt1+tt2*(poisson%rho(ix,iy,iz)-poisson_ref%rho(ix,iy,iz))*bf%weight(ix,iy,iz)
        enddo
        enddo
        enddo
        rhs(ibf)=-tt1*voxel
        !-------------------------------------------------------------
    enddo
    bf%secder(1:bf%nbf,1:bf%nbf)=secder(1:bf%nbf,1:bf%nbf)
    bf%firstder(1:bf%nbf)=-rhs(1:bf%nbf)
    tt1=0.d0
    do iz=1,poisson_ref%ngpz
    do iy=1,poisson_ref%ngpy
    do ix=1,poisson_ref%ngpx
        tt1=tt1+(poisson%rho(ix,iy,iz)-poisson_ref%rho(ix,iy,iz))**2*bf%weight(ix,iy,iz)
    enddo
    enddo
    enddo
    bf%zeroder=tt1*voxel

    do iat=1,atoms%nat
        write(*,'(a,3es14.5)') 'SECDER ',secder(iat,1:3)
    enddo
    do iat=1,atoms%nat
        write(*,'(a,3es14.5)') 'SECDER ',secder(atoms%nat+iat,1:3)
    enddo
    !-----------------------------------------------------------------
    !lwork=max(nbf**2,1000)
    !allocate(secder_t(nbf,nbf))
    !allocate(eval(nbf),work(lwork))
    !secder_t(1:nbf,1:nbf)=secder(1:nbf,1:nbf)
    !call DSYEV('V','U',nbf,secder_t,nbf,eval,work,lwork,info)
    !do ibf=1,nbf
    !    write(*,'(a,i6,es14.5)') 'EVAL ',ibf,eval(ibf)
    !enddo
    !do ibf=1,nbf
    !    write(*,'(a,i6,4es14.5)') 'EVEC ',ibf,secder_t(1:nbf,ibf)
    !enddo
    !-----------------------------------------------------------------
    center(1)=0.d0
    center(2)=0.d0
    center(3)=0.d0
    do iat=1,atoms%nat
        center(1)=center(1)+atoms%ratp(1,iat)
        center(2)=center(2)+atoms%ratp(2,iat)
        center(3)=center(3)+atoms%ratp(3,iat)
    enddo
    center(1)=center(1)/atoms%nat
    center(2)=center(2)/atoms%nat
    center(3)=center(3)/atoms%nat
    secder(nbf+1      ,nbf+1      )=0.d0
    secder(nbf+2      ,nbf+2      )=0.d0
    do iat=1,atoms%nat
        secder(iat          ,nbf+1        )=1.d0
        secder(nbf+1        ,iat          )=1.d0
        secder(iat          ,nbf+2        )=atoms%ratp(1,iat)-center(1)
        secder(nbf+2        ,iat          )=atoms%ratp(1,iat)-center(1)
        secder(atoms%nat+iat,nbf+2        )=1.d0
        secder(nbf+2        ,atoms%nat+iat)=1.d0
    enddo
    rhs(nbf+1)=atoms%qtot-sum(atoms%zat)-sum(bf%qcore) !-sum(atoms%zat)
    cons(1)=0.d0
    cons(2)=0.d0
    cons(3)=0.d0
    do iat=1,atoms%nat
        cons(1)=cons(1)+(atoms%zat(iat)+bf%qcore(iat))*(atoms%ratp(1,iat)-center(1))
        cons(2)=cons(2)+(atoms%zat(iat)+bf%qcore(iat))*(atoms%ratp(2,iat)-center(2))
        cons(3)=cons(3)+(atoms%zat(iat)+bf%qcore(iat))*(atoms%ratp(3,iat)-center(3))
    enddo
    rhs(nbf+2)=atoms%dpm(1)-cons(1)
    allocate(ipiv(nbf+nc))
    call DGETRF(nbf+nc,nbf+nc,secder,nbf+nc,ipiv,info)
    call DGETRS('N',nbf+nc,1,secder,nbf+nc,ipiv,rhs,nbf+nc,info)
    do iat=1,atoms%nat
        c_s(iat)=rhs(iat+0        )
        c_p(iat)=rhs(iat+atoms%nat)
        write(*,'(2f18.8)') c_s(iat),c_p(iat)
    enddo
    bf%csp(1:bf%nbf)=rhs(1:bf%nbf)
    !c_s(1)=-2.d0
    !c_s(2)=-2.d0
    !c_s(3)=-2.d0
    deallocate(secder,rhs,ipiv)
end subroutine get_linearcoeff_rho
!*****************************************************************************************
subroutine get_linearcoeff(parini,fitpar,bf,poisson,poisson_ref,atoms,trial_energy,c_s,c_p)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_cent2, only: get_cost_secder
    use mod_atoms, only: typ_atoms
    use mod_trial_energy, only: typ_trial_energy
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_fitpar), intent(inout):: fitpar
    type(typ_bf), intent(inout):: bf
    type(typ_poisson), intent(inout):: poisson
    type(typ_poisson), intent(in):: poisson_ref
    type(typ_atoms), intent(in):: atoms
    type(typ_trial_energy), intent(inout):: trial_energy
    real(8), intent(out):: c_s(atoms%nat), c_p(atoms%nat)
    !local variables
    real(8):: voxel, center(3), cons(3)
    real(8):: tt1, tt2, tt3, tt4
    integer:: ix, iy, iz, iat, info, nbf, nc, itypat
    integer:: ibf, jbf, lwork
    real(8), allocatable:: secder(:,:), rhs(:)
    real(8), allocatable:: squarefit_raw(:,:)
    real(8), allocatable:: rhs_raw(:)
    real(8), allocatable:: secder_t(:,:), work(:), eval(:)
    real(8), allocatable:: be_s(:), gwe_s(:), gwe_p(:,:)
    real(8), allocatable:: qcore(:), bc(:), gwc(:)
    integer, allocatable:: ipiv(:)
    nbf=2*atoms%nat
    nc=2
    allocate(secder(nbf+nc,nbf+nc),source=0.d0)
    allocate(rhs(nbf+nc),source=0.d0)
    allocate(squarefit_raw(bf%nbf,bf%nbf),rhs_raw(bf%nbf))
    allocate(be_s(atoms%nat),gwe_s(atoms%nat),gwe_p(2,atoms%nat))
    allocate(qcore(atoms%nat),bc(atoms%nat),gwc(atoms%nat))
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        be_s(iat)=fitpar%bv_s1(itypat)
        gwe_s(iat)=fitpar%gwv_s1(itypat)
        gwe_p(1,iat)=fitpar%gwv_p1(itypat)
        gwe_p(2,iat)=fitpar%gwv_p2(itypat)
        qcore(iat)=fitpar%qcore_type(itypat)
        bc(iat)=fitpar%bc_s1(itypat)
        gwc(iat)=fitpar%gwc_s1(itypat)
    enddo
    allocate(poisson%pot_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    poisson%pot_ion=poisson%pot
    call get_cost_secder(parini,trial_energy,poisson,atoms,bf%nbf,bf%iat_list,bf%bt,be_s,gwe_s,gwe_p,qcore,bc,gwc,squarefit_raw,rhs_raw)
    deallocate(poisson%pot_ion)
    secder(1:bf%nbf,1:bf%nbf)=squarefit_raw(1:bf%nbf,1:bf%nbf)
    rhs(1:bf%nbf)=rhs_raw(1:bf%nbf)
    !write(*,*) rhs_raw
    bf%secder(1:bf%nbf,1:bf%nbf)=secder(1:bf%nbf,1:bf%nbf)
    bf%firstder(1:bf%nbf)=-rhs(1:bf%nbf)
    tt1=0.d0
    do iz=1,poisson_ref%ngpz
    do iy=1,poisson_ref%ngpy
    do ix=1,poisson_ref%ngpx
        tt1=tt1+(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))**2
    enddo
    enddo
    enddo
    bf%zeroder=tt1*voxel

    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(*,'(a,3es14.5)') 'SECDER ',secder(iat,1:3)
    enddo
    do iat=1,atoms%nat
        write(*,'(a,3es14.5)') 'SECDER ',secder(atoms%nat+iat,1:3)
    enddo
    endif
    !-----------------------------------------------------------------
    !lwork=max(nbf**2,1000)
    !allocate(secder_t(nbf,nbf))
    !allocate(eval(nbf),work(lwork))
    !secder_t(1:nbf,1:nbf)=secder(1:nbf,1:nbf)
    !call DSYEV('V','U',nbf,secder_t,nbf,eval,work,lwork,info)
    !do ibf=1,nbf
    !    write(*,'(a,i6,es14.5)') 'EVAL ',ibf,eval(ibf)
    !enddo
    !do ibf=1,nbf
    !    write(*,'(a,i6,4es14.5)') 'EVEC ',ibf,secder_t(1:nbf,ibf)
    !enddo
    !-----------------------------------------------------------------
    center(1)=0.d0
    center(2)=0.d0
    center(3)=0.d0
    do iat=1,atoms%nat
        center(1)=center(1)+atoms%ratp(1,iat)
        center(2)=center(2)+atoms%ratp(2,iat)
        center(3)=center(3)+atoms%ratp(3,iat)
    enddo
    center(1)=center(1)/atoms%nat
    center(2)=center(2)/atoms%nat
    center(3)=center(3)/atoms%nat
    secder(nbf+1      ,nbf+1      )=0.d0
    secder(nbf+2      ,nbf+2      )=0.d0
    do iat=1,atoms%nat
        secder(iat          ,nbf+1        )=1.d0
        secder(nbf+1        ,iat          )=1.d0
        secder(iat          ,nbf+2        )=atoms%ratp(1,iat)-center(1)
        secder(nbf+2        ,iat          )=atoms%ratp(1,iat)-center(1)
        secder(atoms%nat+iat,nbf+2        )=1.d0
        secder(nbf+2        ,atoms%nat+iat)=1.d0
    enddo
    rhs(nbf+1)=atoms%qtot-sum(atoms%zat)-sum(bf%qcore) !-sum(atoms%zat)
    cons(1)=0.d0
    cons(2)=0.d0
    cons(3)=0.d0
    do iat=1,atoms%nat
        cons(1)=cons(1)+(atoms%zat(iat)+bf%qcore(iat))*(atoms%ratp(1,iat)-center(1))
        cons(2)=cons(2)+(atoms%zat(iat)+bf%qcore(iat))*(atoms%ratp(2,iat)-center(2))
        cons(3)=cons(3)+(atoms%zat(iat)+bf%qcore(iat))*(atoms%ratp(3,iat)-center(3))
    enddo
    rhs(nbf+2)=atoms%dpm(1)-cons(1)
    allocate(ipiv(nbf+nc))
    call DGETRF(nbf+nc,nbf+nc,secder,nbf+nc,ipiv,info)
    call DGETRS('N',nbf+nc,1,secder,nbf+nc,ipiv,rhs,nbf+nc,info)
    do iat=1,atoms%nat
        c_s(iat)=rhs(iat+0        )
        c_p(iat)=rhs(iat+atoms%nat)
    enddo
    bf%csp(1:bf%nbf)=rhs(1:bf%nbf)
    if(parini%mpi_env%iproc==0) then
        do iat=1,atoms%nat
            write(*,'(2f18.8)') c_s(iat),c_p(iat)
        enddo
    endif
    deallocate(secder,rhs,ipiv)
end subroutine get_linearcoeff
!*****************************************************************************************
subroutine cal_pot_gauss_p(parini,poisson,pot,reset,xyz,gw,p,vgrad)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    real(8), intent(inout):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    logical, intent(in):: reset
    real(8), intent(in):: xyz(3), gw, p(3)
    real(8), intent(inout):: vgrad(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    !type(typ_atoms):: atoms
    real(8):: pi, sf, ss0, ss1, ss2, ss3, ss4
    real(8):: x, y, z, r, ddd, p_dot_r
    real(8):: tt0, tt1, tt2, tt3, tt4
    real(8):: uu1, vv1, qq1, sft
    real(8):: yy1, ff1, hh1
    real(8):: sf_1, sf_2, sf_3, a, b, c
    real(8):: tg0, tg1, tg2, tg3, tg4
    real(8):: xyz111(3), hgrid(3,3)
    integer:: itime1, itime2, icount_rate
    integer:: ix, iy, iz, nx, ny, nz
    xyz111=poisson%xyz111
    hgrid=poisson%hgrid
    nx=poisson%ngpx
    ny=poisson%ngpy
    nz=poisson%ngpz
    pi=4.d0*atan(1.d0)
    if(reset) then
        pot=0.d0
    endif
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    !a=-sf_2/(sf_1-sf_2)
    !b= sf_1/(sf_1-sf_2)
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))
    !call system_clock(itime1)
    !$omp parallel default(shared)
    !$omp do collapse(2) schedule(static,16) &
    !$omp private(ix,iy,iz,p_dot_r,x,y,z,r,ss1,ss2,ss3,tg1,tg2,tg3,tt1,tt2,tt3,ddd) firstprivate(xyz,xyz111,hgrid,p,sf_1,sf_2,sf_3,a,b,c,gw,pi)
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
        x=xyz111(1)+(ix-1)*hgrid(1,1)-xyz(1)
        y=xyz111(2)+(iy-1)*hgrid(2,2)-xyz(2)
        z=xyz111(3)+(iz-1)*hgrid(3,3)-xyz(3)
        r=sqrt(x**2+y**2+z**2)
        p_dot_r=x*p(1)+y*p(2)+z*p(3)
        if(r<1.d-10) then
            ss1=0.d0
            ss2=0.d0
            ss3=0.d0
            tg1=0.d0
            tg2=0.d0
            tg3=0.d0
        else
            tt1=sf_1/sqrt(1.d0+(sf_1*gw)**2)
            tt2=sf_2/sqrt(1.d0+(sf_2*gw)**2)
            tt3=sf_3/sqrt(1.d0+(sf_3*gw)**2)
            if((r*tt1)**2>50.d0) then
                ss1=erf(tt1*r)/r**3
                tg1=0.d0
            else
                ddd=exp(-(r*tt1)**2)
                ss1=erf(tt1*r)/r**3-2.d0*sf_1/(r**2*sqrt(pi*(1.d0+(sf_1*gw)**2)))*ddd
                tg1=-4.d0*sf_1**5*gw/sqrt(pi*(1.d0+(sf_1*gw)**2)**5)*ddd
            endif
            if((r*tt2)**2>50.d0) then
                ss2=erf(tt2*r)/r**3
                tg2=0.d0
            else
                ddd=exp(-(r*tt2)**2)
                ss2=erf(tt2*r)/r**3-2.d0*sf_2/(r**2*sqrt(pi*(1.d0+(sf_2*gw)**2)))*ddd
                tg2=-4.d0*sf_2**5*gw/sqrt(pi*(1.d0+(sf_2*gw)**2)**5)*ddd
            endif
            if((r*tt3)**2>50.d0) then
                ss3=erf(tt3*r)/r**3
                tg3=0.d0
            else
                ddd=exp(-(r*tt3)**2)
                ss3=erf(tt3*r)/r**3-2.d0*sf_3/(r**2*sqrt(pi*(1.d0+(sf_3*gw)**2)))*ddd
                tg3=-4.d0*sf_3**5*gw/sqrt(pi*(1.d0+(sf_3*gw)**2)**5)*ddd
            endif
        endif
        pot(ix,iy,iz)=pot(ix,iy,iz)+p_dot_r*(a*ss1+b*ss2+c*ss3)
        vgrad(ix,iy,iz)=p_dot_r*(a*tg1+b*tg2+c*tg3)
    enddo
    enddo
    enddo
    !$omp end do
    !$omp end parallel
    !call system_clock(itime2)
    !call system_clock(count_rate=icount_rate)
    !time_cal_pot_gauss_p=time_cal_pot_gauss_p+real(itime2-itime1,8)/real(icount_rate,8)
end subroutine cal_pot_gauss_p
!*****************************************************************************************
subroutine cal_pot_gauss_s(parini,poisson,pot,reset,xyz,gw,q,vgrad)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    real(8), intent(inout):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    logical, intent(in):: reset
    real(8), intent(in):: xyz(3), gw, q
    real(8), intent(inout):: vgrad(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    !type(typ_atoms):: atoms
    real(8):: pi, sf_1, sf_2, sf_3, ss0, ss1, ss2, ss3, ss4
    real(8):: x, y, z, r, tt0, tt1, tt2, tt3, tt4, a, b, c, uu1, ww1, qq1, vv1, sft
    real(8):: tg0, tg1, tg2, tg3, tg4
    real(8):: dd0, dd1, dd2, dd3, dd4
    real(8):: xyz111(3), hgrid(3,3)
    integer:: itime1, itime2, icount_rate
    integer:: ix, iy, iz, nx, ny, nz
    xyz111=poisson%xyz111
    hgrid=poisson%hgrid
    nx=poisson%ngpx
    ny=poisson%ngpy
    nz=poisson%ngpz
    pi=4.d0*atan(1.d0)
    if(reset) then
        pot=0.d0
    endif
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    !a=-sf_2/(sf_1-sf_2)
    !b= sf_1/(sf_1-sf_2)
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))

    !call system_clock(itime1)
    !$omp parallel default(shared)
    !$omp do collapse(2) schedule(static,16) &
    !$omp private(ix,iy,iz,x,y,z,r,ss1,ss2,ss3,tg1,tg2,tg3,tt1,tt2,tt3) firstprivate(xyz,xyz111,hgrid,q,sf_1,sf_2,sf_3,a,b,c,gw,pi)
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
        x=xyz111(1)+(ix-1)*hgrid(1,1)-xyz(1)
        y=xyz111(2)+(iy-1)*hgrid(2,2)-xyz(2)
        z=xyz111(3)+(iz-1)*hgrid(3,3)-xyz(3)
        r=sqrt(x**2+y**2+z**2)
        if(r<1.d-10) then
            ss1=2.d0*sf_1/sqrt(pi*(1.d0+(sf_1*gw)**2))
            ss2=2.d0*sf_2/sqrt(pi*(1.d0+(sf_2*gw)**2))
            ss3=2.d0*sf_3/sqrt(pi*(1.d0+(sf_3*gw)**2))
            tg1=-2.d0*sf_1**3*gw/sqrt(pi*(1.d0+(sf_1*gw)**2)**3)
            tg2=-2.d0*sf_2**3*gw/sqrt(pi*(1.d0+(sf_2*gw)**2)**3)
            tg3=-2.d0*sf_3**3*gw/sqrt(pi*(1.d0+(sf_3*gw)**2)**3)
        else
            tt1=sf_1/sqrt(1.d0+(sf_1*gw)**2)
            tt2=sf_2/sqrt(1.d0+(sf_2*gw)**2)
            tt3=sf_3/sqrt(1.d0+(sf_3*gw)**2)
            ss1=erf(tt1*r)/r
            ss2=erf(tt2*r)/r
            ss3=erf(tt3*r)/r
            if((r*tt1)**2>50.d0) then
                tg1=0.d0
            else
                tg1=-2.d0*sf_1**3*gw/sqrt(pi*(1.d0+(sf_1*gw)**2)**3)*exp(-(r*tt1)**2)
            endif
            if((r*tt2)**2>50.d0) then
                tg2=0.d0
            else
                tg2=-2.d0*sf_2**3*gw/sqrt(pi*(1.d0+(sf_2*gw)**2)**3)*exp(-(r*tt2)**2)
            endif
            if((r*tt3)**2>50.d0) then
                tg3=0.d0
            else
                tg3=-2.d0*sf_3**3*gw/sqrt(pi*(1.d0+(sf_3*gw)**2)**3)*exp(-(r*tt3)**2)
            endif
        endif
        pot(ix,iy,iz)=pot(ix,iy,iz)+q*(a*ss1+b*ss2+c*ss3)
        vgrad(ix,iy,iz)=q*(a*tg1+b*tg2+c*tg3)
    enddo
    enddo
    enddo
    !$omp end do
    !$omp end parallel
    !call system_clock(itime2)
    !call system_clock(count_rate=icount_rate)
    !time_cal_pot_gauss_s=time_cal_pot_gauss_s+real(itime2-itime1,8)/real(icount_rate,8)
end subroutine cal_pot_gauss_s
!*****************************************************************************************
subroutine cal_pot_r2gauss_s(parini,poisson,pot,reset,xyz,gw,q,vgrad)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(in):: poisson
    real(8), intent(inout):: pot(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    logical, intent(in):: reset
    real(8), intent(in):: xyz(3), gw, q
    real(8), intent(inout):: vgrad(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    !type(typ_atoms):: atoms
    real(8):: pi, sf, ss0, ss1, ss2, ss3, ss4
    real(8):: x, y, z, r, uu1, ww1, qq1, vv1, sft, ddd
    real(8):: tt0, tt1, tt2, tt3, tt4
    real(8):: sf_1, sf_2, sf_3, a, b, c
    real(8):: tg0, tg1, tg2, tg3, tg4, tg5
    real(8):: sg0, sg1, sg2, sg3, sg4
    real(8):: ee0, ee1
    real(8):: xyz111(3), hgrid(3,3)
    integer:: itime1, itime2, icount_rate
    integer:: ix, iy, iz, nx, ny, nz
    xyz111=poisson%xyz111
    hgrid=poisson%hgrid
    nx=poisson%ngpx
    ny=poisson%ngpy
    nz=poisson%ngpz
    pi=4.d0*atan(1.d0)
    if(reset) then
        pot=0.d0
    endif
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    !a=-sf_2/(sf_1-sf_2)
    !b= sf_1/(sf_1-sf_2)
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))

    !call system_clock(itime1)
    !$omp parallel default(shared)
    !$omp do collapse(2) schedule(static,16) &
    !$omp private(ix,iy,iz,x,y,z,r,ss1,ss2,ss3,tg1,tg2,tg3,tt1,tt2,tt3,ddd) firstprivate(xyz,xyz111,hgrid,q,sf_1,sf_2,sf_3,a,b,c,gw,pi)
    do iz=1,nz
    do iy=1,ny
    do ix=1,nx
        x=xyz111(1)+(ix-1)*hgrid(1,1)-xyz(1)
        y=xyz111(2)+(iy-1)*hgrid(2,2)-xyz(2)
        z=xyz111(3)+(iz-1)*hgrid(3,3)-xyz(3)
        r=sqrt(x**2+y**2+z**2)
        if(r<1.d-10) then
            ss1=2.d0*sf_1*(3.d0+2.d0*(sf_1*gw)**2)/(3.d0*sqrt(pi*(1.d0+(sf_1*gw)**2)**3))
            ss2=2.d0*sf_2*(3.d0+2.d0*(sf_2*gw)**2)/(3.d0*sqrt(pi*(1.d0+(sf_2*gw)**2)**3))
            ss3=2.d0*sf_3*(3.d0+2.d0*(sf_3*gw)**2)/(3.d0*sqrt(pi*(1.d0+(sf_3*gw)**2)**3))
            tg1=-2.d0*sf_1**3*gw*(5.d0+7.d0*(sf_1*gw)**2+2.d0*(sf_1*gw)**4)/(3.d0*sqrt(pi*(1.d0+(sf_1*gw)**2)**7))
            tg2=-2.d0*sf_2**3*gw*(5.d0+7.d0*(sf_2*gw)**2+2.d0*(sf_2*gw)**4)/(3.d0*sqrt(pi*(1.d0+(sf_2*gw)**2)**7))
            tg3=-2.d0*sf_3**3*gw*(5.d0+7.d0*(sf_3*gw)**2+2.d0*(sf_3*gw)**4)/(3.d0*sqrt(pi*(1.d0+(sf_3*gw)**2)**7))
        else
            tt1=sf_1/sqrt(1.d0+(sf_1*gw)**2)
            tt2=sf_2/sqrt(1.d0+(sf_2*gw)**2)
            tt3=sf_3/sqrt(1.d0+(sf_3*gw)**2)
            if((r*tt1)**2>50.d0) then
                ss1=erf(tt1*r)/r
                tg1=0.d0
            else
                ddd=exp(-(r*tt1)**2)
                ss1=erf(tt1*r)/r-2.d0*sf_1**3*gw**2/(3.d0*sqrt(pi*(1.d0+(sf_1*gw)**2)**3))*ddd
                tg1=-2.d0*sf_1**3*gw*(5.d0+2.d0*(sf_1*gw)**4+(sf_1*gw)**2*(7.d0+2.d0*(sf_1*r)**2))/(3.d0*sqrt(pi*(1.d0+(sf_1*gw)**2)**7))*ddd
            endif
            if((r*tt2)**2>50.d0) then
                ss2=erf(tt2*r)/r
                tg2=0.d0
            else
                ddd=exp(-(r*tt2)**2)
                ss2=erf(tt2*r)/r-2.d0*sf_2**3*gw**2/(3.d0*sqrt(pi*(1.d0+(sf_2*gw)**2)**3))*ddd
                tg2=-2.d0*sf_2**3*gw*(5.d0+2.d0*(sf_2*gw)**4+(sf_2*gw)**2*(7.d0+2.d0*(sf_2*r)**2))/(3.d0*sqrt(pi*(1.d0+(sf_2*gw)**2)**7))*ddd
            endif
            if((r*tt3)**2>50.d0) then
                ss3=erf(tt3*r)/r
                tg3=0.d0
            else
                ddd=exp(-(r*tt3)**2)
                ss3=erf(tt3*r)/r-2.d0*sf_3**3*gw**2/(3.d0*sqrt(pi*(1.d0+(sf_3*gw)**2)**3))*ddd
                tg3=-2.d0*sf_3**3*gw*(5.d0+2.d0*(sf_3*gw)**4+(sf_3*gw)**2*(7.d0+2.d0*(sf_3*r)**2))/(3.d0*sqrt(pi*(1.d0+(sf_3*gw)**2)**7))*ddd
            endif
        endif
        pot(ix,iy,iz)=pot(ix,iy,iz)+q*(a*ss1+b*ss2+c*ss3)
        vgrad(ix,iy,iz)=q*(a*tg1+b*tg2+c*tg3)
    enddo
    enddo
    enddo
    !$omp end do
    !$omp end parallel
    !call system_clock(itime2)
    !call system_clock(count_rate=icount_rate)
    !time_cal_pot_r2gauss_s=time_cal_pot_r2gauss_s+real(itime2-itime1,8)/real(icount_rate,8)
end subroutine cal_pot_r2gauss_s
!*****************************************************************************************
subroutine get_dpm_fit_bf_cent2(atoms,dpx,dpy,dpz,dpm_err)
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: dpm_err 
    real(8), intent(out):: dpx, dpy, dpz
    !local variables
    integer:: iat
    real(8):: centroid_x, centroid_y, centroid_z
    dpx=0.d0
    dpy=0.d0
    dpz=0.d0
    centroid_x=sum(atoms%ratp(1,:))/atoms%nat
    centroid_y=sum(atoms%ratp(2,:))/atoms%nat
    centroid_z=sum(atoms%ratp(3,:))/atoms%nat
    write(*,'(a,5f8.2)') 'DPM ',atoms%qat(1),atoms%qat(2),atoms%ratp(1,2)-atoms%ratp(1,1),atoms%zat(1),atoms%zat(2)
    do iat=1, atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
            dpx=dpx+(atoms%ratp(1,iat)-centroid_x)*(atoms%qat(iat)+2.d0)
            dpy=dpy+(atoms%ratp(2,iat)-centroid_y)*(atoms%qat(iat)+2.d0)
            dpz=dpz+(atoms%ratp(3,iat)-centroid_z)*(atoms%qat(iat)+2.d0)
        elseif(trim(atoms%sat(iat))=='O') then
            dpx=dpx+(atoms%ratp(1,iat)-centroid_x)*(atoms%qat(iat)+4.d0)
            dpy=dpy+(atoms%ratp(2,iat)-centroid_y)*(atoms%qat(iat)+4.d0)
            dpz=dpz+(atoms%ratp(3,iat)-centroid_z)*(atoms%qat(iat)+4.d0)
        endif
    enddo
    dpm_err=((dpx-atoms%dpm(1))**2+(dpy-atoms%dpm(2))**2+(dpz-atoms%dpm(3))**2)
end subroutine get_dpm_fit_bf_cent2
!*****************************************************************************************
subroutine get_poisson_ref(parini,gausswidth_ion,poisson,atoms,fitpar)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, update_rat, atom_deallocate_old
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    real(8), intent(in):: gausswidth_ion
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_fitpar), intent(inout):: fitpar
    !local variables
    type(typ_poisson):: poisson_ion
    integer:: igpx, igpy, igpz, iat
    integer:: jgpx, jgpy, jgpz, itypat
    real(8):: rgcut_a, qtot, pi
    real(8):: ehartree_scn_excl, tt1, center(3), q_tmp(1)
    real(8):: dx, dy, dz, r2, coeff, gwt, x, y, z
    !real(8):: xmin, ymin, zmin, xmax, ymax, zmax
    real(8), allocatable::  gausswidth(:)
    real(8), allocatable::  rho(:,:,:)
    real(8), allocatable::  pot(:,:,:)
    integer:: nbgpx, nbgpy, nbgpz, nex, ney, nez, nd
    pi=4.d0*atan(1.d0)
    if(parini%ewald) then
        write(*,*) 'ERROR: ewald=True is wrong when reading from cube file.'
        stop
    endif
    if(trim(parini%psolver)=='kwald') then
        write(*,*) 'ERROR: psolver=kwald is wrong for grid base charge density.'
        stop
    endif
    if(parini%cal_scn .and. parini%screening_factor==0.d0) then
        stop 'ERROR: cal_scn is TRUE and screening_factor is 0, MEANINGLESS!'
    endif
    !-------------------------------------------------------
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Li') atoms%zat(iat)=3.d0
        if(trim(atoms%sat(iat))=='S' ) atoms%zat(iat)=6.d0
    enddo
    !-------------------------------------------------------
    if(.true.) then
    nd=2
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call f_free(poisson%rho)
    poisson%ngpz=poisson%ngpz/nd
    poisson%ngpy=poisson%ngpy/nd
    poisson%ngpx=poisson%ngpx/nd
    poisson%rho=f_malloc0([1.to.(poisson%ngpx),1.to.(poisson%ngpy),1.to.(poisson%ngpz)],id='poisson%rho')
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx,igpy,igpz)=rho(nd*igpx-nd+1,nd*igpy-nd+1,nd*igpz-nd+1)
    enddo
    enddo
    enddo
    deallocate(rho)
    poisson%hgrid=poisson%hgrid*real(nd,kind=8)
    endif
    !-------------------------------------------------------
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    !write(*,*) allocated(poisson%rho)
    call f_free(poisson%rho)
    nex=20
    ney=20
    nez=20
    poisson%rho=f_malloc0([1.to.(poisson%ngpx+2*nex),1.to.(poisson%ngpy+2*ney),1.to.(poisson%ngpz+2*nez)],id='poisson%rho')
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx+nex,igpy+ney,igpz+nez)=rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call update_ratp(atoms)
    do iat=1,atoms%nat
        atoms%ratp(1,iat)=atoms%ratp(1,iat)+real(nex,kind=8)*poisson%hgrid(1,1)
        atoms%ratp(2,iat)=atoms%ratp(2,iat)+real(ney,kind=8)*poisson%hgrid(2,2)
        atoms%ratp(3,iat)=atoms%ratp(3,iat)+real(nez,kind=8)*poisson%hgrid(3,3)
    enddo
    call update_rat(atoms)
    deallocate(rho)
    poisson%ngpx=poisson%ngpx+2*nex
    poisson%ngpy=poisson%ngpy+2*ney
    poisson%ngpz=poisson%ngpz+2*nez
    !-------------------------------------------------------
    atoms%boundcond='free'
    poisson%bc=atoms%boundcond
    allocate(gausswidth(atoms%nat))
    gausswidth=gausswidth_ion !parini%gaussian_width !TO_BE_CORRECTED
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    poisson%task_finit=""
    call init_hartree(parini,atoms,poisson,gausswidth)
    poisson%cell(1)=poisson%hgrid(1,1)*poisson%ngpx
    poisson%cell(2)=poisson%hgrid(2,2)*poisson%ngpy
    poisson%cell(3)=poisson%hgrid(3,3)*poisson%ngpz
    !-------------------------------------------------------
    poisson_ion%alpha=gausswidth_ion !parini%gaussian_width !TO_BE_CORRECTED
    poisson_ion%rgcut=parini%rgcut_ewald*poisson_ion%alpha
    poisson_ion%ngpx=poisson%ngpx
    poisson_ion%ngpy=poisson%ngpy
    poisson_ion%ngpz=poisson%ngpz
    poisson_ion%hgrid(1:3,1:3)=0.d0
    poisson_ion%hgrid(1,1)=poisson%hgrid(1,1)
    poisson_ion%hgrid(2,2)=poisson%hgrid(2,2)
    poisson_ion%hgrid(3,3)=poisson%hgrid(3,3)
    poisson_ion%xyz111=poisson%xyz111
    rgcut_a=parini%rgcut_ewald*maxval(gausswidth) !parini%gaussian_width !3.d0
    nbgpx=int(rgcut_a/poisson_ion%hgrid(1,1))+2
    nbgpy=int(rgcut_a/poisson_ion%hgrid(2,2))+2
    nbgpz=int(rgcut_a/poisson_ion%hgrid(3,3))+2
    poisson_ion%task_finit="alloc_rho"
    call init_hartree(parini,atoms,poisson_ion,gausswidth)
    poisson_ion%reset_rho=.true.
    poisson_ion%nat=atoms%nat
    poisson_ion%cv=atoms%cellvec
    poisson_ion%bc=atoms%boundcond
    poisson_ion%q(1:poisson_ion%nat)=atoms%zat(1:atoms%nat)
    poisson_ion%gw(1:poisson_ion%nat)=gausswidth(1:atoms%nat)
    call get_rat(atoms,poisson_ion%rcart)
    !call put_charge_density(parini,poisson_ion)
    nbgpx=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(1,1))+2
    nbgpy=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(2,2))+2
    nbgpz=int(7.d0*maxval(gausswidth(1:atoms%nat))/poisson_ion%hgrid(3,3))+2

    poisson_ion%rho=0.d0
    do iat=1,atoms%nat
    itypat=atoms%itypat(iat)
    q_tmp(1)=atoms%zat(iat)*fitpar%bz_s(itypat)
    call put_gto_sym_ortho(parini,poisson_ion%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwz_s(itypat), &
        rgcut*fitpar%gwz_s(itypat),poisson_ion%xyz111,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
        poisson_ion%hgrid,poisson_ion%rho)
    q_tmp(1)=atoms%zat(iat)*(1.d0-fitpar%bz_s(itypat))
    call put_r2gto_sym_ortho(parini,poisson_ion%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwz_s(itypat), &
        rgcut*fitpar%gwz_s(itypat),poisson_ion%xyz111,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
        poisson_ion%hgrid,poisson_ion%rho)
    enddo

    !poisson_ion%rho=0.d0
    !do iat=1,atoms%nat
    !gwt=gausswidth(iat)
    !coeff=atoms%zat(iat)/(gwt**3*pi**1.5d0)
    !jgpx=int(poisson_ion%rcart(1,iat)/poisson%hgrid(1,1))
    !jgpy=int(poisson_ion%rcart(2,iat)/poisson%hgrid(2,2))
    !jgpz=int(poisson_ion%rcart(3,iat)/poisson%hgrid(3,3))
    !do igpz=jgpz-nbgpz,jgpz+nbgpz
    !do igpy=jgpy-nbgpy,jgpy+nbgpy
    !do igpx=jgpx-nbgpx,jgpx+nbgpx
    !    dx=(igpx-1)*poisson%hgrid(1,1)-poisson_ion%rcart(1,iat)
    !    dy=(igpy-1)*poisson%hgrid(2,2)-poisson_ion%rcart(2,iat)
    !    dz=(igpz-1)*poisson%hgrid(3,3)-poisson_ion%rcart(3,iat)
    !    r2=dx**2+dy**2+dz**2
    !    if(r2<10.d0**2*gwt**2) then
    !        poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*exp(-r2/gwt**2)
    !    endif
    !enddo
    !enddo
    !enddo
    !enddo
    tt1=0.d0
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        tt1=tt1+poisson_ion%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    tt1=tt1*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    if(parini%mpi_env%iproc==0) then
    write(*,*) 'TT1 ',tt1
    endif
    !-------------------------------------------------------
    !write(*,*) atoms%zat(:)
    qtot=0.d0
    do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
            do igpx=1,poisson%ngpx
                poisson%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)-poisson%rho(igpx,igpy,igpz)
                !poisson%rho(igpx,igpy,igpz)=-poisson%rho(igpx,igpy,igpz)
                qtot=qtot+poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    qtot=qtot*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    atoms%qtot=qtot
    if(parini%mpi_env%iproc==0) then
        write(*,'(a,f20.12)') 'qtot= ',qtot
    endif
    center(1)=0.d0
    center(2)=0.d0
    center(3)=0.d0
    do iat=1,atoms%nat
        center(1)=center(1)+atoms%ratp(1,iat)
        center(2)=center(2)+atoms%ratp(2,iat)
        center(3)=center(3)+atoms%ratp(3,iat)
    enddo
    center(1)=center(1)/atoms%nat
    center(2)=center(2)/atoms%nat
    center(3)=center(3)/atoms%nat
    atoms%dpm(1)=0.d0
    atoms%dpm(2)=0.d0
    atoms%dpm(3)=0.d0
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        x=poisson%xyz111(1)+(igpx-1)*poisson%hgrid(1,1)-center(1)
        y=poisson%xyz111(2)+(igpy-1)*poisson%hgrid(2,2)-center(2)
        z=poisson%xyz111(3)+(igpz-1)*poisson%hgrid(3,3)-center(3)
        atoms%dpm(1)=atoms%dpm(1)+x*poisson%rho(igpx,igpy,igpz)
        atoms%dpm(2)=atoms%dpm(2)+y*poisson%rho(igpx,igpy,igpz)
        atoms%dpm(3)=atoms%dpm(3)+z*poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    atoms%dpm(1)=atoms%dpm(1)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    atoms%dpm(2)=atoms%dpm(2)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    atoms%dpm(3)=atoms%dpm(3)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    write(*,'(a,3f10.5)') 'DPM= ',atoms%dpm(1),atoms%dpm(2),atoms%dpm(3)
    !-------------------------------------------------------
    call update_ratp(atoms)
    call get_hartree(parini,poisson,atoms,gausswidth,ehartree_scn_excl)
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,es24.10,es14.5)') 'ehartree_scn_excl ',ehartree_scn_excl,poisson%screening_factor
    endif
    !-------------------------------------------------------
    if(.true.) then
    nd=1
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(pot(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
        pot(igpx,igpy,igpz)=poisson%pot(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    call f_free(poisson%rho)
    call f_free(poisson%pot)
    poisson%ngpz=poisson%ngpz/nd
    poisson%ngpy=poisson%ngpy/nd
    poisson%ngpx=poisson%ngpx/nd
    poisson%rho=f_malloc0([1.to.(poisson%ngpx),1.to.(poisson%ngpy),1.to.(poisson%ngpz)],id='poisson%rho')
    poisson%pot=f_malloc0([1.to.(poisson%ngpx),1.to.(poisson%ngpy),1.to.(poisson%ngpz)],id='poisson%pot')
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        poisson%rho(igpx,igpy,igpz)=rho(nd*igpx-nd+1,nd*igpy-nd+1,nd*igpz-nd+1)
        poisson%pot(igpx,igpy,igpz)=pot(nd*igpx-nd+1,nd*igpy-nd+1,nd*igpz-nd+1)
    enddo
    enddo
    enddo
    deallocate(rho)
    deallocate(pot)
    poisson%hgrid=poisson%hgrid*real(nd,kind=8)
    endif
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Li') atoms%zat(iat)=0.d0
        if(trim(atoms%sat(iat))=='S' ) atoms%zat(iat)=0.d0
    enddo
    !-------------------------------------------------------
    !atoms%fat=0.d0
    !call force_gto_sym_ortho(parini,poisson_ion%bc,atoms%nat,poisson_ion%rcart, &
    !    poisson_ion%q,gausswidth,6.d0,poisson_ion%xyz111, &
    !    poisson_ion%ngpx,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
    !    poisson_ion%hgrid,poisson%pot,atoms%fat)
    !if(parini%mpi_env%iproc==0) then
    !do iat=1,atoms%nat
    !    write(*,'(a,i4,3es19.10)') 'FAT ',iat,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
    !enddo
    !endif
    !-------------------------------------------------------
    call fini_hartree(parini,atoms,poisson_ion)
end subroutine get_poisson_ref
!*****************************************************************************************
pure function get_coeff_s(gw1,gw2) result(coeff)
    implicit none
    real(8), intent(in):: gw1, gw2
    !local variables
    real(8):: coeff
    coeff=gw1**3/(gw1**3-gw2**3)
end function get_coeff_s
!*****************************************************************************************
pure function get_coeff_p(gw1,gw2) result(coeff)
    implicit none
    real(8), intent(in):: gw1, gw2
    !local variables
    real(8):: coeff
    coeff=gw1**5/(gw1**5-gw2**5)
end function get_coeff_p
!*****************************************************************************************
pure function get_coeff_s_grad(gw1,gw2) result(coeff)
    implicit none
    real(8), intent(in):: gw1, gw2
    !local variables
    real(8):: coeff
    coeff=3.d0*gw1**2*gw2**3/(gw1**3-gw2**3)**2
end function get_coeff_s_grad
!*****************************************************************************************
pure function get_coeff_p_grad(gw1,gw2) result(coeff)
    implicit none
    real(8), intent(in):: gw1, gw2
    !local variables
    real(8):: coeff
    coeff=5.d0*gw1**4*gw2**5/(gw1**5-gw2**5)**2
end function get_coeff_p_grad
!*****************************************************************************************
end module mod_fit_bf_cent2
!*****************************************************************************************
