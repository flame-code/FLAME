!*****************************************************************************************
module mod_fit_bf_cent2
    implicit none
    private
    public:: get_basis_functions_cent2
    type typ_fitpar
        private
        integer:: ntypat
        logical, allocatable:: relaxcore(:)
        real(8), allocatable:: gwc_s1(:)
        real(8), allocatable:: gwc_s2(:)
        real(8), allocatable:: gwv_s1(:)
        real(8), allocatable:: gwv_s2(:)
        real(8), allocatable:: gwv_p1(:)
        real(8), allocatable:: gwv_p2(:)
        real(8), allocatable:: grad_gwc_s1(:)
        real(8), allocatable:: grad_gwc_s2(:)
        real(8), allocatable:: grad_gwv_s1(:)
        real(8), allocatable:: grad_gwv_s2(:)
        real(8), allocatable:: grad_gwv_p1(:)
        real(8), allocatable:: grad_gwv_p2(:)
        real(8), allocatable:: qcore_type(:)
        real(8), allocatable:: bc_s1(:)
        real(8), allocatable:: bc_s2(:)
        real(8), allocatable:: bv_s1(:)
        real(8), allocatable:: bv_s2(:)
        real(8), allocatable:: bv_p1(:)
        real(8), allocatable:: bv_p2(:)
        real(8), allocatable:: grad_bc_s1(:,:)
        real(8), allocatable:: grad_bc_s2(:,:)
        real(8), allocatable:: grad_bv_s1(:,:)
        real(8), allocatable:: grad_bv_s2(:,:)
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
        real(8), allocatable:: bcc_1(:)
        real(8), allocatable:: bcc_2(:)
        real(8), allocatable:: bccg_1(:,:)
        real(8), allocatable:: bccg_2(:,:)
        real(8), allocatable:: bc_1(:)
        real(8), allocatable:: bc_2(:)
        real(8), allocatable:: bcg_1(:,:)
        real(8), allocatable:: bcg_2(:,:)
        character(1), allocatable:: bt(:)
        contains
        procedure, private, pass(self):: init_bf
        procedure, private, pass(self):: fini_bf
        procedure, private, pass(self):: set_bt
        procedure, private, pass(self):: set_bf_bc_bcg
        procedure, private, pass(self):: get_pot_single
        procedure, private, pass(self):: get_pot_single_core
        !procedure, private, pass(self):: get_linearcoeff
    end type typ_bf
contains
!*****************************************************************************************
subroutine init_bf(self,nat,poisson)
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_bf), intent(inout):: self
    integer, intent(in):: nat
    type(typ_poisson), intent(in):: poisson
    !local variables
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
    allocate(self%bcc_1(nat))
    allocate(self%bcc_2(nat))
    allocate(self%bccg_1(2,nat))
    allocate(self%bccg_2(2,nat))
    allocate(self%bc_1(self%nbf))
    allocate(self%bc_2(self%nbf))
    allocate(self%bcg_1(2,self%nbf))
    allocate(self%bcg_2(2,self%nbf))
    allocate(self%bt(self%nbf))
    allocate(self%iat_list(self%nbf))
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
    deallocate(self%bcc_1)
    deallocate(self%bcc_2)
    deallocate(self%bccg_1)
    deallocate(self%bccg_2)
    deallocate(self%bc_1)
    deallocate(self%bc_2)
    deallocate(self%bcg_1)
    deallocate(self%bcg_2)
    deallocate(self%bt)
    deallocate(self%iat_list)
end subroutine fini_bf
!*****************************************************************************************
subroutine set_bt(self)
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_bf), intent(inout):: self
    !local variables
    integer:: ibf
    do ibf=1,self%nbf
        self%iat_list(ibf)=modulo(ibf-1,self%nbf/2)+1
        if((2*(ibf-1))/self%nbf==0) then
            self%bt(ibf)='s'
        elseif((2*(ibf-1))/self%nbf==1) then
            self%bt(ibf)='p'
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
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        self%bcc_1(iat)=fitpar%bc_s1(itypat)
        self%bcc_2(iat)=fitpar%bc_s2(itypat)
        self%bccg_1(1,iat)=fitpar%grad_bc_s1(1,itypat)
        self%bccg_1(2,iat)=fitpar%grad_bc_s1(2,itypat)
        self%bccg_2(1,iat)=fitpar%grad_bc_s2(1,itypat)
        self%bccg_2(2,iat)=fitpar%grad_bc_s2(2,itypat)
    enddo
    do ibf=1,self%nbf
        iat=self%iat_list(ibf)
        itypat=atoms%itypat(iat)
        if(self%bt(ibf)=='s') then
            self%bc_1(ibf)=fitpar%bv_s1(itypat)
            self%bc_2(ibf)=fitpar%bv_s2(itypat)
            self%bcg_1(1,ibf)=fitpar%grad_bv_s1(1,itypat)
            self%bcg_1(2,ibf)=fitpar%grad_bv_s1(2,itypat)
            self%bcg_2(1,ibf)=fitpar%grad_bv_s2(1,itypat)
            self%bcg_2(2,ibf)=fitpar%grad_bv_s2(2,itypat)
        elseif(self%bt(ibf)=='p') then
            self%bc_1(ibf)=fitpar%bv_p1(itypat)
            self%bc_2(ibf)=fitpar%bv_p2(itypat)
            self%bcg_1(1,ibf)=fitpar%grad_bv_p1(1,itypat)
            self%bcg_1(2,ibf)=fitpar%grad_bv_p1(2,itypat)
            self%bcg_2(1,ibf)=fitpar%grad_bv_p2(1,itypat)
            self%bcg_2(2,ibf)=fitpar%grad_bv_p2(2,itypat)
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
    if(self%bt(ibf)=='s') then
        q_tmp(1)=1.d0
        call cal_pot_gauss_s(parini,poisson,self%pot_1,.true.,atoms%ratp(1,iat),fitpar%gwv_s1(itypat),q_tmp(1),self%vgrad_1)
        call cal_pot_gauss_s(parini,poisson,self%pot_2,.true.,atoms%ratp(1,iat),fitpar%gwv_s2(itypat),q_tmp(1),self%vgrad_2)
    elseif(self%bt(ibf)=='p') then
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
    if(self%qcore(iat)<0.d0) then
    if(fitpar%gwc_s2(itypat)>0.d0) then
    q_tmp(1)=1.d0
    call cal_pot_gauss_s(parini,poisson,self%wa1,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),self%wa3)
    q_tmp(1)=1.d0
    call cal_pot_gauss_s(parini,poisson,self%wa2,.true.,atoms%ratp(1,iat),fitpar%gwc_s2(itypat),q_tmp(1),self%wa4)
    else
    q_tmp(1)=1.d0
    call cal_pot_gauss_s(parini,poisson,self%wa1,.true.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),self%wa3)
    self%wa2=0.d0
    self%wa4=0.d0
    endif
    endif
end subroutine get_pot_single_core
!*****************************************************************************************
subroutine init_fit_bf(self,ntypat)
    implicit none
    class(typ_fitpar), intent(inout):: self
    integer, intent(in):: ntypat
    !local variables
    self%ntypat=ntypat
    allocate(self%relaxcore(ntypat))
    allocate(self%gwv_s1(ntypat))
    allocate(self%gwv_s2(ntypat))
    allocate(self%gwc_s1(ntypat))
    allocate(self%gwc_s2(ntypat))
    allocate(self%gwv_p1(ntypat))
    allocate(self%gwv_p2(ntypat))
    allocate(self%grad_gwc_s1(ntypat))
    allocate(self%grad_gwc_s2(ntypat))
    allocate(self%grad_gwv_s1(ntypat))
    allocate(self%grad_gwv_s2(ntypat))
    allocate(self%grad_gwv_p1(ntypat))
    allocate(self%grad_gwv_p2(ntypat))
    allocate(self%qcore_type(ntypat))
    allocate(self%bc_s1(ntypat))
    allocate(self%bc_s2(ntypat))
    allocate(self%bv_s1(ntypat))
    allocate(self%bv_s2(ntypat))
    allocate(self%bv_p1(ntypat))
    allocate(self%bv_p2(ntypat))
    allocate(self%grad_bc_s1(2,ntypat))
    allocate(self%grad_bc_s2(2,ntypat))
    allocate(self%grad_bv_s1(2,ntypat))
    allocate(self%grad_bv_s2(2,ntypat))
    allocate(self%grad_bv_p1(2,ntypat))
    allocate(self%grad_bv_p2(2,ntypat))
end subroutine init_fit_bf
!*****************************************************************************************
subroutine fini_fit_bf(self)
    implicit none
    class(typ_fitpar), intent(inout):: self
    !local variables
    deallocate(self%relaxcore)
    deallocate(self%gwv_s1)
    deallocate(self%gwv_s2)
    deallocate(self%gwc_s1)
    deallocate(self%gwc_s2)
    deallocate(self%gwv_p1)
    deallocate(self%gwv_p2)
    deallocate(self%grad_gwc_s1)
    deallocate(self%grad_gwc_s2)
    deallocate(self%grad_gwv_s1)
    deallocate(self%grad_gwv_s2)
    deallocate(self%grad_gwv_p1)
    deallocate(self%grad_gwv_p2)
    deallocate(self%qcore_type)
    deallocate(self%bc_s1)
    deallocate(self%bc_s2)
    deallocate(self%bv_s1)
    deallocate(self%bv_s2)
    deallocate(self%bv_p1)
    deallocate(self%bv_p2)
    deallocate(self%grad_bc_s1)
    deallocate(self%grad_bc_s2)
    deallocate(self%grad_bv_s1)
    deallocate(self%grad_bv_s2)
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
        self%bv_s1(itypat)=get_coeff_s(self%gwv_s1(itypat),self%gwv_s2(itypat))
        self%bv_s2(itypat)=get_coeff_s(self%gwv_s2(itypat),self%gwv_s1(itypat))
        self%bv_p1(itypat)=get_coeff_p(self%gwv_p1(itypat),self%gwv_p2(itypat))
        self%bv_p2(itypat)=get_coeff_p(self%gwv_p2(itypat),self%gwv_p1(itypat))
        self%grad_bv_s1(1,itypat)=-get_coeff_s_grad(self%gwv_s1(itypat),self%gwv_s2(itypat))
        self%grad_bv_s1(2,itypat)= get_coeff_s_grad(self%gwv_s2(itypat),self%gwv_s1(itypat))
        self%grad_bv_s2(1,itypat)= get_coeff_s_grad(self%gwv_s1(itypat),self%gwv_s2(itypat))
        self%grad_bv_s2(2,itypat)=-get_coeff_s_grad(self%gwv_s2(itypat),self%gwv_s1(itypat))
        self%grad_bv_p1(1,itypat)=-get_coeff_p_grad(self%gwv_p1(itypat),self%gwv_p2(itypat))
        self%grad_bv_p1(2,itypat)= get_coeff_p_grad(self%gwv_p2(itypat),self%gwv_p1(itypat))
        self%grad_bv_p2(1,itypat)= get_coeff_p_grad(self%gwv_p1(itypat),self%gwv_p2(itypat))
        self%grad_bv_p2(2,itypat)=-get_coeff_p_grad(self%gwv_p2(itypat),self%gwv_p1(itypat))
    enddo
    do itypat=1,self%ntypat
        self%bc_s1(itypat)=0.d0
        self%bc_s2(itypat)=0.d0
        self%grad_bc_s1(1,itypat)=0.d0
        self%grad_bc_s1(2,itypat)=0.d0
        self%grad_bc_s2(1,itypat)=0.d0
        self%grad_bc_s2(2,itypat)=0.d0
        if(self%qcore_type(itypat)<0.d0 .and. self%gwc_s2(itypat)>0.d0) then
            self%bc_s1(itypat)=get_coeff_s(self%gwc_s1(itypat),self%gwc_s2(itypat))
            self%bc_s2(itypat)=get_coeff_s(self%gwc_s2(itypat),self%gwc_s1(itypat))
            self%grad_bc_s1(1,itypat)=-get_coeff_s_grad(self%gwc_s1(itypat),self%gwc_s2(itypat))
            self%grad_bc_s1(2,itypat)= get_coeff_s_grad(self%gwc_s2(itypat),self%gwc_s1(itypat))
            self%grad_bc_s2(1,itypat)= get_coeff_s_grad(self%gwc_s1(itypat),self%gwc_s2(itypat))
            self%grad_bc_s2(2,itypat)=-get_coeff_s_grad(self%gwc_s2(itypat),self%gwc_s1(itypat))
        endif
        if(self%qcore_type(itypat)<0.d0 .and. .not.(self%gwc_s2(itypat)>0.d0)) then
            self%bc_s1(itypat)=1.d0
            self%bc_s2(itypat)=0.d0
            self%grad_bc_s1(1,itypat)=0.d0
            self%grad_bc_s1(2,itypat)=0.d0
            self%grad_bc_s2(1,itypat)=0.d0
            self%grad_bc_s2(2,itypat)=0.d0
        endif
    enddo
end subroutine set_bc_bcg
!*****************************************************************************************
subroutine report_fit_bf(self,istep,cost,atoms,c_s,c_p,cost_gw)
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_fitpar), intent(inout):: self
    integer, intent(in):: istep
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: cost, c_s(atoms%nat) ,c_p(atoms%nat), cost_gw
    !local variables
    integer:: itypat, iat
    real(8):: tt1, tt2, tt3, tt4, dpm(3)
    write(*,'(a,i6,es14.5)',advance='no') 'core ',istep,cost
    do itypat=1,self%ntypat
        if(self%qcore_type(itypat)<0.d0) then
            tt1=self%grad_gwc_s1(itypat)
            write(*,'(es10.1)',advance='no') tt1
            if(self%gwc_s2(itypat)>0.d0) then
                tt1=self%grad_gwc_s2(itypat)
                write(*,'(es10.1)',advance='no') tt1
            endif
        endif
    enddo
    write(*,*)
    write(*,'(a,i6,es14.5)',advance='no') 'cost ',istep,cost
    do itypat=1,self%ntypat
        tt1=self%grad_gwv_s1(itypat)
        tt2=self%grad_gwv_s2(itypat)
        tt3=self%grad_gwv_p1(itypat)
        tt4=self%grad_gwv_p2(itypat)
        write(*,'(4es10.1)',advance='no') tt1,tt2,tt3,tt4
    enddo
    write(*,'(f8.3)') cost_gw
    write(*,'(a,i6)',advance='no') 'gwv ',istep
    do itypat=1,self%ntypat
        tt1=self%gwv_s1(itypat)
        tt2=self%gwv_s2(itypat)
        tt3=self%gwv_p1(itypat)
        tt4=self%gwv_p2(itypat)
        write(*,'(4f8.3)',advance='no') tt1,tt2,tt3,tt4
    enddo
    write(*,*)
    write(*,'(a,i6)',advance='no') 'gwc ',istep
    do itypat=1,self%ntypat
        if(self%qcore_type(itypat)<0.d0) then
            tt1=self%gwc_s1(itypat)
            write(*,'(f8.3)',advance='no') tt1
            if(self%gwc_s2(itypat)>0.d0) then
                tt1=self%gwc_s2(itypat)
                write(*,'(f8.3)',advance='no') tt1
            endif
        endif
    enddo
    write(*,*)
    write(*,'(a,i6)',advance='no') 'qp ',istep
    dpm(1)=0.d0
    dpm(2)=0.d0
    dpm(3)=0.d0
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        dpm(1)=dpm(1)+(atoms%zat(iat)+self%qcore_type(itypat)+c_s(iat))*atoms%ratp(1,iat)
        write(*,'(f8.3)',advance='no') c_s(iat)
    enddo
    dpm(1)=dpm(1)+sum(c_p(:))
    write(*,'(2f8.3)') sum(c_p(:)),dpm(1)
end subroutine report_fit_bf
!*****************************************************************************************
subroutine get_basis_functions_cent2(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, typ_atoms_arr, atom_deallocate_old
    !use mod_ann, only: typ_ann_arr
    use mod_ann_io_yaml, only: read_data_yaml
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_opt, only: typ_paropt
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson, poisson_ref
    type(typ_atoms):: atoms
    type(typ_atoms_arr):: atoms_arr
    !type(typ_ann_arr):: ann_arr
    type(typ_fitpar):: fitpar, fitpar_t
    type(typ_bf):: bf
    type(typ_paropt):: paropt
    integer:: iat, iconf, istep, ibf, jbf, nwork
    integer:: ix, iy, iz, itypat, nr, ir
    !real(8):: ehartree_scn_excl
    !real(8):: dpx, dpy, dpz
    !real(8):: time1, time2, time3, time4, time5, time6, time7, time8, time9
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6, tt7, cost, cost_gw
    real(8):: gausswidth(1), q_tmp(1), p_tmp(3), voxel
    real(8):: x, y, z !, r
    real(8):: errmax, rmse, alpha, ener
    real(8), allocatable:: c_s(:)
    real(8), allocatable:: c_p(:)
    real(8), allocatable:: xr(:)
    real(8), allocatable:: fr(:)
    real(8), allocatable:: work(:)
    real(8), allocatable:: pot_ion(:,:,:)
    call f_routine(id='get_basis_functions_cent2')
    !pi=4.d0*atan(1.d0)
    !call read_data_yaml(parini,'list_posinp_cent2.yaml',atoms_arr)
    call get_poisson_ref(parini,poisson_ref,atoms)
    do iat=1,atoms%nat
        write(*,*) atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
    enddo
    write(*,*) poisson_ref%xyz111(1),poisson_ref%xyz111(2),poisson_ref%xyz111(3)
    !do iat=1,atoms%nat
    !    if(trim(atoms%sat(iat))=='Li') atoms%zat(iat)=1.d0
    !    if(trim(atoms%sat(iat))=='S' ) atoms%zat(iat)=4.d0
    !enddo

    alpha=parini%paropt_geopt%alphax


    !-------------------------------------------------------
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

    call fitpar%init_fit_bf(parini%ntypat)
    call fitpar_t%init_fit_bf(parini%ntypat)
    call bf%init_bf(atoms%nat,poisson)
    allocate(c_s(atoms%nat))
    allocate(c_p(atoms%nat))
    allocate(pot_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))

    do iat=1,atoms%nat
        do itypat=1,parini%ntypat
            if(trim(atoms%sat(iat))==trim(parini%stypat(itypat))) then
                atoms%itypat(iat)=parini%ltypat(itypat)
                exit
            endif
        enddo
    enddo

    gausswidth=0.5d0

    !q_tmp=2.d0
    !call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,1),q_tmp,gausswidth, &
    !    6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)

    open(unit=2134,file="gw.inp",status='old')
    do itypat=1,parini%ntypat
        read(2134,*) fitpar%relaxcore(itypat),tt1,tt2,tt3,tt4,tt5,tt6,tt7
        fitpar%qcore_type(itypat)=tt1
        fitpar%gwc_s1(itypat)=tt2
        fitpar%gwc_s2(itypat)=tt3
        fitpar%gwv_s1(itypat)=tt4
        fitpar%gwv_s2(itypat)=tt5
        fitpar%gwv_p1(itypat)=tt6
        fitpar%gwv_p2(itypat)=tt7
    enddo
    close(2134)
    call fitpar%set_bc_bcg()
    call bf%set_bf_bc_bcg(atoms,fitpar)

    !do iat=1,atoms%nat
    !    if(trim(atoms%sat(iat))=='O') atoms%zat(iat)=4.d0
    !enddo

    !gausswidth=0.5d0
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        bf%qcore(iat)=fitpar%qcore_type(itypat)
    enddo
    pot_ion=0.d0
    do iat=1,atoms%nat
    itypat=atoms%itypat(iat)
    call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),gausswidth(1),atoms%zat(iat),bf%vgrad_1)
    enddo

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
        if(fitpar%relaxcore(itypat) .and. fitpar%qcore_type(itypat)<0.d0) then
            if(fitpar%gwc_s2(itypat)>0.d0) then
                nr=nr+2
            else
                nr=nr+1
            endif
        endif
    enddo
    allocate(xr(nr),fr(nr))
    !nwork=3*nr !fire
    nwork=nr*nr+3*nr+3*nr*nr+3*nr !mybfgs
    allocate(work(nwork))
    ir=0
    do itypat=1,parini%ntypat
        if(fitpar%relaxcore(itypat) .and. fitpar%qcore_type(itypat)<0.d0) then
            if(fitpar%gwc_s2(itypat)>0.d0) then
                xr(ir+1)=fitpar%gwc_s1(itypat)
                xr(ir+2)=fitpar%gwc_s2(itypat)
                ir=ir+2
            else
                xr(ir+1)=fitpar%gwc_s1(itypat)
                ir=ir+1
            endif
        endif
        ir=ir+1 ; xr(ir)=fitpar%gwv_s1(itypat)
        ir=ir+1 ; xr(ir)=fitpar%gwv_s2(itypat)
        ir=ir+1 ; xr(ir)=fitpar%gwv_p1(itypat)
        ir=ir+1 ; xr(ir)=fitpar%gwv_p2(itypat)
    enddo
    istep=0
    do
    poisson%pot=pot_ion

    voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)

    do iat=1,atoms%nat
        call bf%get_pot_single_core(parini,poisson,atoms,fitpar,iat)
        tt1=bf%qcore(iat)*bf%bcc_1(iat)
        tt2=bf%qcore(iat)*bf%bcc_2(iat)
        poisson%pot=poisson%pot+tt1*bf%wa1+tt2*bf%wa2
    enddo
    call get_linearcoeff(parini,fitpar,bf,poisson,poisson_ref,atoms,c_s,c_p)
    do ibf=1,bf%nbf
        call bf%get_pot_single(parini,poisson,atoms,fitpar,ibf)
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt1=bf%bc_1(ibf)*bf%pot_1(ix,iy,iz)+bf%bc_2(ibf)*bf%pot_2(ix,iy,iz)
            poisson%pot(ix,iy,iz)=poisson%pot(ix,iy,iz)+bf%csp(ibf)*tt1
        enddo
        enddo
        enddo
    enddo
    fitpar%grad_gwc_s1=0.d0
    fitpar%grad_gwc_s2=0.d0
    if(any(fitpar%relaxcore)) then
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        call bf%get_pot_single_core(parini,poisson,atoms,fitpar,iat)
        tt1=0.d0
        tt2=0.d0
        do iz=1,poisson_ref%ngpz
        do iy=1,poisson_ref%ngpy
        do ix=1,poisson_ref%ngpx
            tt3=bf%bccg_1(1,iat)*bf%wa1(ix,iy,iz)+bf%bccg_2(1,iat)*bf%wa2(ix,iy,iz)+bf%bcc_1(iat)*bf%wa3(ix,iy,iz)
            tt4=bf%bccg_1(2,iat)*bf%wa1(ix,iy,iz)+bf%bccg_2(2,iat)*bf%wa2(ix,iy,iz)+bf%bcc_2(iat)*bf%wa4(ix,iy,iz)
            tt1=tt1+bf%qcore(iat)*tt3*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
            tt2=tt2+bf%qcore(iat)*tt4*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
        enddo
        enddo
        enddo
        fitpar%grad_gwc_s1(itypat)=fitpar%grad_gwc_s1(itypat)+2.d0*tt1*voxel
        fitpar%grad_gwc_s2(itypat)=fitpar%grad_gwc_s2(itypat)+2.d0*tt2*voxel
    enddo
    endif
    !----------------------------------------------------------
    cost=bf%zeroder
    do jbf=1,bf%nbf
        do ibf=1,bf%nbf
            cost=cost+bf%csp(ibf)*bf%secder(ibf,jbf)*bf%csp(jbf)
        enddo
    enddo
    do ibf=1,bf%nbf
        cost=cost+2.d0*bf%csp(ibf)*bf%firstder(ibf)
    enddo
    do itypat=1,parini%ntypat
        fitpar%grad_gwv_s1(itypat)=0.d0
        fitpar%grad_gwv_s2(itypat)=0.d0
        fitpar%grad_gwv_p1(itypat)=0.d0
        fitpar%grad_gwv_p2(itypat)=0.d0
    enddo
    do ibf=1,bf%nbf
        iat=bf%iat_list(ibf)
        itypat=atoms%itypat(iat)
        tt1=bf%csp(ibf)*(bf%bcg_1(1,ibf)*bf%subfirstder(1,ibf)+bf%bcg_2(1,ibf)*bf%subfirstder(2,ibf))
        tt2=bf%csp(ibf)*(bf%bcg_1(2,ibf)*bf%subfirstder(1,ibf)+bf%bcg_2(2,ibf)*bf%subfirstder(2,ibf))
        if(bf%bt(ibf)=='s') then
            fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*tt1
            fitpar%grad_gwv_s2(itypat)=fitpar%grad_gwv_s2(itypat)+2.d0*tt2
        elseif(bf%bt(ibf)=='p') then
            fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*tt1
            fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*tt2
        else
            stop 'ERROR: unknwon br when computing gwgrad'
        endif
        do jbf=1,bf%nbf
            tt1=bf%csp(ibf)*bf%bcg_1(1,ibf)*bf%subsecder(1,1,ibf,jbf)*bf%csp(jbf)*bf%bc_1(jbf)
            tt2=bf%csp(ibf)*bf%bcg_1(1,ibf)*bf%subsecder(1,2,ibf,jbf)*bf%csp(jbf)*bf%bc_2(jbf)
            tt3=bf%csp(ibf)*bf%bcg_2(1,ibf)*bf%subsecder(2,1,ibf,jbf)*bf%csp(jbf)*bf%bc_1(jbf)
            tt4=bf%csp(ibf)*bf%bcg_2(1,ibf)*bf%subsecder(2,2,ibf,jbf)*bf%csp(jbf)*bf%bc_2(jbf)
            if(bf%bt(ibf)=='s') then
                fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            elseif(bf%bt(ibf)=='p') then
                fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            else
                stop 'ERROR: unknwon br when computing gwgrad'
            endif
            tt1=bf%csp(ibf)*bf%bcg_1(2,ibf)*bf%subsecder(1,1,ibf,jbf)*bf%csp(jbf)*bf%bc_1(jbf)
            tt2=bf%csp(ibf)*bf%bcg_1(2,ibf)*bf%subsecder(1,2,ibf,jbf)*bf%csp(jbf)*bf%bc_2(jbf)
            tt3=bf%csp(ibf)*bf%bcg_2(2,ibf)*bf%subsecder(2,1,ibf,jbf)*bf%csp(jbf)*bf%bc_1(jbf)
            tt4=bf%csp(ibf)*bf%bcg_2(2,ibf)*bf%subsecder(2,2,ibf,jbf)*bf%csp(jbf)*bf%bc_2(jbf)
            if(bf%bt(ibf)=='s') then
                fitpar%grad_gwv_s2(itypat)=fitpar%grad_gwv_s2(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            elseif(bf%bt(ibf)=='p') then
                fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*(tt1+tt2+tt3+tt4)
            else
                stop 'ERROR: unknwon br when computing gwgrad'
            endif
        enddo
        tt1=bf%csp(ibf)*(bf%bc_1(ibf)*bf%subfirstder_grad(1,ibf))
        tt2=bf%csp(ibf)*(bf%bc_2(ibf)*bf%subfirstder_grad(2,ibf))
        if(bf%bt(ibf)=='s') then
            fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*tt1
            fitpar%grad_gwv_s2(itypat)=fitpar%grad_gwv_s2(itypat)+2.d0*tt2
        elseif(bf%bt(ibf)=='p') then
            fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*tt1
            fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*tt2
        else
            stop 'ERROR: unknwon br when computing gwgrad'
        endif
        do jbf=1,bf%nbf
            tt1=bf%csp(ibf)*bf%bc_1(ibf)*bf%subsecder_grad(1,1,ibf,jbf)*bf%csp(jbf)*bf%bc_1(jbf)
            tt2=bf%csp(ibf)*bf%bc_1(ibf)*bf%subsecder_grad(1,2,ibf,jbf)*bf%csp(jbf)*bf%bc_2(jbf)
            tt3=bf%csp(ibf)*bf%bc_2(ibf)*bf%subsecder_grad(2,1,ibf,jbf)*bf%csp(jbf)*bf%bc_1(jbf)
            tt4=bf%csp(ibf)*bf%bc_2(ibf)*bf%subsecder_grad(2,2,ibf,jbf)*bf%csp(jbf)*bf%bc_2(jbf)
            if(bf%bt(ibf)=='s') then
                fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)+2.d0*(tt1+tt2)
                fitpar%grad_gwv_s2(itypat)=fitpar%grad_gwv_s2(itypat)+2.d0*(tt3+tt4)
            elseif(bf%bt(ibf)=='p') then
                fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)+2.d0*(tt1+tt2)
                fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+2.d0*(tt3+tt4)
            else
                stop 'ERROR: unknwon br when computing gwgrad'
            endif
        enddo
    enddo

    !call get_cost_gw(parini,fitpar,cost_gw)
    !cost=cost+cost_gw
    call get_cost_gw_new(parini,fitpar,fitpar_t,cost_gw)
    do itypat=1,parini%ntypat
        fitpar%grad_gwc_s1(itypat)=cost_gw*fitpar%grad_gwc_s1(itypat)+cost*fitpar_t%grad_gwc_s1(itypat)
        fitpar%grad_gwc_s2(itypat)=cost_gw*fitpar%grad_gwc_s2(itypat)+cost*fitpar_t%grad_gwc_s2(itypat)
        fitpar%grad_gwv_s1(itypat)=cost_gw*fitpar%grad_gwv_s1(itypat)+cost*fitpar_t%grad_gwv_s1(itypat)
        fitpar%grad_gwv_s2(itypat)=cost_gw*fitpar%grad_gwv_s2(itypat)+cost*fitpar_t%grad_gwv_s2(itypat)
        fitpar%grad_gwv_p1(itypat)=cost_gw*fitpar%grad_gwv_p1(itypat)+cost*fitpar_t%grad_gwv_p1(itypat)
        fitpar%grad_gwv_p2(itypat)=cost_gw*fitpar%grad_gwv_p2(itypat)+cost*fitpar_t%grad_gwv_p2(itypat)
    enddo
    cost=cost*cost_gw

    call fitpar%report_fit_bf(istep,cost,atoms,c_s,c_p,cost_gw)

    if(istep==parini%paropt_geopt%nit) exit
    ir=0
    do itypat=1,parini%ntypat
        if(fitpar%relaxcore(itypat) .and. fitpar%qcore_type(itypat)<0.d0) then
            if(fitpar%gwc_s2(itypat)>0.d0) then
                fr(ir+1)=-fitpar%grad_gwc_s1(itypat)
                fr(ir+2)=-fitpar%grad_gwc_s2(itypat)
                ir=ir+2
            else
                fr(ir+1)=-fitpar%grad_gwc_s1(itypat)
                ir=ir+1
            endif
        endif
        ir=ir+1 ; fr(ir)=-fitpar%grad_gwv_s1(itypat)
        ir=ir+1 ; fr(ir)=-fitpar%grad_gwv_s2(itypat)
        ir=ir+1 ; fr(ir)=-fitpar%grad_gwv_p1(itypat)
        ir=ir+1 ; fr(ir)=-fitpar%grad_gwv_p2(itypat)
    enddo
    !call fire(parini,parini%mpi_env%iproc,nr,xr,cost,fr,work,paropt)
    call mybfgs(parini%mpi_env%iproc,nr,xr,cost,fr,nwork,work,paropt)
    if(paropt%iflag<=0) exit
    ir=0
    do itypat=1,parini%ntypat
        if(fitpar%relaxcore(itypat) .and. fitpar%qcore_type(itypat)<0.d0) then
            if(fitpar%gwc_s2(itypat)>0.d0) then
                fitpar%gwc_s1(itypat)=xr(ir+1)
                fitpar%gwc_s2(itypat)=xr(ir+2)
                ir=ir+2
            else
                fitpar%gwc_s1(itypat)=xr(ir+1)
                ir=ir+1
            endif
        endif
        ir=ir+1 ; fitpar%gwv_s1(itypat)=xr(ir)
        ir=ir+1 ; fitpar%gwv_s2(itypat)=xr(ir)
        ir=ir+1 ; fitpar%gwv_p1(itypat)=xr(ir)
        ir=ir+1 ; fitpar%gwv_p2(itypat)=xr(ir)
    enddo
    call fitpar%set_bc_bcg()
    call bf%set_bf_bc_bcg(atoms,fitpar)
    istep=istep+1
    enddo !end of loop over istep
    deallocate(xr,fr,work)
    call yaml_sequence_close()
    !-----------------------------------------------------------------
    errmax=0.d0
    rmse=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        errmax=max(errmax,abs(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)))
        rmse=rmse+(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))**2
    enddo
    enddo
    enddo
    rmse=sqrt(rmse/real(poisson%ngpx*poisson%ngpy*poisson%ngpz,kind=8))
    write(*,'(a,2es14.5)') 'errmax,rmse= ',errmax,rmse
    open(unit=2134,file="gw.out",status='replace')
    do itypat=1,parini%ntypat
        tt1=fitpar%qcore_type(itypat)
        tt2=fitpar%gwc_s1(itypat)
        tt3=fitpar%gwc_s2(itypat)
        tt4=fitpar%gwv_s1(itypat)
        tt5=fitpar%gwv_s2(itypat)
        tt6=fitpar%gwv_p1(itypat)
        tt7=fitpar%gwv_p2(itypat)
        write(2134,'(l2,7f10.6)') fitpar%relaxcore(itypat),tt1,tt2,tt3,tt4,tt5,tt6,tt7
    enddo
    close(2134)
    !do iz=1,poisson%ngpz
    do iz=poisson%ngpz/2,poisson%ngpz/2
    do iy=poisson%ngpy/2,poisson%ngpy/2
    !do ix=poisson%ngpx/2,poisson%ngpx/2
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
        write(83,'(3f8.3,2es19.10)') x,y,z,poisson_ref%pot(ix,iy,iz),poisson%pot(ix,iy,iz)
    enddo
    enddo
    enddo
    !----------------------------------------------------------
    !poisson%rho=0.d0
    !do iat=1,atoms%nat
    !    itypat=atoms%itypat(iat)
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),atoms%zat(iat),gausswidth(1), &
    !    6.d0*gausswidth(1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !if(bf%qcore(iat)<0.d0) then
    !q_tmp(1)=bf%qcore(iat)*fitpar%bc_s1(itypat)
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwc_s1(itypat), &
    !    6.d0*fitpar%gwc_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !q_tmp(1)=bf%qcore(iat)*fitpar%bc_s2(itypat)
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwc_s2(itypat), &
    !    6.d0*fitpar%gwc_s2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !endif
    !q_tmp(1)=c_s(iat)*fitpar%bv_s1(itypat)
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwv_s1(itypat), &
    !    4.d0*fitpar%gwv_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !q_tmp(1)=c_s(iat)*fitpar%bv_s2(itypat)
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwv_s2(itypat), &
    !    4.d0*fitpar%gwv_s2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !p_tmp=0.d0
    !p_tmp(1)=c_p(iat)*fitpar%bv_p1(itypat)
    !call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),p_tmp,fitpar%gwv_p1(itypat), &
    !    6.d0*fitpar%gwv_p1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    !p_tmp(1)=c_p(iat)*fitpar%bv_p2(itypat)
    !call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),p_tmp,fitpar%gwv_p2(itypat), &
    !    6.d0*fitpar%gwv_p2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    !enddo
    !poisson%rho=poisson%rho-poisson_ref%rho
    !call cube_write('diffrho.cube',atoms,poisson,'rho')
    !call get_hartree(parini,poisson,atoms,gausswidth,ener)
    ener=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        ener=ener+poisson_ref%rho(ix,iy,iz)*poisson%pot(ix,iy,iz)
    enddo
    enddo
    enddo
    ener=ener*0.5d0*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    !if(parini%mpi_env%iproc==0) then
    write(*,'(a,es24.10,es14.5)') 'ehartree_scn_excl ',ener,poisson%screening_factor
    !endif
    !----------------------------------------------------------


    !do iconf=1,atoms_arr%nconf
    !    call atom_deallocate_old(atoms_arr%atoms(iconf))
    !enddo
    !deallocate(atoms_arr%atoms)
    call fini_hartree(parini,atoms,poisson_ref)
    call fini_hartree(parini,atoms,poisson)
    call fitpar%fini_fit_bf()
    call fitpar_t%fini_fit_bf()
    call bf%fini_bf()

    call f_release_routine()
end subroutine get_basis_functions_cent2
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
    do itypat=1,parini%ntypat
        dgw=fitpar%gwv_s1(itypat)-fitpar%gwv_s2(itypat)
        if(abs(dgw)<dgw_min) then
            cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
            tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
            fitpar_t%grad_gwv_s1(itypat)=-tt
            fitpar_t%grad_gwv_s2(itypat)= tt
        endif
        dgw=fitpar%gwv_p1(itypat)-fitpar%gwv_p2(itypat)
        if(abs(dgw)<dgw_min) then
            cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
            tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
            fitpar_t%grad_gwv_p1(itypat)=-tt
            fitpar_t%grad_gwv_p2(itypat)= tt
        endif
        if(fitpar%relaxcore(itypat) .and. fitpar%qcore_type(itypat)<0.d0 .and. fitpar%gwc_s2(itypat)>0.d0) then
            dgw=fitpar%gwc_s1(itypat)-fitpar%gwc_s2(itypat)
            if(abs(dgw)<dgw_min) then
                cost_gw=cost_gw+pref*(1.d0-dgw/dgw_min)**4
                tt=pref*4.d0*(1.d0/dgw_min)*(1.d0-dgw/dgw_min)**3
                fitpar_t%grad_gwc_s1(itypat)=-tt
                fitpar_t%grad_gwc_s2(itypat)= tt
            endif
        endif
    enddo
end subroutine get_cost_gw_new
!*****************************************************************************************
subroutine get_cost_gw(parini,fitpar,cost_gw)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_fitpar), intent(inout):: fitpar
    real(8), intent(out):: cost_gw
    !local variables
    real(8):: dgw, dgw_min, tt, pref
    integer:: itypat
    dgw_min=0.5d0
    pref=1.d0*parini%screening_factor**2
    cost_gw=0.d0
    do itypat=1,parini%ntypat
        dgw=fitpar%gwv_s1(itypat)-fitpar%gwv_s2(itypat)
        if(abs(dgw)<dgw_min) then
            cost_gw=cost_gw+pref*(1.d0-(dgw/dgw_min)**2)**10
            tt=pref*20.d0*(dgw/dgw_min**2)*(1.d0-(dgw/dgw_min)**2)**9
            fitpar%grad_gwv_s1(itypat)=fitpar%grad_gwv_s1(itypat)-tt
            fitpar%grad_gwv_s2(itypat)=fitpar%grad_gwv_s2(itypat)+tt
        endif
        dgw=fitpar%gwv_p1(itypat)-fitpar%gwv_p2(itypat)
        if(abs(dgw)<dgw_min) then
            cost_gw=cost_gw+pref*(1.d0-(dgw/dgw_min)**2)**10
            tt=pref*20.d0*(dgw/dgw_min**2)*(1.d0-(dgw/dgw_min)**2)**9
            fitpar%grad_gwv_p1(itypat)=fitpar%grad_gwv_p1(itypat)-tt
            fitpar%grad_gwv_p2(itypat)=fitpar%grad_gwv_p2(itypat)+tt
        endif
    enddo
end subroutine get_cost_gw
!*****************************************************************************************
subroutine get_linearcoeff(parini,fitpar,bf,poisson,poisson_ref,atoms,c_s,c_p)
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
    real(8):: voxel, dpm(3)
    real(8):: tt1, tt2, tt3, tt4
    integer:: ix, iy, iz, iat, info, nbf, nc
    integer:: ibf, jbf, lwork
    real(8), allocatable:: secder(:,:), rhs(:)
    real(8), allocatable:: secder_t(:,:), work(:), eval(:)
    integer, allocatable:: ipiv(:)
    nbf=2*atoms%nat
    nc=1
    allocate(secder(nbf+nc,nbf+nc),source=0.d0)
    allocate(rhs(nbf+nc),source=0.d0)
    voxel=poisson_ref%hgrid(1,1)*poisson_ref%hgrid(2,2)*poisson_ref%hgrid(3,3)
    do ibf=1,bf%nbf
        call bf%get_pot_single(parini,poisson,atoms,fitpar,ibf)
        tt1=0.d0
        tt2=0.d0
        do iz=1,poisson_ref%ngpz
        do iy=1,poisson_ref%ngpy
        do ix=1,poisson_ref%ngpx
            tt1=tt1+bf%pot_1(ix,iy,iz)*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
            tt2=tt2+bf%pot_2(ix,iy,iz)*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
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
            tt1=tt1+bf%vgrad_1(ix,iy,iz)*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
            tt2=tt2+bf%vgrad_2(ix,iy,iz)*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
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
            call bf%get_pot_single(parini,poisson,atoms,fitpar,jbf)
            tt3=0.d0
            do iz=1,poisson_ref%ngpz
            do iy=1,poisson_ref%ngpy
            do ix=1,poisson_ref%ngpx
                tt1=bf%bc_1(ibf)*bf%wa1(ix,iy,iz)+bf%bc_2(ibf)*bf%wa2(ix,iy,iz)
                tt2=bf%bc_1(jbf)*bf%pot_1(ix,iy,iz)+bf%bc_2(jbf)*bf%pot_2(ix,iy,iz)
                tt3=tt3+tt1*tt2
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
                tt1=tt1+bf%pot_1(ix,iy,iz)*bf%wa1(ix,iy,iz)
                tt2=tt2+bf%pot_2(ix,iy,iz)*bf%wa1(ix,iy,iz)
                tt3=tt3+bf%pot_1(ix,iy,iz)*bf%wa2(ix,iy,iz)
                tt4=tt4+bf%pot_2(ix,iy,iz)*bf%wa2(ix,iy,iz)
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
                tt1=tt1+bf%pot_1(ix,iy,iz)*bf%wa3(ix,iy,iz)
                tt2=tt2+bf%pot_2(ix,iy,iz)*bf%wa3(ix,iy,iz)
                tt3=tt3+bf%pot_1(ix,iy,iz)*bf%wa4(ix,iy,iz)
                tt4=tt4+bf%pot_2(ix,iy,iz)*bf%wa4(ix,iy,iz)
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
                tt1=tt1+bf%vgrad_1(ix,iy,iz)*bf%wa1(ix,iy,iz)
                tt2=tt2+bf%vgrad_2(ix,iy,iz)*bf%wa1(ix,iy,iz)
                tt3=tt3+bf%vgrad_1(ix,iy,iz)*bf%wa2(ix,iy,iz)
                tt4=tt4+bf%vgrad_2(ix,iy,iz)*bf%wa2(ix,iy,iz)
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
            tt2=bf%bc_1(ibf)*bf%wa1(ix,iy,iz)+bf%bc_2(ibf)*bf%wa2(ix,iy,iz)
            tt1=tt1+tt2*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
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
        tt1=tt1+(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))**2
    enddo
    enddo
    enddo
    bf%zeroder=tt1*voxel

    do iat=1,atoms%nat
        write(*,'(a,6es14.5)') 'SECDER ',secder(iat,1:6)
    enddo
    do iat=1,atoms%nat
        write(*,'(a,6es14.5)') 'SECDER ',secder(atoms%nat+iat,1:6)
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
    secder(nbf+1      ,nbf+1      )=0.d0
    !secder(nbf+2      ,nbf+2      )=0.d0
    do iat=1,atoms%nat
        secder(iat          ,nbf+1        )=1.d0
        secder(nbf+1        ,iat          )=1.d0
        !secder(iat          ,nbf+2        )=atoms%ratp(1,iat)
        !secder(nbf+2        ,iat          )=atoms%ratp(1,iat)
        !secder(atoms%nat+iat,nbf+2        )=1.d0
        !secder(nbf+2        ,atoms%nat+iat)=1.d0
    enddo
    rhs(nbf+1)=atoms%qtot-sum(atoms%zat)-sum(bf%qcore) !-sum(atoms%zat)
    dpm(1)=0.d0 ; dpm(2)=0.d0 ; dpm(3)=0.d0
    !rhs(nbf+2)=-2.786975d0-(2.d0*atoms%ratp(1,1)+4.d0*atoms%ratp(1,2))
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
    real(8):: x, y, z, r, tt1, p_dot_r
    real(8):: uu1, vv1, qq1, sft
    real(8):: yy1, ff1, hh1
    real(8):: tg0, tg1, tg2, tg3, tg4
    integer:: ix, iy, iz
    pi=4.d0*atan(1.d0)
    if(reset) then
        pot=0.d0
    endif
    sf=parini%screening_factor
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-xyz(1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-xyz(2)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-xyz(3)
        r=sqrt(x**2+y**2+z**2)
        p_dot_r=x*p(1)+y*p(2)+z*p(3)
        if(r<1.d-10) then
            ss0=4.d0*p_dot_r/(3.d0*sqrt(pi)*gw**3)
            tg0=-3.d0*4.d0*p_dot_r/(3.d0*sqrt(pi)*gw**4)
            sft=sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss1=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg1=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
            sft=2.d0*sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss2=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg2=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
            sft=3.d0*sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss3=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg3=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
            sft=4.d0*sf
            tt1=3.d0*sqrt(pi)*gw**3*(1.d0+sft**2*gw**2)**2
            ss4=4.d0*p_dot_r*(1.d0+3.d0*sft**2*gw**2)/tt1
            yy1=3.d0+10.d0*gw**2*sft**2+15.d0*gw**4*sft**4
            tg4=(-4.d0*p_dot_r*yy1)/(3.d0*sqrt(pi)*gw**4*tt1**3)
        else
            ss0=p_dot_r*(-2.d0*r*exp(-(r/gw)**2)/(sqrt(pi)*gw)+erf(r/gw))/r**3
            tg0=-4.d0*p_dot_r*exp(-(r/gw)**2)/(sqrt(pi)*gw**4)
            sft=sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss1=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg1=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
            sft=2.d0*sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss2=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg2=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
            sft=3.d0*sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss3=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg3=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
            sft=4.d0*sf
            tt1=1.d0+sft**2*gw**2
            uu1=exp(-(sft*r)**2/tt1)*sqrt(pi)*gw
            vv1=(1.d0+2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            qq1=sqrt(pi)*r**3*gw*(tt1)**1.5d0
            ss4=p_dot_r*(-2.d0*r*sqrt(tt1)*exp(-(r/gw)**2)+uu1*vv1)/qq1
            hh1=2.d0*r*tt1*(sft**2*gw**4*tt1-2.d0*r**2*(1.d0+3.d0*sft**2*gw**2*tt1))
            yy1=sft**2*gw**5*sqrt(pi*tt1)*erf(r/(gw*sqrt(tt1)))
            ff1=(-4.d0*r**4*sft**4+tt1**2+4.d0*r**2*sft**2*tt1)*exp(-(sft*r)**2/tt1)
            tg4=p_dot_r*(hh1*exp(-(r/gw)**2)-yy1*ff1)/(sqrt(pi)*r**3*gw**4*tt1**4)
        endif
        !write(*,'(3i4,5es14.5)') ix,iy,iz,ss0,ss1,ss2,ss3,ss4
        pot(ix,iy,iz)=pot(ix,iy,iz)+ss0-4.d0*ss1+6.d0*ss2-4.d0*ss3+ss4
        vgrad(ix,iy,iz)=tg0-4.d0*tg1+6.d0*tg2-4.d0*tg3+tg4
        !poisson%pot(ix,iy,iz)=ss0-2.d0*ss1+ss2
    enddo
    enddo
    enddo
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
    real(8):: pi, sf, ss0, ss1, ss2, ss3, ss4
    real(8):: x, y, z, r, tt1, uu1, ww1, qq1, vv1, sft
    real(8):: tg0, tg1, tg2, tg3, tg4
    real(8):: dd0, dd1, dd2, dd3, dd4
    integer:: ix, iy, iz
    pi=4.d0*atan(1.d0)
    if(reset) then
        pot=0.d0
    endif
    sf=parini%screening_factor
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-xyz(1)
        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-xyz(2)
        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-xyz(3)
        r=sqrt(x**2+y**2+z**2)
        if(r<1.d-10) then
            ss0=2.d0/(sqrt(pi)*gw)
            ss1=2.d0/(sqrt(pi)*gw*(1.d0+sf**2*gw**2))
            ss2=2.d0/(sqrt(pi)*gw*(1.d0+(2.d0*sf)**2*gw**2))
            ss3=2.d0/(sqrt(pi)*gw*(1.d0+(3.d0*sf)**2*gw**2))
            ss4=2.d0/(sqrt(pi)*gw*(1.d0+(4.d0*sf)**2*gw**2))
            tg0=-2.d0/(sqrt(pi)*gw**2)
            tg1=-2.d0*(1.d0+3.d0*sf**2*gw**2)/(sqrt(pi)*(gw*(1.d0+sf**2*gw**2))**2)
            tg2=-2.d0*(1.d0+3.d0*(2.d0*sf)**2*gw**2)/(sqrt(pi)*(gw*(1.d0+(2.d0*sf)**2*gw**2))**2)
            tg3=-2.d0*(1.d0+3.d0*(3.d0*sf)**2*gw**2)/(sqrt(pi)*(gw*(1.d0+(3.d0*sf)**2*gw**2))**2)
            tg4=-2.d0*(1.d0+3.d0*(4.d0*sf)**2*gw**2)/(sqrt(pi)*(gw*(1.d0+(4.d0*sf)**2*gw**2))**2)
        else
            dd0=(r/gw)**2
            if(dd0>50.d0) then
                dd0=0.d0
            else
                dd0=exp(-dd0)
            endif
            dd1=(1.d0*sf*r)**2/(1.d0+(1.d0*sf)**2*gw**2)
            if(dd1>50.d0) then
                dd1=0.d0
            else
                dd1=exp(-dd1)
            endif
            dd2=(2.d0*sf*r)**2/(1.d0+(2.d0*sf)**2*gw**2)
            if(dd2>50.d0) then
                dd2=0.d0
            else
                dd2=exp(-dd2)
            endif
            dd3=(3.d0*sf*r)**2/(1.d0+(3.d0*sf)**2*gw**2)
            if(dd3>50.d0) then
                dd3=0.d0
            else
                dd3=exp(-dd3)
            endif
            dd4=(4.d0*sf*r)**2/(1.d0+(4.d0*sf)**2*gw**2)
            if(dd4>50.d0) then
                dd4=0.d0
            else
                dd4=exp(-dd4)
            endif
            ss0=erf(r/gw)/r
            tt1=1.d0+sf**2*gw**2
            ss1=erf(r/(gw*sqrt(tt1)))*dd1/(r*sqrt(tt1))
            tt1=1.d0+(2.d0*sf)**2*gw**2
            ss2=erf(r/(gw*sqrt(tt1)))*dd2/(r*sqrt(tt1))
            tt1=1.d0+(3.d0*sf)**2*gw**2
            ss3=erf(r/(gw*sqrt(tt1)))*dd3/(r*sqrt(tt1))
            tt1=1.d0+(4.d0*sf)**2*gw**2
            ss4=erf(r/(gw*sqrt(tt1)))*dd4/(r*sqrt(tt1))
            tg0=-2.d0*dd0/(sqrt(pi)*gw**2)
            sft=sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=dd1*sqrt(pi)*gw**3*sft**2
            tg1=(-2.d0*r*sqrt(tt1)*qq1*dd0-vv1*ww1)/uu1
            sft=2.d0*sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=dd2*sqrt(pi)*gw**3*sft**2
            tg2=(-2.d0*r*sqrt(tt1)*qq1*dd0-vv1*ww1)/uu1
            sft=3.d0*sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=dd3*sqrt(pi)*gw**3*sft**2
            tg3=(-2.d0*r*sqrt(tt1)*qq1*dd0-vv1*ww1)/uu1
            sft=4.d0*sf
            tt1=1.d0+sft**2*gw**2
            qq1=1.d0+2.d0*sft**2*gw**2
            ww1=(1.d0-2.d0*r**2*sft**2+gw**2*sft**2)*erf(r/(gw*sqrt(tt1)))
            uu1=sqrt(pi)*r*gw**2*tt1**2.5d0
            vv1=dd4*sqrt(pi)*gw**3*sft**2
            tg4=(-2.d0*r*sqrt(tt1)*qq1*dd0-vv1*ww1)/uu1
        endif
        pot(ix,iy,iz)=pot(ix,iy,iz)+q*(ss0-4.d0*ss1+6.d0*ss2-4.d0*ss3+ss4)
        vgrad(ix,iy,iz)=q*(tg0-4.d0*tg1+6.d0*tg2-4.d0*tg3+tg4)
    enddo
    enddo
    enddo
end subroutine cal_pot_gauss_s
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
subroutine get_poisson_ref(parini,poisson,atoms)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, update_rat, atom_deallocate_old
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    !local variables
    type(typ_poisson):: poisson_ion
    integer:: igpx, igpy, igpz, iat
    integer:: jgpx, jgpy, jgpz
    real(8):: rgcut_a, qtot, pi
    real(8):: ehartree_scn_excl, tt1
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
    call cube_read('rho.cube',atoms,poisson)
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Li') atoms%zat(iat)=3.d0
        if(trim(atoms%sat(iat))=='S' ) atoms%zat(iat)=6.d0
    enddo
    !-------------------------------------------------------
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        rho(igpx,igpy,igpz)=poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    write(*,*) allocated(poisson%rho)
    call f_free(poisson%rho)
    nex=60
    ney=60
    nez=60
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
    gausswidth=0.5d0 !parini%gaussian_width !TO_BE_CORRECTED
    poisson%cal_scn=parini%cal_scn
    poisson%screening_factor=parini%screening_factor
    poisson%task_finit=""
    call init_hartree(parini,atoms,poisson,gausswidth)
    poisson%cell(1)=poisson%hgrid(1,1)*poisson%ngpx
    poisson%cell(2)=poisson%hgrid(2,2)*poisson%ngpy
    poisson%cell(3)=poisson%hgrid(3,3)*poisson%ngpz
    !-------------------------------------------------------
    poisson_ion%alpha=0.5d0 !parini%gaussian_width !TO_BE_CORRECTED
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
    gwt=gausswidth(iat)
    coeff=atoms%zat(iat)/(gwt**3*pi**1.5d0)
    jgpx=int(poisson_ion%rcart(1,iat)/poisson%hgrid(1,1))
    jgpy=int(poisson_ion%rcart(2,iat)/poisson%hgrid(2,2))
    jgpz=int(poisson_ion%rcart(3,iat)/poisson%hgrid(3,3))
    do igpz=jgpz-nbgpz,jgpz+nbgpz
    do igpy=jgpy-nbgpy,jgpy+nbgpy
    do igpx=jgpx-nbgpx,jgpx+nbgpx
        dx=(igpx-1)*poisson%hgrid(1,1)-poisson_ion%rcart(1,iat)
        dy=(igpy-1)*poisson%hgrid(2,2)-poisson_ion%rcart(2,iat)
        dz=(igpz-1)*poisson%hgrid(3,3)-poisson_ion%rcart(3,iat)
        r2=dx**2+dy**2+dz**2
        if(r2<10.d0**2*gwt**2) then
            poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*exp(-r2/gwt**2)
        endif
    enddo
    enddo
    enddo
    enddo
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
    atoms%dpm(1)=0.d0
    atoms%dpm(2)=0.d0
    atoms%dpm(3)=0.d0
    do igpz=1,poisson%ngpz
    do igpy=1,poisson%ngpy
    do igpx=1,poisson%ngpx
        x=poisson%xyz111(1)+(igpx-1)*poisson%hgrid(1,1)
        y=poisson%xyz111(2)+(igpy-1)*poisson%hgrid(2,2)
        z=poisson%xyz111(3)+(igpz-1)*poisson%hgrid(3,3)
        atoms%dpm(1)=atoms%dpm(1)+x*poisson%rho(igpx,igpy,igpz)
        atoms%dpm(2)=atoms%dpm(2)+y*poisson%rho(igpx,igpy,igpz)
        atoms%dpm(3)=atoms%dpm(3)+z*poisson%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    atoms%dpm(1)=atoms%dpm(1)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    atoms%dpm(2)=atoms%dpm(2)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    atoms%dpm(3)=atoms%dpm(3)*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    write(*,'(a,3f8.3)') 'DPM= ',atoms%dpm(1),atoms%dpm(2),atoms%dpm(3)
    !-------------------------------------------------------
    call update_ratp(atoms)
    call get_hartree(parini,poisson,atoms,gausswidth,ehartree_scn_excl)
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,es24.10,es14.5)') 'ehartree_scn_excl ',ehartree_scn_excl,poisson%screening_factor
    endif
    !-------------------------------------------------------
    if(.true.) then
    nd=3
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
subroutine get_proc_stake(mpi_env,n,is,ie)
    use mod_flm_futile
    implicit none
    type(mpi_environment), intent(in):: mpi_env
    integer, intent(in):: n
    integer, intent(out):: is, ie
    !local variables
    integer:: m, mproc
    if(mpi_env%nproc>1) then
        m=n/mpi_env%nproc
        is=mpi_env%iproc*m+1
        mproc=mod(n,mpi_env%nproc)
        is=is+max(0,mpi_env%iproc-mpi_env%nproc+mproc)
        if(mpi_env%iproc>mpi_env%nproc-mproc-1) m=m+1
        ie=is+m-1
    else
        is=1
        ie=n
    endif
end subroutine get_proc_stake
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
