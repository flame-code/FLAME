!*****************************************************************************************
module mod_fit_bf_cent2
    implicit none
    private
    public:: get_basis_functions_cent2
    type typ_fitpar
        private
        integer:: ntypat
        real(8), allocatable:: gwc_s1(:)
        real(8), allocatable:: gwc_s2(:)
        real(8), allocatable:: gw_s1(:)
        real(8), allocatable:: gw_s2(:)
        real(8), allocatable:: gw_p1(:)
        real(8), allocatable:: gw_p2(:)
        real(8), allocatable:: gwgrad_s1(:)
        real(8), allocatable:: gwgrad_s2(:)
        real(8), allocatable:: gwgrad_p1(:)
        real(8), allocatable:: gwgrad_p2(:)
        real(8), allocatable:: qcore_type(:)
        contains
        procedure, private, pass(self):: init_fit_bf
        procedure, private, pass(self):: fini_fit_bf
        !procedure, private, pass(self):: set_param
        procedure, private, pass(self):: report_fit_bf
    end type typ_fitpar
    type typ_bf
        private
        integer:: nbf
        !real(8), allocatable:: 
        integer, allocatable:: iat_list(:)
        real(8), allocatable:: pot_s1(:,:,:)
        real(8), allocatable:: pot_s2(:,:,:)
        real(8), allocatable:: pot_p1(:,:,:)
        real(8), allocatable:: pot_p2(:,:,:)
        character(1), allocatable:: bt(:)
        contains
        procedure, private, pass(self):: init_bf
        procedure, private, pass(self):: fini_bf
        procedure, private, pass(self):: set_bt
        procedure, private, pass(self):: get_pot_iat
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
    allocate(self%pot_s1(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%pot_s2(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%pot_p1(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(self%pot_p2(poisson%ngpx,poisson%ngpy,poisson%ngpz))
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
    deallocate(self%pot_s1)
    deallocate(self%pot_s2)
    deallocate(self%pot_p1)
    deallocate(self%pot_p2)
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
subroutine get_pot_iat(self,parini,poisson,atoms,fitpar,ibf,wa1)
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
    real(8), intent(out):: wa1(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    !local variables
    integer:: iat, itypat
    real(8):: bc_s1, bc_s2, bc_p1, bc_p2
    real(8):: q_tmp(1), p_tmp(3)
    iat=self%iat_list(ibf)
    itypat=atoms%itypat(iat)
    if(self%bt(ibf)=='s') then
        bc_s1= fitpar%gw_s1(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
        bc_s2=-fitpar%gw_s2(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
        q_tmp(1)=1.d0
        call cal_pot_gauss_s(parini,poisson,self%pot_s1,.true.,atoms%ratp(1,iat),fitpar%gw_s1(itypat),q_tmp(1),wa1)
        call cal_pot_gauss_s(parini,poisson,self%pot_s2,.true.,atoms%ratp(1,iat),fitpar%gw_s2(itypat),q_tmp(1),wa1)
        wa1=bc_s1*self%pot_s1+bc_s2*self%pot_s2
    elseif(self%bt(ibf)=='p') then
        bc_p1= fitpar%gw_p1(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
        bc_p2=-fitpar%gw_p2(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
        p_tmp=0.d0
        p_tmp(1)=1.d0
        call cal_pot_gauss_p(parini,poisson,self%pot_p1,.true.,atoms%ratp(1,iat),fitpar%gw_p1(itypat),p_tmp(1),wa1)
        call cal_pot_gauss_p(parini,poisson,self%pot_p2,.true.,atoms%ratp(1,iat),fitpar%gw_p2(itypat),p_tmp(1),wa1)
        wa1=bc_p1*self%pot_p1+bc_p2*self%pot_p2
    else
            stop 'ERROR: unknonbt in get_pot_iat'
    endif
end subroutine get_pot_iat
!*****************************************************************************************
subroutine init_fit_bf(self,ntypat)
    implicit none
    class(typ_fitpar), intent(inout):: self
    integer, intent(in):: ntypat
    !local variables
    self%ntypat=ntypat
    allocate(self%gw_s1(ntypat))
    allocate(self%gw_s2(ntypat))
    allocate(self%gwc_s1(ntypat))
    allocate(self%gwc_s2(ntypat))
    allocate(self%gw_p1(ntypat))
    allocate(self%gw_p2(ntypat))
    allocate(self%gwgrad_s1(ntypat))
    allocate(self%gwgrad_s2(ntypat))
    allocate(self%gwgrad_p1(ntypat))
    allocate(self%gwgrad_p2(ntypat))
    allocate(self%qcore_type(ntypat))
end subroutine init_fit_bf
!*****************************************************************************************
subroutine fini_fit_bf(self)
    implicit none
    class(typ_fitpar), intent(inout):: self
    !local variables
    deallocate(self%gw_s1)
    deallocate(self%gw_s2)
    deallocate(self%gwc_s1)
    deallocate(self%gwc_s2)
    deallocate(self%gw_p1)
    deallocate(self%gw_p2)
    deallocate(self%gwgrad_s1)
    deallocate(self%gwgrad_s2)
    deallocate(self%gwgrad_p1)
    deallocate(self%gwgrad_p2)
    deallocate(self%qcore_type)
end subroutine fini_fit_bf
!*****************************************************************************************
subroutine report_fit_bf(self,istep,cost,nat,c_s,c_p,cgrad_s,cgrad_p)
    implicit none
    class(typ_fitpar), intent(inout):: self
    integer, intent(in):: istep, nat
    real(8), intent(in):: cost, c_s(nat) ,c_p(nat) ,cgrad_s(nat) ,cgrad_p(nat)
    !local variables
    integer:: itypat, iat
    real(8):: tt1, tt2, tt3, tt4
    write(*,'(a,i6,es14.5)',advance='no') 'cost ',istep,cost
    do itypat=1,self%ntypat
        tt1=self%gwgrad_s1(itypat)
        tt2=self%gwgrad_s2(itypat)
        tt3=self%gwgrad_p1(itypat)
        tt4=self%gwgrad_p2(itypat)
        write(*,'(4es10.1)',advance='no') tt1,tt2,tt3,tt4
    enddo
    write(*,'(2es10.1)') maxval(abs(cgrad_s(1:nat))),maxval(abs(cgrad_p(1:nat)))
    write(*,'(a,i6)',advance='no') 'gw ',istep
    do itypat=1,self%ntypat
        tt1=self%gw_s1(itypat)
        tt2=self%gw_s2(itypat)
        tt3=self%gw_p1(itypat)
        tt4=self%gw_p2(itypat)
        write(*,'(4f8.3)',advance='no') tt1,tt2,tt3,tt4
    enddo
    write(*,*)
    write(*,'(a,i6)',advance='no') 'qp ',istep
    do iat=1,nat
        write(*,'(f8.3)',advance='no') c_s(iat)
    enddo
    write(*,'(f8.3)') sum(c_p(:))
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
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_poisson):: poisson, poisson_ref
    type(typ_atoms):: atoms
    type(typ_atoms_arr):: atoms_arr
    !type(typ_ann_arr):: ann_arr
    type(typ_fitpar):: fitpar
    type(typ_bf):: bf
    integer:: iat, iconf, istep
    integer:: ix, iy, iz, itypat
    !real(8):: ehartree_scn_excl
    !real(8):: dpx, dpy, dpz
    !real(8):: time1, time2, time3, time4, time5, time6, time7, time8, time9
    real(8):: tt0, tt1, tt2, tt3, tt4, tt5, tt6, tt7, ts1, ts2, tp1, tp2, cost
    real(8):: tts, ttp
    real(8):: gausswidth(1), q_tmp(1), p_tmp(3), voxel
    real(8):: x, y, z !, r
    real(8):: errmax, rmse, alpha, ener
    real(8):: bc_s1, bc_s2, bc_p1, bc_p2
    real(8):: bcg_s1(2), bcg_s2(2), bcg_p1(2), bcg_p2(2), gradtot
    real(8), allocatable:: c_s(:)
    real(8), allocatable:: c_p(:)
    real(8), allocatable:: cgrad_s(:)
    real(8), allocatable:: cgrad_p(:)
    real(8), allocatable:: vgrad_s1_t(:,:,:)
    real(8), allocatable:: vgrad_s2_t(:,:,:)
    real(8), allocatable:: vgrad_p1_t(:,:,:)
    real(8), allocatable:: vgrad_p2_t(:,:,:)
    real(8), allocatable:: vgrad_s1(:,:,:,:)
    real(8), allocatable:: vgrad_s2(:,:,:,:)
    real(8), allocatable:: vgrad_p1(:,:,:,:)
    real(8), allocatable:: vgrad_p2(:,:,:,:)
    !real(8), allocatable:: pot_iat_s(:,:,:,:)
    !real(8), allocatable:: pot_iat_p(:,:,:,:)
    real(8), allocatable:: pot_ion(:,:,:)
    real(8), allocatable:: qcore(:)
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
    call bf%init_bf(atoms%nat,poisson)
    allocate(qcore(atoms%nat))
    allocate(c_s(atoms%nat))
    allocate(c_p(atoms%nat))
    allocate(cgrad_s(atoms%nat))
    allocate(cgrad_p(atoms%nat))
    allocate(vgrad_s1_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_s2_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_p1_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_p2_t(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(vgrad_s1(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat))
    allocate(vgrad_s2(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat))
    allocate(vgrad_p1(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat))
    allocate(vgrad_p2(poisson%ngpx,poisson%ngpy,poisson%ngpz,parini%ntypat))
    !allocate(pot_iat_s(poisson%ngpx,poisson%ngpy,poisson%ngpz,atoms%nat))
    !allocate(pot_iat_p(poisson%ngpx,poisson%ngpy,poisson%ngpz,atoms%nat))
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
        read(2134,*) tt1,tt2,tt3,tt4,tt5,tt6,tt7
        fitpar%qcore_type(itypat)=tt1
        fitpar%gwc_s1(itypat)=tt2
        fitpar%gwc_s2(itypat)=tt3
        fitpar%gw_s1(itypat)=tt4
        fitpar%gw_s2(itypat)=tt5
        fitpar%gw_p1(itypat)=tt6
        fitpar%gw_p2(itypat)=tt7
    enddo
    close(2134)

    !do iat=1,atoms%nat
    !    if(trim(atoms%sat(iat))=='O') atoms%zat(iat)=4.d0
    !enddo

    !gausswidth=0.5d0
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        qcore(iat)=fitpar%qcore_type(itypat)
    enddo
    pot_ion=0.d0
    do iat=1,atoms%nat
    itypat=atoms%itypat(iat)
    call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),gausswidth(1),atoms%zat(iat),vgrad_s1_t)
    if(qcore(iat)<0.d0) then
    bc_s1= fitpar%gwc_s1(itypat)**3/(fitpar%gwc_s1(itypat)**3-fitpar%gwc_s2(itypat)**3)
    bc_s2=-fitpar%gwc_s2(itypat)**3/(fitpar%gwc_s1(itypat)**3-fitpar%gwc_s2(itypat)**3)
    q_tmp(1)=qcore(iat) !*bc_s1
    call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),fitpar%gwc_s1(itypat),q_tmp(1),vgrad_s1_t)
    !q_tmp(1)=qcore(iat)*bc_s2
    !call cal_pot_gauss_s(parini,poisson,pot_ion,.false.,atoms%ratp(1,iat),fitpar%gwc_s2(itypat),q_tmp(1),vgrad_s1_t)
    endif
    enddo

!    do itypat=1,parini%ntypat
!        if(trim(parini%stypat(itypat))=='Mg') then
!            gw_s1(itypat)=2.464d0
!            gw_s2(itypat)=1.344d0
!            gw_p1(itypat)=1.771d0
!            gw_p2(itypat)=1.450d0
!            !2.683   1.092   3.012   1.117
!            !2.631   1.130   2.746   1.110
!            !2.598   1.254   2.719   1.143
!        elseif(trim(parini%stypat(itypat))=='O') then
!            gw_s1(itypat)=1.346d0
!            gw_s2(itypat)=1.289d0
!            gw_p1(itypat)=1.137d0
!            gw_p2(itypat)=1.016d0
!            !1.707   1.217   1.878   0.995
!            !1.552   1.418   1.834   0.984
!            !1.459   1.273   1.826   0.982
!        else
!            stop 'ERROR: unknown type in get_basis_functions_cent2'
!        endif
!    enddo

    !q_tmp=-2.d0*bc_s1
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),q_tmp,gw_s1, &
    !    6.d0*gw_s1,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !q_tmp=-2.d0*bc_s2
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),q_tmp,gw_s2, &
    !    6.d0*gw_s2,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)
    !q_tmp=-2.d0*b3
    !call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),q_tmp,a3, &
    !    6.d0*a3,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
    !    poisson%hgrid,poisson%rho)

    !p_tmp=0.d0
    !p_tmp(1)=0.3634384d0 !1.d0
    !call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,1),p_tmp,gausswidth, &
    !    6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)

    !call get_hartree(parini,poisson,atoms,gausswidth,ehartree_scn_excl)
    
!    q_tmp(1)=1.d0
!    call cal_pot_gauss_s(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,q_tmp(1),vgrad_s1)
!    pot_s1=poisson%pot
!    gw_s1=gw_s1+1.d-2
!    call cal_pot_gauss_s(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,q_tmp(1),vgrad_s2)
!    pot_s2=poisson%pot
!    !do iz=1,poisson%ngpz
!    do iz=poisson%ngpz/2,poisson%ngpz/2
!    do iy=poisson%ngpy/2,poisson%ngpy/2
!    !do ix=poisson%ngpx/2,poisson%ngpx/2
!    do ix=1,poisson%ngpx
!        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
!        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
!        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
!        write(83,'(3f8.3,2es19.10)') x,y,z,(pot_s2(ix,iy,iz)-pot_s1(ix,iy,iz))/1.d-2,vgrad_s1(ix,iy,iz)
!    enddo
!    enddo
!    enddo
!    stop 'TTTTTTTTTTTTTTTTTTTTT'
    
!    p_tmp=0.d0
!    p_tmp(1)=1.d0
!    call cal_pot_gauss_p(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,p_tmp(1),vgrad_p1)
!    pot_p1=poisson%pot
!    gw_s1=gw_s1+1.d-2
!    call cal_pot_gauss_p(parini,poisson,.true.,atoms%ratp(1,1),gw_s1,p_tmp(1),vgrad_p2)
!    pot_p2=poisson%pot
!    !do iz=1,poisson%ngpz
!    do iz=poisson%ngpz/2,poisson%ngpz/2
!    do iy=poisson%ngpy/2,poisson%ngpy/2
!    !do ix=poisson%ngpx/2,poisson%ngpx/2
!    do ix=1,poisson%ngpx
!        x=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)
!        y=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)
!        z=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)
!        write(83,'(3f8.3,2es19.10)') x,y,z,(pot_p2(ix,iy,iz)-pot_p1(ix,iy,iz))/1.d-2,vgrad_p1(ix,iy,iz)
!    enddo
!    enddo
!    enddo
!    stop 'TTTTTTTTTTTTTTTTTTTTT'




!    write(*,*) 'atoms%nat= ',atoms%nat
!    write(*,*) 'parini%ntypat= ',parini%ntypat
!    write(*,*) 'atoms%itypat= ',atoms%itypat
!    stop 'MMMMMMMMMMMMMMMMMMMMMMM'

    do istep=0,parini%paropt_geopt%nit

    poisson%pot=pot_ion

    voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)

    do iat=1,atoms%nat
    itypat=atoms%itypat(iat)
    bc_s1= fitpar%gw_s1(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
    bc_s2=-fitpar%gw_s2(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
    bc_p1= fitpar%gw_p1(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
    bc_p2=-fitpar%gw_p2(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
    q_tmp(1)=1.d0
    call cal_pot_gauss_s(parini,poisson,bf%pot_s1,.true.,atoms%ratp(1,iat),fitpar%gw_s1(itypat),q_tmp(1),vgrad_s1_t)
    call cal_pot_gauss_s(parini,poisson,bf%pot_s2,.true.,atoms%ratp(1,iat),fitpar%gw_s2(itypat),q_tmp(1),vgrad_s2_t)
    p_tmp=0.d0
    p_tmp(1)=1.d0
    call cal_pot_gauss_p(parini,poisson,bf%pot_p1,.true.,atoms%ratp(1,iat),fitpar%gw_p1(itypat),p_tmp(1),vgrad_p1_t)
    call cal_pot_gauss_p(parini,poisson,bf%pot_p2,.true.,atoms%ratp(1,iat),fitpar%gw_p2(itypat),p_tmp(1),vgrad_p2_t)
    enddo !end of loop over iat
    !write(*,*) 'ALIREZA-1'
    call get_linearcoeff(parini,fitpar,bf,poisson,poisson_ref,atoms, &
        vgrad_s1_t,vgrad_s2_t,qcore,c_s,c_p)
    !write(*,*) 'ALIREZA-2'
    !----------------------------------------------------------
    vgrad_s1=0.d0
    vgrad_s2=0.d0
    vgrad_p1=0.d0
    vgrad_p2=0.d0
    do iat=1,atoms%nat
        !write(*,*) 'ALIREZA-3',iat
        itypat=atoms%itypat(iat)
        bc_s1= fitpar%gw_s1(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
        bc_s2=-fitpar%gw_s2(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
        bc_p1= fitpar%gw_p1(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
        bc_p2=-fitpar%gw_p2(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
        q_tmp(1)=1.d0
        call cal_pot_gauss_s(parini,poisson,bf%pot_s1,.true.,atoms%ratp(1,iat),fitpar%gw_s1(itypat),q_tmp(1),vgrad_s1_t)
        call cal_pot_gauss_s(parini,poisson,bf%pot_s2,.true.,atoms%ratp(1,iat),fitpar%gw_s2(itypat),q_tmp(1),vgrad_s2_t)
        p_tmp=0.d0
        p_tmp(1)=1.d0
        call cal_pot_gauss_p(parini,poisson,bf%pot_p1,.true.,atoms%ratp(1,iat),fitpar%gw_p1(itypat),p_tmp(1),vgrad_p1_t)
        call cal_pot_gauss_p(parini,poisson,bf%pot_p2,.true.,atoms%ratp(1,iat),fitpar%gw_p2(itypat),p_tmp(1),vgrad_p2_t)
        bcg_s1(1)=-3.d0*fitpar%gw_s1(itypat)**2*fitpar%gw_s2(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)**2
        bcg_s1(2)= 3.d0*fitpar%gw_s1(itypat)**3*fitpar%gw_s2(itypat)**2/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)**2
        bcg_s2(1)= 3.d0*fitpar%gw_s2(itypat)**3*fitpar%gw_s1(itypat)**2/(fitpar%gw_s2(itypat)**3-fitpar%gw_s1(itypat)**3)**2
        bcg_s2(2)=-3.d0*fitpar%gw_s2(itypat)**2*fitpar%gw_s1(itypat)**3/(fitpar%gw_s2(itypat)**3-fitpar%gw_s1(itypat)**3)**2
        bcg_p1(1)=-5.d0*fitpar%gw_p1(itypat)**4*fitpar%gw_p2(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)**2
        bcg_p1(2)= 5.d0*fitpar%gw_p1(itypat)**5*fitpar%gw_p2(itypat)**4/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)**2
        bcg_p2(1)= 5.d0*fitpar%gw_p2(itypat)**5*fitpar%gw_p1(itypat)**4/(fitpar%gw_p2(itypat)**5-fitpar%gw_p1(itypat)**5)**2
        bcg_p2(2)=-5.d0*fitpar%gw_p2(itypat)**4*fitpar%gw_p1(itypat)**5/(fitpar%gw_p2(itypat)**5-fitpar%gw_p1(itypat)**5)**2
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tts=bc_s1*bf%pot_s1(ix,iy,iz)+bc_s2*bf%pot_s2(ix,iy,iz)
            ttp=bc_p1*bf%pot_p1(ix,iy,iz)+bc_p2*bf%pot_p2(ix,iy,iz)
            poisson%pot(ix,iy,iz)=poisson%pot(ix,iy,iz)+c_s(iat)*tts+c_p(iat)*ttp
            tt1=bcg_s1(1)*bf%pot_s1(ix,iy,iz)+bcg_s2(1)*bf%pot_s2(ix,iy,iz)+bc_s1*vgrad_s1_t(ix,iy,iz)
            tt2=bcg_s1(2)*bf%pot_s1(ix,iy,iz)+bcg_s2(2)*bf%pot_s2(ix,iy,iz)+bc_s2*vgrad_s2_t(ix,iy,iz)
            tt3=bcg_p1(1)*bf%pot_p1(ix,iy,iz)+bcg_p2(1)*bf%pot_p2(ix,iy,iz)+bc_p1*vgrad_p1_t(ix,iy,iz)
            tt4=bcg_p1(2)*bf%pot_p1(ix,iy,iz)+bcg_p2(2)*bf%pot_p2(ix,iy,iz)+bc_p2*vgrad_p2_t(ix,iy,iz)
            vgrad_s1(ix,iy,iz,itypat)=vgrad_s1(ix,iy,iz,itypat)+tt1*c_s(iat)
            vgrad_s2(ix,iy,iz,itypat)=vgrad_s2(ix,iy,iz,itypat)+tt2*c_s(iat)
            vgrad_p1(ix,iy,iz,itypat)=vgrad_p1(ix,iy,iz,itypat)+tt3*c_p(iat)
            vgrad_p2(ix,iy,iz,itypat)=vgrad_p2(ix,iy,iz,itypat)+tt4*c_p(iat)
        enddo
        enddo
        enddo
    enddo
    cost=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        tt0=poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)
        cost=cost+tt0**2
    enddo
    enddo
    enddo
    cost=cost*voxel
    do itypat=1,parini%ntypat
        ts1=0.d0
        ts2=0.d0
        tp1=0.d0
        tp2=0.d0
        do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
        do ix=1,poisson%ngpx
            tt0=poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz)
            ts1=ts1+vgrad_s1(ix,iy,iz,itypat)*tt0
            ts2=ts2+vgrad_s2(ix,iy,iz,itypat)*tt0
            tp1=tp1+vgrad_p1(ix,iy,iz,itypat)*tt0
            tp2=tp2+vgrad_p2(ix,iy,iz,itypat)*tt0
        enddo
        enddo
        enddo
        fitpar%gwgrad_s1(itypat)=2.d0*ts1*voxel
        fitpar%gwgrad_s2(itypat)=2.d0*ts2*voxel
        fitpar%gwgrad_p1(itypat)=2.d0*tp1*voxel
        fitpar%gwgrad_p2(itypat)=2.d0*tp2*voxel
    enddo

    call fitpar%report_fit_bf(istep,cost,atoms%nat,c_s,c_p,cgrad_s,cgrad_p)


    !alpha=4.d-4
    !if(abs(gwgrad_s1(1))<1.d-1 .and. abs(gwgrad_s2(1))<1.d-1) alpha=4.d-2
    fitpar%gw_s1=fitpar%gw_s1- 1.d0*alpha*fitpar%gwgrad_s1
    fitpar%gw_s2=fitpar%gw_s2- 1.d0*alpha*fitpar%gwgrad_s2
    fitpar%gw_p1=fitpar%gw_p1-10.d0*alpha*fitpar%gwgrad_p1
    fitpar%gw_p2=fitpar%gw_p2-10.d0*alpha*fitpar%gwgrad_p2
    !c_s=c_s-1.d-1*cgrad_s
    !c_p=c_p-1.d-1*cgrad_p
    !write(*,*) 'qtot= ',sum(c_s)

    !bc_s1= gw_s1**3/(gw_s1**3-gw_s2**3)
    !bc_s2=-gw_s2**3/(gw_s1**3-gw_s2**3)
    !bc_p1= gw_p1**5/(gw_p1**5-gw_p2**5)
    !bc_p2=-gw_p2**5/(gw_p1**5-gw_p2**5)


    enddo !end of loop over istep


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
        tt4=fitpar%gw_s1(itypat)
        tt5=fitpar%gw_s2(itypat)
        tt6=fitpar%gw_p1(itypat)
        tt7=fitpar%gw_p2(itypat)
        write(2134,'(7f10.6)') tt1,tt2,tt3,tt4,tt5,tt6,tt7
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
    poisson%rho=0.d0
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),atoms%zat(iat),gausswidth(1), &
        6.d0*gausswidth(1),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    if(qcore(iat)<0.d0) then
    bc_s1= fitpar%gwc_s1(itypat)**3/(fitpar%gwc_s1(itypat)**3-fitpar%gwc_s2(itypat)**3)
    bc_s2=-fitpar%gwc_s2(itypat)**3/(fitpar%gwc_s1(itypat)**3-fitpar%gwc_s2(itypat)**3)
    q_tmp(1)=qcore(iat)*bc_s1
    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwc_s1(itypat), &
        6.d0*fitpar%gwc_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    q_tmp(1)=qcore(iat)*bc_s2
    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gwc_s2(itypat), &
        6.d0*fitpar%gwc_s2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    endif
    bc_s1= fitpar%gw_s1(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
    bc_s2=-fitpar%gw_s2(itypat)**3/(fitpar%gw_s1(itypat)**3-fitpar%gw_s2(itypat)**3)
    q_tmp(1)=c_s(iat)*bc_s1
    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gw_s1(itypat), &
        4.d0*fitpar%gw_s1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    q_tmp(1)=c_s(iat)*bc_s2
    call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,fitpar%gw_s2(itypat), &
        4.d0*fitpar%gw_s2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%rho)
    bc_p1= fitpar%gw_p1(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
    bc_p2=-fitpar%gw_p2(itypat)**5/(fitpar%gw_p1(itypat)**5-fitpar%gw_p2(itypat)**5)
    p_tmp=0.d0
    p_tmp(1)=c_p(iat)*bc_p1
    call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),p_tmp,fitpar%gw_p1(itypat), &
        6.d0*fitpar%gw_p1(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    p_tmp(1)=c_p(iat)*bc_p2
    call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),p_tmp,fitpar%gw_p2(itypat), &
        6.d0*fitpar%gw_p2(itypat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    enddo
    !poisson%rho=poisson%rho-poisson_ref%rho
    !call cube_write('diffrho.cube',atoms,poisson,'rho')
    !call get_hartree(parini,poisson,atoms,gausswidth,ener)
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
    !----------------------------------------------------------


    !do iconf=1,atoms_arr%nconf
    !    call atom_deallocate_old(atoms_arr%atoms(iconf))
    !enddo
    !deallocate(atoms_arr%atoms)
    call fini_hartree(parini,atoms,poisson_ref)
    call fini_hartree(parini,atoms,poisson)
    call fitpar%fini_fit_bf()
    call bf%fini_bf()

    call f_release_routine()
end subroutine get_basis_functions_cent2
!*****************************************************************************************
subroutine get_linearcoeff(parini,fitpar,bf,poisson,poisson_ref,atoms,wa1,wa2,qcore,c_s,c_p)
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_fitpar), intent(inout):: fitpar
    type(typ_bf), intent(inout):: bf
    type(typ_poisson), intent(in):: poisson, poisson_ref
    type(typ_atoms), intent(in):: atoms
    real(8), intent(inout):: wa1(poisson_ref%ngpx,poisson_ref%ngpy,poisson_ref%ngpz)
    real(8), intent(inout):: wa2(poisson_ref%ngpx,poisson_ref%ngpy,poisson_ref%ngpz)
    real(8), intent(in):: qcore(atoms%nat)
    real(8), intent(out):: c_s(atoms%nat), c_p(atoms%nat)
    !local variables
    real(8):: bc_s1, bc_s2, bc_p1, bc_p2
    real(8):: q_tmp(1), p_tmp(3)
    real(8):: voxel, dpm(3)
    real(8):: tt1, tt2, tt3, tt4
    integer:: ix, iy, iz, iat, jat, info, nbf, nc, itypat, jtypat
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
        call bf%get_pot_iat(parini,poisson,atoms,fitpar,ibf,wa1)
        do jbf=ibf,bf%nbf
            call bf%get_pot_iat(parini,poisson,atoms,fitpar,jbf,wa2)
            tt1=0.d0
            do iz=1,poisson_ref%ngpz
            do iy=1,poisson_ref%ngpy
            do ix=1,poisson_ref%ngpx
                tt1=tt1+wa1(ix,iy,iz)*wa2(ix,iy,iz)
            enddo
            enddo
            enddo
            secder(ibf,jbf)=tt1*voxel
            secder(jbf,ibf)=tt1*voxel
        enddo
        !-------------------------------------------------------------
        tt1=0.d0
        do iz=1,poisson_ref%ngpz
        do iy=1,poisson_ref%ngpy
        do ix=1,poisson_ref%ngpx
            tt1=tt1+wa1(ix,iy,iz)*(poisson%pot(ix,iy,iz)-poisson_ref%pot(ix,iy,iz))
        enddo
        enddo
        enddo
        rhs(ibf)=-tt1*voxel
        !-------------------------------------------------------------
    enddo
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
    rhs(nbf+1)=-sum(atoms%zat)-sum(qcore) !-sum(atoms%zat)
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
    real(8):: dx, dy, dz, r2, coeff, gwt
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
    if(parini%mpi_env%iproc==0) then
        write(*,*) 'qtot= ',qtot
    endif
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
end module mod_fit_bf_cent2
!*****************************************************************************************
