!*****************************************************************************************
module mod_cent2
    implicit none
    private
    public:: typ_cent2, cent2_analytic, get_cost_secder
    public:: cal_rho_pot_integral_local
    type typ_bf
        private
        integer:: nbf
        integer:: nbgx, nbgy, nbgz
        real(8), allocatable:: gwn(:)
        real(8), allocatable:: gwz(:)
        real(8), allocatable:: bz(:)
        real(8), allocatable:: gwc(:)
        real(8), allocatable:: bc(:)
        real(8), allocatable:: gwe_s(:)
        real(8), allocatable:: be_s(:)
        real(8), allocatable:: gwe_p(:,:)
        real(8), allocatable:: qcore(:)
        real(8), allocatable:: orb(:)
        real(8), allocatable:: hardness(:)
        real(8), allocatable:: re(:,:)
        real(8), allocatable:: rn(:,:)
        integer, allocatable:: imap(:)
        integer, allocatable:: ibf_list_s(:)
        character(2), allocatable:: bt(:)
        contains
        procedure, private, pass(self):: init_bf
        procedure, private, pass(self):: fini_bf
        procedure, private, pass(self):: set_param
    end type typ_bf
    type typ_cent2
        private
        integer:: nbgx, nbgy, nbgz
        real(8), allocatable:: amat(:,:)
        real(8), allocatable:: amat_t(:,:)
        real(8), allocatable:: cep_rhs(:)
        real(8), allocatable:: rho_n(:,:,:)
        real(8), allocatable:: gausswidth_tmp(:)
        type(typ_bf):: bf
        contains
        !procedure, public, pass(self)::
        procedure, public, pass(self):: init_cent2
        procedure, public, pass(self):: fini_cent2
        procedure, private, pass(self):: init_electrostatic_cent2
        procedure, private, pass(self):: fini_electrostatic_cent2
        procedure, private, pass(self):: get_amat_cent2_analytic
        procedure, private, pass(self):: get_eigenval
        procedure, private, pass(self):: get_expansion_coeff
        procedure, private, pass(self):: get_qat_from_chi_dir_cent2
        procedure, private, pass(self):: get_cep_rhs
        procedure, private, pass(self):: cal_electrostatic_ann_cent2
        procedure, private, pass(self):: cal_cent2_energy
        procedure, private, pass(self):: cent2_ehartree_analytic
        procedure, private, pass(self):: cent2_ehartree_analytic_monopole
        !procedure, private, pass(self):: cent2_ehartree_analytic_test
        procedure, private, pass(self):: reverseCEP
        procedure, private, pass(self):: prefit_cent2
        procedure, public, pass(self):: cal_ann_cent2
        procedure, private, pass(self):: calc_atomic_densities
    end type typ_cent2
contains
!*****************************************************************************************
subroutine cal_ann_cent2(self,parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linkedlists, only: typ_linkedlists
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: get_proc_stake
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_poisson):: poisson
    !local variables
    type(typ_pia_arr):: pia_arr_tmp
    type(typ_linkedlists):: linkedlists
    integer:: iat, i, j, ng
    integer:: iats, iate
    real(8):: epot_c, out_ann
    real(8):: dpx, dpy, dpz
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8, time9
    real(8):: tt1, hinv(3,3), vol
    !integer, save:: icall=0
    call f_routine(id='cal_ann_cent2')
    !icall=icall+1
    !if(parini%mpi_env%iproc==0) then
    !    write(*,*) 'ICALL= ',icall
    !endif
    call update_ratp(atoms)
    call self%init_cent2(parini,ann_arr,atoms,poisson)
    if(parini%iverbose>=2) call cpu_time(time1)
    call symfunc%init_symfunc(parini%mpi_env,parini%iverbose,parini%bondbased_ann,parini%symfunc_type_ann)
    if(ann_arr%compute_symfunc) then
        call symfunc%get_symfunc(ann_arr,atoms,.true.)
    else
        symfunc%linked_lists%rcut=ann_arr%rcut
        symfunc%linked_lists%triplex=.true.
        call linkedlists%calc_linkedlists(atoms,.true.,symfunc%linked_lists,pia_arr_tmp,parini%mpi_env,parini%iverbose,parini%bondbased_ann)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fatpq(1:3,1:symfunc%linked_lists%maxbound_rad))
        allocate(ann_arr%stresspq(1:3,1:3,1:symfunc%linked_lists%maxbound_rad))
    endif
    call get_proc_stake(parini%mpi_env,atoms%nat,iats,iate)
    if(parini%iverbose>=2) call cpu_time(time2)
    over_iat: do iat=iats,iate
        i=atoms%itypat(iat)
        ng=ann_arr%ann(i)%nn(0)
        ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,iat)
        if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
            call cal_architecture(ann_arr%ann(i),out_ann)
            call cal_force_chi_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
        elseif(trim(ann_arr%event)=='train') then
            call cal_architecture_der(ann_arr%ann(i),out_ann)
            ann_arr%chi_i(iat)=out_ann
            tt1=tanh(ann_arr%ann(i)%prefactor_chi*out_ann)
            ann_arr%chi_o(iat)=ann_arr%ann(i)%ampl_chi*tt1+ann_arr%ann(i)%chi0
            call convert_ann_epotd(ann_arr%ann(i),ann_arr%num(i),ann_arr%g_per_atom(1,iat))
            ann_arr%g_per_atom(1:ann_arr%num(1),iat)=ann_arr%g_per_atom(1:ann_arr%num(1),iat)*ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)
        else
            stop 'ERROR: undefined content for ann_arr%event'
        endif
    enddo over_iat
    if(parini%iverbose>=2) call cpu_time(time3)
    call self%calc_atomic_densities(parini,atoms,ann_arr,poisson)
    if(parini%iverbose>=2) call cpu_time(time4)
    if(parini%iverbose>=2) write(*,'(a,f8.2)') 'time: two routines: ',time4-time3
    call self%get_amat_cent2_analytic(parini,ann_arr,atoms)
    self%amat=self%amat_t
    if(parini%prefit_ann) then
        ann_arr%chi_o=0.d0
    endif
    call self%get_cep_rhs(parini,ann_arr,atoms,poisson)
    !do iat=1,(atoms%nat+1)*(atoms%nat+1)
    !    write(91,'(i5,f15.7)') iat,ann_arr%a(iat)
    !enddo
    if(parini%prefit_ann) then
        call self%prefit_cent2(parini,ann_arr,atoms,poisson)
    endif
    call self%get_eigenval(parini,atoms)
    if(parini%iverbose>=2) call cpu_time(time5)
    if(parini%iverbose>=2) write(*,'(a,f8.2)') 'time: prefit_cent2: ',time5-time4
    !This must be here since contribution from coulomb
    !interaction is calculated during the process of charge optimization.
    call self%get_qat_from_chi_dir_cent2(parini,ann_arr,atoms,poisson,self%amat)
    if(parini%iverbose>=2) call cpu_time(time7)
    if(parini%iverbose>=2) write(*,'(a,f8.2)') 'time: get_qat_from_chi_dir_cent2: ',time7-time6
    if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2)  then
        do iat=1,atoms%nat
            write(166,*)iat,ann_arr%chi_o(iat)
            write(99,'(i4,a3,3f6.2)') iat,trim(atoms%sat(iat)),atoms%qat(iat),atoms%zat(iat),atoms%zat(iat)+atoms%qat(iat)
        enddo
        write(99,*) '==='
        write(99,*) sum(atoms%qat),sum(atoms%zat)
        write(99,*) '---'
    endif
    endif
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    call self%cal_cent2_energy(parini,atoms,ann_arr,epot_c,poisson)
    if(parini%iverbose>=2) call cpu_time(time9)
    !if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2) write(*,'(a,f8.2)') 'time: cal_cent2_energy: ',time9-time7
    !endif
    if(parini%iverbose>=2) then
        call yaml_mapping_open('Timing of CENT2')
        call yaml_map('calculation of symfunc',time2-time1)
        call yaml_map('neural network process',time3-time2)
        call yaml_map('init_electrostatic',time4-time3)
        call yaml_map('prefit_cent2',time5-time4)
        call yaml_map('get_amat',time6-time5)
        call yaml_map('get_qat_from_chi_dir',time7-time6)
        call yaml_map('force (SR term)',time9-time7)
        call yaml_map('total time',time9-time1)
        call yaml_mapping_close()
    endif !end of if for printing out timing.
    if(trim(ann_arr%event)=='potential' )then
        atoms%epot=epot_c
    elseif(trim(ann_arr%event)=='evalu') then
        atoms%epot=ann_arr%epot_es
    endif
    call get_dpm(atoms,dpx,dpy,dpz,ann_arr%dpm_err)
    if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2) then
        write(1390,'(3es18.8,a3,3es18.8,a3,es18.8)')dpx,dpy,dpz,' | ',atoms%dpm(1),atoms%dpm(2),atoms%dpm(3),' | ',ann_arr%dpm_err
    endif
    endif
    atoms%dpm(1)=dpx
    atoms%dpm(2)=dpy
    atoms%dpm(3)=dpz
    call self%fini_electrostatic_cent2(parini,ann_arr,atoms,poisson)
    !call repulsive_potential_cent(parini,atoms,ann_arr)
    call getvol_alborz(atoms%cellvec,vol)
    call invertmat_alborz(atoms%cellvec,hinv)
    !The following line is inconsistent with the definition of stress tensor
    atoms%stress(1:3,1:3)=atoms%stress(1:3,1:3)*vol
    do i=1,3
    do j=1,3
        atoms%celldv(i,j)=vol*(atoms%stress(i,1)*hinv(j,1)+atoms%stress(i,2)*hinv(j,2)+atoms%stress(i,3)*hinv(j,3))
    enddo
    enddo

    deallocate(symfunc%linked_lists%prime_bound)
    deallocate(symfunc%linked_lists%bound_rad)
    deallocate(symfunc%linked_lists%bound_ang)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call f_free(symfunc%y)
        call f_free(symfunc%y0d)
        call f_free(symfunc%y0dr)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%chi_i)
        deallocate(ann_arr%chi_o)
        deallocate(ann_arr%chi_d)
        deallocate(ann_arr%fat_chi)
        deallocate(ann_arr%fatpq)
        deallocate(ann_arr%stresspq)
    endif
    deallocate(self%amat)
    deallocate(self%amat_t)
    deallocate(self%cep_rhs)
    deallocate(ann_arr%ipiv)
    call symfunc%fini_symfunc()
    call self%fini_cent2(parini,ann_arr,atoms,poisson)
    call f_release_routine()
end subroutine cal_ann_cent2
!*****************************************************************************************
subroutine init_cent2(self,parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    call self%init_electrostatic_cent2(parini,atoms,ann_arr,poisson)
    if(.not. allocated(ann_arr%ipiv)) then
        allocate(ann_arr%ipiv(1:self%bf%nbf+1))
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fat_chi(1:3,1:atoms%nat))
        allocate(ann_arr%chi_i(self%bf%nbf))
        allocate(ann_arr%chi_o(self%bf%nbf))
        allocate(ann_arr%chi_d(1:atoms%nat))
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
    else
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
    endif
    allocate(self%amat((self%bf%nbf+1),(self%bf%nbf+1)))
    allocate(self%amat_t((self%bf%nbf+1),(self%bf%nbf+1)))
    allocate(self%cep_rhs((self%bf%nbf+1)))
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%qq(1:self%bf%nbf+1))
    endif
end subroutine init_cent2
!*****************************************************************************************
subroutine fini_cent2(self,parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%qq)
    endif
end subroutine fini_cent2
!*****************************************************************************************
subroutine init_bf(self,nbgx,nbgy,nbgz,nat)
    implicit none
    class(typ_bf), intent(inout):: self
    integer, intent(in):: nbgx, nbgy, nbgz, nat
    !local variables
    self%nbf=4*nat
    self%nbgx=nbgx
    self%nbgy=nbgy
    self%nbgz=nbgz
    allocate(self%gwn(nat))
    allocate(self%qcore(nat))
    allocate(self%rn(3,nat))
    allocate(self%hardness(self%nbf))
    allocate(self%imap(self%nbf))
    allocate(self%gwz(nat))
    allocate(self%bz(nat))
    allocate(self%gwc(nat))
    allocate(self%bc(nat))
    allocate(self%gwe_s(nat))
    allocate(self%be_s(nat))
    allocate(self%gwe_p(2,nat))
    allocate(self%orb(self%nbf))
    allocate(self%re(3,self%nbf))
    allocate(self%bt(self%nbf))
    allocate(self%ibf_list_s(nat))
end subroutine init_bf
!*****************************************************************************************
subroutine fini_bf(self)
    implicit none
    class(typ_bf), intent(inout):: self
    !local variables
    deallocate(self%gwn)
    deallocate(self%qcore)
    deallocate(self%gwz)
    deallocate(self%bz)
    deallocate(self%gwc)
    deallocate(self%bc)
    deallocate(self%gwe_s)
    deallocate(self%be_s)
    deallocate(self%gwe_p)
    deallocate(self%orb)
    deallocate(self%hardness)
    deallocate(self%imap)
    deallocate(self%rn)
    deallocate(self%re)
    deallocate(self%bt)
    deallocate(self%ibf_list_s)
end subroutine fini_bf
!*****************************************************************************************
subroutine set_param(self,atoms,ann_arr)
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    implicit none
    class(typ_bf), intent(inout):: self
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(in):: ann_arr
    !local variables
    integer:: iat, ibf
    character(2):: bt_t
    do iat=1,atoms%nat
        self%gwn(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
        self%qcore(iat)=ann_arr%ann(atoms%itypat(iat))%qcore
        self%gwz(iat)=ann_arr%ann(atoms%itypat(iat))%gwz
        self%bz(iat)=ann_arr%ann(atoms%itypat(iat))%bz
        self%gwc(iat)=ann_arr%ann(atoms%itypat(iat))%gwc
        self%bc(iat)=ann_arr%ann(atoms%itypat(iat))%bc
        self%rn(1:3,iat)=atoms%ratp(1:3,iat)
    enddo
    do ibf=1,self%nbf
        iat=modulo(ibf-1,self%nbf/4)+1
        self%imap(ibf)=iat
        if((4*(ibf-1))/self%nbf==0) then
            self%bt(ibf)='s'
            self%ibf_list_s(iat)=ibf
        elseif((4*(ibf-1))/self%nbf==1) then
            self%bt(ibf)='px'
        elseif((4*(ibf-1))/self%nbf==2) then
            self%bt(ibf)='py'
        elseif((4*(ibf-1))/self%nbf==3) then
            self%bt(ibf)='pz'
        else
            stop 'ERROR: unknon ibf in set_bt'
        endif
    enddo
    do ibf=1,self%nbf
        iat=self%imap(ibf)
        bt_t=self%bt(ibf)
        !write(*,*) 'IAT ',iat,ibf
        if(bt_t(1:1)=='s') then
            self%gwe_s(iat)=ann_arr%ann(atoms%itypat(iat))%gwe_s
            self%be_s(iat)=ann_arr%ann(atoms%itypat(iat))%be_s
            self%orb(ibf)=1.d0
            self%hardness(ibf)=ann_arr%ann(atoms%itypat(iat))%hardness
        endif
        if(bt_t(1:1)=='p') then
            self%gwe_p(1,iat)=ann_arr%ann(atoms%itypat(iat))%gwe_p(1)
            self%gwe_p(2,iat)=ann_arr%ann(atoms%itypat(iat))%gwe_p(2)
            self%orb(ibf)=0.d0
            self%hardness(ibf)=1.d-1*ann_arr%ann(atoms%itypat(iat))%hardness
        endif
        self%re(1,ibf)=atoms%ratp(1,iat)
        self%re(2,ibf)=atoms%ratp(2,iat)
        self%re(3,ibf)=atoms%ratp(3,iat)
    enddo
end subroutine set_param
!*****************************************************************************************
subroutine calc_atomic_densities(self,parini,atoms,ann_arr,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linked_lists, only: typ_pia_arr
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, ix, iy, iz, agpx, agpy, agpz, ibf, ii, jbf, itypat
    real(8):: ggw, ggw_t, alpha, beta, q_tmp(1), p_tmp(3), tt, ttx, tty, ttz
    real(8):: a1, a2, a3, b1, b2 !, b3
    !character(2):: bt_t
    poisson%rho(:,:,:)=0.d0
    do iat=1,atoms%nat
        !agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+self%bf%nbgx+0
        !agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+self%bf%nbgy+0
        !agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+self%bf%nbgz+0
!        q_tmp(1)=atoms%zat(iat) !1.d0
!        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,self%bf%gwn(iat), &
!            6.d0*maxval(self%bf%gwn),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
!        if(self%bf%gwc(2,iat)>0.d0) then
            !a1=self%bf%gwc(1,iat)
            !a2=self%bf%gwc(2,iat)
            !b1= a1**3/(a1**3-a2**3)
            !b2=-a2**3/(a1**3-a2**3)
            q_tmp(1)=self%bf%qcore(iat)*self%bf%bc(iat)
            call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,self%bf%gwc(iat), &
                6.d0*self%bf%gwc(iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            q_tmp(1)=self%bf%qcore(iat)*(1.d0-self%bf%bc(iat))
            call put_r2gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,self%bf%gwc(iat), &
                6.d0*self%bf%gwc(iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
!        else
!            a1=self%bf%gwc(1,iat)
!            q_tmp(1)=self%bf%qcore(iat)
!            call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,a1, &
!                6.d0*a1,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
!        endif
    enddo !end of loop over iat
    self%rho_n=poisson%rho !save and keep total ionic charge
end subroutine calc_atomic_densities
!*****************************************************************************************
subroutine init_electrostatic_cent2(self,parini,atoms,ann_arr,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, set_typat
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat
    real(8):: time1, time2
    real(8):: max_cellVec
    !real(8):: time1_t, time2_t, time3_t
    if(trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train') then
        ann_arr%syslinsolver='direct'
    else
        ann_arr%syslinsolver=trim(parini%syslinsolver_ann)
    endif
    ann_arr%ener_ref=0.d0
    if(trim(ann_arr%syslinsolver)=='direct' .or. trim(ann_arr%syslinsolver)=='apply_matrix') then
        if(trim(atoms%boundcond)/='free') then
            write(*,*) 'ERROR: syslinsolver=direct can be used only for free BC.'
        endif
        allocate(self%gausswidth_tmp(1:atoms%nat))
        do iat=1,atoms%nat
            self%gausswidth_tmp(iat)=0.d0
            atoms%zat(iat)=ann_arr%ann(atoms%itypat(iat))%zion
            atoms%qat(iat)=ann_arr%ann(atoms%itypat(iat))%qinit
        enddo
        poisson%nat=atoms%nat
        poisson%bc=atoms%boundcond
        poisson%cv=atoms%cellvec
        poisson%cal_scn=parini%cal_scn
        poisson%screening_factor=parini%screening_factor
        poisson%alpha=parini%alpha_ewald
        poisson%task_finit="alloc_rho:set_ngp"
        if(parini%iverbose>=2) call cpu_time(time1)
        call init_hartree(parini,atoms,poisson,self%gausswidth_tmp)
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'init_hartree_time: ',time2-time1
        max_cellVec=100.d0
        !if(maxval(atoms%cellvec)>max_cellVec) then
        !    stop 'ERROR: atoms%cellvec > max_cellVec'
        !endif
        if(parini%iverbose>=2) call cpu_time(time1)
        if(.not. ann_arr%linear_rho_pot_initiated) then
            ann_arr%linear_rho_pot_initiated=.true. 
        endif
        allocate(self%rho_n(poisson%ngpx,poisson%ngpy,poisson%ngpz))
        self%nbgx=int(poisson%rgcut/poisson%hgrid(1,1))+2
        self%nbgy=int(poisson%rgcut/poisson%hgrid(2,2))+2
        self%nbgz=int(poisson%rgcut/poisson%hgrid(3,3))+2
        call self%bf%init_bf(self%nbgx,self%nbgy,self%nbgz,atoms%nat)
        call self%bf%set_param(atoms,ann_arr)
    else
        write(*,*) 'ERROR: unknown value for syslinsolver',trim(ann_arr%syslinsolver)
        stop
    endif
end subroutine init_electrostatic_cent2
!*****************************************************************************************
subroutine get_amat_cent2_analytic(self,parini,ann_arr,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    !local variables
    integer:: ibf, jbf, iat, jat, ierr
    real(8):: pati(3), patj(3)
    real(8):: e_qr0_qr0, e_qr2_qr2, e_qr0_i_qr2_j, e_qr0_j_qr2_i, e_qr0_i_pr1_j
    real(8):: e_qr2_i_pr1_j, e_qr0_j_pr1_i, e_qr2_j_pr1_i
    real(8):: gwr0_iat, gwr2_iat, gwp1_iat, qr0_iat, qr2_iat, gwr0_jat
    real(8):: gwr2_jat, gwp1_jat, qr0_jat, qr2_jat, gwr0, gwr2, gwp1, qr0, qr2
    real(8):: pi_dot_pi
    real(8):: dx, dy, dz, r, pi, tt
    real(8):: ss1, ss2, ss3, gg1, gg2, gg3, tt1, tt2, tt3
    real(8):: a, b, c, sf_1, sf_2, sf_3, sfs(3)
    character(2):: bt_i, bt_j
    pi=4.d0*atan(1.d0)
    sfs(1)=parini%screening_factor
    sfs(2)=parini%screening_factor*1.1d0
    sfs(3)=parini%screening_factor*1.2d0
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))

    self%amat_t=0.d0
    do ibf=1,self%bf%nbf
        iat=self%bf%imap(ibf)
        !itypat=atoms%itypat(iat)
        bt_i=self%bf%bt(ibf)
        if(trim(self%bf%bt(ibf))=='s') then
        elseif(trim(self%bf%bt(ibf))=='px') then
            pati(1)=1.d0 ; pati(2)=0.d0 ; pati(3)=0.d0
        elseif(trim(self%bf%bt(ibf))=='py') then
            pati(1)=0.d0 ; pati(2)=1.d0 ; pati(3)=0.d0
        elseif(trim(self%bf%bt(ibf))=='pz') then
            pati(1)=0.d0 ; pati(2)=0.d0 ; pati(3)=1.d0
        endif
        gwr0=self%bf%gwe_s(iat)
        gwr2=gwr0
        gwp1=self%bf%gwe_p(1,iat)
        qr0=self%bf%be_s(iat)
        qr2=(1.d0-self%bf%be_s(iat))

        if(bt_i(1:1)=='s') then
        !self-interaction of qr0-qr0 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwr0**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwr0**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwr0**2*sf_3**2)
        tt=0.5d0*2.d0*qr0**2*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of qr2-qr2 interaction
        gg1=gwr2*sf_1
        gg2=gwr2*sf_2
        gg3=gwr2*sf_3
        ss1=(2.d0*sf_1*(3.d0+10.d0*gg1**2+9.d0*gg1**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+10.d0*gg2**2+9.d0*gg2**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+10.d0*gg3**2+9.d0*gg3**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg3**2)**2.5d0)
        tt=tt+0.5d0*qr2**2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qr0-qr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_3**2)**1.5d0)
        tt=tt+1.d0*qr0*qr2*(a*ss1+b*ss2+c*ss3)
        else
        !self-interaction of pr1-pr1 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwp1**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwp1**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwp1**2*sf_3**2)
        pi_dot_pi=pati(1)**2+pati(2)**2+pati(3)**2
        tt=0.5d0*4.d0*pi_dot_pi*(a*tt1**3+b*tt2**3+c*tt3**3)/(3.d0*sqrt(pi))
        endif
        self%amat_t(ibf,ibf)=2.d0*tt
        do jbf=1,ibf-1
            jat=self%bf%imap(jbf)
            bt_j=self%bf%bt(jbf)
            if(trim(self%bf%bt(jbf))=='s') then
            elseif(trim(self%bf%bt(jbf))=='px') then
                patj(1)=1.d0 ; patj(2)=0.d0 ; patj(3)=0.d0
            elseif(trim(self%bf%bt(jbf))=='py') then
                patj(1)=0.d0 ; patj(2)=1.d0 ; patj(3)=0.d0
            elseif(trim(self%bf%bt(jbf))=='pz') then
                patj(1)=0.d0 ; patj(2)=0.d0 ; patj(3)=1.d0
            endif
            dx=atoms%ratp(1,iat)-atoms%ratp(1,jat)
            dy=atoms%ratp(2,iat)-atoms%ratp(2,jat)
            dz=atoms%ratp(3,iat)-atoms%ratp(3,jat)
            r=sqrt(dx**2+dy**2+dz**2)
            gwr0_iat=self%bf%gwe_s(iat)
            gwr2_iat=gwr0_iat
            gwp1_iat=self%bf%gwe_p(1,iat)
            qr0_iat=self%bf%be_s(iat)
            qr2_iat=(1.d0-self%bf%be_s(iat))
            gwr0_jat=self%bf%gwe_s(jat)
            gwr2_jat=gwr0_jat
            gwp1_jat=self%bf%gwe_p(1,jat)
            qr0_jat=self%bf%be_s(jat)
            qr2_jat=(1.d0-self%bf%be_s(jat))
            if(bt_i(1:1)=='s' .and. bt_j(1:1)=='s') then
                !-----------------------------------------------------------------------
                !qr0-qr0 interaction
                e_qr0_qr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwr0_iat,gwr0_jat,qr0_iat,qr0_jat)
                !-----------------------------------------------------------------------
                !qr2-qr2 interaction
                e_qr2_qr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwr2_iat,gwr2_jat,qr2_iat,qr2_jat)
                !-----------------------------------------------------------------------
                !qr0-qr2 interaction
                e_qr0_i_qr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_iat,gwr2_jat,qr0_iat,qr2_jat)
                e_qr0_j_qr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_jat,gwr2_iat,qr0_jat,qr2_iat)
                tt=e_qr0_qr0+e_qr2_qr2+e_qr0_i_qr2_j+e_qr0_j_qr2_i
            elseif(bt_i(1:1)=='p' .and. bt_j(1:1)=='p' .and. iat/=jat) then
                tt=get_ener_pr1_pr1(a,b,c,dx,dy,dz,r,sfs,gwp1_iat,gwp1_jat,pati,patj)
            elseif(iat/=jat) then
                if(bt_i(1:1)=='s') then
                    !-----------------------------------------------------------------------
                    !qr0-pr1 interaction
                    e_qr0_i_pr1_j=get_ener_qr0_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwr0_iat,gwp1_jat,qr0_iat,patj)
                    !-----------------------------------------------------------------------
                    !pr1-qr2 interaction
                    e_qr2_i_pr1_j=get_ener_qr2_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwr2_iat,gwp1_jat,qr2_iat,patj)
                    tt=e_qr0_i_pr1_j+e_qr2_i_pr1_j
                else
                    !-----------------------------------------------------------------------
                    !qr0-pr1 interaction
                    e_qr0_j_pr1_i=get_ener_qr0_pr1(a,b,c, dx, dy, dz,r,sfs,gwr0_jat,gwp1_iat,qr0_jat,pati)
                    !-----------------------------------------------------------------------
                    !pr1-qr2 interaction
                    e_qr2_j_pr1_i=get_ener_qr2_pr1(a,b,c, dx, dy, dz,r,sfs,gwr2_jat,gwp1_iat,qr2_jat,pati)
                    tt=e_qr0_j_pr1_i+e_qr2_j_pr1_i
                    !write(*,*) gwr0_jat,gwp1_iat,qr0_jat,pati
                    !write(*,*) gwr2_jat,gwp1_iat,qr2_jat,pati
                    !write(*,'(a,2i4,2es24.15)') 'PPP ',ibf,jbf,e_qr0_j_pr1_i,e_qr2_j_pr1_i
                endif
            else
                tt=0.d0
            endif
            !write(*,'(a,2i4,es24.15,2(1x,a2))') 'AMAT ',ibf,jbf,tt,trim(self%bf%bt(ibf)),trim(self%bf%bt(jbf))
            self%amat_t(ibf,jbf)=tt
            self%amat_t(jbf,ibf)=self%amat_t(ibf,jbf)
        enddo
        self%amat_t(ibf,ibf)=self%amat_t(ibf,ibf)+self%bf%hardness(ibf)
        self%amat_t(ibf,self%bf%nbf+1)=self%bf%orb(ibf)
        self%amat_t(self%bf%nbf+1,ibf)=self%bf%orb(ibf)
    enddo !end of loop over ibf
    self%amat_t(self%bf%nbf+1,self%bf%nbf+1)=0.d0
    !if(parini%mpi_env%iproc==0) then
    !do ibf=1,self%bf%nbf
    !    do jbf=1,self%bf%nbf
    !        write(27,'(f10.5)',advance='no') self%amat_t(ibf,jbf)
    !    enddo
    !write(27,*)
    !enddo
    !do jbf=1,self%bf%nbf
    !do ibf=1,self%bf%nbf
    !iat=self%bf%imap(ibf)
    !jat=self%bf%imap(jbf)
    !write(*,'(a,2i4,2es24.15,es14.5,2i4)') 'AMAT ',ibf,jbf,self%amat(ibf,jbf), &
    !    self%amat_t(ibf,jbf),abs(self%amat(ibf,jbf)-self%amat_t(ibf,jbf)),iat,jat
    !enddo
    !enddo
    !write(*,'(a,es14.5)') 'AMAT DIFF ',maxval(abs(self%amat-self%amat_t))
    !endif
    !call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)
end subroutine get_amat_cent2_analytic
!*****************************************************************************************
subroutine get_cep_rhs(self,parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: info, iat, jat, ibf, ierr
    real(8):: tt, one, ttt, p(3)
    real(8):: a, b, c, sf_1, sf_2, sf_3, sfs(3)
    real(8):: e_qr0_qcr0, e_qcr0_j_pr1_i, e_qr2_qcr2, e_qcr0_j_qr2_i
    real(8):: e_qr0_i_qcr2_j, e_qcr2_j_pr1_i
    real(8):: gwr0_iat, gwr2_iat, gwp1_iat, qr0_iat, qr2_iat, gwcr0_jat, gwcr2_jat, qcr0_jat, qcr2_jat
    real(8):: gwcr0, gwcr2, qcr0, qcr2, gwr0, gwr2, qr0, qr2
    real(8):: gg1, gg2, gg3, hh1, hh2, hh3, ss1, ss2, ss3, tt1, tt2, tt3
    real(8):: dx, dy, dz, r, pi
    character(2):: bt_t
    pi=4.d0*atan(1.d0)
    sfs(1)=parini%screening_factor
    sfs(2)=parini%screening_factor*1.1d0
    sfs(3)=parini%screening_factor*1.2d0
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))
    do ibf=1,self%bf%nbf
        iat=self%bf%imap(ibf)
        gwr0_iat=self%bf%gwe_s(iat)
        gwr2_iat=gwr0_iat
        gwp1_iat=self%bf%gwe_p(1,iat)
        qr0_iat=1.d0*self%bf%be_s(iat)
        qr2_iat=1.d0*(1.d0-self%bf%be_s(iat))
        p(1)=0.d0 ; p(2)=0.d0 ; p(3)=0.d0
        bt_t=self%bf%bt(ibf)
        if(bt_t=='px') p(1)=1.d0
        if(bt_t=='py') p(2)=1.d0
        if(bt_t=='pz') p(3)=1.d0
        ttt=0.d0
        do jat=1,atoms%nat
        if(jat==iat) cycle
        gwcr0_jat=self%bf%gwc(jat)
        gwcr2_jat=self%bf%gwc(jat)
        qcr0_jat=self%bf%qcore(jat)*self%bf%bc(jat)
        qcr2_jat=self%bf%qcore(jat)*(1.d0-self%bf%bc(jat))
        !-----------------------------------------------------------------------
        dx=atoms%ratp(1,iat)-atoms%ratp(1,jat)
        dy=atoms%ratp(2,iat)-atoms%ratp(2,jat)
        dz=atoms%ratp(3,iat)-atoms%ratp(3,jat)
        r=sqrt(dx**2+dy**2+dz**2)
        if(bt_t(1:1)=='s') then
        !qcr0-qr0 interaction
        e_qr0_qcr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwr0_iat,gwcr0_jat,qr0_iat,qcr0_jat)
        !-----------------------------------------------------------------------
        !qcr2-qr2 interaction
        e_qr2_qcr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwr2_iat,gwcr2_jat,qr2_iat,qcr2_jat)
        !-----------------------------------------------------------------------
        !qcr0-qr2 interaction
        e_qcr0_j_qr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_jat,gwr2_iat,qcr0_jat,qr2_iat)
        !-----------------------------------------------------------------------
        !qr0-qcr2 interaction
        e_qr0_i_qcr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_iat,gwcr2_jat,qr0_iat,qcr2_jat)
        ttt=ttt+e_qr0_qcr0+e_qr2_qcr2+e_qcr0_j_qr2_i+e_qr0_i_qcr2_j
        else
        !qcr0-pr1 interaction
        e_qcr0_j_pr1_i=get_ener_qr0_pr1(a,b,c, dx, dy, dz,r,sfs,gwcr0_jat,gwp1_iat,qcr0_jat,p)
        !-----------------------------------------------------------------------
        !pr1-qcr2 interaction
        e_qcr2_j_pr1_i=get_ener_qr2_pr1(a,b,c, dx, dy, dz,r,sfs,gwcr2_jat,gwp1_iat,qcr2_jat,p)
        ttt=ttt+e_qcr0_j_pr1_i+e_qcr2_j_pr1_i
        endif
        !-----------------------------------------------------------------------
        enddo
        if(bt_t(1:1)=='s') then
        gwcr0=self%bf%gwc(iat)
        gwcr2=self%bf%gwc(iat)
        qcr0=self%bf%qcore(iat)*self%bf%bc(iat)
        qcr2=self%bf%qcore(iat)*(1.d0-self%bf%bc(iat))
        gwr0=self%bf%gwe_s(iat)
        gwr2=gwr0
        qr0=1.d0*self%bf%be_s(iat)
        qr2=1.d0*(1.d0-self%bf%be_s(iat))
        gg1=gwr2*sf_1
        gg2=gwr2*sf_2
        gg3=gwr2*sf_3
        hh1=gwcr2*sf_1
        hh2=gwcr2*sf_2
        hh3=gwcr2*sf_3
        !self-interaction of qcr0-qr0 interaction
        tt1=sf_1/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_1**2)
        tt2=sf_2/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_2**2)
        tt3=sf_3/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_3**2)
        ttt=ttt+1.0d0*2.d0*qcr0*qr0*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of qcr2-qr2 interaction
        ss1=(2.d0*sf_1*(3.d0+5.d0*(hh1**2+gg1**2)+2.0d0*(hh1**4+gg1**4)+5.d0*hh1**2*gg1**2))/(3.d0*sqrt(pi)*(1.d0+hh1**2+gg1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+5.d0*(hh2**2+gg2**2)+2.0d0*(hh2**4+gg2**4)+5.d0*hh2**2*gg2**2))/(3.d0*sqrt(pi)*(1.d0+hh2**2+gg2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+5.d0*(hh3**2+gg3**2)+2.0d0*(hh3**4+gg3**4)+5.d0*hh3**2*gg3**2))/(3.d0*sqrt(pi)*(1.d0+hh3**2+gg3**2)**2.5d0)
        ttt=ttt+1.0d0*qcr2*qr2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qcr0-qr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_3**2)**1.5d0)
        ttt=ttt+1.d0*qcr0*qr2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qr0-qcr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_3**2)**1.5d0)
        ttt=ttt+1.d0*qr0*qcr2*(a*ss1+b*ss2+c*ss3)
        endif
        !if(parini%mpi_env%iproc==0) then
        !write(*,'(i4,2f20.10,es14.5)') ibf,tt,ttt,ttt-tt
        !endif
        bt_t=self%bf%bt(ibf)
        if(bt_t(1:1)=='s') then
            self%cep_rhs(ibf)=-(atoms%zat(iat)+self%bf%qcore(iat))*self%bf%hardness(ibf)-ttt
        else
            self%cep_rhs(ibf)=-ttt
        endif
    enddo
!    stop
end subroutine get_cep_rhs
!*****************************************************************************************
subroutine get_qat_from_chi_dir_cent2(self,parini,ann_arr,atoms,poisson,amat)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: amat(self%bf%nbf+1,self%bf%nbf+1)
    !local variables
    integer:: info, iat, jat, ibf, ierr
    real(8):: tt, one, ttt, p(3)
    real(8):: a, b, c, sf_1, sf_2, sf_3, sfs(3)
    real(8):: e_qr0_qcr0, e_qcr0_j_pr1_i, e_qr2_qcr2, e_qcr0_j_qr2_i
    real(8):: e_qr0_i_qcr2_j, e_qcr2_j_pr1_i
    real(8):: gwr0_iat, gwr2_iat, gwp1_iat, qr0_iat, qr2_iat, gwcr0_jat, gwcr2_jat, qcr0_jat, qcr2_jat
    real(8):: gwcr0, gwcr2, qcr0, qcr2, gwr0, gwr2, qr0, qr2
    real(8):: gg1, gg2, gg3, hh1, hh2, hh3, ss1, ss2, ss3, tt1, tt2, tt3
    real(8):: dx, dy, dz, r, pi
    character(2):: bt_t
    real(8), allocatable:: amat_t(:,:)
    real(8), allocatable:: refvec(:,:)
    real(8), allocatable:: dvec(:,:)
    pi=4.d0*atan(1.d0)
    allocate(amat_t(self%bf%nbf+1,self%bf%nbf+1))
    one=1.d0
    amat_t=amat
    call DGETRF(self%bf%nbf+1,self%bf%nbf+1,amat_t,self%bf%nbf+1,ann_arr%ipiv,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
        stop
    endif
    !call MPI_BARRIER(parini%mpi_env%mpi_comm,ierr)
    !stop 'KKKKKKKKKKKKKKKKKKKKKKKKKK'
    allocate(refvec(3,atoms%nat),source=0.d0)
    !call get_ref_vector(parini,atoms,ann_arr,refvec)
    allocate(dvec(3,atoms%nat))
    do ibf=1,self%bf%nbf
        iat=self%bf%imap(ibf)
        if(trim(self%bf%bt(ibf))=='px') dvec(1,iat)=ann_arr%chi_o(ibf)
        if(trim(self%bf%bt(ibf))=='py') dvec(2,iat)=ann_arr%chi_o(ibf)
        if(trim(self%bf%bt(ibf))=='pz') dvec(3,iat)=ann_arr%chi_o(ibf)
    enddo
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(66,'(i5,2x,a2,2x,6es14.5)') iat,trim(atoms%sat(iat)), &
            refvec(1,iat),refvec(2,iat),refvec(3,iat),dvec(1,iat),dvec(2,iat),dvec(3,iat)
    enddo
    endif
!HERE
    do ibf=1,self%bf%nbf
        if(parini%mpi_env%iproc==0) then
        write(*,*) 'RRRRRRRRRRRRRRRRRR ',ann_arr%chi_o(ibf),self%cep_rhs(ibf)
        endif
        iat=self%bf%imap(ibf)
        !write(*,'(2i5,2x,a2)') ibf,iat,trim(self%bf%bt(ibf))
        !if(.not. trim(self%bf%bt(ibf))=='s') then
        !    if(trim(atoms%sat(iat))=='Li') then
        !        ann_arr%chi_o(ibf)=0.d0
        !    else
        !        ann_arr%chi_o(ibf)=0.d0
        !    endif
        !endif
    enddo
    do ibf=1,self%bf%nbf
!HERE
!        if(parini%mpi_env%iproc==0) then
!        write(*,*) 'rrrrrrrrrrrrrrrrrr ',ann_arr%chi_o(ibf),self%cep_rhs(ibf)
!        endif
        ann_arr%qq(ibf)=self%cep_rhs(ibf)-ann_arr%chi_o(ibf)
    enddo
    ann_arr%qq(self%bf%nbf+1)=atoms%qtot-sum(atoms%zat)-sum(self%bf%qcore)
    call DGETRS('N',self%bf%nbf+1,1,amat_t,self%bf%nbf+1,ann_arr%ipiv,ann_arr%qq,self%bf%nbf+1,info)
    if(parini%mpi_env%iproc==0) then
!HERE
!    do ibf=1,self%bf%nbf
!        write(78,'(i6,f20.10)') ibf,ann_arr%qq(ibf)
!    enddo
    !do ibf=1,self%bf%nbf
    !    write(*,'(a,i6,f20.10)') 'CEPSOLUTION= ',ibf,ann_arr%qq(ibf)
    !enddo
    endif
    !stop 'YYYYYYYYYYYYYYY'
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat(1:atoms%nat)=self%bf%qcore(1:atoms%nat)
    do ibf=1,self%bf%nbf
        if(trim(self%bf%bt(ibf))=='s') then
            iat=self%bf%imap(ibf)
            atoms%qat(iat)=atoms%qat(iat)+ann_arr%qq(ibf)
        endif
    enddo
    do iat=1,atoms%nat
        write(20,'(a3,4es18.6)') atoms%sat(iat),atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat),atoms%qat(iat)+atoms%zat(iat)
    enddo
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        call yaml_map('Lagrange',ann_arr%qq(self%bf%nbf+1))
    endif
    deallocate(amat_t)
end subroutine get_qat_from_chi_dir_cent2
!*****************************************************************************************
subroutine cal_cent2_energy(self,parini,atoms,ann_arr,epot_c,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat
    real(8):: tt1, tt2
    tt1=0.d0
    tt2=0.d0
    do iat=1,atoms%nat
        tt1=tt1+ann_arr%chi_o(iat)*(atoms%qat(iat)+atoms%zat(iat))
        tt2=tt2+(atoms%qat(iat)+atoms%zat(iat))**2*0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness
    enddo
    call self%cal_electrostatic_ann_cent2(parini,atoms,ann_arr,poisson)
    if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2)  then
        write(12,*) ann_arr%epot_es,atoms%epot
    endif
    endif
    epot_c=ann_arr%epot_es+tt1+tt2+ann_arr%ener_ref
    if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2)  then
        write(16,'(5es18.8)') epot_c,ann_arr%epot_es,tt1,tt2,ann_arr%ener_ref
        write(14,*) epot_c,atoms%epot
    endif
    endif
end subroutine cal_cent2_energy
!*****************************************************************************************
subroutine cal_electrostatic_ann_cent2(self,parini,atoms,ann_arr,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_poisson), intent(inout):: poisson
    !local variables
    type(typ_poisson):: poisson_force
    integer:: iat, ix, iy, iz, ibf
    real(8):: tt, alpha, beta, ggw, ggw_t
    real(8),allocatable::fat_t(:,:)
    real(8), allocatable:: gausswidth(:)
    call self%cent2_ehartree_analytic(parini,ann_arr,atoms,ann_arr%epot_es)
    if(trim(ann_arr%event)=='potential') then
        !Force
    endif
end subroutine cal_electrostatic_ann_cent2
!*****************************************************************************************
subroutine fini_electrostatic_cent2(self,parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_poisson), intent(inout):: poisson
    call fini_hartree(parini,atoms,poisson)
    call self%bf%fini_bf()
    deallocate(self%rho_n)
    deallocate(self%gausswidth_tmp)
end subroutine fini_electrostatic_cent2
!*****************************************************************************************
subroutine get_eigenval(self,parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(in):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    !local variables
    integer:: info, nwork
    real(8), allocatable:: amat_t(:,:)
    real(8), allocatable:: eval(:), work(:)
    nwork=max(self%bf%nbf**2,100)
    allocate(amat_t(self%bf%nbf,self%bf%nbf))
    allocate(eval(self%bf%nbf))
    allocate(work(nwork))
    amat_t(1:self%bf%nbf,1:self%bf%nbf)=self%amat(1:self%bf%nbf,1:self%bf%nbf)
    call DSYEV('N','U',self%bf%nbf,amat_t,self%bf%nbf,eval,work,nwork,info)
    if(parini%mpi_env%iproc==0) then
        write(13,'(a6,es18.8,a10,es18.8)') 'MAX: ',maxval(eval),' | MIN: ',minval(eval)
    endif
    deallocate(amat_t)
    deallocate(eval)
    deallocate(work)
end subroutine get_eigenval
!*****************************************************************************************
subroutine get_dpm(atoms,dpx,dpy,dpz,dpm_err)
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
    write(*,'(a,5f8.2)') 'dipole ',atoms%qat(1),atoms%qat(2),atoms%ratp(1,2)-atoms%ratp(1,1),atoms%zat(1),atoms%zat(2)
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
end subroutine get_dpm
!*****************************************************************************************
subroutine prefit_cent2(self,parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, atom_deallocate_old
    use mod_electrostatics, only: typ_poisson
    use mod_trial_energy, only: typ_trial_energy, trial_energy_deallocate
    use mod_trial_energy, only: get_trial_energy, get_rmse
    use mod_processors, only: get_proc_stake
    use mod_flm_futile
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    type(typ_atoms):: atoms_ref
    type(typ_poisson):: poisson_ref
    integer:: itrial, iat
    integer:: nbgx, nbgy, nbgz
    real(8):: xyz(3), tt, rmse, err_U_SRS, U_SRS
    real(8):: ttr0, ttr2, qr0, gwr0, qr2, gwr2, sfs(3), a, b, c, dx, dy, dz, r
    real(8):: p(3), ttp1, ttp2, a1, a2, b1, b2
    real(8):: qavg(10), qvar(10)
    real(8):: cavg(10), cvar(10)
    real(8):: one, dpm(3), ehartree_analytic
    type(typ_trial_energy):: trial_energy
    real(8):: time1, time2
    real(8), allocatable:: amat(:,:)
    real(8), allocatable:: squarefit_raw(:,:)
    real(8), allocatable:: rhs_raw(:)
    integer:: ibf, jbf
    logical, save:: done=.false.
    if(done) return
    one=1.d0
    !call cpu_time(time1)
    call cube_read('rho.cube',atoms_ref,poisson_ref)
    call get_trial_energy(parini,atoms_ref,poisson_ref,1.5d0,self%bf%nbf,self%bf%bz,self%bf%gwz,trial_energy,atoms%qtot,atoms%dpm)
    call fini_hartree(parini,atoms_ref,poisson_ref)
    call atom_deallocate_old(atoms_ref)
    !call cpu_time(time2)
    !if(parini%mpi_env%iproc==0) then
    !    write(*,*) 'time elapsed in get_trial_energy ',time2-time1
    !endif
    nbgx=int(poisson%rgcut/poisson%hgrid(1,1))+2
    nbgy=int(poisson%rgcut/poisson%hgrid(2,2))+2
    nbgz=int(poisson%rgcut/poisson%hgrid(3,3))+2
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,f8.3,3i5)') 'RGCUT ',poisson%rgcut,nbgx,nbgy,nbgz
    write(*,'(a,3i5)') 'ngpx,ngpy,ngpz= ',poisson%ngpx,poisson%ngpy,poisson%ngpz
    write(*,'(a,3f10.6)') 'hgrid= ',poisson%hgrid(1,1),poisson%hgrid(2,2),poisson%hgrid(3,3)
    endif
    !-----------------------------------------------------------------
    call cpu_time(time1)
    !allocate(amat(self%bf%nbf+1,self%bf%nbf+1),source=0.d0)
    !if(parini%mpi_env%nproc>1) then
    !call fmpi_allreduce(amat(1,1),(self%bf%nbf+1)**2,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    !endif
    !amat(self%bf%nbf+1,self%bf%nbf+1)=0.d0
    !do jbf=1,(self%bf%nbf+1)
    !do ibf=1,(self%bf%nbf+1)
    !    self%amat(ibf,jbf)=amat(ibf,jbf)
    !    !write(92,'(i5,f15.7)') (jat-1)*(atoms%nat+1)+iat,amat(iat,jat)
    !enddo
    !enddo
    allocate(squarefit_raw(self%bf%nbf,self%bf%nbf),rhs_raw(self%bf%nbf))
    call get_cost_secder(parini,trial_energy,poisson,atoms,self%bf%nbf,self%bf%imap,self%bf%bt, &
    self%bf%be_s,self%bf%gwe_s,self%bf%gwe_p,self%bf%qcore,self%bf%bc,self%bf%gwc,squarefit_raw,rhs_raw)
    call cpu_time(time2)
    !if(parini%mpi_env%iproc==0) then
    write(*,*) 'time EP ',time2-time1
    !endif
    !-----------------------------------------------------------------
    call self%get_expansion_coeff(parini,ann_arr,atoms,squarefit_raw,rhs_raw)
!HERE
!    if(parini%mpi_env%iproc==0) then
!    do ibf=1,self%bf%nbf
!        write(77,'(i6,f20.10)') ibf,ann_arr%qq(ibf)
!    enddo
!    endif
    !-----------------------------------------------------------------
    atoms%qat(1:atoms%nat)=self%bf%qcore(1:atoms%nat)
    do ibf=1,self%bf%nbf
        if(trim(self%bf%bt(ibf))=='s') then
            iat=self%bf%imap(ibf)
            atoms%qat(iat)=atoms%qat(iat)+ann_arr%qq(ibf)
        endif
    enddo
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(*,'(a,i6,f7.3)') 'QQQ ',iat,atoms%zat(iat)+atoms%qat(iat)
    enddo
    dpm(1)=0.d0
    dpm(2)=0.d0
    dpm(3)=0.d0
    do iat=1,atoms%nat
        dpm(1)=dpm(1)+(atoms%zat(iat)+atoms%qat(iat))*atoms%ratp(1,iat)
        dpm(2)=dpm(2)+(atoms%zat(iat)+atoms%qat(iat))*atoms%ratp(2,iat)
        dpm(3)=dpm(3)+(atoms%zat(iat)+atoms%qat(iat))*atoms%ratp(3,iat)
    enddo
    do ibf=1,self%bf%nbf
        if(trim(self%bf%bt(ibf))=='px') then
            dpm(1)=dpm(1)+ann_arr%qq(ibf)
            write(*,'(a,i4,f8.3,a)') 'dpm= ',self%bf%imap(ibf),ann_arr%qq(ibf),' x'
        endif
        if(trim(self%bf%bt(ibf))=='py') then
            dpm(2)=dpm(2)+ann_arr%qq(ibf)
            write(*,'(a,i4,f8.3,a)') 'dpm= ',self%bf%imap(ibf),ann_arr%qq(ibf),' y'
        endif
        if(trim(self%bf%bt(ibf))=='pz') then
            dpm(3)=dpm(3)+ann_arr%qq(ibf)
            write(*,'(a,i4,f8.3,a)') 'dpm= ',self%bf%imap(ibf),ann_arr%qq(ibf),' z'
        endif
    enddo
    if(parini%mpi_env%iproc==0) then
        write(*,'(a,3f10.5)') 'DPM= ',dpm(1),dpm(2),dpm(3)
    endif
    !-----------------------------------------------------------------
    endif
    call self%reverseCEP(parini,ann_arr,atoms,poisson,self%amat)
    call get_rmse(trial_energy,self%bf%nbf,ann_arr%qq,rmse)
    rmse=1.d3*rmse !convert to mHa
    !call self%cal_etrial_cent2(parini,ann_arr,atoms,poisson,U_SRS)
    !call self%cent2_ehartree_analytic_test(parini,ann_arr,atoms,poisson,U_SRS)
    !call self%cent2_ehartree_analytic_monopole(parini,ann_arr,atoms,ehartree_analytic)
    call self%cent2_ehartree_analytic(parini,ann_arr,atoms,U_SRS)
    !if(parini%mpi_env%iproc==0) then
    !write(*,'(a,2es24.15,es14.5)') 'EH compare: ',U_SRS,ehartree_analytic,U_SRS-ehartree_analytic
    !endif
    !stop
    call prefit_cent2_output(parini,ann_arr,atoms,qavg,qvar,cavg,cvar)
    err_U_SRS=1.d3*(U_SRS-trial_energy%ehartree_scn_excl)/atoms%nat
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,2es24.15)') 'USRS ',U_SRS,trial_energy%ehartree_scn_excl
    write(*,'(a,2f10.3,8f7.3)') 'OPT ',rmse,err_U_SRS, &
        qavg(1),qavg(2),qvar(1),qvar(2),cavg(1),cavg(2),cvar(1),cvar(2)
    do itrial=1,trial_energy%ntrial
        write(*,'(a,i3,2es17.8)') 'ETS ',trial_energy%iat_list(itrial),trial_energy%E_all(itrial),trial_energy%energy(itrial)
    enddo
    endif
    call trial_energy_deallocate(trial_energy)
    done=.true.
end subroutine prefit_cent2
!*****************************************************************************************
subroutine get_cost_secder(parini,trial_energy,poisson,atoms,nbf,imap,bt,be_s,gwe_s,gwe_p,qcore,bc,gwc,squarefit_raw,rhs_raw)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_trial_energy, only: typ_trial_energy, trial_energy_deallocate
    use mod_trial_energy, only: get_trial_energy
    use mod_processors, only: get_proc_stake
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_trial_energy), intent(inout):: trial_energy
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: nbf, imap(nbf)
    character(2), intent(in):: bt(nbf)
    real(8), intent(in):: be_s(atoms%nat), gwe_s(atoms%nat), gwe_p(2,atoms%nat)
    real(8), intent(in):: qcore(atoms%nat), bc(atoms%nat), gwc(atoms%nat)
    real(8), intent(out):: squarefit_raw(nbf,nbf), rhs_raw(nbf)
    !local variables
    integer:: iat, ibf, ibfs, ibfe, jbf, itrial, itrials, itriale
    real(8):: sfs(3), a, b, c, dx, dy, dz, r, a1, a2, b1, b2, gwr0, gwr2
    real(8):: tt, ttr0, ttr2, qr0, qr2, p(3), ttp1, ttp2, xyz(3)
    call get_proc_stake(parini%mpi_env,nbf,ibfs,ibfe)
    trial_energy%EP=0.0
    sfs(1)=parini%screening_factor
    sfs(2)=parini%screening_factor*1.1d0
    sfs(3)=parini%screening_factor*1.2d0
    a=sfs(2)*sfs(3)*(sfs(2)+sfs(3))/((sfs(2)-sfs(1))*(sfs(3)-sfs(1))*(sfs(1)+sfs(2)+sfs(3)))
    b=sfs(1)*sfs(3)*(sfs(1)+sfs(3))/((sfs(3)-sfs(2))*(sfs(1)-sfs(2))*(sfs(1)+sfs(2)+sfs(3)))
    c=sfs(2)*sfs(1)*(sfs(2)+sfs(1))/((sfs(1)-sfs(3))*(sfs(2)-sfs(3))*(sfs(1)+sfs(2)+sfs(3)))
    do ibf=ibfs,ibfe
        iat=imap(ibf)
        do itrial=1,trial_energy%ntrial
            xyz(1:3)=atoms%ratp(1:3,trial_energy%iat_list(itrial))+trial_energy%disp(1:3,itrial)
            dx=atoms%ratp(1,iat)-xyz(1)
            dy=atoms%ratp(2,iat)-xyz(2)
            dz=atoms%ratp(3,iat)-xyz(3)
            r=sqrt(dx**2+dy**2+dz**2)
            if(trim(bt(ibf))=='s') then
                qr0=1.d0*be_s(iat)
                gwr0=gwe_s(iat)
                ttr0=get_ener_qr0_qr0(a,b,c,r,sfs,1.d0,gwr0,1.d0,qr0)
                qr2=1.d0*(1.d0-be_s(iat))
                gwr2=gwr0
                ttr2=get_ener_qr0_qr2(a,b,c,r,sfs,1.d0,gwr2,1.d0,qr2)
                tt=ttr0+ttr2
                !write(67,'(2i7,2es19.10,es14.5)') ibf,itrial,tt,ttr0+ttr2,ttr0+ttr2-tt
            elseif(trim(bt(ibf))=='px') then
                a1=gwe_p(1,iat)
                a2=gwe_p(2,iat)
                b1= a1**5/(a1**5-a2**5)
                b2=-a2**5/(a1**5-a2**5)
                p(1)=1.d0 !*b1
                p(2)=0.d0
                p(3)=0.d0
                ttp1=get_ener_qr0_pr1(a,b,c,dx,dy,dz,r,sfs,1.d0,a1,1.d0,p)
                !p(1)=1.d0*b2
                !ttp2=get_ener_qr0_pr1(a,b,c,dx,dy,dz,r,sfs,1.d0,a2,1.d0,p)
                tt=ttp1 !+ttp2
                !write(67,'(2i7,2es19.10,es14.5)') ibf,itrial,tt,ttp1+ttp2,ttp1+ttp2-tt
            elseif(trim(bt(ibf))=='py') then
                a1=gwe_p(1,iat)
                a2=gwe_p(2,iat)
                b1= a1**5/(a1**5-a2**5)
                b2=-a2**5/(a1**5-a2**5)
                p(1)=0.d0
                p(2)=1.d0 !*b1
                p(3)=0.d0
                ttp1=get_ener_qr0_pr1(a,b,c,dx,dy,dz,r,sfs,1.d0,a1,1.d0,p)
                !p(2)=1.d0*b2
                !ttp2=get_ener_qr0_pr1(a,b,c,dx,dy,dz,r,sfs,1.d0,a2,1.d0,p)
                tt=ttp1 !+ttp2
                !write(67,'(2i7,2es19.10,es14.5)') ibf,itrial,tt,ttp1+ttp2,ttp1+ttp2-tt
            elseif(trim(bt(ibf))=='pz') then
                a1=gwe_p(1,iat)
                a2=gwe_p(2,iat)
                b1= a1**5/(a1**5-a2**5)
                b2=-a2**5/(a1**5-a2**5)
                p(1)=0.d0
                p(2)=0.d0
                p(3)=1.d0 !*b1
                ttp1=get_ener_qr0_pr1(a,b,c,dx,dy,dz,r,sfs,1.d0,a1,1.d0,p)
                !p(3)=1.d0*b2
                !ttp2=get_ener_qr0_pr1(a,b,c,dx,dy,dz,r,sfs,1.d0,a2,1.d0,p)
                tt=ttp1 !+ttp2
                !write(67,'(2i7,2es19.10,es14.5)') ibf,itrial,tt,ttp1+ttp2,ttp1+ttp2-tt
            endif
            trial_energy%EP(ibf,itrial)=tt
        enddo !end of loop over itrial
    enddo !end of loop over ibf
    if(parini%mpi_env%nproc>1) then
    call fmpi_allreduce(trial_energy%EP(1,1),nbf*trial_energy%ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    endif
    call get_proc_stake(parini%mpi_env,trial_energy%ntrial,itrials,itriale)
    trial_energy%EP_n=0.d0
    do itrial=itrials,itriale
        tt=0.d0
        do iat=1,atoms%nat
        xyz(1:3)=atoms%ratp(1:3,trial_energy%iat_list(itrial))+trial_energy%disp(1:3,itrial)
        dx=atoms%ratp(1,iat)-xyz(1)
        dy=atoms%ratp(2,iat)-xyz(2)
        dz=atoms%ratp(3,iat)-xyz(3)
        r=sqrt(dx**2+dy**2+dz**2)
        qr0=qcore(iat)*bc(iat)
        gwr0=gwc(iat)
        ttr0=get_ener_qr0_qr0(a,b,c,r,sfs,1.d0,gwr0,1.d0,qr0)
        qr2=qcore(iat)*(1.d0-bc(iat))
        gwr2=gwr0
        ttr2=get_ener_qr0_qr2(a,b,c,r,sfs,1.d0,gwr2,1.d0,qr2)
        tt=tt+(ttr0+ttr2)
        enddo
        trial_energy%EP_n(itrial)=tt
        !write(68,'(i7,2es19.10,es14.5)') itrial,tt,ttr0+ttr2,ttr0+ttr2-tt
    enddo !end of loop over itrial
    if(parini%mpi_env%nproc>1) then
    call fmpi_allreduce(trial_energy%EP_n(1),trial_energy%ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    endif
    squarefit_raw=0.d0
    do ibf=1,nbf
        do jbf=1,nbf
            tt=0.d0
            do itrial=1,trial_energy%ntrial
                tt=tt+2.d0*trial_energy%EP(ibf,itrial)*trial_energy%EP(jbf,itrial)
            enddo
            squarefit_raw(ibf,jbf)=tt
        enddo
    enddo
    do ibf=1,nbf
        tt=0.d0
        do itrial=1,trial_energy%ntrial
            tt=tt+2.d0*trial_energy%EP(ibf,itrial)*(trial_energy%energy(itrial)-trial_energy%EP_n(itrial))
        enddo
        rhs_raw(ibf)=tt
    enddo
end subroutine get_cost_secder
!*****************************************************************************************
subroutine get_expansion_coeff(self,parini,ann_arr,atoms,squarefit_raw,rhs_raw)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_cent2), intent(in):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: squarefit_raw(self%bf%nbf,self%bf%nbf), rhs_raw(self%bf%nbf)
    !local variables
    integer:: iat, itrial, info, itypat, lwork, iter
    integer:: ibf, jbf
    real(8):: hh, qtarget, tt, frac, ff
    real(8):: ggw, ggw_t, alpha, beta
    real(8), allocatable:: squarefit(:,:), squarefit_t(:,:)
    real(8), allocatable:: real_eigenval(:), work(:)
    integer, allocatable:: nat_type(:)
    allocate(nat_type(parini%ntypat))
    lwork=max(self%bf%nbf**2,100)
    allocate(squarefit(self%bf%nbf,self%bf%nbf))
    allocate(squarefit_t(self%bf%nbf+1,self%bf%nbf+1))
    allocate(real_eigenval(1:self%bf%nbf),work(lwork))
    squarefit=0.d0
    squarefit_t=0.d0
    do ibf=1,self%bf%nbf
        do jbf=1,self%bf%nbf
            squarefit(ibf,jbf)=squarefit_raw(ibf,jbf)
            squarefit_t(ibf,jbf)=squarefit_raw(ibf,jbf)
        enddo
    enddo
    call DSYEV('V','U',self%bf%nbf,squarefit,self%bf%nbf,real_eigenval,work,lwork,info)
    if(parini%mpi_env%iproc==0) then
        !write(33,'(10es10.1)') squarefit(1:10,1)
        !write(33,'(10es10.1)') squarefit(1:10,2)
        !write(33,'(10es10.1)') squarefit(1:10,3)
        !write(33,'(10es10.1)') squarefit(1:10,4)
        !write(33,'(10es10.1)') squarefit(1:10,5)
        !write(33,'(10es10.1)') squarefit(1:10,6)
        !write(33,'(10es10.1)') squarefit(1:10,7)
        !write(33,'(10es10.1)') squarefit(1:10,8)
        !write(33,'(10es10.1)') squarefit(1:10,9)
        !write(33,'(10es10.1)') squarefit(1:10,10)
    do ibf=1,self%bf%nbf
        write(*,'(a,i6,es14.5)') 'EVAL ',ibf,real_eigenval(ibf)
    enddo
    endif
    !-----------------------------------------------------------------
    ff=0.d0
    do ibf=self%bf%nbf-atoms%nat+1,self%bf%nbf
        ff=ff+real_eigenval(ibf)
    enddo
    ff=ff/real(atoms%nat,kind=8)
    if(parini%mpi_env%iproc==0) then
        write(*,'(a,es14.5)') 'average nat largest eval ',ff
    endif
    nat_type=0.d0
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        nat_type(itypat)=nat_type(itypat)+1
    enddo
    if(parini%mpi_env%iproc==0) then
        write(*,*) 'nat_type= ',nat_type
    endif
    frac=0.d0
    do iter=0,0
    do itypat=1,parini%ntypat
        ann_arr%ann(itypat)%qtarget=0.d0
    enddo
    if(iter>0) then
    do ibf=1,self%bf%nbf
        if(trim(self%bf%bt(ibf))=='s') then
            iat=self%bf%imap(ibf)
            itypat=atoms%itypat(iat)
            ann_arr%ann(itypat)%qtarget=ann_arr%ann(itypat)%qtarget+ann_arr%qq(ibf)
        endif
    enddo
    do itypat=1,parini%ntypat
        ann_arr%ann(itypat)%qtarget=ann_arr%ann(itypat)%qtarget/nat_type(itypat)
    enddo
    endif
    squarefit_t=0.d0
    do ibf=1,self%bf%nbf
        do jbf=1,self%bf%nbf
            !squarefit(ibf,jbf)=squarefit_raw(ibf,jbf)
            squarefit_t(ibf,jbf)=squarefit_raw(ibf,jbf)
        enddo
    enddo
    if(iter==0) then
        hh=0.d0
    else
         hh=frac*ff
    endif
    !hh=frac*real_eigenval(self%bf%nbf-atoms%nat+1)
    do ibf=1,self%bf%nbf
        !squarefit(ibf,ibf)=squarefit(ibf,ibf)+hh
        if(trim(self%bf%bt(ibf))=='s') then
            iat=self%bf%imap(ibf)
            itypat=atoms%itypat(iat)
            squarefit_t(ibf,ibf)=squarefit_t(ibf,ibf)+hh
        endif
    enddo
    !call DSYEV('N','U',self%bf%nbf,squarefit,self%bf%nbf,real_eigenval,work,lwork,info)
    !if(parini%mpi_env%iproc==0) then
    !do ibf=1,self%bf%nbf
    !    write(*,'(a,i6,es14.5)') 'EVAL ',ibf,real_eigenval(ibf)
    !enddo
    !endif
    do ibf=1,self%bf%nbf
        squarefit_t(ibf,self%bf%nbf+1)=self%bf%orb(ibf)
        squarefit_t(self%bf%nbf+1,ibf)=self%bf%orb(ibf)
    enddo
    squarefit_t(self%bf%nbf+1,self%bf%nbf+1)=0.d0
    call DGETRF(self%bf%nbf+1,self%bf%nbf+1,squarefit_t,self%bf%nbf+1,ann_arr%ipiv,info)
    ann_arr%qq(self%bf%nbf+1)=atoms%qtot-sum(atoms%zat)-sum(self%bf%qcore)
    do ibf=1,self%bf%nbf
        ann_arr%qq(ibf)=rhs_raw(ibf)
    enddo
    do itypat=1,parini%ntypat
        qtarget=ann_arr%ann(itypat)%qtarget
        if(parini%mpi_env%iproc==0) then
        write(*,'(2a,a1,f8.3)') 'QTARGET_',trim(parini%stypat(itypat)),' ',qtarget
        endif
    enddo
    do ibf=1,self%bf%nbf
        if(trim(self%bf%bt(ibf))=='s') then
        iat=self%bf%imap(ibf)
        itypat=atoms%itypat(iat)
        qtarget=ann_arr%ann(itypat)%qtarget
        ann_arr%qq(ibf)=ann_arr%qq(ibf)+hh*qtarget
        endif
    enddo
    if(parini%mpi_env%iproc==0) then
        write(34,'(es10.1)') ann_arr%qq(self%bf%nbf+1)
    endif
    call DGETRS('N',self%bf%nbf+1,1,squarefit_t,self%bf%nbf+1,ann_arr%ipiv,ann_arr%qq,self%bf%nbf+1,info)
    frac=1.d-5
    enddo
    if(parini%mpi_env%iproc==0) then
        write(34,'(10f10.2)') ann_arr%qq(1:10)
        write(34,'(f6.2)') ann_arr%qq(1)+ann_arr%qq(2)+ann_arr%qq(6)+ann_arr%qq(7)
        write(34,'(f6.2)') atoms%ratp(1,2)-atoms%ratp(1,1)
    endif
    !-----------------------------------------------------------------
    deallocate(squarefit)
    deallocate(squarefit_t)
    deallocate(real_eigenval,work)
end subroutine get_expansion_coeff
!*****************************************************************************************
subroutine cent2_ehartree_analytic_test(self,parini,ann_arr,atoms,poisson,U_SRS)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(out):: U_SRS
    !local variables
    integer:: iat, ix, ibf
    real(8):: tt, a1, q_tmp(1)
    real(8), allocatable:: gausswidth(:)
    real(8), allocatable:: pat(:,:)
    allocate(pat(3,atoms%nat),source=0.d0)
    do ibf=1,self%bf%nbf
        iat=self%bf%imap(ibf)
        if(trim(self%bf%bt(ibf))=='px') pat(1,iat)=ann_arr%qq(ibf)
        if(trim(self%bf%bt(ibf))=='py') pat(2,iat)=ann_arr%qq(ibf)
        if(trim(self%bf%bt(ibf))=='pz') pat(3,iat)=ann_arr%qq(ibf)
    enddo
    poisson%rho=0.d0
    do iat=1,atoms%nat
        ibf=self%bf%ibf_list_s(iat)
        q_tmp(1)=ann_arr%qq(ibf)*self%bf%be_s(iat)
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,self%bf%gwe_s(iat), &
            5.d0*self%bf%gwe_s(iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)

        q_tmp(1)=ann_arr%qq(ibf)*(1.d0-self%bf%be_s(iat))
        call put_r2gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,self%bf%gwe_s(iat), &
            5.d0*self%bf%gwe_s(iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)

        a1=self%bf%gwe_p(1,iat)
        call put_gto_p_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),pat(1,iat),a1, &
            6.d0*a1,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)

        q_tmp(1)=self%bf%qcore(iat)*self%bf%bc(iat)
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,self%bf%gwc(iat), &
            6.d0*self%bf%gwc(iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        q_tmp(1)=self%bf%qcore(iat)*(1.d0-self%bf%bc(iat))
        call put_r2gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,self%bf%gwc(iat), &
            6.d0*self%bf%gwc(iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    enddo
    allocate(gausswidth(atoms%nat),source=1.d0)
    call get_hartree(parini,poisson,atoms,gausswidth,U_SRS)
end subroutine cent2_ehartree_analytic_test
!*****************************************************************************************
subroutine cent2_ehartree_analytic_monopole(self,parini,ann_arr,atoms,ehartree)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
    !local variables
    integer:: iat, jat, ibf, jbf
    real(8):: epot, pi, dx, dy, dz, r, tt1, tt2, tt3
    real(8):: ss1, ss2, ss3
    real(8):: gg1, gg2, gg3, hh1, hh2, hh3
    real(8):: a, b, c, sf_1, sf_2, sf_3, sfs(3)
    real(8):: e_qr0_qr0, e_pr1_pr1, e_qr0_j_pr1_i, e_qr0_i_pr1_j, e_qr2_qr2
    real(8):: e_qr0_i_qr2_j, e_qr0_j_qr2_i, e_qr2_j_pr1_i, e_qr2_i_pr1_j, e_qcr0_qcr0
    real(8):: e_qcr2_qcr2, e_qcr0_i_qcr2_j, e_qcr0_j_qcr2_i, e_qcr0_qr0, e_qr0_qcr0
    real(8):: e_qcr0_j_pr1_i, e_qcr0_i_pr1_j, e_qcr2_qr2, e_qr2_qcr2, e_qcr0_i_qr2_j
    real(8):: e_qr0_j_qcr2_i, e_qcr2_j_pr1_i, e_qcr2_i_pr1_j
    real(8):: e_qr0_i_qcr2_j, e_qcr0_j_qr2_i
    real(8):: gwcr0_iat, gwcr2_iat, qcr0_iat, qcr2_iat, gwcr0_jat, gwcr2_jat, qcr0_jat
    real(8):: qcr2_jat, gwr0_iat, gwr2_iat, gwp1_iat, qr0_iat, qr2_iat, gwr0_jat
    real(8):: gwr2_jat, gwp1_jat, qr0_jat, qr2_jat
    real(8):: pi_dot_pi
    real(8):: gwr0, gwr2, gwp1, qr0, qr2, gwcr0, gwcr2, qcr0, qcr2
    real(8), allocatable:: pat(:,:)
    pi=4.d0*atan(1.d0)
    allocate(pat(3,atoms%nat),source=0.d0)
    do ibf=1,self%bf%nbf
        iat=self%bf%imap(ibf)
        if(trim(self%bf%bt(ibf))=='px') pat(1,iat)=ann_arr%qq(ibf)
        if(trim(self%bf%bt(ibf))=='py') pat(2,iat)=ann_arr%qq(ibf)
        if(trim(self%bf%bt(ibf))=='pz') pat(3,iat)=ann_arr%qq(ibf)
    enddo
    sfs(1)=parini%screening_factor
    sfs(2)=parini%screening_factor*1.1d0
    sfs(3)=parini%screening_factor*1.2d0
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))
    ehartree=0.d0
    do iat=1,atoms%nat
        gwcr0=self%bf%gwc(iat)
        gwcr2=self%bf%gwc(iat)
        qcr0=self%bf%qcore(iat)*self%bf%bc(iat)
        qcr2=self%bf%qcore(iat)*(1.d0-self%bf%bc(iat))
        !self-interaction of qcr0-qcr0 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwcr0**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwcr0**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwcr0**2*sf_3**2)
        ehartree=ehartree+0.5d0*2.d0*qcr0**2*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of qcr2-qcr2 interaction
        hh1=gwcr2*sf_1
        hh2=gwcr2*sf_2
        hh3=gwcr2*sf_3
        ss1=(2.d0*sf_1*(3.d0+10.d0*hh1**2+9.d0*hh1**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*hh1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+10.d0*hh2**2+9.d0*hh2**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*hh2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+10.d0*hh3**2+9.d0*hh3**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*hh3**2)**2.5d0)
        ehartree=ehartree+0.5d0*qcr2**2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qcr0-qcr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwcr0**2+2.d0*gwcr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwcr0**2+2.d0*gwcr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwcr0**2+2.d0*gwcr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwcr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwcr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwcr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qcr0*qcr2*(a*ss1+b*ss2+c*ss3)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ibf=self%bf%ibf_list_s(iat)
        gwr0=self%bf%gwe_s(iat)
        gwr2=gwr0
        gwp1=self%bf%gwe_p(1,iat)
        qr0=ann_arr%qq(ibf)*self%bf%be_s(iat)
        qr2=ann_arr%qq(ibf)*(1.d0-self%bf%be_s(iat))
        !self-interaction of qr0-qr0 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwr0**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwr0**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwr0**2*sf_3**2)
        ehartree=ehartree+0.5d0*2.d0*qr0**2*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of qr2-qr2 interaction
        gg1=gwr2*sf_1
        gg2=gwr2*sf_2
        gg3=gwr2*sf_3
        ss1=(2.d0*sf_1*(3.d0+10.d0*gg1**2+9.d0*gg1**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+10.d0*gg2**2+9.d0*gg2**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+10.d0*gg3**2+9.d0*gg3**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg3**2)**2.5d0)
        ehartree=ehartree+0.5d0*qr2**2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qr0-qr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qr0*qr2*(a*ss1+b*ss2+c*ss3)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !self-interaction of qcr0-qr0 interaction
        tt1=sf_1/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_1**2)
        tt2=sf_2/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_2**2)
        tt3=sf_3/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_3**2)
        ehartree=ehartree+1.0d0*2.d0*qcr0*qr0*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of qcr2-qr2 interaction
        ss1=(2.d0*sf_1*(3.d0+5.d0*(hh1**2+gg1**2)+2.0d0*(hh1**4+gg1**4)+5.d0*hh1**2*gg1**2))/(3.d0*sqrt(pi)*(1.d0+hh1**2+gg1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+5.d0*(hh2**2+gg2**2)+2.0d0*(hh2**4+gg2**4)+5.d0*hh2**2*gg2**2))/(3.d0*sqrt(pi)*(1.d0+hh2**2+gg2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+5.d0*(hh3**2+gg3**2)+2.0d0*(hh3**4+gg3**4)+5.d0*hh3**2*gg3**2))/(3.d0*sqrt(pi)*(1.d0+hh3**2+gg3**2)**2.5d0)
        ehartree=ehartree+1.0d0*qcr2*qr2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qcr0-qr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qcr0*qr2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qr0-qcr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qr0*qcr2*(a*ss1+b*ss2+c*ss3)
    enddo
    do iat=1,atoms%nat
    do jat=iat+1,atoms%nat
        gwcr0_iat=self%bf%gwc(iat)
        gwcr2_iat=self%bf%gwc(iat)
        qcr0_iat=self%bf%qcore(iat)*self%bf%bc(iat)
        qcr2_iat=self%bf%qcore(iat)*(1.d0-self%bf%bc(iat))
        gwcr0_jat=self%bf%gwc(jat)
        gwcr2_jat=self%bf%gwc(jat)
        qcr0_jat=self%bf%qcore(jat)*self%bf%bc(jat)
        qcr2_jat=self%bf%qcore(jat)*(1.d0-self%bf%bc(jat))
        ibf=self%bf%ibf_list_s(iat)
        gwr0_iat=self%bf%gwe_s(iat)
        gwr2_iat=gwr0_iat
        gwp1_iat=self%bf%gwe_p(1,iat)
        qr0_iat=ann_arr%qq(ibf)*self%bf%be_s(iat)
        qr2_iat=ann_arr%qq(ibf)*(1.d0-self%bf%be_s(iat))
        jbf=self%bf%ibf_list_s(jat)
        gwr0_jat=self%bf%gwe_s(jat)
        gwr2_jat=gwr0_jat
        gwp1_jat=self%bf%gwe_p(1,jat)
        qr0_jat=ann_arr%qq(jbf)*self%bf%be_s(jat)
        qr2_jat=ann_arr%qq(jbf)*(1.d0-self%bf%be_s(jat))
        !-----------------------------------------------------------------------
        dx=atoms%ratp(1,iat)-atoms%ratp(1,jat)
        dy=atoms%ratp(2,iat)-atoms%ratp(2,jat)
        dz=atoms%ratp(3,iat)-atoms%ratp(3,jat)
        r=sqrt(dx**2+dy**2+dz**2)
        !-----------------------------------------------------------------------
        !qr0-qr0 interaction
        e_qr0_qr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwr0_iat,gwr0_jat,qr0_iat,qr0_jat)
        !-----------------------------------------------------------------------
        !qr2-qr2 interaction
        e_qr2_qr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwr2_iat,gwr2_jat,qr2_iat,qr2_jat)
        !-----------------------------------------------------------------------
        !qr0-qr2 interaction
        e_qr0_i_qr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_iat,gwr2_jat,qr0_iat,qr2_jat)
        e_qr0_j_qr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_jat,gwr2_iat,qr0_jat,qr2_iat)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !qcr0-qcr0 interaction
        e_qcr0_qcr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwcr0_iat,gwcr0_jat,qcr0_iat,qcr0_jat)
        !-----------------------------------------------------------------------
        !qcr2-qcr2 interaction
        e_qcr2_qcr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwcr2_iat,gwcr2_jat,qcr2_iat,qcr2_jat)
        !-----------------------------------------------------------------------
        !qcr0-qcr2 interaction
        e_qcr0_i_qcr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_iat,gwcr2_jat,qcr0_iat,qcr2_jat)
        e_qcr0_j_qcr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_jat,gwcr2_iat,qcr0_jat,qcr2_iat)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !qcr0-qr0 interaction
        e_qcr0_qr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwcr0_iat,gwr0_jat,qcr0_iat,qr0_jat)
        e_qr0_qcr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwr0_iat,gwcr0_jat,qr0_iat,qcr0_jat)
        !-----------------------------------------------------------------------
        !qcr2-qr2 interaction
        e_qcr2_qr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwcr2_iat,gwr2_jat,qcr2_iat,qr2_jat)
        e_qr2_qcr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwr2_iat,gwcr2_jat,qr2_iat,qcr2_jat)
        !-----------------------------------------------------------------------
        !qcr0-qr2 interaction
        e_qcr0_i_qr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_iat,gwr2_jat,qcr0_iat,qr2_jat)
        e_qcr0_j_qr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_jat,gwr2_iat,qcr0_jat,qr2_iat)
        !-----------------------------------------------------------------------
        !qr0-qcr2 interaction
        e_qr0_j_qcr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_jat,gwcr2_iat,qr0_jat,qcr2_iat)
        e_qr0_i_qcr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_iat,gwcr2_jat,qr0_iat,qcr2_jat)
        ehartree=ehartree+e_qr0_qr0+e_qr2_qr2+ &
            e_qr0_i_qr2_j+e_qr0_j_qr2_i+e_qcr0_qcr0+ &
            e_qcr2_qcr2+e_qcr0_i_qcr2_j+e_qcr0_j_qcr2_i+e_qcr0_qr0+e_qr0_qcr0+ &
            e_qcr2_qr2+e_qr2_qcr2+e_qcr0_i_qr2_j+ &
            e_qr0_i_qcr2_j+e_qcr0_j_qr2_i+ &
            e_qr0_j_qcr2_i
    enddo
    enddo
    deallocate(pat)
end subroutine cent2_ehartree_analytic_monopole
!*****************************************************************************************
subroutine cent2_ehartree_analytic(self,parini,ann_arr,atoms,ehartree)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: ehartree
    !local variables
    integer:: iat, jat, ibf, jbf
    real(8):: epot, pi, dx, dy, dz, r, tt1, tt2, tt3
    real(8):: ss1, ss2, ss3
    real(8):: gg1, gg2, gg3, hh1, hh2, hh3
    real(8):: a, b, c, sf_1, sf_2, sf_3, sfs(3)
    real(8):: e_qr0_qr0, e_pr1_pr1, e_qr0_j_pr1_i, e_qr0_i_pr1_j, e_qr2_qr2
    real(8):: e_qr0_i_qr2_j, e_qr0_j_qr2_i, e_qr2_j_pr1_i, e_qr2_i_pr1_j, e_qcr0_qcr0
    real(8):: e_qcr2_qcr2, e_qcr0_i_qcr2_j, e_qcr0_j_qcr2_i, e_qcr0_qr0, e_qr0_qcr0
    real(8):: e_qcr0_j_pr1_i, e_qcr0_i_pr1_j, e_qcr2_qr2, e_qr2_qcr2, e_qcr0_i_qr2_j
    real(8):: e_qr0_j_qcr2_i, e_qcr2_j_pr1_i, e_qcr2_i_pr1_j
    real(8):: e_qr0_i_qcr2_j, e_qcr0_j_qr2_i
    real(8):: gwcr0_iat, gwcr2_iat, qcr0_iat, qcr2_iat, gwcr0_jat, gwcr2_jat, qcr0_jat
    real(8):: qcr2_jat, gwr0_iat, gwr2_iat, gwp1_iat, qr0_iat, qr2_iat, gwr0_jat
    real(8):: gwr2_jat, gwp1_jat, qr0_jat, qr2_jat
    real(8):: pi_dot_pi
    real(8):: gwr0, gwr2, gwp1, qr0, qr2, gwcr0, gwcr2, qcr0, qcr2
    real(8), allocatable:: pat(:,:)
    pi=4.d0*atan(1.d0)
    allocate(pat(3,atoms%nat),source=0.d0)
    do ibf=1,self%bf%nbf
        iat=self%bf%imap(ibf)
        if(trim(self%bf%bt(ibf))=='px') pat(1,iat)=ann_arr%qq(ibf)
        if(trim(self%bf%bt(ibf))=='py') pat(2,iat)=ann_arr%qq(ibf)
        if(trim(self%bf%bt(ibf))=='pz') pat(3,iat)=ann_arr%qq(ibf)
    enddo
    sfs(1)=parini%screening_factor
    sfs(2)=parini%screening_factor*1.1d0
    sfs(3)=parini%screening_factor*1.2d0
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))
    ehartree=0.d0
    do iat=1,atoms%nat
        gwcr0=self%bf%gwc(iat)
        gwcr2=self%bf%gwc(iat)
        qcr0=self%bf%qcore(iat)*self%bf%bc(iat)
        qcr2=self%bf%qcore(iat)*(1.d0-self%bf%bc(iat))
        !self-interaction of qcr0-qcr0 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwcr0**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwcr0**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwcr0**2*sf_3**2)
        ehartree=ehartree+0.5d0*2.d0*qcr0**2*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of qcr2-qcr2 interaction
        hh1=gwcr2*sf_1
        hh2=gwcr2*sf_2
        hh3=gwcr2*sf_3
        ss1=(2.d0*sf_1*(3.d0+10.d0*hh1**2+9.d0*hh1**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*hh1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+10.d0*hh2**2+9.d0*hh2**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*hh2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+10.d0*hh3**2+9.d0*hh3**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*hh3**2)**2.5d0)
        ehartree=ehartree+0.5d0*qcr2**2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qcr0-qcr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwcr0**2+2.d0*gwcr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwcr0**2+2.d0*gwcr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwcr0**2+2.d0*gwcr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwcr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwcr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwcr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qcr0*qcr2*(a*ss1+b*ss2+c*ss3)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        ibf=self%bf%ibf_list_s(iat)
        gwr0=self%bf%gwe_s(iat)
        gwr2=gwr0
        gwp1=self%bf%gwe_p(1,iat)
        qr0=ann_arr%qq(ibf)*self%bf%be_s(iat)
        qr2=ann_arr%qq(ibf)*(1.d0-self%bf%be_s(iat))
        !self-interaction of qr0-qr0 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwr0**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwr0**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwr0**2*sf_3**2)
        ehartree=ehartree+0.5d0*2.d0*qr0**2*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of pr1-pr1 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwp1**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwp1**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwp1**2*sf_3**2)
        pi_dot_pi=pat(1,iat)**2+pat(2,iat)**2+pat(3,iat)**2
        ehartree=ehartree+0.5d0*4.d0*pi_dot_pi*(a*tt1**3+b*tt2**3+c*tt3**3)/(3.d0*sqrt(pi))
        !self-interaction of qr2-qr2 interaction
        gg1=gwr2*sf_1
        gg2=gwr2*sf_2
        gg3=gwr2*sf_3
        ss1=(2.d0*sf_1*(3.d0+10.d0*gg1**2+9.d0*gg1**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+10.d0*gg2**2+9.d0*gg2**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+10.d0*gg3**2+9.d0*gg3**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg3**2)**2.5d0)
        ehartree=ehartree+0.5d0*qr2**2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qr0-qr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwr0**2+2.d0*gwr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qr0*qr2*(a*ss1+b*ss2+c*ss3)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !self-interaction of qcr0-qr0 interaction
        tt1=sf_1/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_1**2)
        tt2=sf_2/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_2**2)
        tt3=sf_3/sqrt(1.d0+(gwcr0**2+gwr0**2)*sf_3**2)
        ehartree=ehartree+1.0d0*2.d0*qcr0*qr0*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of qcr2-qr2 interaction
        ss1=(2.d0*sf_1*(3.d0+5.d0*(hh1**2+gg1**2)+2.0d0*(hh1**4+gg1**4)+5.d0*hh1**2*gg1**2))/(3.d0*sqrt(pi)*(1.d0+hh1**2+gg1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+5.d0*(hh2**2+gg2**2)+2.0d0*(hh2**4+gg2**4)+5.d0*hh2**2*gg2**2))/(3.d0*sqrt(pi)*(1.d0+hh2**2+gg2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+5.d0*(hh3**2+gg3**2)+2.0d0*(hh3**4+gg3**4)+5.d0*hh3**2*gg3**2))/(3.d0*sqrt(pi)*(1.d0+hh3**2+gg3**2)**2.5d0)
        ehartree=ehartree+1.0d0*qcr2*qr2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qcr0-qr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwcr0**2+2.d0*gwr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwcr0**2+gwr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qcr0*qr2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qr0-qcr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwr0**2+2.d0*gwcr2**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwr0**2+gwcr2**2)*sf_3**2)**1.5d0)
        ehartree=ehartree+1.d0*qr0*qcr2*(a*ss1+b*ss2+c*ss3)
    enddo
    do iat=1,atoms%nat
    do jat=iat+1,atoms%nat
        gwcr0_iat=self%bf%gwc(iat)
        gwcr2_iat=self%bf%gwc(iat)
        qcr0_iat=self%bf%qcore(iat)*self%bf%bc(iat)
        qcr2_iat=self%bf%qcore(iat)*(1.d0-self%bf%bc(iat))
        gwcr0_jat=self%bf%gwc(jat)
        gwcr2_jat=self%bf%gwc(jat)
        qcr0_jat=self%bf%qcore(jat)*self%bf%bc(jat)
        qcr2_jat=self%bf%qcore(jat)*(1.d0-self%bf%bc(jat))
        ibf=self%bf%ibf_list_s(iat)
        gwr0_iat=self%bf%gwe_s(iat)
        gwr2_iat=gwr0_iat
        gwp1_iat=self%bf%gwe_p(1,iat)
        qr0_iat=ann_arr%qq(ibf)*self%bf%be_s(iat)
        qr2_iat=ann_arr%qq(ibf)*(1.d0-self%bf%be_s(iat))
        jbf=self%bf%ibf_list_s(jat)
        gwr0_jat=self%bf%gwe_s(jat)
        gwr2_jat=gwr0_jat
        gwp1_jat=self%bf%gwe_p(1,jat)
        qr0_jat=ann_arr%qq(jbf)*self%bf%be_s(jat)
        qr2_jat=ann_arr%qq(jbf)*(1.d0-self%bf%be_s(jat))
        !-----------------------------------------------------------------------
        dx=atoms%ratp(1,iat)-atoms%ratp(1,jat)
        dy=atoms%ratp(2,iat)-atoms%ratp(2,jat)
        dz=atoms%ratp(3,iat)-atoms%ratp(3,jat)
        r=sqrt(dx**2+dy**2+dz**2)
        !-----------------------------------------------------------------------
        !qr0-qr0 interaction
        e_qr0_qr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwr0_iat,gwr0_jat,qr0_iat,qr0_jat)
        !-----------------------------------------------------------------------
        !pr1-pr1 interaction
        e_pr1_pr1=get_ener_pr1_pr1(a,b,c,dx,dy,dz,r,sfs,gwp1_iat,gwp1_jat,pat(1,iat),pat(1,jat))
        !-----------------------------------------------------------------------
        !qr0-pr1 interaction
        e_qr0_j_pr1_i=get_ener_qr0_pr1(a,b,c, dx, dy, dz,r,sfs,gwr0_jat,gwp1_iat,qr0_jat,pat(1,iat))
        e_qr0_i_pr1_j=get_ener_qr0_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwr0_iat,gwp1_jat,qr0_iat,pat(1,jat))
        !-----------------------------------------------------------------------
        !qr2-qr2 interaction
        e_qr2_qr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwr2_iat,gwr2_jat,qr2_iat,qr2_jat)
        !-----------------------------------------------------------------------
        !qr0-qr2 interaction
        e_qr0_i_qr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_iat,gwr2_jat,qr0_iat,qr2_jat)
        e_qr0_j_qr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_jat,gwr2_iat,qr0_jat,qr2_iat)
        !-----------------------------------------------------------------------
        !pr1-qr2 interaction
        e_qr2_j_pr1_i=get_ener_qr2_pr1(a,b,c, dx, dy, dz,r,sfs,gwr2_jat,gwp1_iat,qr2_jat,pat(1,iat))
        e_qr2_i_pr1_j=get_ener_qr2_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwr2_iat,gwp1_jat,qr2_iat,pat(1,jat))
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !qcr0-qcr0 interaction
        e_qcr0_qcr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwcr0_iat,gwcr0_jat,qcr0_iat,qcr0_jat)
        !-----------------------------------------------------------------------
        !qcr2-qcr2 interaction
        e_qcr2_qcr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwcr2_iat,gwcr2_jat,qcr2_iat,qcr2_jat)
        !-----------------------------------------------------------------------
        !qcr0-qcr2 interaction
        e_qcr0_i_qcr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_iat,gwcr2_jat,qcr0_iat,qcr2_jat)
        e_qcr0_j_qcr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_jat,gwcr2_iat,qcr0_jat,qcr2_iat)
        !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        !qcr0-qr0 interaction
        e_qcr0_qr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwcr0_iat,gwr0_jat,qcr0_iat,qr0_jat)
        e_qr0_qcr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwr0_iat,gwcr0_jat,qr0_iat,qcr0_jat)
        !-----------------------------------------------------------------------
        !qcr0-pr1 interaction
        e_qcr0_j_pr1_i=get_ener_qr0_pr1(a,b,c, dx, dy, dz,r,sfs,gwcr0_jat,gwp1_iat,qcr0_jat,pat(1,iat))
        e_qcr0_i_pr1_j=get_ener_qr0_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwcr0_iat,gwp1_jat,qcr0_iat,pat(1,jat))
        !-----------------------------------------------------------------------
        !qcr2-qr2 interaction
        e_qcr2_qr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwcr2_iat,gwr2_jat,qcr2_iat,qr2_jat)
        e_qr2_qcr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwr2_iat,gwcr2_jat,qr2_iat,qcr2_jat)
        !-----------------------------------------------------------------------
        !qcr0-qr2 interaction
        e_qcr0_i_qr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_iat,gwr2_jat,qcr0_iat,qr2_jat)
        e_qcr0_j_qr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwcr0_jat,gwr2_iat,qcr0_jat,qr2_iat)
        !-----------------------------------------------------------------------
        !qr0-qcr2 interaction
        e_qr0_j_qcr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_jat,gwcr2_iat,qr0_jat,qcr2_iat)
        e_qr0_i_qcr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0_iat,gwcr2_jat,qr0_iat,qcr2_jat)
        !-----------------------------------------------------------------------
        !pr1-qcr2 interaction
        e_qcr2_j_pr1_i=get_ener_qr2_pr1(a,b,c, dx, dy, dz,r,sfs,gwcr2_jat,gwp1_iat,qcr2_jat,pat(1,iat))
        e_qcr2_i_pr1_j=get_ener_qr2_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwcr2_iat,gwp1_jat,qcr2_iat,pat(1,jat))
        ehartree=ehartree+e_qr0_qr0+e_pr1_pr1+e_qr0_j_pr1_i+e_qr0_i_pr1_j+e_qr2_qr2+ &
            e_qr0_i_qr2_j+e_qr0_j_qr2_i+e_qr2_j_pr1_i+e_qr2_i_pr1_j+e_qcr0_qcr0+ &
            e_qcr2_qcr2+e_qcr0_i_qcr2_j+e_qcr0_j_qcr2_i+e_qcr0_qr0+e_qr0_qcr0+ &
            e_qcr0_j_pr1_i+e_qcr0_i_pr1_j+e_qcr2_qr2+e_qr2_qcr2+e_qcr0_i_qr2_j+ &
            e_qr0_i_qcr2_j+e_qcr0_j_qr2_i+ &
            e_qr0_j_qcr2_i+e_qcr2_j_pr1_i+e_qcr2_i_pr1_j
    enddo
    enddo
    deallocate(pat)
end subroutine cent2_ehartree_analytic
!*****************************************************************************************
subroutine prefit_cent2_output(parini,ann_arr,atoms,qavg,qvar,cavg,cvar)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    real(8), intent(inout):: qavg(10), qvar(10)
    real(8), intent(inout):: cavg(10), cvar(10)
    !local variables
    integer:: iat, itypat
    real(8):: q
    real(8):: qmax(10), qmin(10)
    real(8):: cmax(10), cmin(10)
    integer, allocatable:: nat_type(:)
    allocate(nat_type(parini%ntypat),source=0)
    qavg=0.d0
    cavg=0.d0
    qmin=huge(1.d0)
    cmin=huge(1.d0)
    qmax=-huge(1.d0)
    cmax=-huge(1.d0)
    do iat=1,atoms%nat
        itypat=atoms%itypat(iat)
        nat_type(itypat)=nat_type(itypat)+1
        q=atoms%zat(iat)+atoms%qat(iat)
        qavg(itypat)=qavg(itypat)+q
        qmax(itypat)=max(qmax(itypat),q)
        qmin(itypat)=min(qmin(itypat),q)
        cavg(itypat)=cavg(itypat)+ann_arr%chi_o(iat)
        cmax(itypat)=max(cmax(itypat),ann_arr%chi_o(iat))
        cmin(itypat)=min(cmin(itypat),ann_arr%chi_o(iat))
    enddo
    do itypat=1,parini%ntypat
        qavg(itypat)=qavg(itypat)/real(nat_type(itypat),kind=8)
        qvar(itypat)=qmax(itypat)-qmin(itypat)
        cavg(itypat)=cavg(itypat)/real(nat_type(itypat),kind=8)
        cvar(itypat)=cmax(itypat)-cmin(itypat)
    enddo
    deallocate(nat_type)
end subroutine prefit_cent2_output
!*****************************************************************************************
subroutine reverseCEP(self,parini,ann_arr,atoms,poisson,amat)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_electrostatics, only: typ_poisson
    use mod_ann_io_yaml, only: write_yaml_conf_train
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: amat(self%bf%nbf+1,self%bf%nbf+1)
    !local variables
    integer:: iat, ibf, jbf
    real(8):: tt, tt1, tt2
    type(typ_file_info):: file_info
    real(8), allocatable:: ww(:)
    allocate(ww(self%bf%nbf))
    do ibf=1,self%bf%nbf
        ww(ibf)=self%cep_rhs(ibf)
    enddo
    do ibf=1,self%bf%nbf
        tt=0.d0
        do jbf=1,self%bf%nbf
            tt=tt+amat(ibf,jbf)*ann_arr%qq(jbf)
        enddo
        ann_arr%chi_o(ibf)=-tt+ww(ibf)
    enddo
    tt1=0.d0
    tt2=0.d0
    do ibf=1,self%bf%nbf
        iat=self%bf%imap(ibf)
        tt1=tt1+ann_arr%chi_o(ibf)*(ann_arr%qq(ibf)+atoms%zat(iat))
        tt2=tt2+(ann_arr%qq(ibf)+atoms%zat(iat))**2*0.5d0*self%bf%hardness(ibf)
    enddo
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,3f10.3)') 'ECH ',tt1,tt2,tt1+tt2
    do ibf=1,self%bf%nbf
        !write(*,'(a,i4,2f7.3)') 'CHI ',iat,ann_arr%chi_o(iat),atoms%zat(iat)+ann_arr%qq(iat)
        write(51,'(a,i3,a,es19.10,a)') '  - [',ibf,', ',ann_arr%chi_o(ibf),']'
    enddo
    file_info%filename_positions='posout.yaml'
    file_info%print_force=.true.
    file_info%file_position='new'
    call write_yaml_conf_train(file_info,atoms,ann_arr,.true.)
    endif
    deallocate(ww)
end subroutine reverseCEP
!*****************************************************************************************
subroutine cal_rho_pot_integral_local(xyz,xyz111,ngpx,ngpy,ngpz,hgrid,rgcut,rho,pot,ener)
    implicit none
    integer, intent(in):: ngpx, ngpy, ngpz
    real(8), intent(in):: hgrid(3,3), xyz(3), xyz111(3), rgcut
    real(8), intent(in):: rho(ngpx,ngpy,ngpz), pot(ngpx,ngpy,ngpz)
    real(8), intent(out):: ener
    !local variables
    integer:: igpx, igpy, igpz
    integer:: nbgx, nbgy, nbgz, agpx, agpy, agpz
    real(8):: res
    nbgx=int(rgcut/hgrid(1,1))+2
    nbgy=int(rgcut/hgrid(2,2))+2
    nbgz=int(rgcut/hgrid(3,3))+2
    agpx=int((xyz(1)-xyz111(1))/hgrid(1,1))
    agpy=int((xyz(2)-xyz111(2))/hgrid(2,2))
    agpz=int((xyz(3)-xyz111(3))/hgrid(3,3))
    res=0.d0
    do igpz=agpz-nbgz,agpz+nbgz
        do igpy=agpy-nbgy,agpy+nbgy
            do igpx=agpx-nbgx,agpx+nbgx
                res=res+rho(igpx,igpy,igpz)*pot(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    ener=res*(hgrid(1,1)*hgrid(2,2)*hgrid(3,3))
end subroutine cal_rho_pot_integral_local
!*****************************************************************************************
subroutine cent2_analytic(parini,atoms,gwr0,gwp1,gwr2,qr0,pat,qr2)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: gwr0(atoms%nat), gwp1(atoms%nat), gwr2(atoms%nat)
    real(8), intent(in):: qr0(atoms%nat), pat(3,atoms%nat), qr2(atoms%nat)
    !local variables
    integer:: iat, jat
    real(8):: epot, pi, dx, dy, dz, r, tt1, tt2, tt3
    real(8):: ss1, ss2, ss3
    real(8):: gg1, gg2, gg3, ggg
    real(8):: a, b, c, sf_1, sf_2, sf_3, sfs(3)
    real(8):: e_qr0_qr0, e_pr1_pr1, e_qr0_j_pr1_i, e_qr0_i_pr1_j, e_qr2_qr2
    real(8):: e_qr0_i_qr2_j, e_qr0_j_qr2_i, e_qr2_j_pr1_i, e_qr2_i_pr1_j
    real(8):: pi_dot_pi
    pi=4.d0*atan(1.d0)
    sfs(1)=parini%screening_factor
    sfs(2)=parini%screening_factor*1.1d0
    sfs(3)=parini%screening_factor*1.2d0
    sf_1=parini%screening_factor
    sf_2=parini%screening_factor*1.1d0
    sf_3=parini%screening_factor*1.2d0
    a=sf_2*sf_3*(sf_2+sf_3)/((sf_2-sf_1)*(sf_3-sf_1)*(sf_1+sf_2+sf_3))
    b=sf_1*sf_3*(sf_1+sf_3)/((sf_3-sf_2)*(sf_1-sf_2)*(sf_1+sf_2+sf_3))
    c=sf_2*sf_1*(sf_2+sf_1)/((sf_1-sf_3)*(sf_2-sf_3)*(sf_1+sf_2+sf_3))
    epot=0.d0
    do iat=1,atoms%nat
        !self-interaction of qr0-qr0 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwr0(iat)**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwr0(iat)**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwr0(iat)**2*sf_3**2)
        epot=epot+0.5d0*2.d0*qr0(iat)**2*(a*tt1+b*tt2+c*tt3)/sqrt(pi)
        !self-interaction of pr1-pr1 interaction
        tt1=sf_1/sqrt(1.d0+2.d0*gwp1(iat)**2*sf_1**2)
        tt2=sf_2/sqrt(1.d0+2.d0*gwp1(iat)**2*sf_2**2)
        tt3=sf_3/sqrt(1.d0+2.d0*gwp1(iat)**2*sf_3**2)
        pi_dot_pi=pat(1,iat)**2+pat(2,iat)**2+pat(3,iat)**2
        epot=epot+0.5d0*4.d0*pi_dot_pi*(a*tt1**3+b*tt2**3+c*tt3**3)/(3.d0*sqrt(pi))
        !self-interaction of qr2-qr2 interaction
        gg1=gwr2(iat)*sf_1
        gg2=gwr2(iat)*sf_2
        gg3=gwr2(iat)*sf_3
        ss1=(2.d0*sf_1*(3.d0+10.d0*gg1**2+9.d0*gg1**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg1**2)**2.5d0)
        ss2=(2.d0*sf_2*(3.d0+10.d0*gg2**2+9.d0*gg2**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg2**2)**2.5d0)
        ss3=(2.d0*sf_3*(3.d0+10.d0*gg3**2+9.d0*gg3**4))/(3.d0*sqrt(pi)*(1.d0+2.d0*gg3**2)**2.5d0)
        epot=epot+0.5d0*qr2(iat)**2*(a*ss1+b*ss2+c*ss3)
        !self-interaction of qr0-qr2 interaction
        tt1=2.d0*sf_1*(3.d0+(3.d0*gwr0(iat)**2+2.d0*gwr2(iat)**2)*sf_1**2)
        tt2=2.d0*sf_2*(3.d0+(3.d0*gwr0(iat)**2+2.d0*gwr2(iat)**2)*sf_2**2)
        tt3=2.d0*sf_3*(3.d0+(3.d0*gwr0(iat)**2+2.d0*gwr2(iat)**2)*sf_3**2)
        ss1=tt1/(3.d0*sqrt(pi)*(1.d0+(gwr0(iat)**2+gwr2(iat)**2)*sf_1**2)**1.5d0)
        ss2=tt2/(3.d0*sqrt(pi)*(1.d0+(gwr0(iat)**2+gwr2(iat)**2)*sf_2**2)**1.5d0)
        ss3=tt3/(3.d0*sqrt(pi)*(1.d0+(gwr0(iat)**2+gwr2(iat)**2)*sf_3**2)**1.5d0)
        epot=epot+1.d0*qr0(iat)*qr2(iat)*(a*ss1+b*ss2+c*ss3)
    enddo
    do iat=1,atoms%nat
    do jat=iat+1,atoms%nat
        dx=atoms%ratp(1,iat)-atoms%ratp(1,jat)
        dy=atoms%ratp(2,iat)-atoms%ratp(2,jat)
        dz=atoms%ratp(3,iat)-atoms%ratp(3,jat)
        r=sqrt(dx**2+dy**2+dz**2)
        !-----------------------------------------------------------------------
        !qr0-qr0 interaction
        e_qr0_qr0=get_ener_qr0_qr0(a,b,c,r,sfs,gwr0(iat),gwr0(jat),qr0(iat),qr0(jat))
        !-----------------------------------------------------------------------
        !pr1-pr1 interaction
        e_pr1_pr1=get_ener_pr1_pr1(a,b,c,dx,dy,dz,r,sfs,gwp1(iat),gwp1(jat),pat(1,iat),pat(1,jat))
        !-----------------------------------------------------------------------
        !qr0-pr1 interaction
        e_qr0_j_pr1_i=get_ener_qr0_pr1(a,b,c, dx, dy, dz,r,sfs,gwr0(jat),gwp1(iat),qr0(jat),pat(1,iat))
        e_qr0_i_pr1_j=get_ener_qr0_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwr0(iat),gwp1(jat),qr0(iat),pat(1,jat))
        !-----------------------------------------------------------------------
        !qr2-qr2 interaction
        e_qr2_qr2=get_ener_qr2_qr2(a,b,c,r,sfs,gwr2(iat),gwr2(jat),qr2(iat),qr2(jat))
        !-----------------------------------------------------------------------
        !qr0-qr2 interaction
        e_qr0_i_qr2_j=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0(iat),gwr2(jat),qr0(iat),qr2(jat))
        e_qr0_j_qr2_i=get_ener_qr0_qr2(a,b,c,r,sfs,gwr0(jat),gwr2(iat),qr0(jat),qr2(iat))
        !-----------------------------------------------------------------------
        !pr1-qr2 interaction
        e_qr2_j_pr1_i=get_ener_qr2_pr1(a,b,c, dx, dy, dz,r,sfs,gwr2(jat),gwp1(iat),qr2(jat),pat(1,iat))
        e_qr2_i_pr1_j=get_ener_qr2_pr1(a,b,c,-dx,-dy,-dz,r,sfs,gwr2(iat),gwp1(jat),qr2(iat),pat(1,jat))
        epot=epot+e_qr0_qr0+e_pr1_pr1+e_qr0_j_pr1_i+e_qr0_i_pr1_j+e_qr2_qr2+ &
            e_qr0_i_qr2_j+e_qr0_j_qr2_i+e_qr2_j_pr1_i+e_qr2_i_pr1_j
    enddo
    enddo
    atoms%epot=epot
end subroutine cent2_analytic
!*****************************************************************************************
pure function get_ener_qr0_qr0(a,b,c,r,sfs,gw1,gw2,q1,q2) result(ener)
    implicit none
    real(8), intent(in):: a, b, c, r, sfs(3), gw1, gw2, q1, q2
    !local variables
    real(8):: ener, tt1, tt2, tt3
    tt1=sfs(1)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(1)**2)
    tt2=sfs(2)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(2)**2)
    tt3=sfs(3)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(3)**2)
    ener=q1*q2*(a*erf(tt1*r)+b*erf(tt2*r)+c*erf(tt3*r))/r
end function get_ener_qr0_qr0
!*****************************************************************************************
pure function get_ener_pr1_pr1(a,b,c,dx,dy,dz,r,sfs,gw1,gw2,p1,p2) result(ener)
    implicit none
    real(8), intent(in):: a, b, c, dx, dy, dz, r, sfs(3), gw1, gw2, p1(3), p2(3)
    !local variables
    real(8):: ener, tt1, tt2, tt3, ww1, ww2, ww3, vv1, vv2, vv3, vvv, ss1, ss2, ss3
    real(8):: uu1, uu2, uu3, pi
    real(8):: pj_dot_rij, pi_dot_rij, pi_dot_pj
    pi=4.d0*atan(1.d0)
    tt1=sfs(1)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(1)**2)
    tt2=sfs(2)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(2)**2)
    tt3=sfs(3)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(3)**2)
    pj_dot_rij=dx*p2(1)+dy*p2(2)+dz*p2(3)
    pi_dot_rij=dx*p1(1)+dy*p1(2)+dz*p1(3)
    pi_dot_pj=p1(1)*p2(1)+p1(2)*p2(2)+p1(3)*p2(3)
    ss1=exp(-(tt1*r)**2)
    ss2=exp(-(tt2*r)**2)
    ss3=exp(-(tt3*r)**2)
    uu1=erf(tt1*r)
    uu2=erf(tt2*r)
    uu3=erf(tt3*r)
    ww1=-pi_dot_pj*(2.d0*tt1*ss1/(sqrt(pi)*r**2)-uu1/r**3)
    ww2=-pi_dot_pj*(2.d0*tt2*ss2/(sqrt(pi)*r**2)-uu2/r**3)
    ww3=-pi_dot_pj*(2.d0*tt3*ss3/(sqrt(pi)*r**2)-uu3/r**3)
    vvv=-pj_dot_rij*pi_dot_rij
    vv1=vvv*(-4.d0*tt1**3*ss1/(sqrt(pi)*r**2)-6.d0*tt1*ss1/(sqrt(pi)*r**4)+3.d0*uu1/r**5)
    vv2=vvv*(-4.d0*tt2**3*ss2/(sqrt(pi)*r**2)-6.d0*tt2*ss2/(sqrt(pi)*r**4)+3.d0*uu2/r**5)
    vv3=vvv*(-4.d0*tt3**3*ss3/(sqrt(pi)*r**2)-6.d0*tt3*ss3/(sqrt(pi)*r**4)+3.d0*uu3/r**5)
    ener=a*(ww1+vv1)+b*(ww2+vv2)+c*(ww3+vv3)
end function get_ener_pr1_pr1
!*****************************************************************************************
pure function get_ener_qr0_pr1(a,b,c,dx,dy,dz,r,sfs,gwq,gwp,q,p) result(ener)
    implicit none
    real(8), intent(in):: a, b, c, dx, dy, dz, r, sfs(3), gwq, gwp, q, p(3)
    !local variables
    real(8):: ener, tt1, tt2, tt3, vv1, vv2, vv3, vvv, ss1, ss2, ss3
    real(8):: uu1, uu2, uu3, pi
    real(8):: p_dot_rij
    pi=4.d0*atan(1.d0)
    tt1=sfs(1)/sqrt(1.d0+(gwq**2+gwp**2)*sfs(1)**2)
    tt2=sfs(2)/sqrt(1.d0+(gwq**2+gwp**2)*sfs(2)**2)
    tt3=sfs(3)/sqrt(1.d0+(gwq**2+gwp**2)*sfs(3)**2)
    ss1=exp(-(tt1*r)**2)
    ss2=exp(-(tt2*r)**2)
    ss3=exp(-(tt3*r)**2)
    uu1=erf(tt1*r)
    uu2=erf(tt2*r)
    uu3=erf(tt3*r)
    p_dot_rij=dx*p(1)+dy*p(2)+dz*p(3)
    vvv=p_dot_rij*q
    vv1=vvv*(2.d0*tt1*ss1/(sqrt(pi)*r**2)-uu1/r**3)
    vv2=vvv*(2.d0*tt2*ss2/(sqrt(pi)*r**2)-uu2/r**3)
    vv3=vvv*(2.d0*tt3*ss3/(sqrt(pi)*r**2)-uu3/r**3)
    ener=a*vv1+b*vv2+c*vv3
end function get_ener_qr0_pr1
!*****************************************************************************************
pure function get_ener_qr2_qr2(a,b,c,r,sfs,gw1,gw2,q1,q2) result(ener)
    implicit none
    real(8), intent(in):: a, b, c, r, sfs(3), gw1, gw2, q1, q2
    !local variables
    real(8):: ener, pi, tt1, tt2, tt3
    real(8):: bb1, bb2, bb3, ff1, ff2, ff3, ss1, ss2, ss3
    real(8):: uu1, uu2, uu3, dd1, dd2, dd3, gg1, gg2, gg3
    pi=4.d0*atan(1.d0)
    tt1=sfs(1)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(1)**2)
    tt2=sfs(2)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(2)**2)
    tt3=sfs(3)/sqrt(1.d0+(gw1**2+gw2**2)*sfs(3)**2)
    ss1=exp(-(tt1*r)**2)
    ss2=exp(-(tt2*r)**2)
    ss3=exp(-(tt3*r)**2)
    uu1=erf(tt1*r)
    uu2=erf(tt2*r)
    uu3=erf(tt3*r)
    gg1=1.d0+(gw1**2+gw2**2)*sfs(1)**2
    gg2=1.d0+(gw1**2+gw2**2)*sfs(2)**2
    gg3=1.d0+(gw1**2+gw2**2)*sfs(3)**2
    bb1=gw1**2*gw2**2*r*sfs(1)**2*(-4.d0*r**2*sfs(1)**2+6.d0*gg1)
    bb2=gw1**2*gw2**2*r*sfs(2)**2*(-4.d0*r**2*sfs(2)**2+6.d0*gg2)
    bb3=gw1**2*gw2**2*r*sfs(3)**2*(-4.d0*r**2*sfs(3)**2+6.d0*gg3)
    ff1=sfs(1)**3*(bb1-6.d0*r*gg1**2*(gw1**2+gw2**2))
    ff2=sfs(2)**3*(bb2-6.d0*r*gg2**2*(gw1**2+gw2**2))
    ff3=sfs(3)**3*(bb3-6.d0*r*gg3**2*(gw1**2+gw2**2))
    dd1=(ss1*ff1+9.d0*sqrt(pi)*gg1**3.5d0*uu1)/(9.d0*sqrt(pi)*r*gg1**3.5d0)
    dd2=(ss2*ff2+9.d0*sqrt(pi)*gg2**3.5d0*uu2)/(9.d0*sqrt(pi)*r*gg2**3.5d0)
    dd3=(ss3*ff3+9.d0*sqrt(pi)*gg3**3.5d0*uu3)/(9.d0*sqrt(pi)*r*gg3**3.5d0)
    ener=q1*q2*(a*dd1+b*dd2+c*dd3)
end function get_ener_qr2_qr2
!*****************************************************************************************
pure function get_ener_qr0_qr2(a,b,c,r,sfs,gwr0,gwr2,qr0,qr2) result(ener)
    implicit none
    real(8), intent(in):: a, b, c, r, sfs(3), gwr0, gwr2, qr0, qr2
    !local variables
    real(8):: ener, pi, tt1, tt2, tt3
    real(8):: ww1, ww2, ww3, ss1, ss2, ss3
    real(8):: uu1, uu2, uu3, gg1, gg2, gg3
    pi=4.d0*atan(1.d0)
    tt1=sfs(1)/sqrt(1.d0+(gwr0**2+gwr2**2)*sfs(1)**2)
    tt2=sfs(2)/sqrt(1.d0+(gwr0**2+gwr2**2)*sfs(2)**2)
    tt3=sfs(3)/sqrt(1.d0+(gwr0**2+gwr2**2)*sfs(3)**2)
    ss1=exp(-(tt1*r)**2)
    ss2=exp(-(tt2*r)**2)
    ss3=exp(-(tt3*r)**2)
    uu1=erf(tt1*r)
    uu2=erf(tt2*r)
    uu3=erf(tt3*r)
    gg1=1.d0+(gwr2**2+gwr0**2)*sfs(1)**2
    gg2=1.d0+(gwr2**2+gwr0**2)*sfs(2)**2
    gg3=1.d0+(gwr2**2+gwr0**2)*sfs(3)**2
    ww1=-2.d0*sfs(1)**3*(gwr2**2)*ss1/(3.d0*sqrt(pi)*gg1**1.5d0)+1.d0*uu1/r
    ww2=-2.d0*sfs(2)**3*(gwr2**2)*ss2/(3.d0*sqrt(pi)*gg2**1.5d0)+1.d0*uu2/r
    ww3=-2.d0*sfs(3)**3*(gwr2**2)*ss3/(3.d0*sqrt(pi)*gg3**1.5d0)+1.d0*uu3/r
    ener=qr0*qr2*(a*ww1+b*ww2+c*ww3)
end function get_ener_qr0_qr2
!*****************************************************************************************
pure function get_ener_qr2_pr1(a,b,c,dx,dy,dz,r,sfs,gwq,gwp,q,p) result(ener)
    implicit none
    real(8), intent(in):: a, b, c, dx, dy, dz, r, sfs(3), gwq, gwp, q, p(3)
    !local variables
    real(8):: ener, tt1, tt2, tt3, ww1, ww2, ww3, ss1, ss2, ss3
    real(8):: aa1, aa2, aa3, bb1, bb2, bb3, uu1, uu2, uu3, gg1, gg2, gg3, pi
    real(8):: p_dot_rij
    pi=4.d0*atan(1.d0)
    tt1=sfs(1)/sqrt(1.d0+(gwq**2+gwp**2)*sfs(1)**2)
    tt2=sfs(2)/sqrt(1.d0+(gwq**2+gwp**2)*sfs(2)**2)
    tt3=sfs(3)/sqrt(1.d0+(gwq**2+gwp**2)*sfs(3)**2)
    ss1=exp(-(tt1*r)**2)
    ss2=exp(-(tt2*r)**2)
    ss3=exp(-(tt3*r)**2)
    uu1=erf(tt1*r)
    uu2=erf(tt2*r)
    uu3=erf(tt3*r)
    gg1=1.d0+(gwq**2+gwp**2)*sfs(1)**2
    gg2=1.d0+(gwq**2+gwp**2)*sfs(2)**2
    gg3=1.d0+(gwq**2+gwp**2)*sfs(3)**2
    aa1=3.d0*(1.d0+sfs(1)**4*(gwp**4+gwq**4))+2.d0*gwq**2*sfs(1)**2*(3.d0+r**2*sfs(1)**2)
    aa2=3.d0*(1.d0+sfs(2)**4*(gwp**4+gwq**4))+2.d0*gwq**2*sfs(2)**2*(3.d0+r**2*sfs(2)**2)
    aa3=3.d0*(1.d0+sfs(3)**4*(gwp**4+gwq**4))+2.d0*gwq**2*sfs(3)**2*(3.d0+r**2*sfs(3)**2)
    bb1=6.d0*gwp**2*sfs(1)**2*(1.d0+gwq**2*sfs(1)**2)
    bb2=6.d0*gwp**2*sfs(2)**2*(1.d0+gwq**2*sfs(2)**2)
    bb3=6.d0*gwp**2*sfs(3)**2*(1.d0+gwq**2*sfs(3)**2)
    ww1=(-2.d0*ss1*sfs(1)*(aa1+bb1))/(3.d0*sqrt(pi)*r**2*gg1**2.5d0)+uu1/r**3
    ww2=(-2.d0*ss2*sfs(2)*(aa2+bb2))/(3.d0*sqrt(pi)*r**2*gg2**2.5d0)+uu2/r**3
    ww3=(-2.d0*ss3*sfs(3)*(aa3+bb3))/(3.d0*sqrt(pi)*r**2*gg3**2.5d0)+uu3/r**3
    p_dot_rij=dx*p(1)+dy*p(2)+dz*p(3)
    ener=-p_dot_rij*q*(a*ww1+b*ww2+c*ww3)
end function get_ener_qr2_pr1
!*****************************************************************************************
function cutoff_function(r, rc) result(fc)
    implicit none
    real(8), intent(in):: r, rc
    !local variables
    real(8):: fc, pi
    if(r<rc) then
        fc=(1.d0-(r/rc)**2)**3
    else
        fc=0.d0
    endif
    !pi=4.d0*atan(1.d0)
    !if(r<rc*.5d0) then
    !    fc=1.d0
    !elseif(r<rc) then
    !    fc=cos((1.d0-((1.d0-((r-rc*0.5d0)/(rc*0.5d0))**2)**3))*pi*0.5d0)
    !else
    !    fc=0.d0
    !endif
    !Compare it with: PHYSICAL REVIEW B 93, 155203 (2016) -> comparison is done:
    !the second derivative in the cutoff function in PRB paper does not vanish.
end function cutoff_function
!*****************************************************************************************
function cutoff_function_der(r, rc) result(fcd)
    implicit none
    real(8), intent(in):: r, rc
    !local variables
    real(8):: fcd, pi
    if(r<rc) then
        fcd=-6.d0*(r/rc**2)*(1.d0-(r/rc)**2)**2
    else
        fcd=0.d0
    endif
end function cutoff_function_der
!*****************************************************************************************
end module mod_cent2
!*****************************************************************************************
