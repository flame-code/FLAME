!*****************************************************************************************
module mod_cent2
    implicit none
    private
    public:: typ_cent2
    type typ_cent2
        private
        integer:: nbgx, nbgy, nbgz
        real(8), allocatable:: rho_tmp(:,:,:)
        real(8), allocatable:: rho_e_all(:,:,:,:)
        real(8), allocatable:: rho_n_all(:,:,:,:)
        real(8), allocatable:: gausswidth_n(:)
        real(8), allocatable:: gausswidth_tmp(:)
        contains
        !procedure, public, pass(self)::
        !procedure, public, pass(self):: init_cent2
        !procedure, public, pass(self):: fini_cent2
        procedure, private, pass(self):: calc_atomic_densities
        procedure, private, pass(self):: init_electrostatic_cent2
        procedure, private, pass(self):: fini_electrostatic_cent2
        procedure, private, pass(self):: get_amat_cent2
        procedure, public, pass(self):: cal_ann_cent2
    end type typ_cent2
contains
!*****************************************************************************************
subroutine cal_ann_cent2(self,parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    use yaml_output
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_poisson):: poisson
    !local variables
    type(typ_pia_arr):: pia_arr_tmp
    integer:: iat, i, j, ng
    integer:: iats, iate, mat, mproc
    real(8):: epot_c, out_ann
    real(8):: dpx, dpy, dpz
    real(8):: time1, time2, time3, time4, time5, time6, time7, timet1, timet2
    real(8):: tt1, hinv(3,3), vol
    call f_routine(id='cal_ann_cent2')
    call update_ratp(atoms)
    if(.not. allocated(ann_arr%ipiv))allocate(ann_arr%ipiv(1:atoms%nat+1))
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fat_chi(1:3,1:atoms%nat))
        allocate(ann_arr%chi_i(1:atoms%nat))
        allocate(ann_arr%chi_o(1:atoms%nat))
        allocate(ann_arr%chi_d(1:atoms%nat))
        allocate(ann_arr%a(1:(atoms%nat+1)*(atoms%nat+1)))
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
    if(parini%iverbose>=2) call cpu_time(time1)
    call self%init_electrostatic_cent2(parini,atoms,ann_arr,poisson)
    if(parini%iverbose>=2) call cpu_time(time2)
    !if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2) write(*,*) 'init_time: ',time2-time1
    !endif
    call symfunc%init_symfunc(parini%mpi_env)
    if(ann_arr%compute_symfunc) then
        call symfunc%get_symfunc(parini,ann_arr,atoms,.true.)
    else
        symfunc%linked_lists%rcut=ann_arr%rcut
        symfunc%linked_lists%triplex=.true.
        call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr_tmp)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fatpq(1:3,1:symfunc%linked_lists%maxbound_rad))
        allocate(ann_arr%stresspq(1:3,1:3,1:symfunc%linked_lists%maxbound_rad))
    endif
    if(parini%iverbose>=2) call cpu_time(time3)
    if(parini%mpi_env%nproc>1) then
        mat=atoms%nat/parini%mpi_env%nproc
        iats=parini%mpi_env%iproc*mat+1
        mproc=mod(atoms%nat,parini%mpi_env%nproc)
        iats=iats+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
        if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) mat=mat+1
        iate=iats+mat-1
    else
        iats=1
        iate=atoms%nat
        !mat=atoms%nat
    endif
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
    !This must be here since contribution from coulomb
    !interaction is calculated during the process of charge optimization.
    if(parini%iverbose>=2) call cpu_time(time4)
    call get_qat_from_chi_dir_cent2(parini,ann_arr,atoms,poisson,ann_arr%a)
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
    if(parini%iverbose>=2) call cpu_time(timet1)
    !if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2) write(*,*) 'qet_qat_time: ',timet1-time4
    !endif
    if(.not. trim(ann_arr%event)=='potential') then
        call cal_cent2_total_potential(parini,ann_arr,atoms,poisson)
    endif
    if(parini%prefit_ann) then
        call prefit_cent2(parini,ann_arr,atoms,poisson)
    endif
    if(parini%iverbose>=2) call cpu_time(time5)
    !if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2) write(*,*) 'cent2_g_per_time: ',time5-timet1
    !endif
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call cal_cent2_energy(parini,atoms,ann_arr,epot_c,poisson)
    if(parini%iverbose>=2) call cpu_time(timet2)
    !if(parini%mpi_env%iproc==0) then
    if(parini%iverbose>=2) write(*,*) 'get_electrostatic_time: ',timet2-time6
    !endif
    if(parini%iverbose>=2) then
        call cpu_time(time7)
        call yaml_mapping_open('Timing of CENT2')
        call yaml_map('initialize matrix',time2-time1)
        call yaml_map('calculation of symfunc',time3-time2)
        call yaml_map('neural network process',time4-time3)
        call yaml_map('linear equations solver',time5-time4)
        call yaml_map('force (SR term)',time6-time5)
        call yaml_map('energy (SR+LR), force (LR)',time7-time6)
        call yaml_map('total time',time7-time1)
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
        deallocate(ann_arr%a)
        deallocate(ann_arr%fat_chi)
        deallocate(ann_arr%fatpq)
        deallocate(ann_arr%stresspq)
    endif
    deallocate(ann_arr%ipiv)
    call symfunc%fini_symfunc()
    call f_release_routine()
end subroutine cal_ann_cent2
!*****************************************************************************************
subroutine calc_atomic_densities(self,parini,atoms,ann_arr,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    use yaml_output
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_poisson), intent(in):: poisson
    !local variables
    integer:: iat, ix, iy, iz, agpx, agpy, agpz
    real(8):: ggw, ggw_t, alpha, beta, q_tmp(1)
    do iat=1,atoms%nat
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+self%nbgx+0
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+self%nbgy+0
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+self%nbgz+0
        !do iz=agpz-nbgz,agpz+nbgz
        !do iy=agpy-nbgy,agpy+nbgy
        !do ix=agpx-nbgx,agpx+nbgx
        !    dx=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
        !    dy=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
        !    dz=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
        !    dr=sqrt(dx**2+dy**2+dz**2)
        !    linearGridNumber=floor(dr/ann_arr%radpots_cent2%hgp)
        !    rho_e_p1=ann_arr%radpots_cent2%linear_rho_e(linearGridNumber+1,atoms%itypat(iat))
        !    rho_e=ann_arr%radpots_cent2%linear_rho_e(linearGridNumber,atoms%itypat(iat))
        !    tt=(dr/ann_arr%radpots_cent2%hgp-linearGridNumber)*(rho_e_p1-rho_e)+rho_e
        !    rho_all(ix-agpx,iy-agpy,iz-agpz,iat)=tt
        !enddo
        !enddo
        !enddo
        ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        alpha=ggw**3/(ggw**3-ggw_t**3)
        beta=-ggw_t**3/(ggw**3-ggw_t**3)
        q_tmp(1)=alpha
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp,ggw, &
            6.d0*ggw,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        q_tmp(1)=beta
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,ggw_t, &
            6.d0*ggw_t,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        do iz=agpz-self%nbgz,agpz+self%nbgz
        do iy=agpy-self%nbgy,agpy+self%nbgy
        do ix=agpx-self%nbgx,agpx+self%nbgx
            self%rho_e_all(ix-agpx,iy-agpy,iz-agpz,iat)=poisson%rho(ix,iy,iz)
        enddo
        enddo
        enddo
        q_tmp(1)=1.d0
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp,self%gausswidth_n(iat), &
            6.d0*maxval(self%gausswidth_n),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        do iz=agpz-self%nbgz,agpz+self%nbgz
        do iy=agpy-self%nbgy,agpy+self%nbgy
        do ix=agpx-self%nbgx,agpx+self%nbgx
            self%rho_n_all(ix-agpx,iy-agpy,iz-agpz,iat)=poisson%rho(ix,iy,iz)
        enddo
        enddo
        enddo
    enddo !end of loop over iat
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
    integer:: iat, itype
    integer:: agpx, agpy, agpz, ix, iy, iz
    real(8):: time1, time2
    real(8):: max_cellVec, tt
    !real(8):: time1_t, time2_t, time3_t
    real(8), allocatable:: rho(:,:,:)
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
        allocate(self%gausswidth_n(1:atoms%nat))
        do iat=1,atoms%nat
            self%gausswidth_tmp(iat)=0.d0
            self%gausswidth_n(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
            !self%gausswidth(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth
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
        if(maxval(atoms%cellvec)>max_cellVec) then
            stop 'ERROR: atoms%cellvec > max_cellVec'
        endif
        if(parini%iverbose>=2) call cpu_time(time1)
        if(.not. ann_arr%linear_rho_pot_initiated) then
            !call cpu_time(time1_t)
            !call ann_arr%radpots_cent2%init_radpots_cent2(1.d-3,poisson%rgcut,max_cellVec,parini%ntypat)
            !call cpu_time(time2_t)
            !if(parini%mpi_env%iproc==0) then
            !if(parini%iverbose>=2) write(*,*) 'CENT2_NGP: ',ann_arr%radpots_cent2%ngp
            !endif
            ann_arr%linear_rho_pot_initiated=.true. 
            !do itype=1,parini%ntypat
            !    ann_arr%radpots_cent2%gwe(itype)=ann_arr%ann(itype)%gausswidth
            !    ann_arr%radpots_cent2%gwn(itype)=ann_arr%ann(itype)%gausswidth_ion
            !enddo
            !call ann_arr%radpots_cent2%cal_radpots(parini%screening_factor)
            !call cpu_time(time3_t)
            !if(parini%mpi_env%iproc==0) write(*,*) 'time of init_radpots_cent2 ',time2_t-time1_t
            !if(parini%mpi_env%iproc==0) write(*,*) 'time of cal_radpots ',time3_t-time2_t
        endif
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'get_scf_pot_time: ',time2-time1
        if(parini%iverbose>=2) call cpu_time(time1)
        allocate(self%rho_tmp(poisson%ngpx,poisson%ngpy,poisson%ngpz))
        self%nbgx=int(poisson%rgcut/poisson%hgrid(1,1))+3
        self%nbgy=int(poisson%rgcut/poisson%hgrid(2,2))+3
        self%nbgz=int(poisson%rgcut/poisson%hgrid(3,3))+3
        allocate(self%rho_e_all(-self%nbgx:self%nbgx,-self%nbgy:self%nbgy,-self%nbgz:self%nbgz,atoms%nat))
        allocate(self%rho_n_all(-self%nbgx:self%nbgx,-self%nbgy:self%nbgy,-self%nbgz:self%nbgz,atoms%nat))
        call self%calc_atomic_densities(parini,atoms,ann_arr,poisson)
        call self%get_amat_cent2(parini,ann_arr,atoms,poisson,ann_arr%a)
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'get_amt_time: ',time2-time1
        call get_eigenval(parini,atoms,ann_arr%a)
        allocate(poisson%pot_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))
        if(parini%iverbose>=2) call cpu_time(time1)
        allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
        rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
        poisson%pot_ion(:,:,:)=0.d0
        do iat=1,atoms%nat
            !call ann_arr%radpots_cent2%pot_1dto3d('linear_pot_n',atoms%ratp(1,iat), &
            !    atoms%itypat(iat),poisson,.false.,atoms%zat(iat),poisson%pot_ion)
            !call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),atoms%zat(iat),gausswidth(iat), &
            !    6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            poisson%rho=0.d0
            agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+self%nbgx+0
            agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+self%nbgy+0
            agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+self%nbgz+0
            do iz=agpz-self%nbgz,agpz+self%nbgz
            do iy=agpy-self%nbgy,agpy+self%nbgy
            do ix=agpx-self%nbgx,agpx+self%nbgx
                poisson%rho(ix,iy,iz)=atoms%zat(iat)*self%rho_n_all(ix-agpx,iy-agpy,iz-agpz,iat)
            enddo
            enddo
            enddo
            poisson%pot=0.d0
            call get_hartree(parini,poisson,atoms,self%gausswidth_tmp,tt)
            poisson%pot_ion=poisson%pot_ion+poisson%pot
        enddo
        poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
        deallocate(rho)
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'pot_ion_time: ',time2-time1
    else
        write(*,*) 'ERROR: unknown value for syslinsolver',trim(ann_arr%syslinsolver)
        stop
    endif
end subroutine init_electrostatic_cent2
!*****************************************************************************************
subroutine get_amat_cent2(self,parini,ann_arr,atoms,poisson,amat)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_ann, only: typ_ann_arr
    !use mod_radpots_cent2, only: radial_to_3d
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(out):: amat(1:atoms%nat+1,1:atoms%nat+1)
    !local variables
    integer:: iat, jat, ix, iy, iz
    integer:: agpx, agpy, agpz 
    real(8):: tt
    real(8), allocatable:: gausswidth(:)
    allocate(gausswidth(atoms%nat))
    self%rho_tmp(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    do iat=1,atoms%nat
        !call ann_arr%radpots_cent2%pot_1dto3d('linear_pot_e',atoms%ratp(1,iat),atoms%itypat(iat),poisson,.true.,one,pot)
        poisson%rho=0.d0
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+self%nbgx+0
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+self%nbgy+0
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+self%nbgz+0
        do iz=agpz-self%nbgz,agpz+self%nbgz
        do iy=agpy-self%nbgy,agpy+self%nbgy
        do ix=agpx-self%nbgx,agpx+self%nbgx
            poisson%rho(ix,iy,iz)=self%rho_e_all(ix-agpx,iy-agpy,iz-agpz,iat)
        enddo
        enddo
        enddo
        poisson%pot=0.d0
        call get_hartree(parini,poisson,atoms,gausswidth,tt)
        do jat=1,iat
            agpx=int(atoms%ratp(1,jat)/poisson%hgrid(1,1))+self%nbgx+0
            agpy=int(atoms%ratp(2,jat)/poisson%hgrid(2,2))+self%nbgy+0
            agpz=int(atoms%ratp(3,jat)/poisson%hgrid(3,3))+self%nbgz+0
            tt=0.d0
            do iz=agpz-self%nbgz,agpz+self%nbgz
            do iy=agpy-self%nbgy,agpy+self%nbgy
            do ix=agpx-self%nbgx,agpx+self%nbgx
                tt=tt+self%rho_e_all(ix-agpx,iy-agpy,iz-agpz,jat)*poisson%pot(ix,iy,iz)
            enddo
            enddo
            enddo
            amat(iat,jat)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
            amat(jat,iat)=amat(iat,jat)
        enddo
        amat(iat,iat)=amat(iat,iat)+ann_arr%ann(atoms%itypat(iat))%hardness
        amat(iat,atoms%nat+1)=1.d0
        amat(atoms%nat+1,iat)=1.d0
    enddo !end of loop over iat
    amat(atoms%nat+1,atoms%nat+1)=0.d0
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=self%rho_tmp(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    deallocate(gausswidth)
end subroutine get_amat_cent2
!*****************************************************************************************
subroutine get_qat_from_chi_dir_cent2(parini,ann_arr,atoms,poisson,amat)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: info, iat
    real(8):: tt, one
    real(8), allocatable:: a(:,:)
    real(8), allocatable:: rho(:,:,:)
    real(8):: ggw, ggw_t, alpha, beta, q_tmp(1)
    allocate(a(atoms%nat+1,atoms%nat+1))
    one=1.d0
    a=amat
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%qq(1:atoms%nat+1))
    endif
    call DGETRF(atoms%nat+1,atoms%nat+1,a,atoms%nat+1,ann_arr%ipiv,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
        stop
    endif
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    do iat=1,atoms%nat
        !call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_e',atoms%ratp(1,iat), &
        !    atoms%itypat(iat),poisson,one,poisson%pot_ion,tt)
        ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        alpha=ggw**3/(ggw**3-ggw_t**3)
        beta=-ggw_t**3/(ggw**3-ggw_t**3)
        q_tmp(1)=alpha
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp,ggw, &
            6.d0*ggw,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        q_tmp(1)=beta
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,ggw_t, &
            6.d0*ggw_t,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        call cal_rho_pot_integral_local(atoms%ratp(1,iat),poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rgcut, &
            poisson%rho,poisson%pot_ion,tt)
        ann_arr%qq(iat)=-ann_arr%chi_o(iat)-atoms%zat(iat)*ann_arr%ann(atoms%itypat(iat))%hardness-tt
    enddo !iat
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    deallocate(rho)
    ann_arr%qq(atoms%nat+1)=-1.d0*sum(atoms%zat) !atoms%qtot !ASK Dr this should be -ztot or qat+zat or atoms%qtot
    call DGETRS('N',atoms%nat+1,1,a,atoms%nat+1,ann_arr%ipiv,ann_arr%qq,atoms%nat+1,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat(1:atoms%nat)=ann_arr%qq(1:atoms%nat)
    do iat=1,atoms%nat
        write(20,'(a3,4es18.6)') atoms%sat(iat),atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat),atoms%qat(iat)+atoms%zat(iat)
    enddo
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        call yaml_map('Lagrange',ann_arr%qq(atoms%nat+1))
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%qq)
    endif
    deallocate(a)
end subroutine get_qat_from_chi_dir_cent2
!*****************************************************************************************
subroutine cal_cent2_total_potential(parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    !use mod_radpots_cent2, only: radial_to_3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat
    real(8), allocatable:: gausswidth(:)
    real(8), allocatable:: rho(:,:,:)
    real(8):: ggw, ggw_t, alpha, beta, tt, q_tmp(1)
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(gausswidth(atoms%nat))
    rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    !poisson%pot(:,:,:)=0.d0
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=0.d0
    do iat=1,atoms%nat
        !call ann_arr%radpots_cent2%pot_1dto3d('linear_pot_e',atoms%ratp(1,iat), &
        !    atoms%itypat(iat),poisson,.false.,atoms%qat(iat),poisson%pot)
        ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        alpha=ggw**3/(ggw**3-ggw_t**3)
        beta=-ggw_t**3/(ggw**3-ggw_t**3)
        q_tmp(1)=alpha*atoms%qat(iat)
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp(1),ggw, &
            6.d0*ggw,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        q_tmp(1)=beta*atoms%qat(iat)
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp(1),ggw_t, &
            6.d0*ggw_t,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    enddo
    poisson%pot=0.d0
    call get_hartree(parini,poisson,atoms,gausswidth,tt)
    poisson%pot=poisson%pot+poisson%pot_ion
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    deallocate(rho)
    deallocate(gausswidth)
end subroutine cal_cent2_total_potential
!*****************************************************************************************
subroutine cal_cent2_energy(parini,atoms,ann_arr,epot_c,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
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
    call cal_electrostatic_ann_cent2(parini,atoms,ann_arr,poisson)
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
subroutine cal_electrostatic_ann_cent2(parini,atoms,ann_arr,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_poisson), intent(inout):: poisson
    !local variables
    type(typ_poisson):: poisson_force
    integer:: iat
    real(8):: tt, tte, ttn
    real(8):: alpha, beta, ggw, ggw_t
    real(8),allocatable::fat_t(:,:)
    real(8), allocatable:: gausswidth(:)
    real(8), allocatable:: rho(:,:,:)
    real(8):: q_tmp(1)
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(gausswidth(atoms%nat))
    rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    do iat=1,atoms%nat
        gausswidth(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
    enddo
    if(trim(ann_arr%event)=='potential') then
        !poisson%pot(:,:,:)=0.d0
        poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=0.d0
        do iat=1,atoms%nat
            !call ann_arr%radpots_cent2%pot_1dto3d('linear_pot_e',atoms%ratp(1,iat), &
            !    atoms%itypat(iat),poisson,.false.,atoms%qat(iat),poisson%pot)
            ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
            ggw_t=0.95d0*ggw
            alpha=ggw**3/(ggw**3-ggw_t**3)
            beta=-ggw_t**3/(ggw**3-ggw_t**3)
            q_tmp(1)=alpha*atoms%qat(iat)
            call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp(1),ggw, &
                6.d0*ggw,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            q_tmp(1)=beta*atoms%qat(iat)
            call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp(1),ggw_t, &
                6.d0*ggw_t,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            !call ann_arr%radpots_cent2%pot_1dto3d('linear_pot_n',atoms%ratp(1,iat), &
            !    atoms%itypat(iat),poisson,.false.,atoms%zat(iat),poisson%pot)
            call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),atoms%zat(iat),gausswidth(iat), &
                6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        enddo
        poisson%pot=0.d0
        call get_hartree(parini,poisson,atoms,gausswidth,tt)
    endif
    tt=0.d0
    do iat=1,atoms%nat
        !call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_e',atoms%ratp(1,iat), &
        !    atoms%itypat(iat),poisson,atoms%qat(iat),poisson%pot,tte)
        ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        alpha=ggw**3/(ggw**3-ggw_t**3)
        beta=-ggw_t**3/(ggw**3-ggw_t**3)
        q_tmp(1)=alpha*atoms%qat(iat)
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp(1),ggw, &
            6.d0*ggw,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        q_tmp(1)=beta*atoms%qat(iat)
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp(1),ggw_t, &
            6.d0*ggw_t,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        call cal_rho_pot_integral_local(atoms%ratp(1,iat),poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rgcut, &
            poisson%rho,poisson%pot,tte)
        !call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_n',atoms%ratp(1,iat), &
        !    atoms%itypat(iat),poisson,atoms%zat(iat),poisson%pot,ttn)
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),atoms%zat(iat),gausswidth(iat), &
            6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        call cal_rho_pot_integral_local(atoms%ratp(1,iat),poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rgcut, &
            poisson%rho,poisson%pot,ttn)
        tt=tt+ttn+tte
    enddo
    ann_arr%epot_es=0.5d0*tt
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    deallocate(rho)
    deallocate(gausswidth)
    if(trim(ann_arr%event)=='potential') then
        !Force
        poisson%rcart=atoms%ratp
        poisson_force=poisson
        poisson_force%q(:)=atoms%qat(:)
        poisson_force%gw(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
        poisson_force%gw_ewald(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
        !atoms%fat=0.d0
        allocate(fat_t(1:3,atoms%nat))
        fat_t(1:3,1:atoms%nat)=0.d0
        call force_gto_sym_ortho(parini,poisson_force%bc,poisson_force%nat,poisson_force%rcart, &
            poisson_force%q,poisson_force%gw_ewald,poisson_force%rgcut,poisson_force%xyz111, &
            poisson_force%lda,poisson_force%ngpx, &
            poisson_force%ngpy,poisson_force%ngpz,poisson_force%hgrid,poisson_force%pot,fat_t)
        do iat=1,atoms%nat
            ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
            ggw_t=0.95d0*ggw
            alpha=ggw**3/(ggw**3-ggw_t**3)
            atoms%fat(1:3,iat)=atoms%fat(1:3,iat)+alpha*fat_t(1:3,iat)
        enddo
        poisson%rcart=atoms%ratp
        poisson_force=poisson
        poisson_force%q(:)=atoms%qat(:)
        poisson_force%gw(:)=0.95d0*ann_arr%ann(atoms%itypat(:))%gausswidth
        poisson_force%gw_ewald(:)=0.95d0*ann_arr%ann(atoms%itypat(:))%gausswidth
        fat_t(1:3,1:atoms%nat)=0.d0
        call force_gto_sym_ortho(parini,poisson_force%bc,poisson_force%nat,poisson_force%rcart, &
            poisson_force%q,poisson_force%gw_ewald,poisson_force%rgcut,poisson_force%xyz111, &
            poisson_force%lda,poisson_force%ngpx, &
            poisson_force%ngpy,poisson_force%ngpz,poisson_force%hgrid,poisson_force%pot,fat_t)
        do iat=1,atoms%nat
            ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
            ggw_t=0.95d0*ggw
            beta=-ggw_t**3/(ggw**3-ggw_t**3)
            atoms%fat(1:3,iat)=atoms%fat(1:3,iat)+beta*fat_t(1:3,iat)
        enddo
        deallocate(fat_t)
        poisson_force%q(:)=atoms%zat(:)
        poisson_force%gw(:)=ann_arr%ann(atoms%itypat(:))%gausswidth_ion
        poisson_force%gw_ewald(:)=ann_arr%ann(atoms%itypat(:))%gausswidth_ion
        call force_gto_sym_ortho(parini,poisson_force%bc,poisson_force%nat,poisson_force%rcart, &
            poisson_force%q,poisson_force%gw_ewald,poisson_force%rgcut,poisson_force%xyz111, &
            poisson_force%lda,poisson_force%ngpx, &
            poisson_force%ngpy,poisson_force%ngpz,poisson_force%hgrid,poisson_force%pot,atoms%fat)
    endif
end subroutine cal_electrostatic_ann_cent2
!*****************************************************************************************
subroutine fini_electrostatic_cent2(self,parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    class(typ_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_poisson), intent(inout):: poisson
    call fini_hartree(parini,atoms,poisson)
    deallocate(self%rho_e_all)
    deallocate(self%rho_n_all)
    deallocate(self%rho_tmp)
    deallocate(self%gausswidth_n)
    deallocate(self%gausswidth_tmp)
end subroutine fini_electrostatic_cent2
!*****************************************************************************************
subroutine get_eigenval(parini,atoms,amat)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: info, lwork
    real(8), allocatable:: eigen_val_amat(:,:)
    real(8), allocatable:: real_eigenval(:), work(:)
    lwork=max(atoms%nat**2,100)
    allocate(eigen_val_amat(1:atoms%nat,1:atoms%nat))
    allocate(real_eigenval(1:atoms%nat))
    allocate(work(1:4*atoms%nat))
    eigen_val_amat(1:atoms%nat,1:atoms%nat)=amat(1:atoms%nat,1:atoms%nat)
    call DSYEV('N','U',atoms%nat,eigen_val_amat,atoms%nat,real_eigenval,work,lwork,info)
    if(parini%mpi_env%iproc==0) then
    write(13,'(a6,es18.8,a10,es18.8)') 'MAX: ',maxval(real_eigenval),' | MIN: ',minval(real_eigenval)
    endif
    deallocate(eigen_val_amat)
    deallocate(real_eigenval)
    deallocate(work)
end subroutine get_eigenval
!*****************************************************************************************
subroutine get_dpm(atoms,dpx,dpy,dpz,dpm_err)
    use mod_atoms, only: typ_atoms
    use dynamic_memory
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
    do iat=1, atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
            dpx=dpx+(atoms%ratp(1,iat)-centroid_x)*(atoms%qat(iat)+2.d0)
            dpy=dpy+(atoms%ratp(2,iat)-centroid_y)*(atoms%qat(iat)+2.d0)
            dpz=dpz+(atoms%ratp(3,iat)-centroid_z)*(atoms%qat(iat)+2.d0)
        elseif(trim(atoms%sat(iat))=='O') then
            dpx=dpx+(atoms%ratp(1,iat)-centroid_x)*(atoms%qat(iat)+6.d0)
            dpy=dpy+(atoms%ratp(2,iat)-centroid_y)*(atoms%qat(iat)+6.d0)
            dpz=dpz+(atoms%ratp(3,iat)-centroid_z)*(atoms%qat(iat)+6.d0)
        endif
    enddo
    dpm_err=((dpx-atoms%dpm(1))**2+(dpy-atoms%dpm(2))**2+(dpz-atoms%dpm(3))**2)
end subroutine get_dpm
!*****************************************************************************************
subroutine prefit_cent2(parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use mod_trial_energy, only: typ_trial_energy, trial_energy_deallocate
    use wrapper_MPI, only: fmpi_allreduce, FMPI_SUM, fmpi_wait
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: itrial, iat, jat, ii
    integer:: nbgx, nbgy, nbgz
    real(8):: xyz(3), tt, rmse, err_U_SRS, U_SRS
    real(8):: qavg_Mg, qavg_O, qvar_Mg, qvar_O
    real(8):: cavg_Mg, cavg_O, cvar_Mg, cvar_O
    real(8):: pi
    real(8), allocatable:: EP(:,:)
    real(8), allocatable:: E_all(:)
    real(8), allocatable:: linear_rho_t(:)
    real(8):: ggw, ggw_t, alpha, beta
    real(8), allocatable:: rho(:,:,:)
    real(8), allocatable:: gausswidth(:)
    real(8), allocatable:: EP_n(:)
    real(8):: one
    type(typ_trial_energy), pointer:: trial_energy=>null()
    real(8):: time1, time2
    integer:: mtrial, itrials, itriale, mproc, mat, iats, iate, ireq
    logical, save:: done=.false.
    if(done) return
    one=1.d0
    !call cpu_time(time1)
    call cal_trial_from_cube(parini,trial_energy)
    !call cpu_time(time2)
    !if(parini%mpi_env%iproc==0) then
    !    write(*,*) 'time elapsed in cal_trial_from_cube ',time2-time1
    !endif
    allocate(EP_n(trial_energy%ntrial))
    allocate(E_all(trial_energy%ntrial))
    allocate(EP(atoms%nat,trial_energy%ntrial))
    nbgx=int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy=int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz=int(poisson%rgcut/poisson%hgrid(3,3))+3
    if(parini%mpi_env%iproc==0) then
    write(*,*) 'RGCUT ',poisson%rgcut,nbgx,nbgy,nbgz
    endif
    !-----------------------------------------------------------------
    call cpu_time(time1)
    EP=0.0
    if(parini%mpi_env%nproc>1) then
        mat=atoms%nat/parini%mpi_env%nproc
        iats=parini%mpi_env%iproc*mat+1
        mproc=mod(atoms%nat,parini%mpi_env%nproc)
        iats=iats+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
        if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) mat=mat+1
        iate=iats+mat-1
    else
        iats=1
        iate=atoms%nat
    endif

    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    allocate(gausswidth(atoms%nat))
    rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    do iat=iats,iate
    !do iat=1,atoms%nat
        !call ann_arr%radpots_cent2%pot_1dto3d('linear_pot_e',atoms%ratp(1,iat),atoms%itypat(iat),poisson,.true.,one,poisson%pot)
        ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        alpha=ggw**3/(ggw**3-ggw_t**3)
        beta=-ggw_t**3/(ggw**3-ggw_t**3)
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),alpha,ggw, &
            6.d0*ggw,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),beta,ggw_t, &
            6.d0*ggw_t,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        poisson%pot=0.d0
        call get_hartree(parini,poisson,atoms,gausswidth,tt)
        !do iz=1,poisson%ngpz
        !do iy=1,poisson%ngpy
        !do ix=1,poisson%ngpx
        !write(80+iat,'(3i6,f10.5)') ix,iy,iz,poisson%pot(ix,iy,iz)
        !enddo
        !enddo
        !enddo
        do itrial=1,trial_energy%ntrial
            xyz(1:3)=atoms%ratp(1:3,trial_energy%iat_list(itrial))+trial_energy%disp(1:3,itrial)
            !call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_t',xyz, &
            !    0,poisson,one,poisson%pot,tt)
            call put_gto_sym_ortho(parini,poisson%bc,.true.,1,xyz,1.d0,1.d0, &
                6.d0*1.d0,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
            call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
                poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rgcut, &
                poisson%rho,poisson%pot,tt)
            EP(iat,itrial)=tt
        enddo !end of loop over itrial
    enddo !end of loop over iat
    if(parini%mpi_env%nproc>1) then
    call fmpi_allreduce(EP(1,1),atoms%nat*trial_energy%ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm) !,request=ireq)
    endif
    if(parini%mpi_env%nproc>1) then
        mtrial=trial_energy%ntrial/parini%mpi_env%nproc
        itrials=parini%mpi_env%iproc*mtrial+1
        mproc=mod(trial_energy%ntrial,parini%mpi_env%nproc)
        itrials=itrials+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
        if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) mtrial=mtrial+1
        itriale=itrials+mtrial-1
    else
        itrials=1
        itriale=trial_energy%ntrial
    endif
    EP_n=0.d0
    do itrial=itrials,itriale
    !do itrial=1,trial_energy%ntrial
        xyz(1:3)=atoms%ratp(1:3,trial_energy%iat_list(itrial))+trial_energy%disp(1:3,itrial)
        !call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_t',xyz, &
        !    0,poisson,one,poisson%pot_ion,tt)
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,xyz,1.d0,1.d0, &
            6.d0*1.d0,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        call cal_rho_pot_integral_local(xyz,poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rgcut, &
            poisson%rho,poisson%pot_ion,tt)
        EP_n(itrial)=tt
    enddo !end of loop over itrial
    if(parini%mpi_env%nproc>1) then
    !call fmpi_wait(ireq)
    call fmpi_allreduce(EP_n(1),trial_energy%ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    endif
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    deallocate(rho)
    deallocate(gausswidth)
    call cpu_time(time2)
    !if(parini%mpi_env%iproc==0) then
    write(*,*) 'time EP ',time2-time1
    !endif
    !-----------------------------------------------------------------
    call get_expansion_coeff(parini,ann_arr,atoms,trial_energy%ntrial, &
        trial_energy%energy,EP,EP_n)
    !-----------------------------------------------------------------
    call reverseCEP(parini,ann_arr,atoms,poisson,ann_arr%a)
    atoms%qat(1:atoms%nat)=ann_arr%qq(1:atoms%nat)
    rmse=0.d0
    do itrial=1,trial_energy%ntrial
        tt=0.d0
        do iat=1,atoms%nat
            tt=tt+atoms%qat(iat)*EP(iat,itrial)
        enddo
        E_all(itrial)=EP_n(itrial)+tt
        rmse=rmse+(E_all(itrial)-trial_energy%energy(itrial))**2
    enddo
    rmse=1.d3*sqrt(rmse/trial_energy%ntrial)
    call cal_etrial_cent2(parini,ann_arr,atoms,poisson,U_SRS)
    call prefit_cent2_output(ann_arr,atoms,qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O)
    err_U_SRS=1.d3*(U_SRS-trial_energy%ehartree_scn_excl)/atoms%nat
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,2es24.15)') 'USRS ',U_SRS,trial_energy%ehartree_scn_excl
    write(*,'(a,2f10.3,8f7.3)') 'OPT ',rmse,err_U_SRS, &
        qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O
    do itrial=1,trial_energy%ntrial
        write(*,'(a,i3,2es24.15)') 'ETS ',trial_energy%iat_list(itrial),E_all(itrial),trial_energy%energy(itrial)
    enddo
    endif
    deallocate(EP_n)
    deallocate(E_all)
    call trial_energy_deallocate(trial_energy)
    done=.true.
end subroutine prefit_cent2
!*****************************************************************************************
subroutine get_expansion_coeff(parini,ann_arr,atoms,ntrial,energy,EP,EP_n)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: ntrial
    real(8), intent(in):: energy(ntrial), EP(atoms%nat,ntrial), EP_n(ntrial)
    !local variables
    integer:: iat, jat, itrial, info, itypat, iter, lwork
    real(8):: hh_Mg, hh_O, hh, qtarget_Mg, qtarget_O, qtarget
    real(8):: tt, fract
    real(8), allocatable:: squarefit(:,:), squarefit_t(:,:)
    real(8), allocatable:: real_eigenval(:), work(:)
    lwork=max(atoms%nat**2,100)
    allocate(squarefit(atoms%nat,atoms%nat))
    allocate(squarefit_t(atoms%nat+1,atoms%nat+1))
    allocate(real_eigenval(1:atoms%nat),work(lwork))
    squarefit=0.d0
    squarefit_t=0.d0
    do iat=1,atoms%nat
        do jat=1,atoms%nat
            tt=0.d0
            do itrial=1,ntrial
                tt=tt+2.d0*EP(iat,itrial)*EP(jat,itrial)
            enddo
            squarefit(iat,jat)=tt
            squarefit_t(iat,jat)=tt
        enddo
    enddo
    call DSYEV('N','U',atoms%nat,squarefit,atoms%nat,real_eigenval,work,lwork,info)
    !if(parini%mpi_env%iproc==0) then
    !do iat=1,atoms%nat
    !    write(*,'(a,i6,es14.5)') 'EVAL ',iat,real_eigenval(iat)
    !enddo
    !endif
    do iat=1,atoms%nat
        do jat=1,atoms%nat
            squarefit(iat,jat)=squarefit_t(iat,jat)
        enddo
    enddo
    !hh_Mg=40.d-2
    !hh_O=40.d-2
    !-----------------------------------------------------------------
    !do iter=0,0
    hh_Mg=1.d-3*real_eigenval(atoms%nat)
    hh_O=1.d-3*real_eigenval(atoms%nat)
    !fract=min(1.d-3,1.d-5*real(iter,8))
    !if(iter>0) then
    !hh_Mg=fract*real_eigenval(atoms%nat)
    !hh_O=fract*real_eigenval(atoms%nat)
    !else
    !hh_Mg=0.d0
    !hh_O=0.d0
    !endif
    !do jat=1,atoms%nat
    !    do iat=1,atoms%nat
    !        squarefit_t(iat,jat)=squarefit(iat,jat)
    !    enddo
    !enddo
    !if(iter>0) then
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') hh=hh_Mg
        if(trim(atoms%sat(iat))=='O' ) hh=hh_O
        squarefit(iat,iat)=squarefit(iat,iat)+hh
        squarefit_t(iat,iat)=squarefit_t(iat,iat)+hh
    enddo
    !endif
    call DSYEV('N','U',atoms%nat,squarefit,atoms%nat,real_eigenval,work,lwork,info)
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(*,'(a,i6,es14.5)') 'EVAL ',iat,real_eigenval(iat)
    enddo
    endif
    squarefit_t(1:atoms%nat,atoms%nat+1)=1.d0
    squarefit_t(atoms%nat+1,1:atoms%nat)=1.d0
    squarefit_t(atoms%nat+1,atoms%nat+1)=0.d0
    call DGETRF(atoms%nat+1,atoms%nat+1,squarefit_t,atoms%nat+1,ann_arr%ipiv,info)
    ann_arr%qq(atoms%nat+1)=-sum(atoms%zat)
    do iat=1,atoms%nat
        tt=0.d0
        do itrial=1,ntrial
            tt=tt+2.d0*EP(iat,itrial)*(energy(itrial)-EP_n(itrial))
        enddo
        ann_arr%qq(iat)=tt
    enddo
    do itypat=1,parini%ntypat
        if(trim(parini%stypat(itypat))=='Mg') then
            qtarget_Mg=-ann_arr%ann(itypat)%zion+ann_arr%ann(itypat)%qtarget
            if(parini%mpi_env%iproc==0) then
            write(*,'(a,f8.3)') 'QTARGET_Mg ',qtarget_Mg
            endif
        endif
        if(trim(parini%stypat(itypat))=='O') then
            qtarget_O=-ann_arr%ann(itypat)%zion-ann_arr%ann(itypat)%qtarget
            if(parini%mpi_env%iproc==0) then
            write(*,'(a,f8.3)') 'QTARGET_O  ',qtarget_O
            endif
        endif
    enddo
    !if(iter==0) then
    !    qtarget_Mg=0.d0
    !    qtarget_O=0.d0
    !endif
    !qtarget_Mg=-2.d0+ann_arr%ann(1)%qtarget
    !qtarget_O =-4.d0-ann_arr%ann(1)%qtarget
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') hh=hh_Mg
        if(trim(atoms%sat(iat))=='O' ) hh=hh_O
        if(trim(atoms%sat(iat))=='Mg') qtarget=qtarget_Mg
        if(trim(atoms%sat(iat))=='O' ) qtarget=qtarget_O
        ann_arr%qq(iat)=ann_arr%qq(iat)+hh*qtarget
    enddo
    call DGETRS('N',atoms%nat+1,1,squarefit_t,atoms%nat+1,ann_arr%ipiv,ann_arr%qq,atoms%nat+1,info)
    !qtarget_Mg=0.d0
    !qtarget_O=0.d0
    !do iat=1,atoms%nat
    !    if(trim(atoms%sat(iat))=='Mg') then
    !        !qtarget_Mg=-ann_arr%ann(itypat)%zion+ann_arr%ann(itypat)%qtarget
    !        qtarget_Mg=qtarget_Mg+ann_arr%qq(iat)
    !    endif
    !    if(trim(atoms%sat(iat))=='O') then
    !        !qtarget_O=-ann_arr%ann(itypat)%zion-ann_arr%ann(itypat)%qtarget
    !        qtarget_O=qtarget_O+ann_arr%qq(iat)
    !    endif
    !enddo
    !qtarget_Mg=qtarget_Mg/(atoms%nat/2)
    !qtarget_O=qtarget_O/(atoms%nat/2)
    !if(parini%mpi_env%iproc==0) then
    !do itypat=1,parini%ntypat
    !    if(trim(parini%stypat(itypat))=='Mg') then
    !        write(*,'(a,f8.3,es14.5)') 'QTARGET_Mg ',ann_arr%ann(itypat)%zion+qtarget_Mg,fract
    !    endif
    !    if(trim(parini%stypat(itypat))=='O') then
    !        write(*,'(a,f8.3)') 'QTARGET_O  ',ann_arr%ann(itypat)%zion+qtarget_O
    !    endif
    !enddo
    !endif
    !enddo !end of loop over iter
    !-----------------------------------------------------------------

    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(*,'(a,i6,f7.3)') 'QQQ ',iat,atoms%zat(iat)+ann_arr%qq(iat)
    enddo
    endif
    deallocate(squarefit)
    deallocate(squarefit_t)
    deallocate(real_eigenval,work)
end subroutine get_expansion_coeff
!*****************************************************************************************
subroutine cal_etrial_cent2(parini,ann_arr,atoms,poisson,U_SRS)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    !use mod_radpots_cent2, only: radial_to_3d
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(out):: U_SRS
    !local variables
    integer:: iat, ix, iy, iz
    real(8):: one, tt
    real(8):: ggw, ggw_t, alpha, beta
    real(8), allocatable:: gausswidth(:)
    real(8), allocatable:: qat_tmp(:)
    real(8), allocatable:: rho(:,:,:)
    one=1.d0
    allocate(gausswidth(atoms%nat))
    allocate(qat_tmp(atoms%nat))
    !poisson%pot(:,:,:)=0.d0
    !do iat=1,atoms%nat
    !    call ann_arr%radpots_cent2%pot_1dto3d('linear_pot_e',atoms%ratp(1,iat), &
    !        atoms%itypat(iat),poisson,.false.,atoms%qat(iat),poisson%pot)
    !enddo
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    do iat=1,atoms%nat
    ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
    ggw_t=0.95d0*ggw
    alpha=ggw**3/(ggw**3-ggw_t**3)
    gausswidth(iat)=ggw
    qat_tmp(iat)=atoms%qat(iat)*alpha
    enddo
    call put_gto_sym_ortho(parini,poisson%bc,.true.,atoms%nat,atoms%ratp,qat_tmp,gausswidth, &
        6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    do iat=1,atoms%nat
    ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
    ggw_t=0.95d0*ggw
    beta=-ggw_t**3/(ggw**3-ggw_t**3)
    gausswidth(iat)=ggw_t
    qat_tmp(iat)=atoms%qat(iat)*beta
    enddo
    call put_gto_sym_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,qat_tmp,gausswidth, &
        6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    poisson%pot=0.d0
    call get_hartree(parini,poisson,atoms,gausswidth,tt)

    poisson%pot=poisson%pot+poisson%pot_ion
    gausswidth=0.5d0 !TO_BE_CORRECTED
    atoms%fat=0.d0
    call force_gto_sym_ortho(parini,atoms%boundcond,atoms%nat,atoms%ratp, &
        atoms%zat,gausswidth,6.d0,poisson%xyz111, &
        poisson%ngpx,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%pot,atoms%fat)
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        !write(*,'(a,i4,3es19.10)') 'FAT ',iat,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
        write(61,'(a,i3,3(a2,es24.15),a)') '  - [',iat,', ',atoms%fat(1,iat),', ',atoms%fat(2,iat),', ',atoms%fat(3,iat),']'
    enddo
    endif
    !tt=0.d0
    !do iat=1,atoms%nat
    !    !call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_e',atoms%ratp(1,iat), &
    !    !    atoms%itypat(iat),poisson,atoms%qat(iat),poisson%pot,tte)
    !    call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_n',atoms%ratp(1,iat), &
    !        atoms%itypat(iat),poisson,atoms%zat(iat),poisson%pot,ttn)
    !    !tt=tt+tte+ttn
    !    tt=tt+ttn
    !enddo
    call put_gto_sym_ortho(parini,poisson%bc,.false.,atoms%nat,atoms%ratp,atoms%zat,gausswidth, &
        6.d0*maxval(gausswidth),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
    tt=0.d0
    do iz=1,poisson%ngpz
    do iy=1,poisson%ngpy
    do ix=1,poisson%ngpx
        tt=tt+poisson%rho(ix,iy,iz)*poisson%pot(ix,iy,iz)
    enddo
    enddo
    enddo
    U_SRS=0.5d0*tt*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    deallocate(rho)
    deallocate(qat_tmp)
    deallocate(gausswidth)
end subroutine cal_etrial_cent2
!*****************************************************************************************
subroutine prefit_cent2_output(ann_arr,atoms,qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    real(8), intent(inout):: qavg_Mg, qavg_O, qvar_Mg, qvar_O
    real(8), intent(inout):: cavg_Mg, cavg_O, cvar_Mg, cvar_O
    !local variables
    integer:: iat
    real(8):: q
    real(8):: qmax_Mg, qmin_Mg, qmax_O, qmin_O
    real(8):: cmax_Mg, cmin_Mg, cmax_O, cmin_O
    qavg_Mg=0.d0
    qavg_O=0.d0
    qmax_Mg=-huge(1.d0)
    qmax_O=-huge(1.d0)
    qmin_Mg=huge(1.d0)
    qmin_O=huge(1.d0)
    cavg_Mg=0.d0
    cavg_O=0.d0
    cmax_Mg=-huge(1.d0)
    cmax_O=-huge(1.d0)
    cmin_Mg=huge(1.d0)
    cmin_O=huge(1.d0)
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') then
            q=atoms%zat(iat)+atoms%qat(iat)
            qavg_Mg=qavg_Mg+q
            qmax_Mg=max(qmax_Mg,q)
            qmin_Mg=min(qmin_Mg,q)
            cavg_Mg=cavg_Mg+ann_arr%chi_o(iat)
            cmax_Mg=max(cmax_Mg,ann_arr%chi_o(iat))
            cmin_Mg=min(cmin_Mg,ann_arr%chi_o(iat))

        endif
        if(trim(atoms%sat(iat))=='O') then
            q=atoms%zat(iat)+atoms%qat(iat)
            qavg_O=qavg_O+q
            qmax_O=max(qmax_O,q)
            qmin_O=min(qmin_O,q)
            cavg_O=cavg_O+ann_arr%chi_o(iat)
            cmax_O=max(cmax_O,ann_arr%chi_o(iat))
            cmin_O=min(cmin_O,ann_arr%chi_o(iat))
        endif
    enddo
    qavg_Mg=qavg_Mg/(atoms%nat/2.d0)
    qavg_O=qavg_O/(atoms%nat/2.d0)
    qvar_Mg=qmax_Mg-qmin_Mg
    qvar_O=qmax_O-qmin_O
    cavg_Mg=cavg_Mg/(atoms%nat/2.d0)
    cavg_O=cavg_O/(atoms%nat/2.d0)
    cvar_Mg=cmax_Mg-cmin_Mg
    cvar_O=cmax_O-cmin_O
end subroutine prefit_cent2_output
!*****************************************************************************************
subroutine reverseCEP(parini,ann_arr,atoms,poisson,amat)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_electrostatics, only: typ_poisson
    use mod_ann_io_yaml, only: write_yaml_conf_train
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: iat, jat
    real(8):: tt, tt1, tt2, one
    real(8):: ggw, ggw_t, alpha, beta, q_tmp(1)
    type(typ_file_info):: file_info
    real(8), allocatable:: ww(:)
    real(8), allocatable:: rho(:,:,:)
    one=1.d0
    allocate(ww(atoms%nat))
    allocate(rho(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    do iat=1,atoms%nat
        !call ann_arr%radpots_cent2%energy_1dto3d('linear_rho_e',atoms%ratp(1,iat), &
        !    atoms%itypat(iat),poisson,one,poisson%pot_ion,tt)
        ggw=ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        alpha=ggw**3/(ggw**3-ggw_t**3)
        beta=-ggw_t**3/(ggw**3-ggw_t**3)
        q_tmp(1)=alpha
        call put_gto_sym_ortho(parini,poisson%bc,.true.,1,atoms%ratp(1,iat),q_tmp,ggw, &
            6.d0*ggw,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        q_tmp(1)=beta
        call put_gto_sym_ortho(parini,poisson%bc,.false.,1,atoms%ratp(1,iat),q_tmp,ggw_t, &
            6.d0*ggw_t,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rho)
        call cal_rho_pot_integral_local(atoms%ratp(1,iat),poisson%xyz111, &
            poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,poisson%rgcut, &
            poisson%rho,poisson%pot_ion,tt)
        ww(iat)=-atoms%zat(iat)*ann_arr%ann(atoms%itypat(iat))%hardness-tt
    enddo
    do iat=1,atoms%nat
        tt=0.d0
        do jat=1,atoms%nat
            tt=tt+amat(iat,jat)*ann_arr%qq(jat)
        enddo
        ann_arr%chi_o(iat)=-tt+ww(iat)
    enddo
    tt1=0.d0
    tt2=0.d0
    do iat=1,atoms%nat
        tt1=tt1+ann_arr%chi_o(iat)*(ann_arr%qq(iat)+atoms%zat(iat))
        tt2=tt2+(ann_arr%qq(iat)+atoms%zat(iat))**2*0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness
    enddo
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,3f10.3)') 'ECH ',tt1,tt2,tt1+tt2
    do iat=1,atoms%nat
        !write(*,'(a,i4,2f7.3)') 'CHI ',iat,ann_arr%chi_o(iat),atoms%zat(iat)+ann_arr%qq(iat)
        write(51,'(a,i3,a,es19.10,a)') '  - [',iat,', ',ann_arr%chi_o(iat),']'
    enddo
    file_info%filename_positions='posout.yaml'
    file_info%print_force=.true.
    file_info%file_position='new'
    call write_yaml_conf_train(file_info,atoms,ann_arr,.true.)
    endif
    deallocate(ww)
    poisson%rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    deallocate(rho)
end subroutine reverseCEP
!*****************************************************************************************
subroutine cal_trial_from_cube(parini,trial_energy)
    use mod_parini, only: typ_parini
    use mod_trial_energy, only: typ_trial_energy, trial_energy_allocate
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, atom_deallocate_old
    use wrapper_MPI, only: fmpi_allreduce, FMPI_SUM
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_trial_energy), pointer, intent(inout):: trial_energy
    !local variables
    type(typ_poisson):: poisson
    type(typ_poisson):: poisson_scn
    type(typ_poisson):: poisson_ion
    type(typ_atoms):: atoms
    integer:: istat, igpx, igpy, igpz, iat, ntrial, itrial, nsegx, nsegy, nsegz
    integer:: jgpx, jgpy, jgpz
    real(8):: rgcut_a, qtot, pi!, qtot_e, qtot_i
    real(8):: epot_scn, ehartree_scn_excl, tt1, tt2
    real(8):: xyz(3), dxyz(3), epot_trial, gwt
    real(8):: dx, dy, dz, r2, coeff, rloc, c1, c2
    real(8):: xmin, ymin, zmin, xmax, ymax, zmax
    real(8):: q_one(1), gw_one(1)
    real(8), allocatable::  gausswidth(:)
    real(8), allocatable::  rat_trial(:,:)
    integer:: nbgpx, nbgpy, nbgpz, ix, iy, iz
    integer:: mtrial, itrials, itriale, mproc
    !real(8) :: xyz(3)
    !integer:: ny,nz
    !real(8):: max_ion, max_ele
    pi=4.d0*atan(1.d0)
    !if(.not. parini%gaussian_width>0.d0) then
    !    stop 'ERROR: gaussian_width must be set.'
    !endif
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
    atoms%boundcond='free'
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
    !if(trim(atoms%sat(iat))=='Mg') gwt=0.6d0
    !if(trim(atoms%sat(iat))=='O' ) gwt=0.3d0
    gwt=gausswidth(iat)
    coeff=atoms%zat(iat)/(gwt**3*pi**1.5d0)
    !coeff=atoms%zat(iat)*2.d0/(3.d0*gwt**5*pi**1.5d0)
    !coeff=atoms%zat(iat)*4.d0/(15.d0*gwt**7*pi**1.5d0)
    if(trim(atoms%sat(iat))=='Mg') then
        !0.65406138674  2 -5.223929095  0.913704167481045 rloc nloc c1 .. cnloc
        rloc=0.65406138674d0
        c1=-5.223929095d0
        c2=0.913704167481045d0
    endif
    if(trim(atoms%sat(iat))=='O') then
        !0.3454999999    2 -11.7435870154  1.90653967947 rloc nloc c1 .. cnloc
        rloc=0.3454999999d0
        c1=-11.7435870154d0
        c2=1.90653967947d0
    endif
    !c1=-c1
    !c2=-c2
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
        !if(r2<10.d0**2*rloc**2) then
            poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*exp(-r2/gwt**2)
            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*r2*exp(-r2/gwt**2)
            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+coeff*r2**2*exp(-r2/gwt**2)
            !tt1=exp(-r2/(2.d0*rloc**2))
            !tt2=-3.d0*rloc**4*(c1-2.d0*c2)+rloc**2*(c1-7.d0*c2)*r2+c2*r2**2+rloc**3*sqrt(2.d0/pi)*atoms%zat(iat)
            !poisson_ion%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)+tt1*tt2/(4.d0*rloc**6*pi)
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
        !write(33,'(3i5,f8.5)') igpx,igpy,igpz,poisson_ion%rho(igpx,igpy,igpz)
    enddo
    enddo
    enddo
    tt1=tt1*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
    if(parini%mpi_env%iproc==0) then
    write(*,*) 'TT1 ',tt1
    endif
    !stop 'WWWWWWWWWWWWWWW'
    !-------------------------------------------------------
    !call get_rat_iat(atoms,1,xyz)
    !ny=int(xyz(2)/poisson%hgrid(2,2)) 
    !nz=int(xyz(3)/poisson%hgrid(3,3)) 
    !max_ion=maxval(poisson_ion%rho(:,:,:))
    !max_ele=maxval(poisson%rho(:,:,:))
    !do igpz=1,poisson%ngpz
    !    do igpy=1,poisson%ngpy
    !        do igpx=1,poisson%ngpx
    !            if (poisson_ion%rho(igpx,igpy,igpz)==max_ion) write(*,*) 'ind max_ion', igpx,igpy,igpz
    !            if (poisson%rho(igpx,igpy,igpz)==max_ele) then
    !                write(*,*) 'ind max_ele', igpx,igpy,igpz
    !                ny = igpy
    !                nz = igpz
    !            end if
    !        end do
    !    end do
    !end do
    !do igpx=1,poisson%ngpx
    !    write(1370,'(3es14.6)') igpx*poisson%hgrid(1,1),poisson%rho(igpx,ny,nz),poisson_ion%rho(igpx,ny,nz)
    !end do
    !-------------------------------------------------------
    qtot=0.d0
    !qtot_e=0.d0
    !qtot_i=0.d0
    do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
            do igpx=1,poisson%ngpx
    !            qtot_e=qtot_e+poisson%rho(igpx,igpy,igpz)
    !            qtot_i=qtot_i+poisson_ion%rho(igpx,igpy,igpz)
                poisson%rho(igpx,igpy,igpz)=poisson_ion%rho(igpx,igpy,igpz)-poisson%rho(igpx,igpy,igpz)
                !poisson%rho(igpx,igpy,igpz)=-poisson%rho(igpx,igpy,igpz)
                qtot=qtot+poisson%rho(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    !do igpx=1,poisson%ngpx
    !    write(1371,'(2es14.6)') igpx*poisson%hgrid(1,1),poisson%rho(igpx,ny,nz)
    !end do
    qtot=qtot*poisson_ion%hgrid(1,1)*poisson_ion%hgrid(2,2)*poisson_ion%hgrid(3,3)
    if(parini%mpi_env%iproc==0) then
    write(*,*) 'qtot= ',qtot
    endif
    !write(*,*) 'qtot_e= ',qtot_e
    !write(*,*) 'qtot_i= ',qtot_i
    !-------------------------------------------------------
    call update_ratp(atoms)
    call get_hartree(parini,poisson,atoms,gausswidth,ehartree_scn_excl)
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,es24.15,es14.5)') 'ehartree_scn_excl ',ehartree_scn_excl,poisson%screening_factor
    endif
    atoms%fat=0.d0
    call force_gto_sym_ortho(parini,poisson_ion%bc,atoms%nat,poisson_ion%rcart, &
        poisson_ion%q,gausswidth,6.d0,poisson_ion%xyz111, &
        poisson_ion%ngpx,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz, &
        poisson_ion%hgrid,poisson%pot,atoms%fat)
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(*,'(a,i4,3es19.10)') 'FAT ',iat,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
    enddo
    endif
    !-------------------------------------------------------
    poisson_ion%gw(1:poisson_ion%nat)=1.d0
    xmin= huge(1.d0)
    ymin= huge(1.d0)
    zmin= huge(1.d0)
    xmax=-huge(1.d0)
    ymax=-huge(1.d0)
    zmax=-huge(1.d0)
    do iat=1,atoms%nat
        if(poisson_ion%rcart(1,iat)<xmin) xmin=poisson_ion%rcart(1,iat)
        if(poisson_ion%rcart(2,iat)<ymin) ymin=poisson_ion%rcart(2,iat)
        if(poisson_ion%rcart(3,iat)<zmin) zmin=poisson_ion%rcart(3,iat)
        if(poisson_ion%rcart(1,iat)>xmax) xmax=poisson_ion%rcart(1,iat)
        if(poisson_ion%rcart(2,iat)>ymax) ymax=poisson_ion%rcart(2,iat)
        if(poisson_ion%rcart(3,iat)>zmax) zmax=poisson_ion%rcart(3,iat)
    enddo
    dxyz(1)=1.5d0
    dxyz(2)=1.5d0
    dxyz(3)=1.5d0
    nsegx=int((poisson%hgrid(1,1)*poisson%ngpx-16.d0)/dxyz(1))+1
    nsegy=int((poisson%hgrid(2,2)*poisson%ngpy-16.d0)/dxyz(2))+1
    nsegz=int((poisson%hgrid(3,3)*poisson%ngpz-16.d0)/dxyz(3))+1
    dxyz(1)=(poisson%hgrid(1,1)*poisson%ngpx-16.d0)/real(nsegx,kind=8)
    dxyz(2)=(poisson%hgrid(2,2)*poisson%ngpy-16.d0)/real(nsegy,kind=8)
    dxyz(3)=(poisson%hgrid(3,3)*poisson%ngpz-16.d0)/real(nsegz,kind=8)
    !nseg=10
    !ntrial=(nseg+1)**3
    ntrial=(nsegx+1)*(nsegy+1)*(nsegz+1)
    call trial_energy_allocate(ntrial,trial_energy)
    trial_energy%ehartree_scn_excl=ehartree_scn_excl
    if(parini%mpi_env%iproc==0) then
    write(*,'(a,3i3,i6)') 'nsegx,nsegy,nsegz,ntrial ',nsegx,nsegy,nsegz,ntrial
    endif
    !ntrial=atoms%nat
    allocate(rat_trial(3,ntrial))
    itrial=0
    do iz=0,nsegz
    do iy=0,nsegy
    do ix=0,nsegx
        itrial=itrial+1
        rat_trial(1,itrial)=8.d0+dxyz(1)*ix
        rat_trial(2,itrial)=8.d0+dxyz(2)*iy
        rat_trial(3,itrial)=8.d0+dxyz(3)*iz
    enddo
    enddo
    enddo
    !itrial=0
    !do iat=1,atoms%nat
    !    itrial=itrial+1
    !    rat_trial(1,itrial)=poisson_ion%rcart(1,iat)
    !    rat_trial(2,itrial)=poisson_ion%rcart(2,iat)
    !    rat_trial(3,itrial)=poisson_ion%rcart(3,iat)
    !enddo
    if(parini%mpi_env%nproc>1) then
        mtrial=ntrial/parini%mpi_env%nproc
        itrials=parini%mpi_env%iproc*mtrial+1
        mproc=mod(ntrial,parini%mpi_env%nproc)
        itrials=itrials+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
        if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) mtrial=mtrial+1
        itriale=itrials+mtrial-1
    else
        itrials=1
        itriale=ntrial
    endif
    !write(*,'(a,4i8)') 'iproc,itrials,itriale,ntrial ',parini%mpi_env%iproc,itrials,itriale,ntrial
    trial_energy%energy=0.d0
    trial_energy%disp=0.d0
    trial_energy%iat_list=0
    !do itrial=1,ntrial
    do itrial=itrials,itriale
        xyz(1)=rat_trial(1,itrial)-poisson_ion%rcart(1,1)
        xyz(2)=rat_trial(2,itrial)-poisson_ion%rcart(2,1)
        xyz(3)=rat_trial(3,itrial)-poisson_ion%rcart(3,1)
        !    put_gto_sym_ortho(parini,bc,reset,nat,rxyz,qat,gw,rgcut,xyz111,ngx,ngy,ngz,hgrid,rho)
        q_one(1)=1.d0
        gw_one(1)=1.d0
        call put_gto_sym_ortho(parini,poisson_ion%bc,.true.,1,rat_trial(1,itrial),q_one,gw_one, &
            6.d0,poisson_ion%xyz111,poisson_ion%ngpx,poisson_ion%ngpy,poisson_ion%ngpz,poisson_ion%hgrid,poisson_ion%rho)
        epot_trial=0.d0
        do igpz=1,poisson%ngpz
        do igpy=1,poisson%ngpy
        do igpx=1,poisson%ngpx
            epot_trial=epot_trial+poisson_ion%rho(igpx,igpy,igpz)*poisson%pot(igpx,igpy,igpz)
        enddo
        enddo
        enddo
        epot_trial=epot_trial*(poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3))
        trial_energy%energy(itrial)=epot_trial
        trial_energy%disp(1,itrial)=xyz(1)
        trial_energy%disp(2,itrial)=xyz(2)
        trial_energy%disp(3,itrial)=xyz(3)
        trial_energy%iat_list(itrial)=1
        !if(ix==0 .and. iy==0 .and. iz==0) then
        !    epot_trial0=epot_trial
        !    !write(*,'(a,i5,es24.15,es14.5)') 'iat,epot_trial ',iat,epot_trial,poisson%screening_factor
        !else
        !    write(*,'(a,i5,4es24.15,es14.5)') 'iat,epot_trial ',iat,dxyz(1),dxyz(2),dxyz(3),epot_trial-epot_trial0,poisson%screening_factor
            !if(parini%mpi_env%iproc==0) then
            !write(*,'(a,i5,4es24.15,es14.5)') 'iat,epot_trial ',1,xyz(1),xyz(2),xyz(3),epot_trial,poisson%screening_factor
            !write(71,'(a,i3,4(a2,es24.15),a,es14.5,a)') '  - [',1,', ',xyz(1),', ',xyz(2),', ',xyz(3),', ',epot_trial,', ',poisson%screening_factor,']'
            !endif
        !endif
    enddo
    if(parini%mpi_env%nproc>1) then
    call fmpi_allreduce(trial_energy%energy(1),ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(trial_energy%disp(1,1),3*ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    call fmpi_allreduce(trial_energy%iat_list(1),ntrial,op=FMPI_SUM,comm=parini%mpi_env%mpi_comm)
    endif

    if(parini%mpi_env%iproc==0) then
    write(*,'(a,6f8.1)') 'MINMAX ',xmin,ymin,zmin,xmax,ymax,zmax
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(1,1)*poisson%ngpx
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(2,2)*poisson%ngpy
    write(*,'(a,1f8.1)') 'BBBBBB ',poisson%hgrid(3,3)*poisson%ngpz
    endif
    !-------------------------------------------------------
    !call cube_write('total_rho.cube',atoms,poisson,'rho')
    !call cube_write('total_pot.cube',atoms,poisson,'pot')
    call fini_hartree(parini,atoms,poisson)
    call fini_hartree(parini,atoms,poisson_ion)
    call atom_deallocate_old(atoms)
end subroutine cal_trial_from_cube
!*****************************************************************************************
subroutine cal_rho_pot_integral_local(xyz,xyz111,ngpx,ngpy,ngpz,hgrid,rgcut,rho,pot,ener)
    implicit none
    integer, intent(in):: ngpx, ngpy, ngpz
    real(8), intent(in):: hgrid(3,3), xyz(3), xyz111(3), rgcut
    real(8), intent(in):: rho(ngpx,ngpy,ngpz), pot(ngpx,ngpy,ngpz)
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
                res=res+rho(igpx,igpy,igpz)*pot(igpx,igpy,igpz)
            enddo
        enddo
    enddo
    ener=res*(hgrid(1,1)*hgrid(2,2)*hgrid(3,3))
end subroutine cal_rho_pot_integral_local
!*****************************************************************************************
end module mod_cent2
!*****************************************************************************************
