!*****************************************************************************************
module mod_cent2
    implicit none
    private
    public:: typ_cent2
    type typ_cent2
        contains
        !procedure, public, pass(self)::
        procedure, public, nopass:: cal_ann_cent2
    end type typ_cent2
contains
!*****************************************************************************************
subroutine cal_ann_cent2(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_electrostatics, only: typ_poisson
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_poisson):: poisson
    !local variables
    type(typ_pia_arr):: pia_arr_tmp
    integer:: iat, i, j, ng
    real(8):: epot_c, out_ann
    real(8):: dpx, dpy, dpz
    real(8):: time1, time2, time3, time4, time5, time6, time7, timet1, timet2
    real(8):: tt1, hinv(3,3), vol
    call f_routine(id='cal_ann_cent2')
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
    call init_electrostatic_cent2(parini,atoms,ann_arr,ann_arr%a,poisson)
    if(parini%iverbose>=2) call cpu_time(time2)
    if(parini%iverbose>=2) write(*,*) 'init_time: ',time2-time1
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
    over_iat: do iat=1,atoms%nat
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
    if(parini%iverbose>=2)  then
        do iat=1,atoms%nat
            write(166,*)iat,ann_arr%chi_o(iat)
            write(99,'(i4,a3,3f6.2)') iat,trim(atoms%sat(iat)),atoms%qat(iat),atoms%zat(iat),atoms%zat(iat)+atoms%qat(iat)
        enddo
        write(99,*) '==='
        write(99,*) sum(atoms%qat),sum(atoms%zat)
        write(99,*) '---'
    endif
    if(parini%iverbose>=2) call cpu_time(timet1)
    if(parini%iverbose>=2) write(*,*) 'qet_qat_time: ',timet1-time4
    if(.not. trim(ann_arr%event)=='potential') then
        call cal_cent2_total_potential(parini,ann_arr,atoms,poisson)
    endif
    if(parini%prefit_ann) then
        call prefit_cent2(parini,ann_arr,atoms,poisson)
    endif
    if(parini%iverbose>=2) call cpu_time(time5)
    if(parini%iverbose>=2) write(*,*) 'cent2_g_per_time: ',time5-timet1
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call cal_cent2_energy(parini,atoms,ann_arr,epot_c,poisson)
    if(parini%iverbose>=2) call cpu_time(timet2)
    if(parini%iverbose>=2) write(*,*) 'get_electrostatic_time: ',timet2-time6
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
    if(parini%iverbose>=2) then
        write(1390,'(3es18.8,a3,3es18.8,a3,es18.8)')dpx,dpy,dpz,' | ',atoms%dpm(1),atoms%dpm(2),atoms%dpm(3),' | ',ann_arr%dpm_err
    endif
    atoms%dpm(1)=dpx
    atoms%dpm(2)=dpy
    atoms%dpm(3)=dpz
    call fini_electrostatic_cent2(parini,ann_arr,atoms,poisson)
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
    call f_release_routine()
end subroutine cal_ann_cent2
!*****************************************************************************************
subroutine init_electrostatic_cent2(parini,atoms,ann_arr,amat,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, set_typat
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: amat(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, itype
    real(8):: hgp
    real(8):: time1, time2
    real(8):: max_cellVec
    real(8),allocatable:: gausswidth(:)
    associate(epot_es=>ann_arr%epot_es)
    !pi=4.d0*atan(1.d0)
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
        allocate(gausswidth(1:atoms%nat))
        do iat=1,atoms%nat
            gausswidth(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth
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
        call init_hartree(parini,atoms,poisson,gausswidth)
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'init_hartree_time: ',time2-time1
        hgp=1.d-3
        !poisson%ngp=ceiling(2.d0*(poisson%rgcut+parini%max_cellVec)/hgp)
        max_cellVec=100.d0
        if(maxval(atoms%cellvec)>max_cellVec) then
            stop 'ERROR: atoms%cellvec > max_cellVec'
        endif
        poisson%ngp=ceiling(2.d0*(poisson%rgcut+max_cellVec)/hgp)
        if(.not. poisson%linear_allocated) then
            poisson%linear_allocated=.true.
            allocate(poisson%linear_rho_e(0:poisson%ngp,1:parini%ntypat))
            allocate(poisson%linear_pot_e(0:poisson%ngp,1:parini%ntypat))
            allocate(poisson%linear_rho_n(0:poisson%ngp,1:parini%ntypat))
            allocate(poisson%linear_pot_n(0:poisson%ngp,1:parini%ntypat))
            poisson%linear_rho_e(0:poisson%ngp,1:parini%ntypat)=0.d0
            poisson%linear_pot_e(0:poisson%ngp,1:parini%ntypat)=0.d0
            poisson%linear_rho_n(0:poisson%ngp,1:parini%ntypat)=0.d0
            poisson%linear_pot_n(0:poisson%ngp,1:parini%ntypat)=0.d0
        else
            poisson%linear_rho_e(0:poisson%ngp,1:parini%ntypat)=0.d0
            poisson%linear_pot_e(0:poisson%ngp,1:parini%ntypat)=0.d0
            poisson%linear_rho_n(0:poisson%ngp,1:parini%ntypat)=0.d0
            poisson%linear_pot_n(0:poisson%ngp,1:parini%ntypat)=0.d0
        endif
        if(parini%iverbose>=2) call cpu_time(time1)
        if(.not. ann_arr%linear_rho_pot_initiated) then
            if(parini%iverbose>=2) write(*,*) 'CENT2_NGP: ',poisson%ngp
            ann_arr%linear_rho_pot_initiated=.true. 
            allocate(ann_arr%linear_rho_e(0:poisson%ngp,parini%ntypat))
            allocate(ann_arr%linear_rho_n(0:poisson%ngp,parini%ntypat))
            allocate(ann_arr%linear_pot_e(0:poisson%ngp,parini%ntypat))
            allocate(ann_arr%linear_pot_n(0:poisson%ngp,parini%ntypat))
            do itype=1,parini%ntypat
                call get_scf_pot_cent2_twogauss(poisson%ngp,ann_arr%ann(itype)%gausswidth,&
                                       parini%screening_factor,poisson%linear_rho_e(0,itype),poisson%linear_pot_e(0,itype))
                call get_scf_pot_cent2_onegauss(poisson%ngp,ann_arr%ann(itype)%gausswidth_ion,&
                                       parini%screening_factor,poisson%linear_rho_n(0,itype),poisson%linear_pot_n(0,itype))
                ann_arr%linear_rho_e(0:poisson%ngp,itype)=poisson%linear_rho_e(0:poisson%ngp,itype)
                ann_arr%linear_rho_n(0:poisson%ngp,itype)=poisson%linear_rho_n(0:poisson%ngp,itype)
                ann_arr%linear_pot_e(0:poisson%ngp,itype)=poisson%linear_pot_e(0:poisson%ngp,itype)
                ann_arr%linear_pot_n(0:poisson%ngp,itype)=poisson%linear_pot_n(0:poisson%ngp,itype)
            enddo
        else
            do itype=1,parini%ntypat
                poisson%linear_rho_e(0:poisson%ngp,itype)=ann_arr%linear_rho_e(0:poisson%ngp,itype)
                poisson%linear_rho_n(0:poisson%ngp,itype)=ann_arr%linear_rho_n(0:poisson%ngp,itype)
                poisson%linear_pot_e(0:poisson%ngp,itype)=ann_arr%linear_pot_e(0:poisson%ngp,itype)
                poisson%linear_pot_n(0:poisson%ngp,itype)=ann_arr%linear_pot_n(0:poisson%ngp,itype)
            enddo
        endif
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'get_scf_pot_time: ',time2-time1
        call update_ratp(atoms)
        if(parini%iverbose>=2) call cpu_time(time1)
        call get_amat_cent2(ann_arr,atoms,poisson,amat)
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'get_amt_time: ',time2-time1
        call get_eigenval(atoms,amat)
        allocate(poisson%pot_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))
        hgp=1.d-3
        if(parini%iverbose>=2) call cpu_time(time1)
        poisson%pot_ion(:,:,:)=0.d0
        do iat=1,atoms%nat
            call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_n(0,atoms%itypat(iat)), &
                atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,.false.,atoms%zat(iat),poisson%pot_ion)
        enddo
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'pot_ion_time: ',time2-time1
    else
        write(*,*) 'ERROR: unknown value for syslinsolver',trim(ann_arr%syslinsolver)
        stop
    endif
    deallocate(gausswidth)
    end associate
end subroutine init_electrostatic_cent2
!*****************************************************************************************
subroutine get_amat_cent2(ann_arr,atoms,poisson,a)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_poisson
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(out):: a(1:atoms%nat+1,1:atoms%nat+1)
    !local variables
    integer:: iat, jat, ix, iy, iz
    integer:: nbgx, nbgy, nbgz, linearGridNumber
    integer:: agpx, agpy, agpz 
    real(8):: dx, dy, dz, dr, hgp, tt
    real(8):: rho_e, rho_e_p1, one
    real(8), allocatable:: rho_all(:,:,:,:), pot(:,:,:)
    one=1.d0
    hgp=1.d-3
    nbgx=int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy=int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz=int(poisson%rgcut/poisson%hgrid(3,3))+3
    allocate(rho_all(-nbgx:nbgx,-nbgy:nbgy,-nbgz:nbgz,atoms%nat))
    allocate(pot(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do iat=1,atoms%nat
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+nbgx+0
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+nbgy+0
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+nbgz+0
        do iz=agpz-nbgz,agpz+nbgz
        do iy=agpy-nbgy,agpy+nbgy
        do ix=agpx-nbgx,agpx+nbgx
            dx=poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
            dy=poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
            dz=poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
            dr=sqrt(dx**2+dy**2+dz**2)
            linearGridNumber=floor(dr/hgp)
            rho_e_p1=poisson%linear_rho_e(linearGridNumber+1,atoms%itypat(iat))
            rho_e=poisson%linear_rho_e(linearGridNumber,atoms%itypat(iat))
            rho_all(ix-agpx,iy-agpy,iz-agpz,iat)=(dr/hgp-linearGridNumber)*(rho_e_p1-rho_e)+rho_e
        enddo
        enddo
        enddo
    enddo !end of loop over iat
    do iat=1,atoms%nat
        call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,.true.,one,pot)
        do jat=1,iat
            agpx=int(atoms%ratp(1,jat)/poisson%hgrid(1,1))+nbgx+0
            agpy=int(atoms%ratp(2,jat)/poisson%hgrid(2,2))+nbgy+0
            agpz=int(atoms%ratp(3,jat)/poisson%hgrid(3,3))+nbgz+0
            tt=0.d0
            do iz=agpz-nbgz,agpz+nbgz
            do iy=agpy-nbgy,agpy+nbgy
            do ix=agpx-nbgx,agpx+nbgx
                tt=tt+rho_all(ix-agpx,iy-agpy,iz-agpz,jat)*pot(ix,iy,iz)
            enddo
            enddo
            enddo
            a(iat,jat)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
            a(jat,iat)=a(iat,jat)
        enddo
        a(iat,iat)=a(iat,iat)+ann_arr%ann(atoms%itypat(iat))%hardness
        a(iat,atoms%nat+1)=1.d0
        a(atoms%nat+1,iat)=1.d0
    enddo !end of loop over iat
    a(atoms%nat+1,atoms%nat+1)=0.d0
    deallocate(rho_all)
    deallocate(pot)
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
    type(typ_poisson), intent(in):: poisson
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: info, iat
    real(8):: hgp, tt, one
    real(8):: a(atoms%nat+1,atoms%nat+1)
    associate(nat=>atoms%nat)
    one=1.d0
    a=amat
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%qq(1:nat+1))
    endif
    call DGETRF(nat+1,nat+1,a,nat+1,ann_arr%ipiv,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
        stop
    endif
    hgp=1.d-3
    do iat=1,atoms%nat
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,one,poisson%pot_ion,tt)
        ann_arr%qq(iat)=-ann_arr%chi_o(iat)-atoms%zat(iat)*ann_arr%ann(atoms%itypat(iat))%hardness-tt
    enddo !iat
    ann_arr%qq(nat+1)=-1.d0*sum(atoms%zat) !atoms%qtot !ASK Dr this should be -ztot or qat+zat or atoms%qtot
    call DGETRS('N',nat+1,1,a,nat+1,ann_arr%ipiv,ann_arr%qq,nat+1,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat(1:nat)=ann_arr%qq(1:nat)
    do iat=1,nat
        write(20,'(a3,4es18.6)') atoms%sat(iat),atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat),atoms%qat(iat)+atoms%zat(iat)
    enddo
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        call yaml_map('Lagrange',ann_arr%qq(nat+1))
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%qq)
    endif
    end associate
end subroutine get_qat_from_chi_dir_cent2
!*****************************************************************************************
subroutine cal_cent2_total_potential(parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat
    real(8):: hgp
    hgp=1.d-3
    poisson%pot(:,:,:)=0.d0
    do iat=1,atoms%nat
        call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,.false.,atoms%qat(iat),poisson%pot)
    enddo
    poisson%pot=poisson%pot+poisson%pot_ion
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
    if(parini%iverbose>=2)  then
        write(12,*) ann_arr%epot_es,atoms%epot
    endif
    epot_c=ann_arr%epot_es+tt1+tt2+ann_arr%ener_ref
    if(parini%iverbose>=2)  then
        write(16,'(5es18.8)') epot_c,ann_arr%epot_es,tt1,tt2,ann_arr%ener_ref
        write(14,*) epot_c,atoms%epot
    endif
end subroutine cal_cent2_energy
!*****************************************************************************************
subroutine cal_electrostatic_ann_cent2(parini,atoms,ann_arr,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
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
    real(8):: hgp
    real(8):: alpha,beta,ggw,ggw_t
    real(8),allocatable::fat_t(:,:)
    hgp=1.d-3
    if(trim(ann_arr%event)=='potential') then
        poisson%pot(:,:,:)=0.d0
        do iat=1,atoms%nat
            call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
                atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,.false.,atoms%qat(iat),poisson%pot)
            call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_n(0,atoms%itypat(iat)), &
                atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,.false.,atoms%zat(iat),poisson%pot)
        enddo
    endif
    tt=0.d0
    do iat=1,atoms%nat
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,atoms%qat(iat),poisson%pot,tte)
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_n(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,atoms%zat(iat),poisson%pot,ttn)
        tt=tt+ttn+tte
    enddo
    ann_arr%epot_es=0.5d0*tt
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
subroutine fini_electrostatic_cent2(parini,ann_arr,atoms,poisson)
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
    call fini_hartree(parini,atoms,poisson)
end subroutine fini_electrostatic_cent2
!*****************************************************************************************
subroutine get_scf_pot_cent2_onegauss(ngp,gw,scf,rho,pot)
    use mod_qat_target, only: cal_powern_screened_poisson_gaussian, cal_pot_hartree
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: gw, scf
    real(8), intent(out):: rho(0:ngp)
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp
    real(8):: den_coeff, hgp
    real(8):: rrad(0:ngp), pot_scn(0:ngp),weight (0:ngp)
    real(8):: pi
    pi=4.d0*atan(1.d0)
    den_coeff=1.d0/((pi**1.5d0)*(gw**3))
    hgp=1.d-3
    rho=0.d0
    pot=0.d0
    do igp=0,ngp
        rrad(igp)= igp*hgp
        if(rrad(igp)>(10.d0*gw)) then
            rho(igp)=0.d0
        else
            rho(igp)=den_coeff*Exp(-1.d0*(rrad(igp)**2/gw**2)) ! Charge density with Q=1
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
subroutine get_scf_pot_cent2_twogauss(ngp,gw,scf,rho,pot)
    use mod_qat_target, only: cal_powern_screened_poisson_gaussian, cal_pot_hartree
    implicit none
    integer, intent(in):: ngp
    real(8), intent(in):: gw, scf
    real(8), intent(out):: rho(0:ngp)
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer:: igp
    real(8):: den_coeff, hgp
    real(8):: rrad(0:ngp), weight(0:ngp)
    real(8):: pi, gw_t, alpha, beta !, tt
    real(8), allocatable:: pot_scn(:), pot_scn_t(:)
    pi=4.d0*atan(1.d0)
    allocate(pot_scn(0:ngp),pot_scn_t(0:ngp))
    den_coeff=1.d0/((pi**1.5d0)*(gw**3))
    hgp=1.d-3
    rho=0.d0
    pot=0.d0
    do igp=0,ngp
        rrad(igp)= igp*hgp
        if(rrad(igp)>(10.d0*gw)) then
            rho(igp)=0.d0
        else
            rho(igp)=den_coeff*exp(-1.d0*(rrad(igp)**2/gw**2)) ! Charge density with Q=1
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
subroutine get_eigenval(atoms,amat)
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: info
    real(8), allocatable:: eigen_val_amat(:,:)
    real(8), allocatable:: real_eigenval(:), work(:)
    allocate(eigen_val_amat(1:atoms%nat,1:atoms%nat))
    allocate(real_eigenval(1:atoms%nat))
    allocate(work(1:4*atoms%nat))
    eigen_val_amat(1:atoms%nat,1:atoms%nat)=amat(1:atoms%nat,1:atoms%nat)
    call DSYEV('N','U',atoms%nat,eigen_val_amat,atoms%nat,real_eigenval,work,4*atoms%nat,info)
    write(13,'(a6,es18.8,a10,es18.8)') 'MAX: ',maxval(real_eigenval),' | MIN: ',minval(real_eigenval)
    deallocate(eigen_val_amat)
    deallocate(real_eigenval)
    deallocate(work)
end subroutine
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
subroutine prefit_cent2(parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: itrial, iat, jat, ii
    integer:: nbgx, nbgy, nbgz, info, itypat
    real(8):: xyz(3), hgp, tt, rmse, err_U_SRS, U_SRS
    real(8):: qavg_Mg, qavg_O, qvar_Mg, qvar_O
    real(8):: cavg_Mg, cavg_O, cvar_Mg, cvar_O
    real(8):: pi
    real(8), allocatable:: EP(:,:)
    real(8), allocatable:: E_all(:)
    real(8), allocatable:: linear_rho_t(:)
    real(8), allocatable:: squarefit(:,:), squarefit_t(:,:)
    real(8), allocatable:: real_eigenval(:), work(:)
    real(8), allocatable:: EP_n(:)
    real(8):: hh_Mg, hh_O, hh, qtarget_Mg, qtarget_O, qtarget
    real(8):: one
    logical, save:: done=.false.
    if(done) return
    one=1.d0
    hgp=1.d-3
    pi=4.d0*atan(1.d0)
    associate(nat=>atoms%nat)
    if(.not. associated(atoms%trial_energy)) then
        stop 'ERROR: this part of CENT2 requires atoms%trial_energy being associated'
    endif
    allocate(EP_n(atoms%trial_energy%ntrial))
    allocate(linear_rho_t(0:poisson%ngp))
    do ii=0,poisson%ngp
        tt= ii*hgp
        if(tt>(10.d0*1.d0)) then
            linear_rho_t(ii)=0.d0
        else
            linear_rho_t(ii)=1.d0/((pi**1.5d0)*(1.d0**3))*Exp(-1.d0*(tt**2/1.d0**2))
        endif
    enddo
    allocate(E_all(atoms%trial_energy%ntrial))
    allocate(EP(atoms%nat,atoms%trial_energy%ntrial))
    nbgx=int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy=int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz=int(poisson%rgcut/poisson%hgrid(3,3))+3
    write(*,*) 'RGCUT ',poisson%rgcut,nbgx,nbgy,nbgz
    !-----------------------------------------------------------------
    do iat=1,nat
        call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,.true.,one,poisson%pot)
        do itrial=1,atoms%trial_energy%ntrial
            xyz(1:3)=atoms%ratp(1:3,atoms%trial_energy%iat_list(itrial))+atoms%trial_energy%disp(1:3,itrial)
            !call put_gto_sym_ortho(parini,poisson%bc,.true.,1,xyz,one,one, &
            !    poisson%rgcut,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,trial_rho)
            call radial_to_3d_energy(poisson%ngp,hgp,linear_rho_t,xyz, &
                poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid, &
                poisson%rgcut,one,poisson%pot,tt)
            EP(iat,itrial)=tt
        enddo !end of loop over itrial
    enddo !end of loop over iat
    do itrial=1,atoms%trial_energy%ntrial
        xyz(1:3)=atoms%ratp(1:3,atoms%trial_energy%iat_list(itrial))+atoms%trial_energy%disp(1:3,itrial)
        call radial_to_3d_energy(poisson%ngp,hgp,linear_rho_t,xyz, &
            poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid, &
            poisson%rgcut,one,poisson%pot_ion,tt)
        EP_n(itrial)=tt
    enddo !end of loop over itrial
    !-----------------------------------------------------------------
    allocate(squarefit(atoms%nat,atoms%nat))
    allocate(squarefit_t(atoms%nat+1,atoms%nat+1))
    allocate(real_eigenval(1:atoms%nat),work(atoms%nat*atoms%nat))
    squarefit=0.d0
    squarefit_t=0.d0
    do iat=1,atoms%nat
        do jat=1,atoms%nat
            tt=0.d0
            do itrial=1,atoms%trial_energy%ntrial
                tt=tt+2.d0*EP(iat,itrial)*EP(jat,itrial)
            enddo
            squarefit(iat,jat)=tt
            squarefit_t(iat,jat)=tt
        enddo
    enddo
    call DSYEV('N','U',atoms%nat,squarefit,atoms%nat,real_eigenval,work,atoms%nat**2,info)
    do iat=1,atoms%nat
        do jat=1,atoms%nat
            squarefit(iat,jat)=squarefit_t(iat,jat)
        enddo
    enddo
    !hh_Mg=40.d-2
    !hh_O=40.d-2
    hh_Mg=1.d-3*real_eigenval(atoms%nat)
    hh_O=1.d-3*real_eigenval(atoms%nat)
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') hh=hh_Mg
        if(trim(atoms%sat(iat))=='O' ) hh=hh_O
        squarefit(iat,iat)=squarefit(iat,iat)+hh
        squarefit_t(iat,iat)=squarefit_t(iat,iat)+hh
    enddo
    call DSYEV('N','U',atoms%nat,squarefit,atoms%nat,real_eigenval,work,atoms%nat**2,info)
    do iat=1,atoms%nat
        write(*,'(a,i6,es14.5)') 'EVAL ',iat,real_eigenval(iat)
    enddo
    squarefit_t(1:atoms%nat,atoms%nat+1)=1.d0
    squarefit_t(atoms%nat+1,1:atoms%nat)=1.d0
    squarefit_t(atoms%nat+1,atoms%nat+1)=0.d0
    call DGETRF(nat+1,nat+1,squarefit_t,nat+1,ann_arr%ipiv,info)
    ann_arr%qq(nat+1)=-sum(atoms%zat)
    do iat=1,nat
        tt=0.d0
        do itrial=1,atoms%trial_energy%ntrial
            tt=tt+2.d0*EP(iat,itrial)*(atoms%trial_energy%energy(itrial)-EP_n(itrial))
        enddo
        ann_arr%qq(iat)=tt
    enddo
    do itypat=1,parini%ntypat
        if(trim(parini%stypat(itypat))=='Mg') then
            qtarget_Mg=-ann_arr%ann(itypat)%zion+ann_arr%ann(1)%spring_const
            write(*,'(a,f8.3)') 'QTARGET_Mg ',qtarget_Mg
        endif
        if(trim(parini%stypat(itypat))=='O') then
            qtarget_O=-ann_arr%ann(itypat)%zion-ann_arr%ann(1)%spring_const
            write(*,'(a,f8.3)') 'QTARGET_O  ',qtarget_O
        endif
    enddo
    !qtarget_Mg=-2.d0+ann_arr%ann(1)%spring_const
    !qtarget_O =-4.d0-ann_arr%ann(1)%spring_const
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') hh=hh_Mg
        if(trim(atoms%sat(iat))=='O' ) hh=hh_O
        if(trim(atoms%sat(iat))=='Mg') qtarget=qtarget_Mg
        if(trim(atoms%sat(iat))=='O' ) qtarget=qtarget_O
        ann_arr%qq(iat)=ann_arr%qq(iat)+hh*qtarget
    enddo
    call DGETRS('N',nat+1,1,squarefit_t,nat+1,ann_arr%ipiv,ann_arr%qq,nat+1,info)
    do iat=1,atoms%nat
        write(*,'(a,i6,f7.3)') 'QQQ ',iat,atoms%zat(iat)+ann_arr%qq(iat)
    enddo
    deallocate(squarefit)
    deallocate(squarefit_t)
    deallocate(real_eigenval,work)
    !-----------------------------------------------------------------
    call reverseCEP(ann_arr,atoms,poisson,ann_arr%a)
    atoms%qat(1:atoms%nat)=ann_arr%qq(1:atoms%nat)
    rmse=0.d0
    do itrial=1,atoms%trial_energy%ntrial
        tt=0.d0
        do iat=1,nat
            tt=tt+atoms%qat(iat)*EP(iat,itrial)
        enddo
        E_all(itrial)=EP_n(itrial)+tt
        rmse=rmse+(E_all(itrial)-atoms%trial_energy%energy(itrial))**2
    enddo
    rmse=1.d3*sqrt(rmse/atoms%trial_energy%ntrial)
    call cal_etrial_cent2(parini,ann_arr,atoms,poisson,hgp,U_SRS)
    call prefit_cent2_output(ann_arr,atoms,qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O)
    err_U_SRS=1.d3*(U_SRS-atoms%epot)/nat
    write(*,'(a,2es24.15)') 'USRS ',U_SRS,atoms%epot
    write(*,'(a,2f10.3,8f7.3)') 'OPT ',rmse,err_U_SRS, &
        qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O
    do itrial=1,atoms%trial_energy%ntrial
        write(*,'(a,i3,2es24.15)') 'ETS ',atoms%trial_energy%iat_list(itrial),E_all(itrial),atoms%trial_energy%energy(itrial)
    enddo
    deallocate(EP_n)
    deallocate(E_all)
    done=.true.
    end associate
end subroutine prefit_cent2
!*****************************************************************************************
subroutine cal_etrial_cent2(parini,ann_arr,atoms,poisson,hgp,U_SRS)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: hgp
    real(8), intent(out):: U_SRS
    !local variables
    integer:: iat
    real(8):: one, tt, tte, ttn
    real(8), allocatable:: gausswidth(:)
    one=1.d0
    allocate(gausswidth(atoms%nat))
    poisson%pot(:,:,:)=0.d0
    do iat=1,atoms%nat
        call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,.false.,atoms%qat(iat),poisson%pot)
    enddo
    poisson%pot=poisson%pot+poisson%pot_ion
    gausswidth=0.5d0
    atoms%fat=0.d0
    call force_gto_sym_ortho(parini,atoms%boundcond,atoms%nat,atoms%ratp, &
        atoms%zat,gausswidth,6.d0,poisson%xyz111, &
        poisson%ngpx,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
        poisson%hgrid,poisson%pot,atoms%fat)
    do iat=1,atoms%nat
        !write(*,'(a,i4,3es19.10)') 'FAT ',iat,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
        write(61,'(a,i3,3(a2,es24.15),a)') '  - [',iat,', ',atoms%fat(1,iat),', ',atoms%fat(2,iat),', ',atoms%fat(3,iat),']'
    enddo
    tt=0.d0
    do iat=1,atoms%nat
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,atoms%qat(iat),poisson%pot,tte)
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_n(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,atoms%zat(iat),poisson%pot,ttn)
        tt=tt+tte+ttn
    enddo
    U_SRS=0.5d0*tt
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
subroutine reverseCEP(ann_arr,atoms,poisson,amat)
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use yaml_output
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(in):: poisson
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: iat, jat
    real(8):: hgp, tt, tt1, tt2, one
    real(8), allocatable:: ww(:)
    one=1.d0
    hgp=1.d-3
    allocate(ww(atoms%nat))
    do iat=1,atoms%nat
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,one,poisson%pot_ion,tt)
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
    write(*,'(a,3f10.3)') 'ECH ',tt1,tt2,tt1+tt2
    do iat=1,atoms%nat
        !write(*,'(a,i4,2f7.3)') 'CHI ',iat,ann_arr%chi_o(iat),atoms%zat(iat)+ann_arr%qq(iat)
        write(51,'(a,i3,a,es19.10,a)') '  - [',iat,', ',ann_arr%chi_o(iat),']'
    enddo
    deallocate(ww)
end subroutine reverseCEP
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
end module mod_cent2
!*****************************************************************************************
