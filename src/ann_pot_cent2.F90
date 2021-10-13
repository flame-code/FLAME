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
    integer:: iat, jat, i, j, ng
    real(8):: epot_c, out_ann
    real(8):: dpx, dpy, dpz
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8, timet1, timet2
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol
    call f_routine(id='cal_ann_cent2')
    if (.not. allocated(ann_arr%ipiv))allocate(ann_arr%ipiv(1:atoms%nat+1))
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%fat_chi(1:3,1:atoms%nat))
        allocate(ann_arr%chi_i(1:atoms%nat))
        allocate(ann_arr%chi_o(1:atoms%nat))
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
    if(parini%iverbose>=2) call cpu_time(time1)
    call init_electrostatic_cent2(parini,atoms,ann_arr,ann_arr%a,poisson)
    if(parini%iverbose>=2) call cpu_time(time2)
    if(parini%iverbose>=2) write(*,*) 'init_time: ' , time2-time1
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
        end do
        write(99,*) '==='
        write(99,*) sum(atoms%qat),sum(atoms%zat)
        write(99,*) '---'
    end if
    if(parini%iverbose>=2) call cpu_time(timet1)
    if(parini%iverbose>=2) write(*,*) 'qet_qat_time: ' , timet1-time4
    if(.not. trim(ann_arr%event)=='potential') call cent2_g_per_atom(parini,ann_arr,atoms,poisson,ann_arr%a)
    if(parini%prefit_ann) then
        call prefit_cent2(parini,ann_arr,atoms,poisson)
    endif
    if(parini%iverbose>=2) call cpu_time(time5)
    if(parini%iverbose>=2) write(*,*) 'cent2_g_per_time: ' , time5-timet1
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call get_electrostatic_cent2(parini,atoms,ann_arr,epot_c,ann_arr%a,poisson)
    if(parini%iverbose>=2) call cpu_time(timet2)
    if(parini%iverbose>=2) write(*,*) 'get_electrostatic_time: ' , timet2-time6
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
    elseif(trim(ann_arr%event)=='train') then
        atoms%epot=ann_arr%epot_trial
    elseif(trim(ann_arr%event)=='evalu') then
        atoms%epot=ann_arr%epot_es
    end if
    !if(trim(ann_arr%event)=='evalu') then
    !    tt1=0.d0
    !    tt2=0.d0
    !    tt3=0.d0
    !    do iat=1,atoms%nat
    !        fx_es=atoms%fat(1,iat)-ann_arr%fat_chi(1,iat)
    !        fy_es=atoms%fat(2,iat)-ann_arr%fat_chi(2,iat)
    !        fz_es=atoms%fat(3,iat)-ann_arr%fat_chi(3,iat)
    !        tt1=tt1+fx_es**2+fy_es**2+fz_es**2
    !        tt2=tt2+ann_arr%fat_chi(1,iat)**2+ann_arr%fat_chi(2,iat)**2+ann_arr%fat_chi(3,iat)**2
    !        tt3=tt3+fx_es*ann_arr%fat_chi(1,iat)+fy_es*ann_arr%fat_chi(2,iat)+fz_es*ann_arr%fat_chi(3,iat)
    !    enddo
    !    tt1=sqrt(tt1)
    !    tt2=sqrt(tt2)
    !    ann_arr%fchi_angle=tt3/(tt1*tt2)
    !    ann_arr%fchi_norm=tt2/tt1
    !endif
    call get_dpm(atoms,dpx,dpy,dpz,ann_arr%dpm_err)
    if(parini%iverbose>=2) then
        write(1390,'(3es18.8,a3,3es18.8,a3,es18.8)')dpx,dpy,dpz,' | ',atoms%dpm(1),atoms%dpm(2),atoms%dpm(3),' | ',ann_arr%dpm_err
    end if
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
        deallocate(ann_arr%fat_chi)
        deallocate(ann_arr%fatpq)
        deallocate(ann_arr%stresspq)
    endif
    deallocate(ann_arr%ipiv)
    call f_release_routine()
end subroutine cal_ann_cent2
!*****************************************************************************************
subroutine init_electrostatic_cent2(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, set_typat
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, jat, igx, igy, igz, itype
    integer:: nbgx, nbgy, nbgz
    integer:: linearGridNumber
    real(8):: vol, c, hgp
    real(8):: dx, dy, dz, dr, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8, timet1, timet2
    real(8):: max_cellVec
    real(8),allocatable :: gausswidth(:), rho(:), pot(:)
    associate(epot_es=>ann_arr%epot_es)
    pi=4.d0*atan(1.d0)
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
        end do
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
        if(parini%iverbose>=2) write(*,*) 'init_hartree_time: ' , time2-time1
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
        end if
        if(parini%iverbose>=2) call cpu_time(time1)
        if( .not. ann_arr%linear_rho_pot_initiated) then 
            if(parini%iverbose>=2) write(*,*) 'CENT2_NGP: ',poisson%ngp
            ann_arr%linear_rho_pot_initiated=.true. 
            allocate(ann_arr%linear_rho_e(0:poisson%ngp,parini%ntypat))
            allocate(ann_arr%linear_rho_n(0:poisson%ngp,parini%ntypat))
            allocate(ann_arr%linear_pot_e(0:poisson%ngp,parini%ntypat))
            allocate(ann_arr%linear_pot_n(0:poisson%ngp,parini%ntypat))
            do itype=1,parini%ntypat
                !subroutine get_scf_pot_cent2_onegauss(cv,ngp,rgcut,gw,scf,rho,pot)
                call get_scf_pot_cent2_twogauss(atoms%cellvec,poisson%ngp,poisson%rgcut,ann_arr%ann(itype)%gausswidth,&
                                       parini%screening_factor,poisson%linear_rho_e(0,itype),poisson%linear_pot_e(0,itype))
                call get_scf_pot_cent2_onegauss(atoms%cellvec,poisson%ngp,poisson%rgcut,ann_arr%ann(itype)%gausswidth_ion,&
                                       parini%screening_factor,poisson%linear_rho_n(0,itype),poisson%linear_pot_n(0,itype))
                ann_arr%linear_rho_e(0:poisson%ngp,itype)=poisson%linear_rho_e(0:poisson%ngp,itype)
                ann_arr%linear_rho_n(0:poisson%ngp,itype)=poisson%linear_rho_n(0:poisson%ngp,itype)
                ann_arr%linear_pot_e(0:poisson%ngp,itype)=poisson%linear_pot_e(0:poisson%ngp,itype)
                ann_arr%linear_pot_n(0:poisson%ngp,itype)=poisson%linear_pot_n(0:poisson%ngp,itype)
            end do
        else
            do itype=1,parini%ntypat
                poisson%linear_rho_e(0:poisson%ngp,itype)=ann_arr%linear_rho_e(0:poisson%ngp,itype)
                poisson%linear_rho_n(0:poisson%ngp,itype)=ann_arr%linear_rho_n(0:poisson%ngp,itype)
                poisson%linear_pot_e(0:poisson%ngp,itype)=ann_arr%linear_pot_e(0:poisson%ngp,itype)
                poisson%linear_pot_n(0:poisson%ngp,itype)=ann_arr%linear_pot_n(0:poisson%ngp,itype)
            end do
        endif
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'get_scf_pot_time: ' , time2-time1
        call update_ratp(atoms)
        if(parini%iverbose>=2) call cpu_time(time1)
        if( .not. ann_arr%amat_initiated) then 
            call get_amat_cent2(ann_arr,atoms,poisson,a)
        end if
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'get_amt_time: ' , time2-time1
        call get_eigenval(atoms,a)
        allocate(poisson%pot_ion(poisson%ngpx,poisson%ngpy,poisson%ngpz))
        nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
        nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
        nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
        !hgp = 2.d0*(poisson%rgcut+maxval(atoms%cellvec))/poisson%ngp
        !hgp = 100.d0/poisson%ngp
        hgp=1.d-3
        if(parini%iverbose>=2) call cpu_time(time1)
        poisson%pot_ion(:,:,:)=0.d0
        do iat = 1 , atoms%nat
            do igx = 1 , poisson%ngpx
                do igy = 1 , poisson%ngpy
                    do igz = 1 , poisson%ngpz
                        dx = poisson%xyz111(1)+(igx-1)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                        dy = poisson%xyz111(2)+(igy-1)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                        dz = poisson%xyz111(3)+(igz-1)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                        dr = sqrt(dx**2+dy**2+dz**2)
                        linearGridNumber=floor(dr/hgp)
                        poisson%pot_ion(igx,igy,igz)=poisson%pot_ion(igx,igy,igz)+atoms%zat(iat)*((dr/hgp-linearGridNumber)*&
                            (poisson%linear_pot_n(linearGridNumber+1,atoms%itypat(iat))&
                            -poisson%linear_pot_n(linearGridNumber,atoms%itypat(iat)))&
                            +poisson%linear_pot_n(linearGridNumber,atoms%itypat(iat)))
                    end do
                end do
            end do
        end do
        if(parini%iverbose>=2) call cpu_time(time2)
        if(parini%iverbose>=2) write(*,*) 'pot_ion_time: ' , time2-time1
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
    real(8) ,intent(out):: a(1:atoms%nat+1,1:atoms%nat+1)
    !local variables
    integer :: iat,jat,igx,igy,igz
    integer :: nbgx, nbgy, nbgz, linearGridNumber
    integer :: agpx, agpy, agpz 
    real(8) :: dx ,dy ,dz ,dr, hgp, tt
    real(8) :: rho_e,rho_e_p1,pot_e,pot_e_p1
    real(8), allocatable:: grid_rho_new(:,:,:,:), grid_pot_new(:,:,:)
    !real(8) :: time1, time2
    hgp=1.d-3
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    !call cpu_time(time1)
    allocate(grid_rho_new(-nbgx:nbgx,-nbgy:nbgy,-nbgz:nbgz,atoms%nat))
    allocate(grid_pot_new(poisson%ngpx,poisson%ngpy,poisson%ngpz))
    do iat=1,atoms%nat
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+nbgx+0
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+nbgy+0
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+nbgz+0
        do igz=agpz-nbgz,agpz+nbgz
        do igy=agpy-nbgy,agpy+nbgy
        do igx=agpx-nbgx,agpx+nbgx
            dx=poisson%xyz111(1)+(igx-1)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
            dy=poisson%xyz111(2)+(igy-1)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
            dz=poisson%xyz111(3)+(igz-1)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
            dr=sqrt(dx**2+dy**2+dz**2)
            linearGridNumber=floor(dr/hgp)
            rho_e_p1=poisson%linear_rho_e(linearGridNumber+1,atoms%itypat(iat))
            rho_e=poisson%linear_rho_e(linearGridNumber,atoms%itypat(iat))
            grid_rho_new(igx-agpx,igy-agpy,igz-agpz,iat)=(dr/hgp-linearGridNumber)*(rho_e_p1-rho_e)+rho_e
        enddo
        enddo
        enddo
    enddo
    do iat=1,atoms%nat
        do igz=1,poisson%ngpz
        do igy=1,poisson%ngpy
        do igx=1,poisson%ngpx
            dx=poisson%xyz111(1)+(igx-1)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
            dy=poisson%xyz111(2)+(igy-1)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
            dz=poisson%xyz111(3)+(igz-1)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
            dr=sqrt(dx**2+dy**2+dz**2)
            linearGridNumber=floor(dr/hgp)
            pot_e_p1=poisson%linear_pot_e(linearGridNumber+1,atoms%itypat(iat))
            pot_e=poisson%linear_pot_e(linearGridNumber,atoms%itypat(iat))
            grid_pot_new(igx,igy,igz)=(dr/hgp-linearGridNumber)*(pot_e_p1-pot_e)+pot_e
        enddo
        enddo
        enddo
        do jat=1,iat
            agpx=int(atoms%ratp(1,jat)/poisson%hgrid(1,1))+nbgx+0
            agpy=int(atoms%ratp(2,jat)/poisson%hgrid(2,2))+nbgy+0
            agpz=int(atoms%ratp(3,jat)/poisson%hgrid(3,3))+nbgz+0
            tt=0.d0
            do igz=agpz-nbgz,agpz+nbgz
            do igy=agpy-nbgy,agpy+nbgy
            do igx=agpx-nbgx,agpx+nbgx
                tt=tt+grid_rho_new(igx-agpx,igy-agpy,igz-agpz,jat)*grid_pot_new(igx,igy,igz)
                !write(1000+jat*80+iat,'(3i4,2f20.15)') igx,igy,igz,grid_rho_new(igx-agpx,igy-agpy,igz-agpz,jat),grid_pot_new(igx,igy,igz)
            enddo
            enddo
            enddo
            !close(1000+jat*80+iat)
            a(iat,jat)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
            a(jat,iat)=a(iat,jat)
        enddo
        a(iat,iat)=a(iat,iat)+ann_arr%ann(atoms%itypat(iat))%hardness
        a(iat,atoms%nat+1)=1.d0
        a(atoms%nat+1,iat)=1.d0
    enddo !end of loop over iat
    a(atoms%nat+1,atoms%nat+1)=0.d0
    deallocate(grid_rho_new)
    deallocate(grid_pot_new)
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
    integer:: info , iat
    integer:: linearGridNumber
    integer:: nbgx, nbgy, nbgz
    integer:: igx, igy, igz
    integer:: agpx, agpy, agpz
    real(8):: dx, dy, dz, dr
    real(8):: hgp, tt, one
    real(8) ::rho_val 
    real(8) :: a(atoms%nat+1,atoms%nat+1)
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
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    do iat=1, atoms%nat
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,one,poisson%pot_ion,tt)
        ann_arr%qq(iat)=-ann_arr%chi_o(iat)-atoms%zat(iat)*ann_arr%ann(atoms%itypat(iat))%hardness-tt
    end do !iat
    ann_arr%qq(nat+1)=-1.d0*sum(atoms%zat) !atoms%qtot !ASK Dr this should be -ztot or qat+zat or atoms%qtot
    call DGETRS('N',nat+1,1,a,nat+1,ann_arr%ipiv,ann_arr%qq,nat+1,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat(1:nat)=ann_arr%qq(1:nat)
    do iat=1,nat
        write(20,'(a3,4es18.6)') atoms%sat(iat),atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat),atoms%qat(iat)+atoms%zat(iat)
    end do
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
subroutine cent2_g_per_atom(parini,ann_arr,atoms,poisson,amat)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: info, iat, jat, kat, igx, igy, igz, linearGridNumber
    integer:: nbgx, nbgy, nbgz
    integer:: agpx, agpy, agpz
    real(8):: dx, dy, dz, dr, hgp, tt
    real(8):: E_par(1:atoms%nat), rhs(1:atoms%nat+1)
    real(8):: pot_val, one
    real(8):: a(atoms%nat+1,atoms%nat+1)
    real(8), allocatable:: trial_rho(:,:,:),trial_gw(:),trial_qat(:)
    real(8), allocatable:: pot_single(:,:,:)
    associate(nat=>atoms%nat)
    one=1.d0
    a=amat
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    hgp=1.d-3
    poisson%pot(:,:,:)=0.d0
    do iat=1,nat
        call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,.false.,atoms%qat(iat),poisson%pot)
    enddo
    poisson%pot=poisson%pot+poisson%pot_ion
    if(.not. ann_arr%EPar_initiated) then
        allocate(trial_rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz))
        allocate(pot_single(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz))
        allocate(trial_qat(1:atoms%nat))
        allocate(trial_gw(1:atoms%nat))
        trial_qat(:)=0.d0
        trial_qat(atoms%trial_ref_nat(1))=1.d0
        trial_gw=1.d0
        call put_gto_sym_ortho(parini,poisson%bc,.true.,atoms%nat,atoms%ratp,trial_qat,trial_gw,&
                poisson%rgcut,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,trial_rho)
        !write(*,*) 'trial_rho',maxval(trial_rho),maxloc(trial_rho)
        do iat=1, nat
            call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
                atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,.true.,one,pot_single)
            tt=0.d0
            do igz=1,poisson%ngpz
                do igy=1,poisson%ngpy
                    do igx=1,poisson%ngpx
                        tt=tt+pot_single(igx,igy,igz)*trial_rho(igx,igy,igz)
                    enddo
                enddo
            enddo
            ann_arr%EP(iat)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
        enddo
        deallocate(trial_rho)
        deallocate(pot_single)
        deallocate(trial_qat)
        deallocate(trial_gw)
    endif
    if(.not. ann_arr%chiQPar_initiated) then
        call DGETRF(nat+1,nat+1,a,nat+1,ann_arr%ipiv,info)
        if(info/=0) then
            write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
            stop
        endif
        do jat=1, nat
            rhs(1:nat+1)=0.d0
            rhs(jat)=-1.d0
            call DGETRS('N',nat+1,1,a,nat+1,ann_arr%ipiv,rhs,nat+1,info)
            if(info/=0) then
                write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
                stop
            endif
            do iat=1, nat
                ann_arr%Xq(iat,jat)=rhs(iat)
            end do
        end do
    end if
    do kat=1,ann_arr%num(1)
        do jat=1, nat
            tt=0.d0
            do iat=1, nat
                tt=tt+ann_arr%EP(iat)*ann_arr%Xq(iat,jat)*ann_arr%g_per_atom(kat,jat)
            end do
            ann_arr%g_per_atom(kat,jat)=tt
        end do
    end do
    !deallocate(ann_arr%ipiv)
    end associate
end subroutine cent2_g_per_atom
!*****************************************************************************************
subroutine get_electrostatic_cent2(parini,atoms,ann_arr,epot_c,amat,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    real(8):: a(atoms%nat+1,atoms%nat+1)
    associate(epot_es=>ann_arr%epot_es)
    a=amat
    tt1=0.d0
    tt2=0.d0
    do iat=1,atoms%nat
        tt1=tt1+ann_arr%chi_o(iat)*(atoms%qat(iat)+atoms%zat(iat))
        tt2=tt2+(atoms%qat(iat)+atoms%zat(iat))**2*0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness
    enddo
    call cal_electrostatic_ann_cent2(parini,atoms,ann_arr,a,poisson)
    if(parini%iverbose>=2)  then
        write(12,*) epot_es, atoms%epot
    end if
    epot_c=epot_es+tt1+tt2+ann_arr%ener_ref
    if(parini%iverbose>=2)  then
        write(16,'(5es18.8)') epot_c, epot_es, tt1, tt2, ann_arr%ener_ref
        write(14,*) epot_c, atoms%epot
    end if
    end associate
end subroutine get_electrostatic_cent2
!*****************************************************************************************
subroutine cal_electrostatic_ann_cent2(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    type(typ_poisson):: poisson_local
    type(typ_poisson):: poisson_force
    integer:: iat, jat
    integer:: igx, igy, igz
    integer:: nbgx, nbgy, nbgz, linearGridNumber
    integer :: agpx, agpy, agpz!, tt1(1:3)
    real(8):: tt,ehartree,rho_val, tte, ttn
    real(8):: dx, dy, dz, dr, hgp
    real(8):: ehartree_2
    real(8):: alpha,beta,ggw,ggw_t
    real(8),allocatable::fd_fat(:,:),fat_t(:,:)
    real(8),allocatable::trial_rho(:,:,:),trial_gw(:),trial_qat(:)
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    hgp=1.d-3
    if(trim(ann_arr%event)=='potential') then
        nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
        nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
        nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
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
    if(trim(ann_arr%event)/='potential' ) then
    allocate(trial_rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz))
    allocate(trial_qat(1:atoms%nat),trial_gw(1:atoms%nat))
    trial_qat(:)=0.d0
    trial_qat(atoms%trial_ref_nat(1))=1.d0
    trial_gw=1.d0
    call put_gto_sym_ortho(parini,poisson%bc,.true.,atoms%nat,atoms%ratp,trial_qat,trial_gw,&
            poisson%rgcut,poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,trial_rho)
    tt=0.d0
    !do iat=1, atoms%nat
        do igx =1,poisson%ngpx 
            do igy = 1,poisson%ngpy
                do igz = 1,poisson%ngpz
                    tt=tt+trial_rho(igx,igy,igz)*poisson%pot(igx,igy,igz)
                end do
            end do
        end do
    !end do !iat
    ann_arr%epot_trial=0.5d0*tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    endif
    if(trim(ann_arr%event)=='potential' )then
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
    do iat = 1, atoms%nat
        ggw = ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        alpha=ggw**3/(ggw**3-ggw_t**3)
        atoms%fat(1:3,iat)=atoms%fat(1:3,iat)+alpha*fat_t(1:3,iat)
    end do
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
    do iat = 1,atoms%nat
        ggw = ann_arr%ann(atoms%itypat(iat))%gausswidth
        ggw_t=0.95d0*ggw
        beta=-ggw_t**3/(ggw**3-ggw_t**3)
        atoms%fat(1:3,iat)=atoms%fat(1:3,iat)+beta*fat_t(1:3,iat)
    end do
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
subroutine get_scf_pot_cent2_onegauss(cv,ngp,rgcut,gw,scf,rho,pot)
    use mod_qat_target, only: cal_powern_screened_poisson_gaussian, cal_pot_hartree
    implicit none
    integer, intent(in) :: ngp
    real(8), intent(in) :: cv(1:3,1:3)
    real(8), intent(in) :: rgcut,gw,scf
    real(8), intent(out):: rho(0:ngp)
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer :: igp
    real(8) :: den_coeff, hgp
    real(8) :: rrad(0:ngp), pot_scn(0:ngp),weight (0:ngp)
    real(8) :: pi!, tt1
    pi = 4.d0*atan(1.d0)
    den_coeff = 1.d0/((pi**1.5d0)*(gw**3))
    hgp=1.d-3
    rho=0.d0
    pot=0.d0
    do igp = 0 , ngp
        rrad(igp)= igp*hgp
        if (rrad(igp)>(10.d0*gw)) then
            rho(igp) = 0.d0
        else
            rho(igp) = den_coeff*Exp(-1.d0*(rrad(igp)**2/gw**2)) ! Charge density with Q=1
        endif
    enddo
    if(scf>0.d0) then
        call set_weight(ngp,rrad,weight)
        call cal_powern_screened_poisson_gaussian(ngp,rrad,weight,1.d0,gw,scf,4,pot_scn)
    end if
    call cal_pot_hartree(ngp,rrad,rho,pot)
    do igp = 0 , ngp
        pot(igp)=pot(igp)-pot_scn(igp)
    end do
end subroutine get_scf_pot_cent2_onegauss
!*****************************************************************************************
subroutine get_scf_pot_cent2_twogauss(cv,ngp,rgcut,gw,scf,rho,pot)
    use mod_qat_target, only: cal_powern_screened_poisson_gaussian, cal_pot_hartree
    implicit none
    integer, intent(in) :: ngp
    real(8), intent(in) :: cv(1:3,1:3)
    real(8), intent(in) :: rgcut,gw,scf
    real(8), intent(out):: rho(0:ngp)
    real(8), intent(out):: pot(0:ngp)
    !local variables
    integer :: igp
    real(8) :: den_coeff, hgp
    real(8) :: rrad(0:ngp), weight(0:ngp)
    real(8) :: pi, gw_t, alpha, beta !, tt
    real(8), allocatable:: pot_scn(:), pot_scn_t(:)
    pi = 4.d0*atan(1.d0)
    allocate(pot_scn(0:ngp),pot_scn_t(0:ngp))
    den_coeff = 1.d0/((pi**1.5d0)*(gw**3))
    hgp=1.d-3
    rho=0.d0
    pot=0.d0
    do igp = 0 , ngp
        rrad(igp)= igp*hgp
        if (rrad(igp)>(10.d0*gw)) then
            rho(igp) = 0.d0
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
    den_coeff = 1.d0/((pi**1.5d0)*(gw_t**3))
    do igp = 0 , ngp
        !rrad(igp)= igp*hgp
        if (rrad(igp)>(10.d0*gw_t)) then
            !rho(igp) = 0.d0
        else
            rho(igp)=rho(igp)+beta*den_coeff*exp(-1.d0*(rrad(igp)**2/gw_t**2)) ! Charge density with Q=1
        endif
    enddo
    if(scf>0.d0) then
        call set_weight(ngp,rrad,weight)
        call cal_powern_screened_poisson_gaussian(ngp,rrad,weight,1.d0,gw,scf,4,pot_scn)
        call cal_powern_screened_poisson_gaussian(ngp,rrad,weight,1.d0,gw_t,scf,4,pot_scn_t)
    end if
    !tt=0.d0
    !do igp=0,ngp
    !    tt=tt+rrad(igp)**2*rho(igp)*weight(igp)
    !enddo
    !tt=tt*(4.d0*pi)
    !write(*,*) 'RHO ',tt
    call cal_pot_hartree(ngp,rrad,rho,pot)
    do igp = 0 , ngp
        pot(igp)=pot(igp)-(alpha*pot_scn(igp)+beta*pot_scn_t(igp))
    end do
end subroutine get_scf_pot_cent2_twogauss
!*****************************************************************************************
subroutine get_eigenval(atoms,amat)
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !local variables
    integer :: info
    real(8), allocatable :: eigen_val_amat(:,:),a(:,:)
    real(8), allocatable :: real_eigenval(:),work(:)
    allocate(eigen_val_amat(1:atoms%nat,1:atoms%nat))
    allocate(real_eigenval(1:atoms%nat),work(1:4*atoms%nat))
    eigen_val_amat(1:atoms%nat,1:atoms%nat)=amat(1:atoms%nat,1:atoms%nat)
    call DSYEV('N','U',atoms%nat,eigen_val_amat,atoms%nat,real_eigenval,work,4*atoms%nat,info)
    write(13,'(a6,es18.8,a10,es18.8)') 'MAX: ',maxval(real_eigenval),' | MIN: ',minval(real_eigenval)
    deallocate(eigen_val_amat)
    deallocate(real_eigenval)
end subroutine
!*****************************************************************************************
subroutine get_dpm(atoms,dpx,dpy,dpz,dpm_err)
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(in):: atoms
    real(8), intent(out):: dpm_err 
    real(8), intent(out) :: dpx, dpy, dpz
    !local variables
    integer :: iat
    real(8) :: centroid_x, centroid_y, centroid_z
        dpx = 0.d0
        dpy = 0.d0
        dpz = 0.d0
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

        end do

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
    integer:: itrial, iat, jat, kat, ix, iy, iz, iq, istep, nstep, ii
    integer:: agpx, agpy, agpz
    integer:: nbgx, nbgy, nbgz, linearGridNumber, info, itypat
    real(8):: xyz(3), dx, dy, dz, dr, hgp, tt, rho_val, q, cf, rmse, err_U_SRS
    real(8):: qavg_Mg, qavg_O, qvar_Mg, qvar_O
    real(8):: cavg_Mg, cavg_O, cvar_Mg, cvar_O
    real(8):: g_Mg, g_O
    real(8):: alpha, alphax, alphat, t1, t2, t3, y0, y1, DDOT, rlambda, pi
    external:: DDOT
    real(8), allocatable:: EP(:,:), rhs(:), trial_rho(:,:,:)
    real(8), allocatable:: E_all(:), g(:), gt(:), h(:), chi_old(:)
    real(8), allocatable:: rho_ion(:,:,:)
    real(8), allocatable:: linear_rho_t(:)
    real(8), allocatable:: amat(:,:), amat_t(:,:)
    real(8), allocatable :: real_eigenval(:), work(:)
    real(8), allocatable :: EP_n(:)
    real(8):: gausswidth(400)
    real(8):: hh_Mg, hh_O, hh, qtarget_Mg, qtarget_O, qtarget
    real(8):: one
    logical, save:: done=.false.
    if(done) return
    one=1.d0
    hgp=1.d-3
    pi=4.d0*atan(1.d0)
    associate(nat=>atoms%nat)
    allocate(EP_n(atoms%ntrial))
    allocate(linear_rho_t(0:poisson%ngp))
    do ii=0,poisson%ngp
        tt= ii*hgp
        if (tt>(10.d0*1.d0)) then
            linear_rho_t(ii) = 0.d0
        else
            linear_rho_t(ii) = 1.d0/((pi**1.5d0)*(1.d0**3))*Exp(-1.d0*(tt**2/1.d0**2))
        endif
    enddo
    allocate(E_all(atoms%ntrial+1))
    allocate(g(atoms%nat),gt(atoms%nat),h(atoms%nat),chi_old(atoms%nat))
    allocate(EP(atoms%nat,atoms%ntrial))
    allocate(rhs(atoms%nat+1))
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    write(*,*) 'RGCUT ',poisson%rgcut,nbgx,nbgy,nbgz
    !-----------------------------------------------------------------
    allocate(trial_rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz))
        do iat=1, nat
            call radial_to_3d(poisson%ngp,hgp,poisson%linear_pot_e(0,atoms%itypat(iat)), &
                atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
                poisson%hgrid,.true.,one,poisson%pot)
            do itrial=1,atoms%ntrial
                xyz(1:3)=atoms%ratp(1:3,atoms%trial_ref_nat(itrial))+atoms%trial_ref_disp(1:3,itrial)
        agpx=int(xyz(1)/poisson%hgrid(1,1))+nbgx
        agpy=int(xyz(2)/poisson%hgrid(2,2))+nbgy
        agpz=int(xyz(3)/poisson%hgrid(3,3))+nbgz
        trial_rho=0.d0
        do ix = agpx-nbgx,agpx+nbgx
            do iy = agpy-nbgy,agpy+nbgy
                do iz = agpz-nbgz,agpz+nbgz
                    dx = poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-xyz(1)
                    dy = poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-xyz(2)
                    dz = poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-xyz(3)
                    dr = sqrt(dx**2+dy**2+dz**2)
                    linearGridNumber=floor(dr/hgp)
                    trial_rho(ix,iy,iz)=1.d0*((dr/hgp-linearGridNumber)*&
                        (linear_rho_t(linearGridNumber+1)&
                        -linear_rho_t(linearGridNumber))&
                        +linear_rho_t(linearGridNumber))
                end do
            end do
        end do
                tt=0.d0
                do ix=1,poisson%ngpx
                do iy=1,poisson%ngpy
                do iz=1,poisson%ngpz
                    tt=tt+poisson%pot(ix,iy,iz)*trial_rho(ix,iy,iz)
                enddo
                enddo
                enddo
                EP(iat,itrial)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
                if(iat==1) then
                    tt=0.d0
                    do ix=1,poisson%ngpx
                    do iy=1,poisson%ngpy
                    do iz=1,poisson%ngpz
                        tt=tt+poisson%pot_ion(ix,iy,iz)*trial_rho(ix,iy,iz)
                    enddo
                    enddo
                    enddo
                    EP_n(itrial)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
                endif
            enddo
        enddo
    !-----------------------------------------------------------------
    allocate(amat(atoms%nat,atoms%nat))
    allocate(amat_t(atoms%nat+1,atoms%nat+1))
    allocate(real_eigenval(1:atoms%nat),work(atoms%nat*atoms%nat))
    amat=0.d0
    amat_t=0.d0
    do iat=1,atoms%nat
        do jat=1,atoms%nat
            tt=0.d0
            do itrial=1,atoms%ntrial
                tt=tt+2.d0*EP(iat,itrial)*EP(jat,itrial)
            enddo
            amat(iat,jat)=tt
            amat_t(iat,jat)=tt
        enddo
    enddo
    call DSYEV('N','U',atoms%nat,amat,atoms%nat,real_eigenval,work,atoms%nat**2,info)
    do iat=1,atoms%nat
        do jat=1,atoms%nat
            amat(iat,jat)=amat_t(iat,jat)
        enddo
    enddo
    !hh_Mg=40.d-2
    !hh_O=40.d-2
    hh_Mg=1.d-3*real_eigenval(atoms%nat)
    hh_O=1.d-3*real_eigenval(atoms%nat)
    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Mg') hh=hh_Mg
        if(trim(atoms%sat(iat))=='O' ) hh=hh_O
        amat(iat,iat)=amat(iat,iat)+hh
        amat_t(iat,iat)=amat_t(iat,iat)+hh
    enddo
    call DSYEV('N','U',atoms%nat,amat,atoms%nat,real_eigenval,work,atoms%nat**2,info)
    do iat=1,atoms%nat
        write(*,'(a,i6,es14.5)') 'EVAL ',iat,real_eigenval(iat)
    enddo
    amat_t(1:atoms%nat,atoms%nat+1)=1.d0
    amat_t(atoms%nat+1,1:atoms%nat)=1.d0
    amat_t(atoms%nat+1,atoms%nat+1)=0.d0
    !do jat=1,atoms%nat+1
    !    do iat=1,atoms%nat+1
    !        write(49,'(es19.10)') amat_t(iat,jat)
    !    enddo
    !enddo
    call DGETRF(nat+1,nat+1,amat_t,nat+1,ann_arr%ipiv,info)
    ann_arr%qq(nat+1)=-sum(atoms%zat)
    do iat=1,nat
        tt=0.d0
        do itrial=1,atoms%ntrial
            tt=tt+2.d0*EP(iat,itrial)*(atoms%trial_ref_energy(itrial)-EP_n(itrial))
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
    !do iat=1,atoms%nat
    !    write(49,'(es19.10,a5)') ann_arr%qq(iat),trim(atoms%sat(iat))
    !enddo
    !write(49,'(es19.10)') ann_arr%qq(nat+1)
    call DGETRS('N',nat+1,1,amat_t,nat+1,ann_arr%ipiv,ann_arr%qq,nat+1,info)
    do iat=1,atoms%nat
        write(*,'(a,i6,f7.3)') 'QQQ ',iat,atoms%zat(iat)+ann_arr%qq(iat)
    enddo
    deallocate(EP_n)
    deallocate(amat)
    deallocate(amat_t)
    deallocate(real_eigenval,work)
    !-----------------------------------------------------------------
    do iat=1,nat
        call random_number(tt)
        ann_arr%chi_o(iat)=ann_arr%ann(atoms%itypat(iat))%chi0+2.d-2*(tt-0.5d0)
    enddo
    !-----------------------------------------------------------------
    nstep=100
    do istep=0,nstep
        !do iq=1,1
        !q=0.5d0+0.1d0*iq
        !do iat=1,nat
        !    if(trim(atoms%sat(iat))=='Mg') atoms%qat(iat)=-atoms%zat(iat)+q
        !    if(trim(atoms%sat(iat))=='O' ) atoms%qat(iat)=-atoms%zat(iat)-q
        !enddo
        call prefit_cent2_gradient(parini,ann_arr,atoms,poisson,nbgx,nbgy,nbgz,linear_rho_t,hgp,.true.,EP,cf,rmse,E_all,g)
        call prefit_cent2_output(ann_arr,atoms,qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O)
        err_U_SRS=1.d3*(E_all(atoms%ntrial+1)-atoms%epot)/nat
        write(*,'(a,2es24.15)') 'USRS ',E_all(atoms%ntrial+1),atoms%epot
        !write(*,'(a,i6,2f10.5,es14.5,8f6.2)') 'OPT ',istep,rmse,err_U_SRS,sqrt(sum(g(1:nat)**2)), &
        write(*,'(a,i6,2f10.3,es14.5,8f7.3)') 'OPT ',istep,rmse,err_U_SRS,sqrt(sum(g(1:nat)**2)), &
            qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O
        if(istep==0) then
            do itrial=1,atoms%ntrial
                write(*,'(a,i3,2es24.15)') 'ETS ',atoms%trial_ref_nat(itrial),E_all(itrial),atoms%trial_ref_energy(itrial)
            enddo
            stop 'YYYYYYYYYYYYYYYYYYYYYYYYYYYYYY'
        endif
        if((rmse<10.d0) .or. istep==nstep) exit
        do iat=1,nat
            tt=abs(1.d-3*g(iat))
            ann_arr%chi_o(iat)=ann_arr%chi_o(iat)-sign(min(1.d-2,tt),g(iat))
        enddo
    enddo !end of loop over istep SD
    alphax=2.d-3
    alphat=2.d0*alphax
    alpha=0.d0
    nstep=1000
    do istep=0,nstep
        call prefit_cent2_gradient(parini,ann_arr,atoms,poisson,nbgx,nbgy,nbgz,linear_rho_t,hgp,.false.,EP,cf,rmse,E_all,g)
        call prefit_cent2_output(ann_arr,atoms,qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O)
        err_U_SRS=1.d3*(E_all(atoms%ntrial+1)-atoms%epot)/nat
        write(*,'(a,i6,2f10.3,es14.5,9f6.2)') 'OPT ',istep,rmse,err_U_SRS,sqrt(sum(g(1:nat)**2)), &
            qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O,alpha/alphax
        if((rmse<0.01d0 .and. abs(err_U_SRS)<0.1d0) .or. istep==nstep) then
            write(*,'(a,i6,2f10.3,es14.5,9f6.2)') 'FIN ',istep,rmse,err_U_SRS,sqrt(sum(g(1:nat)**2)), &
                qavg_Mg,qavg_O,qvar_Mg,qvar_O,cavg_Mg,cavg_O,cvar_Mg,cvar_O,alpha/alphax
            do itrial=1,atoms%ntrial
                write(*,'(a,i3,2es24.15)') 'ETE ',atoms%trial_ref_nat(itrial),E_all(itrial),atoms%trial_ref_energy(itrial)
            enddo
            exit
        endif
        if(istep==0) then
            h=-g
        else
            rlambda=(t1-t2)/t3
            t1=DDOT(nat,g,1,g,1)
            t2=DDOT(nat,gt,1,g,1)
            t3=DDOT(nat,gt,1,gt,1)
            rlambda=(t1-t2)/t3
            h=-g+rlambda*h

        endif
        do iat=1,nat
            chi_old(iat)=ann_arr%chi_o(iat)
            ann_arr%chi_o(iat)=ann_arr%chi_o(iat)+alphat*h(iat)
        enddo
        call prefit_cent2_gradient(parini,ann_arr,atoms,poisson,nbgx,nbgy,nbgz,linear_rho_t,hgp,.false.,EP,cf,rmse,E_all,gt)
        y0=DDOT(nat,g,1,h,1)
        y1=DDOT(nat,gt,1,h,1)
        tt=y0/(y0-y1)
        !write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
        alpha=alphat*max(min(tt,2.d0),-0.1d0)
        do iat=1,nat
            ann_arr%chi_o(iat)=chi_old(iat)+alpha*h(iat)
        enddo
        !hold=h
        gt=g
    enddo !end of loop over istep CG
    do iat=1,nat
        write(*,'(a,i4,2f7.3)') 'CHI ',iat,ann_arr%chi_o(iat),atoms%zat(iat)+atoms%qat(iat)
    enddo
    done=.true.
    end associate
end subroutine prefit_cent2
!*****************************************************************************************
subroutine prefit_cent2_gradient(parini,ann_arr,atoms,poisson,nbgx,nbgy,nbgz,linear_rho_t,hgp,applychi,EP,cf,rmse,E_all,g)
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
    integer, intent(in):: nbgx, nbgy, nbgz
    real(8), intent(in):: hgp
    logical, intent(in):: applychi
    real(8), intent(in):: linear_rho_t(0:poisson%ngp)
    real(8), intent(in):: EP(atoms%nat,atoms%ntrial)
    real(8), intent(inout):: cf, rmse, g(atoms%nat), E_all(atoms%ntrial+1)
    !local variables
    integer:: iat, jat, itrial, ix, iy, iz
    integer:: agpx, agpy, agpz
    integer:: linearGridNumber
    real(8):: tt, U_SRS, rho_val
    real(8):: xyz(3), dx, dy, dz, dr, coeff, one, tte, ttn
    real(8):: gausswidth(400)
    real(8), allocatable:: E_par(:), trial_rho(:,:,:)
    one=1.d0
    allocate(E_par(atoms%nat))
    allocate(trial_rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz))
    call reverseCEP(parini,ann_arr,atoms,poisson,ann_arr%a)
    !call get_qat_from_chi_dir_cent2(parini,ann_arr,atoms,poisson,ann_arr%a)
    !do iat=1,atoms%nat
    !    write(*,'(a,i6,f7.3)') 'qqq ',iat,atoms%zat(iat)+ann_arr%qq(iat)
    !enddo
    atoms%qat(1:atoms%nat)=ann_arr%qq(1:atoms%nat)
    poisson%pot(:,:,:)=0.d0
    do iat=1, atoms%nat
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

    do iat=1,atoms%nat
        call radial_to_3d_energy(poisson%ngp,hgp,poisson%linear_rho_e(0,atoms%itypat(iat)), &
            atoms%ratp(1,iat),poisson%xyz111,poisson%ngpx,poisson%ngpy,poisson%ngpz, &
            poisson%hgrid,poisson%rgcut,one,poisson%pot,tt)
        E_par(iat)=tt
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
    E_all(atoms%ntrial+1)=U_SRS
    !write(*,'(a,2es24.15)') 'USRS ',U_SRS,atoms%epot
    do itrial=1,atoms%ntrial
        xyz(1:3)=atoms%ratp(1:3,atoms%trial_ref_nat(itrial))+atoms%trial_ref_disp(1:3,itrial)
        !call put_gto_sym_ortho(parini,poisson%bc,.true.,1,xyz,1.d0,1.d0, &
        !    poisson%rgcut,poisson%ngpx,poisson%ngpy,poisson%ngpz,poisson%hgrid,trial_rho)
        agpx=int(xyz(1)/poisson%hgrid(1,1))+nbgx
        agpy=int(xyz(2)/poisson%hgrid(2,2))+nbgy
        agpz=int(xyz(3)/poisson%hgrid(3,3))+nbgz
        trial_rho=0.d0
        do ix = agpx-nbgx,agpx+nbgx
            do iy = agpy-nbgy,agpy+nbgy
                do iz = agpz-nbgz,agpz+nbgz
                    dx = poisson%xyz111(1)+(ix-1)*poisson%hgrid(1,1)-xyz(1)
                    dy = poisson%xyz111(2)+(iy-1)*poisson%hgrid(2,2)-xyz(2)
                    dz = poisson%xyz111(3)+(iz-1)*poisson%hgrid(3,3)-xyz(3)
                    dr = sqrt(dx**2+dy**2+dz**2)
                    linearGridNumber=floor(dr/hgp)
                    trial_rho(ix,iy,iz)=1.d0*((dr/hgp-linearGridNumber)*&
                        (linear_rho_t(linearGridNumber+1)&
                        -linear_rho_t(linearGridNumber))&
                        +linear_rho_t(linearGridNumber))
                end do
            end do
        end do
        tt=0.d0
        do ix=1,poisson%ngpx
        do iy=1,poisson%ngpy
        do iz=1,poisson%ngpz
            tt=tt+poisson%pot(ix,iy,iz)*trial_rho(ix,iy,iz)
        enddo
        enddo
        enddo
        tt=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
        !write(*,'(i3,2es24.15)') atoms%trial_ref_nat(itrial),tt,atoms%trial_ref_energy(itrial)
        E_all(itrial)=tt
    enddo
    cf=0.d0
    do itrial=1,atoms%ntrial
        cf=cf+(E_all(itrial)-atoms%trial_ref_energy(itrial))**2
    enddo
    rmse=1.d3*sqrt(cf/atoms%ntrial)
    coeff=0.5d0
    !cf=cf+coeff*(E_all(atoms%ntrial+1)-atoms%epot)**2
    if(applychi) then
    do iat=1,atoms%nat
        cf=cf+1.d-2*(ann_arr%chi_o(iat)-ann_arr%ann(atoms%itypat(iat))%chi0)**2
    enddo
    endif
    g=0.d0
    do iat=1,atoms%nat
        do itrial=1,atoms%ntrial
            tt=0.d0
            do jat=1,atoms%nat
                tt=tt+EP(jat,itrial)*ann_arr%Xq(jat,iat)
            enddo
            g(iat)=g(iat)+2.d0*tt*(E_all(itrial)-atoms%trial_ref_energy(itrial))
        enddo
    enddo
    !do iat=1,atoms%nat
    !    tt=0.d0
    !    do jat=1,atoms%nat
    !        tt=tt+E_par(jat)*ann_arr%Xq(jat,iat)
    !    enddo
    !    g(iat)=g(iat)+coeff*2.d0*tt*(E_all(atoms%ntrial+1)-atoms%epot)
    !enddo
    if(applychi) then
    do iat=1,atoms%nat
        g(iat)=g(iat)+1.d-2*2.d0*(ann_arr%chi_o(iat)-ann_arr%ann(atoms%itypat(iat))%chi0)
    enddo
    endif
end subroutine prefit_cent2_gradient
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
    integer:: itrial, iat, jat, kat, ix, iy, iz, iq, istep, nstep
    integer:: agpx, agpy, agpz
    integer:: nbgx, nbgy, nbgz, linearGridNumber
    real(8):: xyz(3), dx, dy, dz, dr, hgp, tt, rho_val, q, cf, rmse, err_U_SRS
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
    integer:: info , iat, jat
    integer:: linearGridNumber
    integer:: nbgx, nbgy, nbgz
    integer:: igx, igy, igz
    integer:: agpx, agpy, agpz
    real(8):: dx, dy, dz, dr
    real(8):: hgp, tt, tt1, tt2, one
    !real(8) :: grid_rho(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    real(8) ::rho_val 
    real(8) :: ww(400)
    associate(nat=>atoms%nat)
    one=1.d0
    hgp=1.d-3
    do iat=1, atoms%nat
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
    do iat=1,nat
        !write(*,'(a,i4,2f7.3)') 'CHI ',iat,ann_arr%chi_o(iat),atoms%zat(iat)+ann_arr%qq(iat)
        write(51,'(a,i3,a,es19.10,a)') '  - [',iat,', ',ann_arr%chi_o(iat),']'
    enddo
    end associate
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
