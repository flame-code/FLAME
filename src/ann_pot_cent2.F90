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
    call f_routine(id='cal_ann_cent1')
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
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
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
            write(99,*) iat,atoms%sat(iat),atoms%qat(iat),atoms%zat(iat)
        end do
        write(99,*) '==='
        write(99,*) sum(atoms%qat),sum(atoms%zat)
        write(99,*) '---'
    end if
    if(parini%iverbose>=2) call cpu_time(timet1)
    if(parini%iverbose>=2) write(*,*) 'qet_qat_time: ' , timet1-time4
    if(.not. trim(ann_arr%event)=='potential') call cent2_g_per_atom(parini,ann_arr,atoms,poisson,ann_arr%a)
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
    if(trim(ann_arr%event)=='potential' )then!
        atoms%epot=epot_c
    elseif(trim(ann_arr%event)=='train'.or. trim(ann_arr%event)=='evalu') then
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
            STOP('ERROR: atoms%cellvec > max_cellVec')
        endif
        poisson%ngp=ceiling(2.d0*(poisson%rgcut+max_cellVec)/hgp)
        if(.not. poisson%linear_allocated) then
            poisson%linear_allocated=.true.
            allocate(poisson%linear_rho_e(1:parini%ntypat,0:poisson%ngp))
            allocate(poisson%linear_pot_e(1:parini%ntypat,0:poisson%ngp))
            allocate(poisson%linear_rho_n(1:parini%ntypat,0:poisson%ngp))
            allocate(poisson%linear_pot_n(1:parini%ntypat,0:poisson%ngp))
            poisson%linear_rho_e(1:parini%ntypat,0:poisson%ngp)=0.d0
            poisson%linear_pot_e(1:parini%ntypat,0:poisson%ngp)=0.d0
            poisson%linear_rho_n(1:parini%ntypat,0:poisson%ngp)=0.d0
            poisson%linear_pot_n(1:parini%ntypat,0:poisson%ngp)=0.d0
        else
            poisson%linear_rho_e(1:parini%ntypat,0:poisson%ngp)=0.d0
            poisson%linear_pot_e(1:parini%ntypat,0:poisson%ngp)=0.d0
            poisson%linear_rho_n(1:parini%ntypat,0:poisson%ngp)=0.d0
            poisson%linear_pot_n(1:parini%ntypat,0:poisson%ngp)=0.d0
        end if
        if(parini%iverbose>=2) call cpu_time(time1)
        if( .not. ann_arr%linear_rho_pot_initiated) then 
            if(parini%iverbose>=2) write(*,*) 'CENT2_NGP: ',poisson%ngp
            ann_arr%linear_rho_pot_initiated=.true. 
            allocate(ann_arr%linear_rho_e(parini%ntypat,0:poisson%ngp))
            allocate(ann_arr%linear_rho_n(parini%ntypat,0:poisson%ngp))
            allocate(ann_arr%linear_pot_e(parini%ntypat,0:poisson%ngp))
            allocate(ann_arr%linear_pot_n(parini%ntypat,0:poisson%ngp))
            do itype=1,parini%ntypat
                !subroutine get_scf_pot_cent2(cv,ngp,rgcut,gw,scf,rho,pot)
                call get_scf_pot_cent2(atoms%cellvec,poisson%ngp,poisson%rgcut,ann_arr%ann(itype)%gausswidth,&
                                       parini%screening_factor,poisson%linear_rho_e(itype,0:poisson%ngp),poisson%linear_pot_e(itype,0:poisson%ngp))
                call get_scf_pot_cent2(atoms%cellvec,poisson%ngp,poisson%rgcut,ann_arr%ann(itype)%gausswidth_ion,&
                                       parini%screening_factor,poisson%linear_rho_n(itype,0:poisson%ngp),poisson%linear_pot_n(itype,0:poisson%ngp))
                ann_arr%linear_rho_e(itype,0:poisson%ngp)=poisson%linear_rho_e(itype,0:poisson%ngp)
                ann_arr%linear_rho_n(itype,0:poisson%ngp)=poisson%linear_rho_n(itype,0:poisson%ngp)
                ann_arr%linear_pot_e(itype,0:poisson%ngp)=poisson%linear_pot_e(itype,0:poisson%ngp)
                ann_arr%linear_pot_n(itype,0:poisson%ngp)=poisson%linear_pot_n(itype,0:poisson%ngp)
            end do
        else
            do itype=1,parini%ntypat
                poisson%linear_rho_e(itype,0:poisson%ngp)=ann_arr%linear_rho_e(itype,0:poisson%ngp)
                poisson%linear_rho_n(itype,0:poisson%ngp)=ann_arr%linear_rho_n(itype,0:poisson%ngp)
                poisson%linear_pot_e(itype,0:poisson%ngp)=ann_arr%linear_pot_e(itype,0:poisson%ngp)
                poisson%linear_pot_n(itype,0:poisson%ngp)=ann_arr%linear_pot_n(itype,0:poisson%ngp)
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
                        dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                        dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                        dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                        dr = sqrt(dx**2+dy**2+dz**2)
                        linearGridNumber=floor(dr/hgp)
                        poisson%pot_ion(igx,igy,igz)=poisson%pot_ion(igx,igy,igz)+atoms%zat(iat)*((dr/hgp-linearGridNumber)*&
                            (poisson%linear_pot_n(atoms%itypat(iat),linearGridNumber+1)&
                            -poisson%linear_pot_n(atoms%itypat(iat),linearGridNumber))&
                            +poisson%linear_pot_n(atoms%itypat(iat),linearGridNumber))
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
    !LOCAL Variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    integer :: iat,jat,igx,igy,igz
    integer :: nbgx, nbgy, nbgz, linearGridNumber
    integer :: agpx, agpy, agpz 
    real(8) :: dx ,dy ,dz ,dr, hgp, tt
    real(8) :: grid_pot(atoms%nat,poisson%ngpx,poisson%ngpy,poisson%ngpz)
    real(8) :: grid_rho(atoms%nat,poisson%ngpx,poisson%ngpy,poisson%ngpz)
    real(8) :: rho_e,rho_e_p1,pot_e,pot_e_p1
    hgp=1.d-3
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    do igz = 1 , poisson%ngpz
        do igy = 1 , poisson%ngpy
            do igx = 1 , poisson%ngpx
                do iat = 1 , atoms%nat
                    dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                    dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                    dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                    dr = sqrt(dx**2+dy**2+dz**2)
                    linearGridNumber=floor(dr/hgp)
                    rho_e_p1=poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber+1)
                    rho_e=poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber)
                    pot_e_p1=poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber+1)
                    pot_e=poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber)
                    grid_rho(iat,igx,igy,igz)=(dr/hgp-linearGridNumber)*(rho_e_p1-rho_e)+rho_e
                    grid_pot(iat,igx,igy,igz)=(dr/hgp-linearGridNumber)*(pot_e_p1-pot_e)+pot_e
                end do
            end do
        end do ! MOVE THIS LOOP TO INSIDE JAT
    end do
    do iat=1, atoms%nat
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+nbgx
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+nbgy
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+nbgz
        do jat=iat, atoms%nat
            tt=0.d0
            do igx = agpx-nbgx,agpx+nbgx
                do igy = agpy-nbgy,agpy+nbgy
                    do igz = agpz-nbgz,agpz+nbgz
                        tt=tt+grid_rho(iat,igx,igy,igz)*grid_pot(jat,igx,igy,igz)
                    end do
                end do
            end do
            a(iat,jat)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
            a(jat,iat)=a(iat,jat)
        end do !jat
        a(iat,iat)=a(iat,iat)+ann_arr%ann(atoms%itypat(iat))%hardness
        a(iat,atoms%nat+1)=1.d0
        a(atoms%nat+1,iat)=1.d0
    end do !iat
    a(atoms%nat+1,atoms%nat+1)=0.d0
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
    real(8):: hgp, tt
    !real(8) :: grid_rho(poisson%ngpx,poisson%ngpy,poisson%ngpz)
    real(8) ::rho_val 
    real(8) :: a(atoms%nat+1,atoms%nat+1)
    associate(nat=>atoms%nat)
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
        tt=0.d0
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+nbgx
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+nbgy
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+nbgz
        do igx = agpx-nbgx,agpx+nbgx
            do igy = agpy-nbgy,agpy+nbgy
                do igz = agpz-nbgz,agpz+nbgz
                    dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                    dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                    dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                    dr = sqrt(dx**2+dy**2+dz**2)
                    linearGridNumber=floor(dr/hgp)
                    rho_val=(dr/hgp-linearGridNumber)*&
                        (poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber+1)&
                        -poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber))&
                        +poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber)
                    tt=tt+rho_val*poisson%pot_ion(igx,igy,igz)
                end do
            end do
        end do
        ann_arr%qq(iat)=-ann_arr%chi_o(iat)-atoms%zat(iat)*ann_arr%ann(atoms%itypat(iat))%hardness-tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
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
    real(8):: q_chi_par(1:atoms%nat,1:atoms%nat)
    real(8):: rho_val 
    real(8):: a(atoms%nat+1,atoms%nat+1)
    associate(nat=>atoms%nat)
    a=amat
    nbgx = int(poisson%rgcut/poisson%hgrid(1,1))+3
    nbgy = int(poisson%rgcut/poisson%hgrid(2,2))+3
    nbgz = int(poisson%rgcut/poisson%hgrid(3,3))+3
    hgp=1.d-3
    poisson%pot(:,:,:)=0.d0
    do iat=1, nat
        do igx = 1 , poisson%ngpx
            do igy = 1 , poisson%ngpy
                do igz = 1 , poisson%ngpz
                    dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                    dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                    dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                    dr = sqrt(dx**2+dy**2+dz**2)
                    linearGridNumber=floor(dr/hgp)
                    poisson%pot(igx,igy,igz)=poisson%pot(igx,igy,igz)+atoms%qat(iat)*((dr/hgp-linearGridNumber)*&
                        (poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber+1)&
                        -poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber))&
                        +poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber))
                end do
            end do
        end do
    end do !iat
    poisson%pot=poisson%pot+poisson%pot_ion
    do iat=1, nat
        tt=0.d0
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+nbgx
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+nbgy
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+nbgz
        do igx = agpx-nbgx,agpx+nbgx
            do igy = agpy-nbgy,agpy+nbgy
                do igz = agpz-nbgz,agpz+nbgz
                    dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                    dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                    dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                    dr = sqrt(dx**2+dy**2+dz**2)
                    linearGridNumber=floor(dr/hgp)
                    rho_val=(dr/hgp-linearGridNumber)*&
                        (poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber+1)&
                        -poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber))&
                        +poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber)
                    tt=tt+rho_val*poisson%pot(igx,igy,igz)
                end do
            end do
        end do
        E_par(iat)=tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    end do !iat
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
            q_chi_par(iat,jat)=rhs(iat)
        end do
    end do
    do kat=1,ann_arr%num(1)
        do jat=1, nat
            tt=0.d0
            do iat=1, nat
                tt=tt+E_par(iat)*q_chi_par(iat,jat)*ann_arr%g_per_atom(kat,jat)
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
    real(8):: tt,ehartree,rho_val
    real(8):: dx, dy, dz, dr, hgp
!    real(8):: grid_rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)
    real(8):: ehartree_2
    real(8),allocatable::fd_fat(:,:)
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
            do igx = 1 , poisson%ngpx
                do igy = 1 , poisson%ngpy
                    do igz = 1 , poisson%ngpz
                        dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                        dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                        dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                        dr = sqrt(dx**2+dy**2+dz**2)
                        linearGridNumber=floor(dr/hgp)
                        poisson%pot(igx,igy,igz)=poisson%pot(igx,igy,igz)+atoms%qat(iat)*((dr/hgp-linearGridNumber)*&
                            (poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber+1)&
                            -poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber))&
                            +poisson%linear_pot_e(atoms%itypat(iat),linearGridNumber))&
                            +atoms%zat(iat)*((dr/hgp-linearGridNumber)*&
                            (poisson%linear_pot_n(atoms%itypat(iat),linearGridNumber+1)&
                             -poisson%linear_pot_n(atoms%itypat(iat),linearGridNumber))&
                            +poisson%linear_pot_n(atoms%itypat(iat),linearGridNumber))
                    end do
                end do
            end do
        end do !iat
    end if
    !grid_rho(1:poisson%ngpx,1:poisson%ngpy,1:poisson%ngpz)=0.d0
    tt=0.d0
    do iat=1, atoms%nat
        agpx=int(atoms%ratp(1,iat)/poisson%hgrid(1,1))+nbgx
        agpy=int(atoms%ratp(2,iat)/poisson%hgrid(2,2))+nbgy
        agpz=int(atoms%ratp(3,iat)/poisson%hgrid(3,3))+nbgz
        do igx = agpx-nbgx,agpx+nbgx
            do igy = agpy-nbgy,agpy+nbgy
                do igz = agpz-nbgz,agpz+nbgz
                    dx = (-nbgx+igx)*poisson%hgrid(1,1)-atoms%ratp(1,iat)
                    dy = (-nbgy+igy)*poisson%hgrid(2,2)-atoms%ratp(2,iat)
                    dz = (-nbgz+igz)*poisson%hgrid(3,3)-atoms%ratp(3,iat)
                    dr = sqrt(dx**2+dy**2+dz**2)
                    linearGridNumber=floor(dr/hgp)
                    !grid_rho(igx,igy,igz)=grid_rho(igx,igy,igz)+atoms%qat(iat)*((dr/hgp-linearGridNumber)*&
                    rho_val=atoms%qat(iat)*((dr/hgp-linearGridNumber)*&
                        (poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber+1)&
                        -poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber))&
                        +poisson%linear_rho_e(atoms%itypat(iat),linearGridNumber))&
                        +atoms%zat(iat)*((dr/hgp-linearGridNumber)*&
                        (poisson%linear_rho_n(atoms%itypat(iat),linearGridNumber+1)&
                        -poisson%linear_rho_n(atoms%itypat(iat),linearGridNumber))&
                        +poisson%linear_rho_n(atoms%itypat(iat),linearGridNumber))
                    tt=tt+rho_val*poisson%pot(igx,igy,igz)
                end do
            end do
        end do
    end do !iat
    !tt=0.d0
    !do igx = 1 , poisson%ngpx
    !    do igy = 1 , poisson%ngpy
    !        do igz = 1 , poisson%ngpz
    !            tt=tt+grid_rho(igx,igy,igz)*poisson%pot(igx,igy,igz)
    !        end do
    !    end do
    !end do
    ann_arr%epot_es=0.5d0*tt*poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    !Force
    poisson%rcart=atoms%ratp
    poisson_force=poisson
    poisson_force%q(:)=atoms%qat(:)
    poisson_force%gw(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
    poisson_force%gw_ewald(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
    !atoms%fat=0.d0
    call force_gto_sym_ortho(parini,poisson_force%bc,poisson_force%nat,poisson_force%rcart, &
        poisson_force%q,poisson_force%gw_ewald,poisson_force%rgcut,poisson_force%lda,poisson_force%ngpx, &
        poisson_force%ngpy,poisson_force%ngpz,poisson_force%hgrid,poisson_force%pot,atoms%fat)
    poisson_force%q(:)=atoms%zat(:)
    poisson_force%gw(:)=ann_arr%ann(atoms%itypat(:))%gausswidth_ion
    poisson_force%gw_ewald(:)=ann_arr%ann(atoms%itypat(:))%gausswidth_ion
    call force_gto_sym_ortho(parini,poisson_force%bc,poisson_force%nat,poisson_force%rcart, &
        poisson_force%q,poisson_force%gw_ewald,poisson_force%rgcut,poisson_force%lda,poisson_force%ngpx, &
        poisson_force%ngpy,poisson_force%ngpz,poisson_force%hgrid,poisson_force%pot,atoms%fat)
    !write(*,*) maxval(poisson_force%pot),minval(poisson_force%pot)
    !write(*,*) poisson_force%q
    !write(*,*) atoms%fat
    !stop
    !poisson_local=poisson
    !poisson_local%q=0.d0
    !poisson_local%q(:)=atoms%qat(:)
    !poisson_local%reset_rho=.true.
    !poisson_local%gw(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
    !poisson_local%gw_ewald(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
    !poisson_local%rcart=atoms%ratp
    !call put_gto_sym_ortho(parini,poisson_local%bc,poisson_local%reset_rho,poisson_local%nat, &
    !    poisson_local%rcart,poisson_local%q,poisson_local%gw,poisson_local%rgcut, &
    !    poisson_local%ngpx,poisson_local%ngpy,poisson_local%ngpz,poisson_local%hgrid,poisson_local%rho)
    !poisson_local%q=0.d0
    !poisson_local%q(:)=atoms%zat(:)
    !poisson_local%reset_rho=.false.
    !poisson_local%gw(:)=ann_arr%ann(atoms%itypat(:))%gausswidth_ion
    !call put_gto_sym_ortho(parini,poisson_local%bc,poisson_local%reset_rho,poisson_local%nat, &
    !    poisson_local%rcart,poisson_local%q,poisson_local%gw,poisson_local%rgcut, &
    !    poisson_local%ngpx,poisson_local%ngpy,poisson_local%ngpz,poisson_local%hgrid,poisson_local%rho)
    !call get_hartree(parini,poisson_local,atoms,poisson_local%gw,ehartree)

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
subroutine get_scf_pot_cent2(cv,ngp,rgcut,gw,scf,rho,pot)
    integer, intent(in) :: ngp
    real(8), intent(in) :: cv(1:3,1:3)
    real(8), intent(in) :: rgcut,gw,scf
    real(8), intent(out):: rho(0:ngp)
    real(8), intent(out):: pot(0:ngp)
    !Local Variables!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
end subroutine
!*****************************************************************************************
subroutine get_eigenval(atoms,amat)
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: amat(atoms%nat+1,atoms%nat+1)
    !Local Variables
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
    !Local Variables
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
