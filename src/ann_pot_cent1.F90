!*****************************************************************************************
subroutine cal_ann_cent1(parini,atoms,symfunc,ann_arr)
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
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol
    real(8):: dpm_err,dpx,dpy,dpz
    call f_routine(id='cal_ann_cent1')
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
        ann_arr%a=0.d0
    else
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%a=0.d0
    endif
    if(parini%iverbose>=2) call cpu_time(time1)
    call init_electrostatic_cent1(parini,atoms,ann_arr,ann_arr%a,poisson)
    if(parini%iverbose>=2) call cpu_time(time2)
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
    call get_qat_from_chi_cent1(parini,ann_arr,atoms,poisson,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time5)
    atoms%stress(1:3,1:3)=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call get_electrostatic_cent1(parini,atoms,ann_arr,epot_c,ann_arr%a,poisson)
    if(parini%iverbose>=2) then
        call cpu_time(time7)
        call yaml_mapping_open('Timing of CENT1')
        call yaml_map('initialize matrix',time2-time1)
        call yaml_map('calculation of symfunc',time3-time2)
        call yaml_map('neural network process',time4-time3)
        call yaml_map('linear equations solver',time5-time4)
        call yaml_map('force (SR term)',time6-time5)
        call yaml_map('energy (SR+LR), force (LR)',time7-time6)
        call yaml_map('total time',time7-time1)
        call yaml_mapping_close()
        !write(*,'(a,f8.3)') 'Timing:cent1: initialize matrix          ',time2-time1
        !write(*,'(a,f8.3)') 'Timing:cent1: calculation of symfunc     ',time3-time2
        !write(*,'(a,f8.3)') 'Timing:cent1: neural network process     ',time4-time3
        !write(*,'(a,f8.3)') 'Timing:cent1: linear equations solver    ',time5-time4
        !write(*,'(a,f8.3)') 'Timing:cent1: force (SR term)            ',time6-time5
        !write(*,'(a,f8.3)') 'Timing:cent1: energy (SR+LR), force (LR) ',time7-time6
        !write(*,'(a,f8.3)') 'Timing:cent1: total time                 ',time7-time1
    endif !end of if for printing out timing.
    atoms%epot=epot_c
    if(trim(ann_arr%event)=='evalu') then
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        do iat=1,atoms%nat
            fx_es=atoms%fat(1,iat)-ann_arr%fat_chi(1,iat)
            fy_es=atoms%fat(2,iat)-ann_arr%fat_chi(2,iat)
            fz_es=atoms%fat(3,iat)-ann_arr%fat_chi(3,iat)
            tt1=tt1+fx_es**2+fy_es**2+fz_es**2
            tt2=tt2+ann_arr%fat_chi(1,iat)**2+ann_arr%fat_chi(2,iat)**2+ann_arr%fat_chi(3,iat)**2
            tt3=tt3+fx_es*ann_arr%fat_chi(1,iat)+fy_es*ann_arr%fat_chi(2,iat)+fz_es*ann_arr%fat_chi(3,iat)
        enddo
        tt1=sqrt(tt1)
        tt2=sqrt(tt2)
        ann_arr%fchi_angle=tt3/(tt1*tt2)
        ann_arr%fchi_norm=tt2/tt1
    endif
    call get_dpm_cent1(atoms,dpx,dpy,dpz,ann_arr%dpm_err)
    if(parini%iverbose>=2) then
        write(1390,'(3es18.8,a3,3es18.8,a3,es18.8)')dpx,dpy,dpz,' | ',atoms%dpm(1),atoms%dpm(2),atoms%dpm(3),' | ',ann_arr%dpm_err
    end if
    atoms%dpm(1)=dpx
    atoms%dpm(2)=dpy
    atoms%dpm(3)=dpz
    call fini_electrostatic_cent1(parini,ann_arr,atoms,poisson)
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
    call f_release_routine()
end subroutine cal_ann_cent1
!*****************************************************************************************
subroutine get_qat_from_chi_cent1(parini,ann_arr,atoms,poisson,a)
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
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: iat
    real(8):: qtot
    character(200):: smsg
    if(trim(ann_arr%syslinsolver)=='direct') then
        if(trim(atoms%boundcond)=='free') then
            call get_qat_from_chi_dir(parini,ann_arr,atoms,a)
            !qtot=0.d0
            !do iat=1,atoms%nat
            !    qtot=qtot+atoms%qat(iat)
            !    write(*,'(i6,a5,2f8.3)') iat,trim(atoms%stypat(atoms%itypat(iat))),atoms%qat(iat),ann_arr%chi_o(iat)
            !enddo
            !write(*,*) 'QTOT= ',qtot
        else
            smsg='ERROR: solving linear system of equations with direct approach is'
            smsg=trim(smsg)//' possible only for free BC, atoms%boundcond= '
            write(*,'(2a)') trim(smsg),trim(atoms%boundcond)
            stop
        endif
    elseif(trim(ann_arr%syslinsolver)=='apply_matrix') then
        if(trim(atoms%boundcond)=='free') then
            call get_qat_from_chi_iter(parini,ann_arr,atoms,a)
        else
            smsg='ERROR: solving linear system of equations with explicitly'
            smsg=trim(smsg)//' applying matrix is possible only for free BC, atoms%boundcond= '
            write(*,'(2a)') trim(smsg),trim(atoms%boundcond)
            stop
        endif
    elseif(trim(ann_arr%syslinsolver)=='operator') then
        call get_qat_from_chi_operator(parini,poisson,ann_arr,atoms)
    else
        stop 'ERROR: unknown syslinsolver'
    endif
    if(parini%iverbose>=2) then
        call yaml_map('charge on atoms',atoms%qat(1:atoms%nat),fmt='(f10.5)')
        !call yaml_sequence_open('charge on atom')
        !do iat=1,atoms%nat
        !    call yaml_sequence(advance='no')
        !    call yaml_scalar(iat)
        !    call yaml_scalar(atoms%qat(iat))
        !    !write(*,*) 'charge on atom ',iat,atoms%qat(iat) 
        !enddo
        !call yaml_sequence_close()
        !do iat=1,atoms%nat
        !    write(*,*) 'charge on atom ',iat,atoms%qat(iat) 
        !enddo
    endif
end subroutine get_qat_from_chi_cent1
!*****************************************************************************************
subroutine get_qat_from_chi_dir(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: info , iat
    associate(nat=>atoms%nat)
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%ipiv(1:nat+1))
        allocate(ann_arr%qq(1:nat+1))
    endif
    call DGETRF(nat+1,nat+1,a,nat+1,ann_arr%ipiv,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
        stop
    endif
    !do iat=1,nat
    !if(trim(atoms%sat(iat))=='Na') then
    !    chi(iat)=chi(iat)/0.93d0
    !elseif(trim(atoms%sat(iat))=='Cl') then
    !    chi(iat)=chi(iat)/3.16d0
    !endif
    !enddo
    ann_arr%qq(1:nat)=-ann_arr%chi_o(1:nat)
    ann_arr%qq(nat+1)=atoms%qtot
    call DGETRS('N',nat+1,1,a,nat+1,ann_arr%ipiv,ann_arr%qq,nat+1,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat(1:nat)=ann_arr%qq(1:nat)
    do iat=1,nat
        write(20,'(a3,4es18.6)') atoms%sat(iat),atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat),atoms%qat(iat)
    end do
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        call yaml_map('Lagrange',ann_arr%qq(nat+1))
        !write(*,*) 'Lagrange ',ann_arr%qq(nat+1)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%ipiv)
        deallocate(ann_arr%qq)
    endif
    end associate
end subroutine get_qat_from_chi_dir
!*****************************************************************************************
subroutine init_electrostatic_cent1(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    real(8),allocatable :: gausswidth(:)
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    associate(epot_es=>ann_arr%epot_es)
    pi=4.d0*atan(1.d0)
    if(trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train') then
        if(trim(atoms%boundcond)=='free' .and. parini%free_bc_direct) then
            ann_arr%syslinsolver='direct'
        else
            ann_arr%syslinsolver=trim(parini%syslinsolver_ann)
        endif
    else
        ann_arr%syslinsolver=trim(parini%syslinsolver_ann)
    endif
    ann_arr%ener_ref=0.d0
    do iat=1,atoms%nat
        ann_arr%ener_ref=ann_arr%ener_ref+ann_arr%ann(atoms%itypat(iat))%ener_ref
    enddo
    if (.not. parini%ewald) then 
        poisson%alpha = maxval(ann_arr%ann(:)%gausswidth)
    else 
        if (parini%alpha_ewald<0.d0) then
            call getvol_alborz(atoms%cellvec,vol)
            c=2.2d0
            poisson%alpha = 1.d0/(c*sqrt(pi)*(atoms%nat/vol**2)**(1.d0/6.d0))
            write(*,*)"optimized alpha = ", poisson%alpha
        else
            poisson%alpha=parini%alpha_ewald
        endif
    end if
    if(trim(ann_arr%syslinsolver)=='direct' .or. trim(ann_arr%syslinsolver)=='apply_matrix') then
        if(trim(atoms%boundcond)/='free') then
            write(*,*) 'ERROR: syslinsolver=direct can be used only for free BC.'
        endif
        call get_amat_cent1(atoms,ann_arr,a)
    elseif(trim(ann_arr%syslinsolver)=='operator') then
        if(trim(atoms%boundcond)=='bulk' .or. trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='free') then
            allocate(gausswidth(atoms%nat))
            gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
            poisson%task_finit="alloc_rho:set_ngp"
            call init_hartree(parini,atoms,poisson,gausswidth)
            deallocate(gausswidth)
        else
            write(*,'(a)',advance='no') 'ERROR: currently syslinsolver=operator only '
            write(*,'(2a)') 'for BC=bulk/slab, but now BC=',trim(trim(atoms%boundcond))
            stop
        endif
    else
        write(*,*) 'ERROR: unknown value for syslinsolver',trim(ann_arr%syslinsolver)
        stop
    endif
    end associate
end subroutine init_electrostatic_cent1
!*****************************************************************************************
subroutine get_amat_cent1(atoms,ann_arr,a)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: iat, jat
    real(8):: dx, dy, dz, r, pi, beta_iat, beta_jat, gama
    pi=4.d0*atan(1.d0)
    call update_ratp(atoms)
    do iat=1,atoms%nat
        a(iat,atoms%nat+1)=1.d0
        a(atoms%nat+1,iat)=1.d0
        beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
        gama=1.d0/sqrt(beta_iat**2+beta_iat**2)
        a(iat,iat)=gama*2.d0/sqrt(pi)+ann_arr%ann(atoms%itypat(iat))%hardness
        do jat=iat+1,atoms%nat
            dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
            dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
            dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth
            gama=1.d0/sqrt(beta_iat**2+beta_jat**2)
            a(iat,jat)=erf(gama*r)/r
            a(jat,iat)=a(iat,jat)
        enddo
    enddo
    a(atoms%nat+1,atoms%nat+1)=0.d0
end subroutine get_amat_cent1
!*****************************************************************************************
subroutine fini_electrostatic_cent1(parini,ann_arr,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    !use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    if(trim(ann_arr%syslinsolver)=='operator') then
        call fini_hartree(parini,atoms,poisson)
    endif
end subroutine fini_electrostatic_cent1
!*****************************************************************************************
subroutine get_electrostatic_cent1(parini,atoms,ann_arr,epot_c,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    associate(epot_es=>ann_arr%epot_es)
    tt1=0.d0
    tt2=0.d0
    do iat=1,atoms%nat
        tt1=tt1+ann_arr%chi_o(iat)*atoms%qat(iat)
        tt2=tt2+atoms%qat(iat)**2*0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness
    enddo
    call cal_electrostatic_ann(parini,atoms,ann_arr,a,poisson)
    epot_c=epot_es+tt1+tt2+ann_arr%ener_ref
    end associate
end subroutine get_electrostatic_cent1
!*****************************************************************************************
subroutine cal_electrostatic_ann(parini,atoms,ann_arr,a,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_poisson
    use hartree_mod, only: method_samare
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: a(atoms%nat+1,atoms%nat+1)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, jat
    real(8):: tt2, tt3, ttf, gama, pi
    real(8):: dx, dy, dz, r, beta_iat, beta_jat, ehartree_t
    real(8), allocatable:: gausswidth(:)
    allocate(gausswidth(1:atoms%nat))
    if(trim(atoms%boundcond)=='free' .and. trim(ann_arr%syslinsolver)/='operator') then
        pi=4.d0*atan(1.d0)
        tt2=0.d0
        tt3=0.d0
        call update_ratp(atoms)
        do iat=1,atoms%nat
            beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
            gama=1.d0/sqrt(beta_iat**2+beta_iat**2)
            tt2=tt2+atoms%qat(iat)**2*gama/sqrt(pi)
            do jat=iat+1,atoms%nat
                dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
                dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
                dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                beta_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth
                gama=1.d0/sqrt(beta_iat**2+beta_jat**2)
                tt3=tt3+atoms%qat(iat)*atoms%qat(jat)*erf(gama*r)/r
                ttf=(2.d0/sqrt(pi)*gama*exp(-(gama*r)**2)/r**2-erf(gama*r)/r**3)*atoms%qat(iat)*atoms%qat(jat)
                atoms%fat(1,iat)=atoms%fat(1,iat)+ttf*dx
                atoms%fat(2,iat)=atoms%fat(2,iat)+ttf*dy
                atoms%fat(3,iat)=atoms%fat(3,iat)+ttf*dz
                atoms%fat(1,jat)=atoms%fat(1,jat)-ttf*dx
                atoms%fat(2,jat)=atoms%fat(2,jat)-ttf*dy
                atoms%fat(3,jat)=atoms%fat(3,jat)-ttf*dz
            enddo
        enddo
        ann_arr%epot_es=tt2+tt3
    elseif(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk' .or. trim(atoms%boundcond)=='free') then
        if (parini%bigdft_kwald) then 
            method_samare = "kwald"
        endif
        gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
        call get_hartree(parini,poisson,atoms,gausswidth,ehartree_t)
        poisson%gw(1:poisson%nat)=poisson%gw_ewald(1:poisson%nat)
        call get_hartree_force(parini,poisson,atoms)
        if (parini%bigdft_kwald) then 
            method_samare = "defau"
        endif
    else
        stop 'ERROR: the requested BCs is not yet implemented.'
    endif
    deallocate(gausswidth)
end subroutine cal_electrostatic_ann
!*****************************************************************************************
subroutine charge_analysis(parini,atoms,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: iat, i !, ii1, ii2
    real(8):: q, c !, tt1, tt2, ss1, ss2, tt1min, tt1max, tt2min, tt2max, ss1min, ss1max, ss2min, ss2max
    real(8):: chi_min_per_conf(10), chi_max_per_conf(10)
    !real(8):: dipole(3)
    !dipole(1)=0.d0 ; dipole(2)=0.d0 ; dipole(3)=0.d0
    !do iat=1,atoms%nat
    !    dipole(1)=dipole(1)+atoms%qat(iat)*atoms%rat(1,iat)
    !    dipole(2)=dipole(2)+atoms%qat(iat)*atoms%rat(2,iat)
    !    dipole(3)=dipole(3)+atoms%qat(iat)*atoms%rat(3,iat)
    !enddo
    !write(91,'(a,3es14.5)') 'dipole moment ',dipole(1),dipole(2),dipole(3)
    !----------------------------------------------------
    chi_min_per_conf(1:10)= 1.d20
    chi_max_per_conf(1:10)=-1.d20
    do iat=1,atoms%nat
        q=atoms%zat(iat)+atoms%qat(iat)
        c=ann_arr%chi_o(iat)
        i=atoms%itypat(iat)
        !write(81,*) i,atoms%stypat(i),parini%stypat(i)
        ann_arr%natsum(i)=ann_arr%natsum(i)+1
        ann_arr%qmin(i)=min(q,ann_arr%qmin(i))
        ann_arr%qmax(i)=max(q,ann_arr%qmax(i))
        ann_arr%qsum(i)=ann_arr%qsum(i)+q
        ann_arr%chi_min(i)=min(c,ann_arr%chi_min(i))
        ann_arr%chi_max(i)=max(c,ann_arr%chi_max(i))
        chi_min_per_conf(i)=min(c,chi_min_per_conf(i))
        chi_max_per_conf(i)=max(c,chi_max_per_conf(i))
        ann_arr%chi_sum(i)=ann_arr%chi_sum(i)+c
    enddo
    do i=1,ann_arr%nann
        ann_arr%chi_delta(i)=max(ann_arr%chi_delta(i),chi_max_per_conf(i)-chi_min_per_conf(i))
    enddo
    !if(parini%iverbose>=1) then
    !write(61,'(a,2(a,4f8.3),es11.2,i5)') trim(str),' Na=',tt1/ii1,tt1min,tt1max,tt1max-tt1min,' Cl=',tt2/ii2,tt2min,tt2max,tt2max-tt2min,tt1+tt2,atoms%nat
    !write(71,'(a,2(a,4f8.3),es11.2,i5)') trim(str),' Na=',ss1/ii1,ss1min,ss1max,ss1max-ss1min,' Cl=',ss2/ii2,ss2min,ss2max,ss2max-ss2min,ss1+ss2,atoms%nat
    !endif
end subroutine charge_analysis
!*****************************************************************************************
subroutine get_qat_from_chi_iter(parini,ann_arr,atoms,a)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: info, iter, i, j, niter
    integer, allocatable:: ipiv(:)
    real(8), allocatable:: qq(:), h(:), Ah(:), g(:)
    real(8):: hAh, ddot, alpha, beta, dot_gold, dot_gnew 
    real(8):: resnormtol, resnorm
    allocate(ipiv(1:atoms%nat+1))
    allocate(qq(1:atoms%nat+1))
    allocate(h(1:atoms%nat+1))
    allocate(Ah(1:atoms%nat+1))
    allocate(g(1:atoms%nat+1))

    qq(1:atoms%nat+1)=0.d0  !??????????????????
    resnormtol=1.d-10
    niter=1000        !?????????????????? 
    iter=0
    g(1:atoms%nat)=-ann_arr%chi_o(1:atoms%nat)
    g(atoms%nat+1)=atoms%qtot
    call DSYMV ('U',atoms%nat+1,1.d0,a,atoms%nat+1,qq,1,-1.d0,g,1)
    h(1:atoms%nat+1)=-g(1:atoms%nat+1)
    dot_gold = DDOT(atoms%nat+1,g,1,g,1)
    do j=1,niter
        iter = iter+1
        call DSYMV ('U',atoms%nat+1,1.d0,a,atoms%nat+1,h,1,0.d0,Ah,1)
        hAh= DDOT(atoms%nat+1,h,1,Ah,1)
        alpha = dot_gold/hAh

        do i=1,atoms%nat+1
            qq(i) = qq(i)+alpha*h(i)
            g(i) = g(i)+alpha*Ah(i)
        enddo
        dot_gnew = ddot(atoms%nat+1,g,1,g,1)
        beta = dot_gnew/dot_gold
        h(1:atoms%nat+1) = -g(1:atoms%nat+1)+beta*h(1:atoms%nat+1)
        dot_gold = dot_gnew
        resnorm = sqrt(dot_gold)
        if (resnorm<resnormtol) exit
    enddo

    atoms%qat(1:atoms%nat)=qq(1:atoms%nat)
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        call yaml_map('Lagrange',qq(atoms%nat+1))
        !write(*,*) 'Lagrange ',qq(atoms%nat+1),iter
    endif
    deallocate(ipiv)
    deallocate(qq)
    deallocate(h)
    deallocate(Ah)
    deallocate(g)
end subroutine get_qat_from_chi_iter
!*****************************************************************************************
subroutine get_ener_gradient_cent1(parini,poisson,ann_arr,atoms,g,qtot)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_poisson),intent(inout):: poisson
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: g(atoms%nat), qtot
    !local variables
    real(8):: gtot
    integer:: iat, igpx, igpy, igpz
    real(8), allocatable:: gausswidth(:)
    allocate(gausswidth(1:atoms%nat))
    gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
    if(.not. poisson%initialized) then
        stop 'ERROR: get_ener_gradient_cent1: poisson is not initialized!'
    endif
    poisson%reset_rho=.true.
    poisson%nat=atoms%nat
    poisson%cv=atoms%cellvec
    poisson%bc=atoms%boundcond
    poisson%q(1:poisson%nat)=atoms%qat(1:atoms%nat)
    poisson%gw(1:poisson%nat)=gausswidth(1:atoms%nat)
    call update_ratp(atoms)
    poisson%rcart(1:3,1:poisson%nat)=atoms%ratp(1:3,1:atoms%nat)
    call put_charge_density(parini,poisson)
    poisson%qgrad(1:atoms%nat)=0.d0
    call get_hartree(parini,poisson,atoms,gausswidth,ann_arr%epot_es)
    poisson%gw=poisson%gw_ewald
    call get_hartree_grad_rho(parini,poisson,atoms,ann_arr%epot_es)
    g=poisson%qgrad
    qtot=0.d0
    gtot=0.d0
    do iat=1,atoms%nat
        g(iat)=g(iat)+ann_arr%chi_o(iat)+atoms%qat(iat)*ann_arr%ann(atoms%itypat(iat))%hardness
        qtot= qtot+ atoms%qat(iat)
        gtot=gtot+g(iat)
    enddo
    do iat=1,atoms%nat
        g(iat)=g(iat)-gtot/atoms%nat
    enddo
    deallocate(gausswidth)
end subroutine get_ener_gradient_cent1
!*****************************************************************************************
subroutine get_qat_from_chi_operator(parini,poisson,ann_arr,atoms)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms, set_qat, update_ratp
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_poisson),intent(inout):: poisson
    !local variables
    integer:: info , iat , iter, igpx, igpy, igpz, niter_sd
    real(8), allocatable:: qq(:),g(:),h(:),gt(:)
    real(8) :: beta, dpm
    real(8) :: DDOT, aa, bb, gnrmtol, tt, qtot, y0, y1, gnrm, gnrm2, gnrm2old
    real(8) :: alpha, alphax, alpha0, de, epotlong_old, sss, qtot_tmp, dipole(3)
    gnrmtol=parini%gnrmtol_eem
    call set_qat(atoms)
    !atoms%qat=atoms%qat*0.8d0
    !do iat=1,atoms%nat
    !    write(*,'(i5,f10.3)') iat,atoms%qat(iat)
    !enddo
    !stop
    !open(unit=1358,file='charges.dat',status='old')
    !do iat=1,atoms%atoms%nat
    !    read(1358,*) tt,atoms%qat(iat)
    !enddo
    !close(1358)
    if(abs(atoms%qtot)>1.d-6 .and. trim(atoms%boundcond)/='free') then
        write(*,'(a)') 'ERROR: Operator approach not ready for ionized systems.'
        write(*,'(a)') '       Correct the initial guess of atomic charge assignment.'
        stop
    else
        sss=sum(atoms%qat(1:atoms%nat))
        sss=sss/real(atoms%nat,8)
        do iat=1,atoms%nat
            atoms%qat(iat)=atoms%qat(iat)-sss
        enddo
    endif
    allocate(qq(1:atoms%nat))
    allocate(g(1:atoms%nat))
    allocate(h(1:atoms%nat))
    allocate(gt(1:atoms%nat))
    !---------------------------------------------------------------
    qtot_tmp=sum(atoms%qat(1:atoms%nat))
    !if(parini%iverbose>=2) then
    !    write(*,'(a,es14.5)') 'cep begin ',qtot_tmp
    !endif
    if(trim(atoms%boundcond)=='slab') then
        alphax=0.4d0*parini%alphax_q
    else
        alphax=1.d0*parini%alphax_q
    endif
    if(parini%iverbose>=2) then
        call yaml_sequence_open('Charge equilibration process')
    endif
    alpha=1.d-1*alphax
    do iter=0,parini%nstep_cep
        if(parini%iverbose>=2) then
            call yaml_sequence(advance='no')
        endif
        call get_ener_gradient_cent1(parini,poisson,ann_arr,atoms,g,qtot)
        gnrm2=DDOT(atoms%nat,g,1,g,1)
        gnrm=sqrt(gnrm2)
        if(iter==0) epotlong_old=ann_arr%epot_es
        de=ann_arr%epot_es-epotlong_old
        if(parini%iverbose>=2) then
            dipole(1)=0.d0 ; dipole(2)=0.d0 ; dipole(3)=0.d0
            call update_ratp(atoms)
            do iat=1,atoms%nat
                !write(30+iter,'(i5,2f20.10)') iat,atoms%qat(iat),g(iat)
                dipole(1)=dipole(1)+atoms%qat(iat)*atoms%ratp(1,iat)
                dipole(2)=dipole(2)+atoms%qat(iat)*atoms%ratp(2,iat)
                dipole(3)=dipole(3)+atoms%qat(iat)*atoms%ratp(3,iat)
            enddo
            !write(*,'(a,4f10.3)') 'qtot,dipole ',qtot,dipole(1),dipole(2),dipole(3)
            !write(*,'(a,i5,es24.15,3es14.5)') 'cep: ',iter,ann_arr%epot_es,de,gnrm,alpha/alphax
            !call yaml_sequence(advance='no')
            call yaml_mapping_open('cep',flow=.true.)
            call yaml_map('iter',iter,fmt='(i5)')
            call yaml_map('epot_es',ann_arr%epot_es,fmt='(es22.13)')
            call yaml_map('de',de,fmt='(es11.2)')
            call yaml_map('gnrm',gnrm,fmt='(es12.3)')
            call yaml_map('stepsize',alpha/alphax,fmt='(es12.3)')
            call yaml_map('qtot',qtot,fmt='(f10.3)')
            call yaml_map('dpx',dipole(1),fmt='(f10.3)')
            call yaml_map('dpy',dipole(2),fmt='(f10.3)')
            call yaml_map('dpz',dipole(3),fmt='(f10.3)')
            call yaml_mapping_close()
        endif
        if(gnrm<1.d-7) then
            !write(*,'(a,i5,es24.15,3es14.5)') 'CEP converged: ', &
            !    iter,ann_arr%epot_es,de,gnrm,alpha/alphax
            if(parini%iverbose>=2) then
                call yaml_sequence_close()
            endif
            call yaml_mapping_open('CEP',flow=.true.)
            call yaml_map('iter',iter,fmt='(i5)')
            call yaml_map('epot_es',ann_arr%epot_es,fmt='(es22.13)')
            call yaml_map('de',de,fmt='(es11.2)')
            call yaml_map('gnrm',gnrm,fmt='(es12.3)')
            call yaml_map('stepsize',alpha/alphax,fmt='(es12.3)')
            call yaml_mapping_close()
            exit
        endif
        if(iter==0) then
            gt=g
            qq=atoms%qat
        endif
        tt=DDOT(atoms%nat,g,1,gt,1)/sqrt(DDOT(atoms%nat,g,1,g,1)*DDOT(atoms%nat,gt,1,gt,1))
        !write(*,'(a,i4,f10.3,2es14.5)') 'angle  ',iter,tt,sum((g-gt)**2),sum((qq-atoms%qat)**2)
        if(tt>0.5d0) then
            alpha=min(alpha*2.0d0,1.d0*alphax)
        else
            alpha=alpha*0.5d0
            ann_arr%epot_es=epotlong_old
            g=gt
            atoms%qat=qq
        endif
        epotlong_old=ann_arr%epot_es
        gt=g
        qq=atoms%qat
        do iat=1,atoms%nat
            atoms%qat(iat)=atoms%qat(iat)-alpha*g(iat)
        enddo
    enddo
    if(.not. (gnrm<1.d-7)) then
        if(parini%iverbose>=2) then
            call yaml_sequence_close()
        endif
    endif
    !do iat=1,atoms%nat
    !    write(*,'(i5,4f10.3)') iat,atoms%qat(iat),atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
    !stop
    niter_sd=iter
    !---------------------------------------------------------------
    iter=0
    alphax=0.8d0
    alpha0=1.d0*alphax
    alpha=0.d0
    gnrm2old=0.d0 !it is not supposed to be used for iter=0
    do
        exit !no CG
        call get_ener_gradient_cent1(parini,poisson,ann_arr,atoms,g,qtot)
        gnrm2=DDOT(atoms%nat,g,1,g,1)
        gnrm=sqrt(gnrm2)
        !dpm=0.d0
        !do iat=1,atoms%nat
        !    dpm=dpm+atoms%qat(iat)*atoms%rat(3,iat)
        !enddo
        !write(33,'(i4,es14.5)') iter,dpm
        !write(*,'(a,i6,2es14.5)') 'get_charge_slab ',iter,gnrm,qtot
        if(iter==0) epotlong_old=ann_arr%epot_es
        de=ann_arr%epot_es-epotlong_old
        qtot_tmp=sum(atoms%qat(1:atoms%nat))
        write(*,'(a,i5,es24.15,4es14.5)') 'iter,gnrm ',iter,ann_arr%epot_es,de,gnrm,alpha/alphax,qtot_tmp
        if(gnrm<gnrmtol) then
            if(parini%iverbose>=1) then
                write(*,'(a,2i5,2es14.5)') 'CEP converged ',iter,iter+niter_sd,gnrm,qtot_tmp
            endif
            exit
        endif
        if(iter>1000) then
            write(*,'(a)') 'ERROR: exceeds maximum number of iterations in CG eem1'
            stop
        endif
        if(iter==0) then
            h(1:atoms%nat)=-g(1:atoms%nat)
        else
            beta=gnrm2/gnrm2old
            h(1:atoms%nat)=-g(1:atoms%nat)+beta*h(1:atoms%nat)
        endif
        !obtaining inforrmation from a trial point
!        qq(1:atoms%nat)=atoms%qat(1:atoms%nat) !saving charges in a temporary array
!        do iat=1,atoms%nat
!            atoms%qat(iat)=atoms%qat(iat)+alpha0*h(iat)
!        enddo
        epotlong_old=ann_arr%epot_es !This must be here to avoid saving eh of trial point
!        call get_ener_gradient_cent1(parini,poisson,ann_arr,atoms,gt,qtot)
!        y0=DDOT(atoms%nat,gt,1,h,1)
!        y1=DDOT(atoms%nat,g,1,h,1)
!        tt=y0/(y0-y1)
!        !write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
!        !alpha=alpha0*tt
!        alpha=alpha0*max(min(tt,2.0d0),-0.2d0)
!        atoms%qat(1:atoms%nat)=qq(1:atoms%nat) !putting back the charges
        !updating charges
        alpha=alphax
        do iat=1,atoms%nat
            atoms%qat(iat)=atoms%qat(iat)+alpha*h(iat)
        enddo
        gnrm2old=gnrm2
        iter=iter+1
    enddo
    call charge_analysis(parini,atoms,ann_arr)
    deallocate(qq)
    deallocate(g)
    deallocate(h)
    deallocate(gt)
end subroutine get_qat_from_chi_operator
!*****************************************************************************************
subroutine get_dpm_cent1(atoms,dpx,dpy,dpz,dpm_err)
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
            dpx=dpx+(atoms%ratp(1,iat)-centroid_x)*(atoms%qat(iat))
            dpy=dpy+(atoms%ratp(2,iat)-centroid_y)*(atoms%qat(iat))
            dpz=dpz+(atoms%ratp(3,iat)-centroid_z)*(atoms%qat(iat))
        end do
        dpm_err=((dpx-atoms%dpm(1))**2+(dpy-atoms%dpm(2))**2+(dpz-atoms%dpm(3))**2)
end subroutine get_dpm_cent1
!*****************************************************************************************
