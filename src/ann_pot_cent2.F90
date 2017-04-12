!*****************************************************************************************
subroutine cal_ann_eem2(parini,atoms,symfunc,ann_arr,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    use mod_electrostatics, only: typ_ewald_p3d
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_ekf), intent(inout):: ekf
    type(typ_ewald_p3d):: ewald_p3d
    !local variables
    integer:: iat, i, j, ng
    real(8):: epot_c, out_ann
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es
    call f_routine(id='cal_ann_eem2')
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        ann_arr%fat_chi=f_malloc0([1.to.3,1.to.atoms%nat],id='fat_chi')
        ann_arr%chi_i=f_malloc0([1.to.atoms%nat],id='ann_arr%chi_i')
        ann_arr%chi_o=f_malloc0([1.to.atoms%nat],id='ann_arr%chi_o')
        ann_arr%chi_d=f_malloc0([1.to.atoms%nat],id='ann_arr%chi_d')
        ann_arr%a=f_malloc0([1.to.(atoms%nat+1)*(atoms%nat+1)],id='a: aq=-chi')
    else
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
        ann_arr%a=0.d0
    endif
    if(trim(ann_arr%event)=='train') then
        !The following is allocated with ekf%num(1), this means number of
        !nodes in the input layer is the same for all atom types.
        !Therefore, it must be fixed later.
        !g_per_atom=f_malloc([1.to.ekf%num(1),1.to.atoms%nat],id='g_per_atom') !HERE
        do i=1,ann_arr%n
            call convert_x_ann(ekf%num(i),ekf%x(ekf%loc(i)),ann_arr%ann(i))
        enddo
    endif
    if(parini%iverbose>=2) call cpu_time(time1)
    !call cal_electrostatic_eem2(parini,'init',atoms,ann_arr,epot_c,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time2)
    !if(ann_arr%compute_symfunc) then !CORRECT_IT
    !    call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    !endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        ann_arr%fatpq=f_malloc([1.to.3,1.to.symfunc%linked_lists%maxbound_rad],id='fatpq')
        ann_arr%stresspq=f_malloc([1.to.3,1.to.3,1.to.symfunc%linked_lists%maxbound_rad],id='stresspq')
    endif
    if(parini%iverbose>=2) call cpu_time(time3)
    goto 1000 !CORRECT_IT
    over_iat: do iat=1,atoms%nat
        i=atoms%itypat(iat)
        ng=ann_arr%ann(i)%nn(0)
        !ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,iat) !CORRECT_IT
        if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
            call cal_architecture(ann_arr%ann(i),out_ann)
            call cal_force_chi_part1(parini,symfunc,iat,atoms,out_ann,ann_arr)
        elseif(trim(ann_arr%event)=='train') then
            !call cal_architecture_der(ann_arr%ann(i),out_ann) !CORRECT_IT
            if(trim(atoms%sat(iat))=='Na') out_ann=-0.25 !CORRECT_IT
            if(trim(atoms%sat(iat))=='Cl') out_ann= 0.25 !CORRECT_IT
            ann_arr%chi_i(iat)=out_ann
            tt1=tanh(ann_arr%ann(i)%prefactor_chi*out_ann)
            ann_arr%chi_o(iat)=ann_arr%ann(i)%ampl_chi*tt1+ann_arr%ann(i)%chi0
            call convert_ann_epotd(ann_arr%ann(i),ekf%num(i),ann_arr%g_per_atom(1,iat))
            ann_arr%g_per_atom(1:ekf%num(1),iat)=ann_arr%g_per_atom(1:ekf%num(1),iat)*ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)
        else
            stop 'ERROR: undefined content for ann_arr%event'
        endif
    enddo over_iat
1000  continue !CORRECT_IT
    do iat=1,atoms%nat
        !write(*,*) iat,trim(atoms%sat(iat))
        if(trim(atoms%sat(iat))=='Na') ann_arr%chi_i(iat)=-0.25d0
        if(trim(atoms%sat(iat))=='Cl') ann_arr%chi_i(iat)= 0.25d0
    enddo
    write(*,*) ann_arr%chi_i(1:atoms%nat)
    !stop 'AFTER LOOP OVER NN'
    !This must be here since contribution from coulomb
    !interaction is calculated during the process of charge optimization.
    atoms%stress(1:3,1:3)=0.d0
    !This msut be here otherwise it will zero forces which were calculated by kwald.
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(parini%iverbose>=2) call cpu_time(time4)
    call get_qat_from_chi2(parini,ann_arr,atoms,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time5)
    !if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then !CORRECT_IT
    !    call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    !endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    !call cal_electrostatic_eem2(parini,'calculate',atoms,ann_arr,epot_c,ann_arr%a)
    if(parini%iverbose>=2) then
        call cpu_time(time7)
        write(*,'(a,f8.3)') 'Timing:cent2: initialize matrix          ',time2-time1
        write(*,'(a,f8.3)') 'Timing:cent2: calculation of symfunc     ',time3-time2
        write(*,'(a,f8.3)') 'Timing:cent2: neural network process     ',time4-time3
        write(*,'(a,f8.3)') 'Timing:cent2: linear equations solver    ',time5-time4
        write(*,'(a,f8.3)') 'Timing:cent2: force (SR term)            ',time6-time5
        write(*,'(a,f8.3)') 'Timing:cent2: energy (SR+LR), force (LR) ',time7-time6
        write(*,'(a,f8.3)') 'Timing:cent2: total time                 ',time7-time1
    endif !end of if for printing out timing.
    atoms%epot=epot_c
    if(trim(ann_arr%event)=='evalu' .and. atoms%nat<=parini%nat_force) then
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
        ann_arr%fchi_norm=tt2/max(tt1,1.d-3)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train' .and. trim(parini%symfunc)/='do_not_save')) then
        call f_free(symfunc%linked_lists%prime_bound)
        call f_free(symfunc%linked_lists%bound_rad)
        call f_free(symfunc%linked_lists%bound_ang)
    endif
    if(trim(ann_arr%event)=='potential' .or. trim(parini%symfunc)=='do_not_save') then
        call f_free(symfunc%y)
        call f_free(symfunc%y0d)
        call f_free(symfunc%y0dr)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        call f_free(ann_arr%chi_i)
        call f_free(ann_arr%chi_o)
        call f_free(ann_arr%chi_d)
        call f_free(ann_arr%a)
        call f_free(ann_arr%fat_chi)
        call f_free(ann_arr%fatpq)
        call f_free(ann_arr%stresspq)
    endif
    if(trim(ann_arr%event)=='train') then
        ekf%g(1:ekf%n)=0.d0
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ekf%num(1)
                ekf%g(ekf%loc(i)+j-1)=ekf%g(ekf%loc(i)+j-1)+atoms%qat(iat)*ann_arr%g_per_atom(j,iat)
            enddo
        enddo
    endif
    call f_release_routine()
end subroutine cal_ann_eem2
!*****************************************************************************************
subroutine get_qat_from_chi2(parini,ann_arr,atoms,a)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    !local variables
end subroutine get_qat_from_chi2
!*****************************************************************************************
