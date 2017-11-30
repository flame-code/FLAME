!*****************************************************************************************
subroutine cal_ann_cent1(parini,atoms,symfunc,ann_arr,ekf)
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
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol
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
    call cal_electrostatic_cent1(parini,'init',atoms,ann_arr,epot_c,ann_arr%a,ewald_p3d)
    if(parini%iverbose>=2) call cpu_time(time2)
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
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
            call convert_ann_epotd(ann_arr%ann(i),ekf%num(i),ann_arr%g_per_atom(1,iat))
            ann_arr%g_per_atom(1:ekf%num(1),iat)=ann_arr%g_per_atom(1:ekf%num(1),iat)*ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)
        else
            stop 'ERROR: undefined content for ann_arr%event'
        endif
    enddo over_iat
    !This must be here since contribution from coulomb
    !interaction is calculated during the process of charge optimization.
    atoms%stress(1:3,1:3)=0.d0
    !This msut be here otherwise it will zero forces which were calculated by kwald.
    atoms%fat(1:3,1:atoms%nat)=0.d0
    if(parini%iverbose>=2) call cpu_time(time4)
    call get_qat_from_chi(parini,ann_arr,atoms,ewald_p3d,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time5)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call cal_electrostatic_cent1(parini,'calculate',atoms,ann_arr,epot_c,ann_arr%a,ewald_p3d)
    if(parini%iverbose>=2) then
        call cpu_time(time7)
        write(*,'(a,f8.3)') 'Timing:cent1: initialize matrix          ',time2-time1
        write(*,'(a,f8.3)') 'Timing:cent1: calculation of symfunc     ',time3-time2
        write(*,'(a,f8.3)') 'Timing:cent1: neural network process     ',time4-time3
        write(*,'(a,f8.3)') 'Timing:cent1: linear equations solver    ',time5-time4
        write(*,'(a,f8.3)') 'Timing:cent1: force (SR term)            ',time6-time5
        write(*,'(a,f8.3)') 'Timing:cent1: energy (SR+LR), force (LR) ',time7-time6
        write(*,'(a,f8.3)') 'Timing:cent1: total time                 ',time7-time1
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
        ann_arr%fchi_norm=tt2/tt1
    endif
    if(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk') then
        call destruct_ewald_p3d(parini,atoms,ewald_p3d)
    endif
    !call repulsive_potential_cent(parini,atoms,ann_arr)
    call getvol_alborz(atoms%cellvec,vol)
    call invertmat_alborz(atoms%cellvec,hinv)
    do i=1,3
    do j=1,3
        atoms%celldv(i,j)=vol*(atoms%stress(i,1)*hinv(j,1)+atoms%stress(i,2)*hinv(j,2)+atoms%stress(i,3)*hinv(j,3))
    enddo
    enddo

    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train' .and. trim(parini%symfunc)/='do_not_save')) then
        deallocate(symfunc%linked_lists%prime_bound)
        deallocate(symfunc%linked_lists%bound_rad)
        deallocate(symfunc%linked_lists%bound_ang)
    endif
    if(trim(ann_arr%event)=='potential' .or. trim(parini%symfunc)=='do_not_save') then
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
    if(trim(ann_arr%event)=='train') then
        ekf%g(1:ekf%n)=0.d0
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ekf%num(1)
                ekf%g(ekf%loc(i)+j-1)=ekf%g(ekf%loc(i)+j-1)+atoms%qat(iat)*ann_arr%g_per_atom(j,iat)
            enddo
        enddo
        do i=1,ann_arr%n
    !        ekf%g(ekf%loc(i)+ekf%num(1)-1)=ekf%g(ekf%loc(i)+ekf%num(1)-1)*1.d-4
            !write(*,*) 'GGG ',ia,ekf%loc(ia)+ekf%num(1)-1
        enddo
    endif
    call f_release_routine()
end subroutine cal_ann_cent1
!*****************************************************************************************
subroutine get_qat_from_chi(parini,ann_arr,atoms,ewald_p3d,a)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: iat
    real(8):: qtot
    character(200):: smsg
    if(trim(parini%syslinsolver_ann)=='direct') then
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
    elseif(trim(parini%syslinsolver_ann)=='apply_matrix') then
        if(trim(atoms%boundcond)=='free') then
            call get_qat_from_chi_iter(parini,ann_arr,atoms,a)
        else
            smsg='ERROR: solving linear system of equations with explicitly'
            smsg=trim(smsg)//' applying matrix is possible only for free BC, atoms%boundcond= '
            write(*,'(2a)') trim(smsg),trim(atoms%boundcond)
            stop
        endif
    elseif(trim(parini%syslinsolver_ann)=='operator') then
        call get_qat_from_chi_operator(parini,ewald_p3d,ann_arr,atoms)
    else
        stop 'ERROR: unknown syslinsolver'
    endif
end subroutine get_qat_from_chi
!*****************************************************************************************
subroutine get_qat_from_chi_dir(parini,ann_arr,atoms,a)
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
    integer:: info !, iat
    associate(nat=>atoms%nat)
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        allocate(ann_arr%ipiv(1:nat+1))
        allocate(ann_arr%qq(1:nat+1))
    endif
    call DGETRF(nat+1,nat+1,a,nat+1,ann_arr%ipiv,info)
    if(info/=0) then
        write(*,'(a,i)') 'ERROR: DGETRF info=',info
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
        write(*,'(a,i)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat(1:nat)=ann_arr%qq(1:nat)
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        write(*,*) 'Lagrange ',ann_arr%qq(nat+1)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        deallocate(ann_arr%ipiv)
        deallocate(ann_arr%qq)
    endif
    end associate
end subroutine get_qat_from_chi_dir
!*****************************************************************************************
subroutine cal_electrostatic_cent1(parini,str_job,atoms,ann_arr,epot_c,a,ewald_p3d)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: str_job
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_c
    real(8), intent(inout):: a(atoms%nat+1,atoms%nat+1)
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    integer:: iat, jat
    real(8):: vol, c
    real(8):: dx, dy, dz, r, tt1, tt2, pi, beta_iat, beta_jat, gama, ttf
    associate(epot_es=>ann_arr%epot_es)
    if(trim(str_job)=='init') then
        pi=4.d0*atan(1.d0)
        ann_arr%ener_ref=0.d0
        do iat=1,atoms%nat
            ann_arr%ener_ref=ann_arr%ener_ref+ann_arr%ann(atoms%itypat(iat))%ener_ref
        enddo
        !ann_arr%ener_ref=ann_arr%ener_ref-1.d-2*atoms%nat**(-1.0/3.0)
        !ann_arr%ener_ref=0.0144273d0*atoms%nat**(-1.0/3.0)+0.081049d0*atoms%nat**(-1.3/2.0)-0.161332d0+0.142238365422159d0  +  atoms%nat*1.d-3/27.2d0
        if (.not. parini%ewald) then 
            ewald_p3d%alpha = maxval(ann_arr%ann(:)%gausswidth)
        else 
            if (parini%alpha_ewald<0.d0) then
                call getvol_alborz(atoms%cellvec,vol)
                c=2.2d0
                ewald_p3d%alpha = 1.d0/(c*sqrt(pi)*(atoms%nat/vol**2)**(1.d0/6.d0))
                write(*,*)"optimized alpha = ", ewald_p3d%alpha
            else
                ewald_p3d%alpha=parini%alpha_ewald
            endif
        end if
        if(trim(atoms%boundcond)=='bulk' .or. trim(atoms%boundcond)=='slab') then
            call construct_ewald_p3d(parini,atoms,ewald_p3d)
        else
            do iat=1,atoms%nat
                a(iat,atoms%nat+1)=1.d0
                a(atoms%nat+1,iat)=1.d0
                beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
                gama=1.d0/sqrt(beta_iat**2+beta_iat**2)
                a(iat,iat)=gama*2.d0/sqrt(pi)+ann_arr%ann(atoms%itypat(iat))%hardness
                do jat=iat+1,atoms%nat
                    dx=atoms%rat(1,jat)-atoms%rat(1,iat)
                    dy=atoms%rat(2,jat)-atoms%rat(2,iat)
                    dz=atoms%rat(3,jat)-atoms%rat(3,iat)
                    r=sqrt(dx*dx+dy*dy+dz*dz)
                    beta_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth
                    gama=1.d0/sqrt(beta_iat**2+beta_jat**2)
                    a(iat,jat)=erf(gama*r)/r
                    a(jat,iat)=a(iat,jat)
                enddo
            enddo
            a(atoms%nat+1,atoms%nat+1)=0.d0
        endif
    elseif(trim(str_job)=='calculate') then
        tt1=0.d0
        tt2=0.d0
        do iat=1,atoms%nat
            tt1=tt1+ann_arr%chi_o(iat)*atoms%qat(iat)
            tt2=tt2+atoms%qat(iat)**2*0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness
        enddo
        call cal_electrostatic_ann(parini,atoms,ann_arr,a,ewald_p3d)
        epot_c=epot_es+tt1+tt2+ann_arr%ener_ref
        if(parini%iverbose>=1) then
        !write(81,'(i6,4es14.5,3f7.1)') atoms%nat,epot_c,tt1,tt2,epot_es,1.d2*tt1/epot_c,1.d2*tt2/epot_c,1.d2*epot_es/epot_c
        endif
    else
        stop 'ERROR: unknown job in cal_electrostatic_eem1'
    endif
    end associate
end subroutine cal_electrostatic_cent1
!*****************************************************************************************
subroutine cal_electrostatic_ann(parini,atoms,ann_arr,a,ewald_p3d)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    use mod_electrostatics, only: typ_ewald_p3d
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: a(atoms%nat+1,atoms%nat+1)
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    !local variables
    integer:: iat, jat
    real(8):: tt2, tt3, ttf, gama, pi
    real(8):: dx, dy, dz, r, beta_iat, beta_jat
    real(8), allocatable:: gausswidth(:)
    allocate(gausswidth(1:atoms%nat))
    if(trim(atoms%boundcond)=='free') then
        pi=4.d0*atan(1.d0)
        tt2=0.d0
        tt3=0.d0
        do iat=1,atoms%nat
            beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
            gama=1.d0/sqrt(beta_iat**2+beta_iat**2)
            tt2=tt2+atoms%qat(iat)**2*gama/sqrt(pi)
            do jat=iat+1,atoms%nat
                dx=atoms%rat(1,jat)-atoms%rat(1,iat)
                dy=atoms%rat(2,jat)-atoms%rat(2,iat)
                dz=atoms%rat(3,jat)-atoms%rat(3,iat)
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
    elseif(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk') then
        if(trim(parini%psolver_ann)/='kwald') then
            if(parini%ewald) then 
                gausswidth(:)=ewald_p3d%alpha
            else
                gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
            endif
            call longerange_forces(atoms,ewald_p3d,gausswidth)
        endif
    else
        stop 'ERROR: the requested BCs is not yet implemented.'
    endif
    deallocate(gausswidth)
end subroutine cal_electrostatic_ann
!*****************************************************************************************
subroutine charge_analysis(parini,atoms,ann_arr)
    use mod_interface
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
    do i=1,ann_arr%n
        ann_arr%chi_delta(i)=max(ann_arr%chi_delta(i),chi_max_per_conf(i)-chi_min_per_conf(i))
    enddo
    !if(parini%iverbose>=1) then
    !write(61,'(a,2(a,4f8.3),es11.2,i5)') trim(str),' Na=',tt1/ii1,tt1min,tt1max,tt1max-tt1min,' Cl=',tt2/ii2,tt2min,tt2max,tt2max-tt2min,tt1+tt2,atoms%nat
    !write(71,'(a,2(a,4f8.3),es11.2,i5)') trim(str),' Na=',ss1/ii1,ss1min,ss1max,ss1max-ss1min,' Cl=',ss2/ii2,ss2min,ss2max,ss2max-ss2min,ss1+ss2,atoms%nat
    !endif
end subroutine charge_analysis
!*****************************************************************************************
subroutine get_qat_from_chi_iter(parini,ann_arr,atoms,a)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use dynamic_memory
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
    associate(nat=>atoms%nat)
    allocate(ipiv(1:nat+1))
    allocate(qq(1:nat+1))
    allocate(h(1:nat+1))
    allocate(Ah(1:nat+1))
    allocate(g(1:nat+1))

    qq(1:nat+1)=0.d0  !??????????????????
    resnormtol=1.d-10
    niter=1000        !?????????????????? 
    iter=0
    g(1:nat)=-ann_arr%chi_o(1:nat)
    g(nat+1)=atoms%qtot
    call DSYMV ('U',nat+1,1.d0,a,nat+1,qq,1,-1.d0,g,1)
    h(1:nat+1)=-g(1:nat+1)
    dot_gold = DDOT(nat+1,g,1,g,1)
    do j=1,niter
        iter = iter+1
        call DSYMV ('U',nat+1,1.d0,a,nat+1,h,1,0.d0,Ah,1)
        hAh= DDOT(nat+1,h,1,Ah,1)
        alpha = dot_gold/hAh

        do i=1,nat+1
            qq(i) = qq(i)+alpha*h(i)
            g(i) = g(i)+alpha*Ah(i)
        enddo
        dot_gnew = ddot(nat+1,g,1,g,1)
        beta = dot_gnew/dot_gold
        h(1:nat+1) = -g(1:nat+1)+beta*h(1:nat+1)
        dot_gold = dot_gnew
        resnorm = sqrt(dot_gold)
        if (resnorm<resnormtol) exit
    enddo

    atoms%qat(1:nat)=qq(1:nat)
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>1) then
        write(*,*) 'Lagrange ',qq(nat+1),iter
    endif
    deallocate(ipiv)
    deallocate(qq)
    deallocate(h)
    deallocate(Ah)
    deallocate(g)
    end associate
end subroutine get_qat_from_chi_iter
!*****************************************************************************************
subroutine cal_ugradient(parini,ewald_p3d,ann_arr,atoms,g,qtot)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(out):: g(atoms%nat), qtot
    !local variables
    real(8):: dpm, pi, gtot
    integer:: iat, igpx, igpy, igpz
    real(8), allocatable:: gausswidth(:)
    pi=4.d0*atan(1.d0)
    allocate(gausswidth(1:atoms%nat))
    gausswidth(:)=ann_arr%ann(atoms%itypat(:))%gausswidth
    call get_hartree(parini,ewald_p3d,atoms,gausswidth,ann_arr%epot_es,g)
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
end subroutine cal_ugradient
!*****************************************************************************************
subroutine get_qat_from_chi_operator(parini,ewald_p3d,ann_arr,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_ewald_p3d),intent(inout):: ewald_p3d
    !local variables
    integer:: info , iat , iter, igpx, igpy, igpz, niter_sd
    real(8), allocatable:: qq(:),g(:),h(:),gt(:)
    real(8) :: beta, dpm
    real(8) :: DDOT, aa, bb, gnrmtol, tt, qtot, y0, y1, gnrm, gnrm2, gnrm2old
    real(8) :: alpha, alphax, alpha0, de, epotlong_old, sss, qtot_tmp, dipole(3)
    associate(nat=>atoms%nat)
    gnrmtol=parini%gnrmtol_eem
    call set_qat(atoms)
    !atoms%qat=atoms%qat*0.8d0
    !do iat=1,nat
    !    write(*,'(i5,f10.3)') iat,atoms%qat(iat)
    !enddo
    !stop
    !open(unit=1358,file='charges.dat',status='old')
    !do iat=1,atoms%nat
    !    read(1358,*) tt,atoms%qat(iat)
    !enddo
    !close(1358)
    if(abs(atoms%qtot)>1.d-10) then
        write(*,'(a)') 'ERROR: Operator approach not ready for ionized systems.'
        write(*,'(a)') '       Correct the initial guess of atomic charge assignment.'
        stop
    else
        sss=sum(atoms%qat(1:nat))
        sss=sss/real(nat,8)
        do iat=1,nat
            atoms%qat(iat)=atoms%qat(iat)-sss
        enddo
    endif
    allocate(qq(1:nat))
    allocate(g(1:nat))
    allocate(h(1:nat))
    allocate(gt(1:nat))
    !---------------------------------------------------------------
    qtot_tmp=sum(atoms%qat(1:atoms%nat))
    if(parini%iverbose>=2) then
        write(*,'(a,es14.5)') 'cep begin ',qtot_tmp
    endif
    if(trim(atoms%boundcond)=='slab') then
        alphax=4.d-1
    else
        alphax=1.d0
    endif
    alpha=1.d-1*alphax
    do iter=0,1000
        call cal_ugradient(parini,ewald_p3d,ann_arr,atoms,g,qtot)
        if(parini%iverbose>=2) then
            dipole(1)=0.d0 ; dipole(2)=0.d0 ; dipole(3)=0.d0
            do iat=1,nat
                !write(30+iter,'(i5,2f20.10)') iat,atoms%qat(iat),g(iat)
                dipole(1)=dipole(1)+atoms%qat(iat)*atoms%rat(1,iat)
                dipole(2)=dipole(2)+atoms%qat(iat)*atoms%rat(2,iat)
                dipole(3)=dipole(3)+atoms%qat(iat)*atoms%rat(3,iat)
            enddo
            write(*,'(a,4f10.3)') 'qtot,dipole ',qtot,dipole(1),dipole(2),dipole(3)
        endif
        gnrm2=DDOT(nat,g,1,g,1)
        gnrm=sqrt(gnrm2)
        if(iter==0) epotlong_old=ann_arr%epot_es
        de=ann_arr%epot_es-epotlong_old
        if(parini%iverbose>=2) then
            write(*,'(a,i5,es24.15,3es14.5)') 'iter,gnrm ',iter,ann_arr%epot_es,de,gnrm,alpha/alphax
        endif
        if(gnrm<1.d-7) exit
        if(iter==0) then
            gt=g
            qq=atoms%qat
        endif
        tt=DDOT(nat,g,1,gt,1)/sqrt(DDOT(nat,g,1,g,1)*DDOT(nat,gt,1,gt,1))
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
        do iat=1,nat
            atoms%qat(iat)=atoms%qat(iat)-alpha*g(iat)
        enddo
    enddo
    !do iat=1,nat
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
        call cal_ugradient(parini,ewald_p3d,ann_arr,atoms,g,qtot)
        gnrm2=DDOT(nat,g,1,g,1)
        gnrm=sqrt(gnrm2)
        dpm=0.d0
        do iat=1,nat
            dpm=dpm+atoms%qat(iat)*atoms%rat(3,iat)
        enddo
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
            h(1:nat)=-g(1:nat)
        else
            beta=gnrm2/gnrm2old
            h(1:nat)=-g(1:nat)+beta*h(1:nat)
        endif
        !obtaining inforrmation from a trial point
!        qq(1:nat)=atoms%qat(1:nat) !saving charges in a temporary array
!        do iat=1,nat
!            atoms%qat(iat)=atoms%qat(iat)+alpha0*h(iat)
!        enddo
        epotlong_old=ann_arr%epot_es !This must be here to avoid saving eh of trial point
!        call cal_ugradient(parini,ewald_p3d,ann_arr,atoms,gt,qtot)
!        y0=DDOT(nat,gt,1,h,1)
!        y1=DDOT(nat,g,1,h,1)
!        tt=y0/(y0-y1)
!        !write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
!        !alpha=alpha0*tt
!        alpha=alpha0*max(min(tt,2.0d0),-0.2d0)
!        atoms%qat(1:nat)=qq(1:nat) !putting back the charges
        !updating charges
        alpha=alphax
        do iat=1,nat
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
    if(parini%iverbose>=2) then
        do iat=1,atoms%nat
            write(*,*) 'charge on atom ',iat,atoms%qat(iat) 
        enddo
    endif
    end associate
end subroutine get_qat_from_chi_operator
!******************************************************************************************************************
