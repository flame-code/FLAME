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
    call cal_electrostatic_eem2(parini,'init',atoms,ann_arr,epot_c,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time2)
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    endif
    if(.not. (trim(parini%task)=='ann' .and. trim(parini%subtask_ann)=='train')) then
        ann_arr%fatpq=f_malloc([1.to.3,1.to.symfunc%linked_lists%maxbound_rad],id='fatpq')
        ann_arr%stresspq=f_malloc([1.to.3,1.to.3,1.to.symfunc%linked_lists%maxbound_rad],id='stresspq')
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
    call get_qat_from_chi2(parini,ann_arr,atoms,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time5)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
    if(parini%iverbose>=2) call cpu_time(time6)
    call cal_electrostatic_eem2(parini,'calculate',atoms,ann_arr,epot_c,ann_arr%a)
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
        ann_arr%fchi_norm=tt2/tt1
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
    integer:: info, iat, jat, ztot
    integer, allocatable:: ipiv(:)
    real(8), allocatable:: qq(:)
    real(8):: pi, sqrtpiinv, beta_iat, beta_jat
    real(8):: dx, dy, dz, r
    real(8):: gamau, alpha_iat, alpha_jat, gamau_ijat, gamau_jiat
    associate(nat=>atoms%nat)
    pi=4.d0*atan(1.d0)
    allocate(ipiv(nat+1))
    allocate(qq(nat+1))
    call DGETRF(nat+1,nat+1,a,nat+1,ipiv,info)
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
    sqrtpiinv=1.d0/sqrt(pi)
    qq(1:nat)=0.d0
    do iat=1,nat
        beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
        alpha_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
        gamau=1.d0/sqrt(beta_iat**2+alpha_iat**2)
        qq(iat)=qq(iat)-ann_arr%chi_o(iat)-sqrtpiinv*gamau*atoms%zat(iat)
        qq(iat)=qq(iat)-ann_arr%ann(atoms%itypat(iat))%hardness*atoms%zat(iat)
        !qq(iat)=qq(iat)-ann_arr%chi_o(iat)-sqrtpiinv*beta_iatinv*atoms%zat(iat)
        do jat=iat+1,atoms%nat
            dx=atoms%rat(1,jat)-atoms%rat(1,iat)
            dy=atoms%rat(2,jat)-atoms%rat(2,iat)
            dz=atoms%rat(3,jat)-atoms%rat(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            beta_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth
            alpha_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth_ion
            gamau_ijat=1.d0/sqrt(beta_iat**2+alpha_jat**2)
            gamau_jiat=1.d0/sqrt(beta_jat**2+alpha_iat**2)
            qq(iat)=qq(iat)-atoms%zat(jat)*erf(gamau_ijat*r)/r
            qq(jat)=qq(jat)-atoms%zat(iat)*erf(gamau_jiat*r)/r
        enddo
    enddo
    qq(nat+1)=atoms%qtot-atoms%ztot
    call DGETRS('N',nat+1,1,a,nat+1,ipiv,qq,nat+1,info)
    if(info/=0) then
        write(*,'(a,i)') 'ERROR: DGETRS info=',info
        stop
    endif
    atoms%qat(1:nat)=qq(1:nat)
    call charge_analysis(parini,atoms,ann_arr)
    if(parini%iverbose>=1) then
        write(*,*) 'Lagrange ',qq(nat+1)
    endif
    deallocate(ipiv,qq)
    end associate
end subroutine get_qat_from_chi2
!*****************************************************************************************
subroutine cal_electrostatic_eem2(parini,str_job,atoms,ann_arr,epot_es,a)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: str_job
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(out):: epot_es
    real(8), intent(out):: a(atoms%nat+1,atoms%nat+1)
    !local variables
    integer:: iat, jat
    real(8):: dx, dy, dz, r, tt1, tt2, tt3, tt4, tt5, tt6 ,tt7, ttf
    real(8):: pi, beta_iat, beta_jat, gama
    real(8):: gamau, alpha_iat, alpha_jat, gamau_ijat, gamau_jiat
    real(8):: sqrt2, sqrt2inv, beta_iatinv, sqrtpiinv
    real(8):: rinv, r2inv, r3inv
    pi=4.d0*atan(1.d0)
    sqrt2=sqrt(2.d0)
    sqrt2inv=1.d0/sqrt2
    sqrtpiinv=1.d0/sqrt(pi)
    atoms%ztot=0.d0
    do iat=1,atoms%nat
        atoms%zat(iat)=ann_arr%ann(atoms%itypat(iat))%zion
        atoms%ztot=atoms%ztot+atoms%zat(iat)
    enddo
    if(trim(str_job)=='init') then
        ann_arr%ener_ref=0.d0
        do iat=1,atoms%nat
            a(iat,atoms%nat+1)=1.d0
            a(atoms%nat+1,iat)=1.d0
            beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
            ann_arr%ener_ref=ann_arr%ener_ref+ann_arr%ann(atoms%itypat(iat))%ener_ref
            gama=1.d0/beta_iat *sqrt2inv
            a(iat,iat)=gama*2.d0*sqrtpiinv+ann_arr%ann(atoms%itypat(iat))%hardness
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
    elseif(trim(str_job)=='calculate') then
        tt1=0.d0
        tt2=0.d0
        tt3=0.d0
        tt4=0.d0
        tt5=0.d0
        tt6=0.d0
        tt7=0.d0
        do iat=1,atoms%nat
            !tt1=tt1+ann_arr%chi_o(iat)*atoms%qat(iat)
            tt1=tt1+ann_arr%chi_o(iat)*(atoms%qat(iat)+atoms%zat(iat))
            tt1=tt1+0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness*(atoms%qat(iat)+atoms%zat(iat))**2
            beta_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth
            alpha_iat=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
            gama=1.d0/beta_iat*sqrt2inv
            gamau=1.d0/sqrt(beta_iat**2+alpha_iat**2)
            !tt2=tt2+atoms%qat(iat)**2*(gama*sqrtpiinv+0.5d0*ann_arr%ann(atoms%itypat(iat))%hardness)
            tt2=tt2+atoms%qat(iat)**2*gama*sqrtpiinv
            tt3=tt3+gamau*sqrtpiinv*atoms%qat(iat)*atoms%zat(iat)
            do jat=iat+1,atoms%nat
                dx=atoms%rat(1,jat)-atoms%rat(1,iat)
                dy=atoms%rat(2,jat)-atoms%rat(2,iat)
                dz=atoms%rat(3,jat)-atoms%rat(3,iat)
                r=sqrt(dx*dx+dy*dy+dz*dz)
                beta_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth
                alpha_jat=ann_arr%ann(atoms%itypat(jat))%gausswidth_ion
                gama=1.d0/sqrt(beta_iat**2+beta_jat**2)
                gamau=1.d0/sqrt(alpha_iat**2+alpha_jat**2)
                gamau_ijat=1.d0/sqrt(beta_iat**2+alpha_jat**2)
                gamau_jiat=1.d0/sqrt(beta_jat**2+alpha_iat**2)
                tt4=tt4+atoms%qat(iat)*atoms%qat(jat)*erf(gama*r)/r
                tt5=tt5+atoms%qat(iat)*atoms%zat(jat)*erf(gamau_ijat*r)/r
                tt6=tt6+atoms%qat(jat)*atoms%zat(iat)*erf(gamau_jiat*r)/r
                tt7=tt7+atoms%zat(jat)*atoms%zat(iat)*erf(gamau*r)/r
                rinv=1.d0/r
                r2inv=rinv*rinv
                r3inv=r2inv*rinv
                ttf=(2.d0*sqrtpiinv*gama*exp(-(gama*r)**2)*r2inv-erf(gama*r)*r3inv)*atoms%qat(iat)*atoms%qat(jat)
                ttf=ttf+(2.d0*sqrtpiinv*gamau_ijat*exp(-(gamau_ijat*r)**2)*r2inv-erf(gamau_ijat*r)*r3inv)*atoms%qat(iat)*atoms%zat(jat)
                ttf=ttf+(2.d0*sqrtpiinv*gamau_jiat*exp(-(gamau_jiat*r)**2)*r2inv-erf(gamau_jiat*r)*r3inv)*atoms%qat(jat)*atoms%zat(iat)
                ttf=ttf+(2.d0*sqrtpiinv*gamau*exp(-(gamau*r)**2)*r2inv-erf(gamau*r)*r3inv)*atoms%zat(iat)*atoms%zat(jat)
                atoms%fat(1,iat)=atoms%fat(1,iat)+ttf*dx
                atoms%fat(2,iat)=atoms%fat(2,iat)+ttf*dy
                atoms%fat(3,iat)=atoms%fat(3,iat)+ttf*dz
                atoms%fat(1,jat)=atoms%fat(1,jat)-ttf*dx
                atoms%fat(2,jat)=atoms%fat(2,jat)-ttf*dy
                atoms%fat(3,jat)=atoms%fat(3,jat)-ttf*dz
            enddo
        enddo
        epot_es=tt1+tt2+tt3+tt4+tt5+tt6+tt7+ann_arr%ener_ref
        if(parini%iverbose>=1) then
        !write(81,'(i6,4es14.5,3f7.1)') atoms%nat,epot_es,tt1,tt2,tt3,1.d2*tt1/epot_es,1.d2*tt2/epot_es,1.d2*tt3/epot_es
        endif
    else
        stop 'ERROR: unknown job in cal_electrostatic_eem2'
    endif
end subroutine cal_electrostatic_eem2
!*****************************************************************************************
