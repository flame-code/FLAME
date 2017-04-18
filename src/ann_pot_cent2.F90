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
    else
        ann_arr%fat_chi=0.d0
        ann_arr%chi_i=0.d0
        ann_arr%chi_o=0.d0
        ann_arr%chi_d=0.d0
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
    call get_qat_from_chi2(parini,ann_arr,atoms)
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
subroutine get_qat_from_chi2(parini,ann_arr,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iter, iat, niter, jat, igw
    real(8):: gnrm, epot_old, de, qtot, gnrm2, gtot, q1
    real(8), allocatable:: grad1(:)
    real(8), allocatable:: grad2(:)
    real(8), allocatable:: rel(:,:)
    real(8):: tt1, tt2, tt3
    allocate(grad1(3*atoms%nat))
    allocate(grad2(atoms%nat))
    allocate(rel(3,atoms%nat))

    rel(1:3,1:atoms%nat)=atoms%rat(1:3,1:atoms%nat) !+2.d-1

    do iat=1,atoms%nat
        if(trim(atoms%sat(iat))=='Na') atoms%qat(iat)= 0.6d0-atoms%zat(1)
        if(trim(atoms%sat(iat))=='Cl') atoms%qat(iat)=-0.6d0-atoms%zat(atoms%nat)
        if(trim(atoms%sat(iat))=='W' ) atoms%qat(iat)= 0.6d0-atoms%zat(1)
        if(trim(atoms%sat(iat))=='S' ) atoms%qat(iat)=-0.6d0-atoms%zat(atoms%nat)
    enddo

    niter=200
    do iter=0,niter
        call cal_potential_cent2(ann_arr,atoms,rel,grad1,grad2)
        if(iter==0) epot_old=atoms%epot
        de=atoms%epot-epot_old
        gnrm=sqrt(sum(grad1**2))
        gtot=sum(grad2(1:atoms%nat))
        grad2(1:atoms%nat)=grad2(1:atoms%nat)-gtot/atoms%nat
        gnrm2=sqrt(sum(grad2(1:atoms%nat)**2))
        qtot=sum(atoms%zat(1:atoms%nat))+sum(atoms%qat(1:atoms%nat))
        q1=atoms%zat(1)+atoms%qat(1)
        write(*,'(a,i5,es24.15,3es11.2,2f8.3)') 'iter,epot,gnrm ',iter,atoms%epot,de,gnrm,gnrm2,q1,qtot
        !write(41,'(a,i5,6f17.10,es14.5)') 'iter,dis ',iter,rel(1:3,1)-atoms%rat(1:3,1),rel(1:3,2)-atoms%rat(1:3,2),gnrm
        write(51,'(i5,2f8.3)') iter,rel(1,1)-atoms%rat(1,1),rel(1,2)-atoms%rat(1,2)
        if(gnrm<5.d-4 .and. gnrm2<1.d-3) exit
        if(iter==niter) exit
        do iat=1,atoms%nat
            rel(1,iat)=rel(1,iat)-2.d-1*grad1(3*iat-2)
            rel(2,iat)=rel(2,iat)-2.d-1*grad1(3*iat-1)
            rel(3,iat)=rel(3,iat)-2.d-1*grad1(3*iat-0)
        enddo
        do iat=1,atoms%nat
            atoms%qat(iat)=atoms%qat(iat)-0.5d0*grad2(iat)
        enddo
        epot_old=atoms%epot
    enddo
    write(*,'(a,i5,2f8.3)') 'DISP ',iter,rel(1,1)-atoms%rat(1,1),rel(1,2)-atoms%rat(1,2)

    deallocate(grad1)
    deallocate(grad2)
end subroutine get_qat_from_chi2
!*****************************************************************************************
subroutine cal_potential_cent2(ann_arr,atoms,rel,grad1,grad2)
    use mod_interface
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: rel(3,atoms%nat)
    real(8), intent(out):: grad1(3,atoms%nat), grad2(atoms%nat)
    !local variables
    integer:: iat
    real(8):: epot_es, dx, dy, dz, epot_es_bps, tt1, tt2
    real(8):: hardness, spring_const
    real(8), allocatable:: grad1_t(:,:)
    real(8), allocatable:: grad2_t(:)
    allocate(grad1_t(3,atoms%nat),source=0.d0)
    allocate(grad2_t(atoms%nat),source=0.d0)
    grad1=0.d0
    grad2=0.d0
    atoms%fat=0.d0
    atoms%epot=0.d0
    do iat=1,atoms%nat !summation over ions/electrons
        atoms%epot=atoms%epot+ann_arr%chi_o(iat)*(atoms%zat(iat)+atoms%qat(iat))
        hardness=ann_arr%ann(atoms%itypat(iat))%hardness
        atoms%epot=atoms%epot+0.5d0*hardness*(atoms%zat(iat)+atoms%qat(iat))**2
        grad2(iat)=grad2(iat)+ann_arr%chi_o(iat)+(atoms%zat(iat)+atoms%qat(iat))*hardness
        dx=rel(1,iat)-atoms%rat(1,iat)
        dy=rel(2,iat)-atoms%rat(2,iat)
        dz=rel(3,iat)-atoms%rat(3,iat)
        spring_const=ann_arr%ann(atoms%itypat(iat))%spring_const
        atoms%epot=atoms%epot+0.5d0*spring_const*(dx**2+dy**2+dz**2)
        grad1(1,iat)=grad1(1,iat)+spring_const*dx
        grad1(2,iat)=grad1(2,iat)+spring_const*dy
        grad1(3,iat)=grad1(3,iat)+spring_const*dz
    enddo
    grad1_t(1:3,1:atoms%nat)=grad1(1:3,1:atoms%nat)
    !call cal_electrostatic_cent2(nat,atoms%rat,rel,atoms%qat,epot_es,atoms%fat,grad1,grad2,epot_atom)
    !call cal_pot_with_bps(atoms%nat,atoms%rat,rel,atoms%qat,atoms%cellvec,epot_es_bps,grad1_t,grad2_t)
    !tt1=sum((grad2-grad2_t)**2)**0.5
    !tt2=sum((grad1-grad1_t)**2)**0.5
    !write(*,'(a,2es24.15,3es14.5)') 'DIFF ',epot_es,epot_es_bps,epot_es-epot_es_bps,tt1,tt2
    !write(*,'(a,2es24.15,3es14.5)') 'DIFF ',epot_es,epot_es_bps,epot_es-epot_es_bps,grad2(1),grad2_t(1)
    !write(*,'(a,2es24.15,3es14.5)') 'DIFF ',epot_es,epot_es_bps,epot_es-epot_es_bps,grad1(1,1),grad1_t(1,1)
    !stop
    epot_es=epot_es_bps
    grad1(1:3,1:atoms%nat)=grad1_t(1:3,1:atoms%nat)
    grad2(1:atoms%nat)=grad2_t(1:atoms%nat)
    atoms%epot=atoms%epot+epot_es
    !epot=epot+ener_ref
    deallocate(grad1_t)
    deallocate(grad2_t)
end subroutine cal_potential_cent2
!*****************************************************************************************
