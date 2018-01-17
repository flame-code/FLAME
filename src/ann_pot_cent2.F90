!*****************************************************************************************
subroutine cal_ann_cent2(parini,atoms,symfunc,ann_arr,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf, typ_cent
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_ekf), intent(inout):: ekf
    !local variables
    type(typ_pia_arr):: pia_arr_tmp
    type(typ_cent):: cent
    integer:: iat, i, j, ng
    real(8):: epot_c, out_ann
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol, fnet(3)
    call f_routine(id='cal_ann_cent2')
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
    call construct_ewald_p3d(parini,atoms,cent%ewald_p3d)
    !call cal_electrostatic_eem2(parini,'init',atoms,ann_arr,epot_c,ann_arr%a)
    if(parini%iverbose>=2) call cpu_time(time2)
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    else
        symfunc%linked_lists%rcut=ann_arr%rcut
        symfunc%linked_lists%triplex=.true.
        call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr_tmp)
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
    call get_qat_from_chi2(parini,ann_arr,atoms,cent)
    if(parini%iverbose>=3) then
        do iat=1,atoms%nat
            write(82,'(i5,1x,a,1x,2f8.3)') iat,trim(atoms%sat(iat)), &
                ann_arr%chi_o(iat),atoms%zat(iat)+atoms%qat(iat)
        enddo
    endif
    !write(*,'(a,4f8.3)') 'TEST ',ann_arr%chi_o(71),ann_arr%chi_o(84), &
    !    atoms%zat(71)+atoms%qat(71),atoms%zat(84)+atoms%qat(84)
    if(parini%iverbose>=2) call cpu_time(time5)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_force_chi_part2(parini,symfunc,atoms,ann_arr)
    endif !end of if for potential
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
    !atoms%epot=epot_c
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
        ann_arr%fchi_norm=tt2/max(tt1,1.d-3)
    endif
    if(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk') then
        call destruct_ewald_p3d(parini,atoms,cent%ewald_p3d)
    endif

    call getvol_alborz(atoms%cellvec,vol)
    call invertmat_alborz(atoms%cellvec,hinv)
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
        !do i=1,ann_arr%n
        !    ekf%g(ekf%loc(i)+ekf%num(1)-1)=ekf%g(ekf%loc(i)+ekf%num(1)-1)*1.d-4
        !    !write(*,*) 'GGG ',ia,ekf%loc(ia)+ekf%num(1)-1
        !enddo
    endif
    if(parini%iverbose>=3) then
        fnet=0.d0
        do iat=1,atoms%nat
            fnet(1)=fnet(1)+atoms%fat(1,iat)
            fnet(2)=fnet(2)+atoms%fat(2,iat)
            fnet(3)=fnet(3)+atoms%fat(3,iat)
        enddo
        write(*,'(a,3es14.5)') 'NET FORCE ',fnet(1),fnet(2),fnet(3)
    endif
    call f_release_routine()
end subroutine cal_ann_cent2
!*****************************************************************************************
subroutine get_qat_from_chi2(parini,ann_arr,atoms,cent)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    !local variables
    integer:: istep, iat
    real(8):: gnrm, epot_old, de, gnrm2, gtot, q1
    real(8):: ttrand(3), qtot
    real(8):: alpha_q, alpha_r, cos_angle_r, cos_angle_q, DDOT, tt1, tt2, tt3
    real(8), allocatable:: qgrad_old(:)
    real(8), allocatable:: rgrad_old(:,:)
    real(8), allocatable:: qat_old(:)
    real(8), allocatable:: rel_old(:,:)
    call init_cent2(parini,ann_arr,atoms,cent)
    allocate(qgrad_old(atoms%nat),rgrad_old(3,atoms%nat))
    allocate(qat_old(atoms%nat),rel_old(3,atoms%nat))
    alpha_q=2.d-1
    alpha_r=2.d-1
    do istep=0,parini%nstep_cep
        call cal_potential_cent2(parini,ann_arr,atoms,cent)
        gnrm=sqrt(sum(cent%rgrad**2))
        gtot=sum(cent%qgrad(1:atoms%nat))
        cent%qgrad(1:atoms%nat)=cent%qgrad(1:atoms%nat)-gtot/atoms%nat
        gnrm2=sqrt(sum(cent%qgrad(1:atoms%nat)**2))
        qtot=sum(atoms%zat(1:atoms%nat))+sum(atoms%qat(1:atoms%nat))
        q1=atoms%zat(1)+atoms%qat(1)
        if(istep==0) then
            epot_old=atoms%epot
            rgrad_old=cent%rgrad
            qgrad_old=cent%qgrad
            rel_old=cent%rel
            qat_old=atoms%qat
        endif
        de=atoms%epot-epot_old
        if(parini%iverbose>=2) then
            write(*,'(a,i5,es24.15,3es11.2,4f8.3)') 'cep: ', &
                istep,atoms%epot,de,gnrm,gnrm2,q1,qtot,alpha_r,alpha_q
            write(51,'(i5,2f8.3)') istep,cent%rel(1,1)-atoms%rat(1,1),cent%rel(1,2)-atoms%rat(1,2)
        endif
        if(gnrm<parini%rgnrmtol .and. gnrm2<parini%qgnrmtol) then
            write(*,'(a,i5,es24.15,2es11.2,2f8.3)') 'CEP converged: ', &
                istep,atoms%epot,gnrm,gnrm2,q1,qtot
            exit
        endif
        if(istep==parini%nstep_cep) then
            write(*,'(a)') 'CEP did not converge, FLAME will stop.'
            stop
        endif
        tt1=DDOT(3*atoms%nat,cent%rgrad,1,rgrad_old,1)
        tt2=DDOT(3*atoms%nat,cent%rgrad,1,cent%rgrad,1)
        tt3=DDOT(3*atoms%nat,rgrad_old,1,rgrad_old,1)
        cos_angle_r=tt1/sqrt(tt2*tt3)
        tt1=DDOT(atoms%nat,cent%qgrad,1,qgrad_old,1)
        tt2=DDOT(atoms%nat,cent%qgrad,1,cent%qgrad,1)
        tt3=DDOT(atoms%nat,qgrad_old,1,qgrad_old,1)
        cos_angle_q=tt1/sqrt(tt2*tt3)
        !write(*,*) 'ANG ',cos_angle_r,cos_angle_q
        if(cos_angle_r>0.5d0 .and. cos_angle_q>0.5d0) then
            alpha_q=alpha_q*1.05d0
            alpha_r=alpha_r*1.05d0
        else
            alpha_q=alpha_q*0.5d0
            alpha_r=alpha_r*0.5d0
        endif
        do iat=1,atoms%nat
            cent%rel(1,iat)=cent%rel(1,iat)-alpha_r*cent%rgrad(1,iat)
            cent%rel(2,iat)=cent%rel(2,iat)-alpha_r*cent%rgrad(2,iat)
            cent%rel(3,iat)=cent%rel(3,iat)-alpha_r*cent%rgrad(3,iat)
            atoms%qat(iat)=atoms%qat(iat)-alpha_q*cent%qgrad(iat)
        enddo
        epot_old=atoms%epot
        rgrad_old=cent%rgrad
        qgrad_old=cent%qgrad
        rel_old=cent%rel
        qat_old=atoms%qat
    enddo
    !write(*,'(a,i5,2f8.3)') 'DISP ',istep,cent%rel(1,1)-atoms%rat(1,1),cent%rel(1,2)-atoms%rat(1,2)
    !write(*,'(a,9f10.2)') 'WHAT ', &
    !                    cent%rel(1:3,70)-atoms%rat(1:3,70), &
    !                    cent%rel(1:3,51)-atoms%rat(1:3,70), &
    !                    atoms%rat(1:3,51)-atoms%rat(1:3,70)

    call cent2_force(parini,ann_arr,atoms,cent)

    call charge_analysis(parini,atoms,ann_arr)
    call final_cent2(cent)
end subroutine get_qat_from_chi2
!*****************************************************************************************
subroutine init_cent2(parini,ann_arr,atoms,cent)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    !local variables
    integer:: iat
    real(8):: qtot, qtot_ion, qtot_ele
    real(8):: ttrand(3)
    real(8):: zion, alpha_elec, alpha_largest, error_erfc_over_r, rcut

    allocate(cent%rgrad(1:3,1:atoms%nat))
    allocate(cent%qgrad(1:atoms%nat))
    allocate(cent%rel(1:3,1:atoms%nat))
    allocate(cent%gwe(1:atoms%nat))
    allocate(cent%gwi(1:atoms%nat))
    allocate(cent%gwit(1:atoms%nat))

    cent%rgrad= 0.d0
    cent%qgrad= 0.d0
    cent%rel= 0.d0
    cent%gwe= 0.d0
    cent%gwi= 0.d0
    cent%gwit= 0.d0

    do iat=1,atoms%nat
        call random_number(ttrand)
        cent%rel(1:3,iat)=atoms%rat(1:3,iat) !+(ttrand(1:3)-0.5d0)*2.d0*1.d-2
    enddo

    cent%ewald_p3d%linked_lists%rcut=parini%rcut_ewald
    !This linked list is used for the short range part of the Ewald.
    call call_linkedlist(parini,atoms,.false.,cent%ewald_p3d%linked_lists,cent%ewald_p3d%pia_arr)

    do iat=1,atoms%nat
        zion=ann_arr%ann(atoms%itypat(iat))%zion
        atoms%zat(iat)=zion
        atoms%qat(iat)=ann_arr%ann(atoms%itypat(iat))%qinit-zion
        !write(*,'(i,3f8.2)') iat,atoms%zat(iat),atoms%qat(iat),atoms%zat(iat)+atoms%qat(iat)
    enddo
    qtot_ion=sum(atoms%zat(1:atoms%nat))
    qtot_ele=sum(atoms%qat(1:atoms%nat))
    qtot=qtot_ion+qtot_ele
    write(*,'(a,3es14.5)') 'Initial total charges: ionic,electronic,net ', &
        qtot_ion,qtot_ele,qtot
    do iat=1,atoms%nat
        cent%gwi(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
        cent%gwe(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth
        cent%gwit(iat)=parini%alpha_ewald
    enddo
    alpha_elec=maxval(cent%gwe(1:atoms%nat))
    alpha_largest=sqrt(alpha_elec**2+parini%alpha_ewald**2)
    alpha_largest=max(alpha_largest,parini%alpha_ewald*sqrt(2.d0))
    rcut=cent%ewald_p3d%linked_lists%rcut
    error_erfc_over_r=erfc(rcut/alpha_largest)/rcut
    write(*,*) 'short range at cut-off: ',error_erfc_over_r !CORRECT_IT
    !cent%ewald_p3d%rgcut=8.d0/0.529d0 !parini%rgcut_ewald*ewald_p3d%alpha !CORRECT_IT
    cent%ewald_p3d%rgcut=sqrt(-log(1.d-8))*alpha_largest
end subroutine init_cent2
!*****************************************************************************************
subroutine final_cent2(cent)
    use mod_interface
    use mod_ann, only: typ_cent
    use dynamic_memory
    implicit none
    type(typ_cent), intent(inout):: cent
    !local variables
    deallocate(cent%ewald_p3d%linked_lists%prime_bound)
    deallocate(cent%ewald_p3d%linked_lists%bound_rad)
    deallocate(cent%ewald_p3d%linked_lists%bound_ang)
    deallocate(cent%ewald_p3d%pia_arr%pia)
    deallocate(cent%rgrad)
    deallocate(cent%qgrad)
    deallocate(cent%rel)
    deallocate(cent%gwi)
    deallocate(cent%gwe)
    deallocate(cent%gwit)
end subroutine final_cent2
!*****************************************************************************************
subroutine cent2_force(parini,ann_arr,atoms,cent)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    !local variables
    integer:: iat
    real(8):: spring_const, dx, dy, dz
    associate(nx=>cent%ewald_p3d%poisson_p3d%ngpx)
    associate(ny=>cent%ewald_p3d%poisson_p3d%ngpy)
    associate(nz=>cent%ewald_p3d%poisson_p3d%ngpz)
    call gauss_force(parini,'bulk',atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,cent%gwit, &
        cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,atoms%fat)
    call cal_shortrange_ewald_force(parini,ann_arr,atoms,cent)
    do iat=1,atoms%nat !summation over ions/electrons
        dx=cent%rel(1,iat)-atoms%rat(1,iat)
        dy=cent%rel(2,iat)-atoms%rat(2,iat)
        dz=cent%rel(3,iat)-atoms%rat(3,iat)
        spring_const=ann_arr%ann(atoms%itypat(iat))%spring_const
        atoms%fat(1,iat)=atoms%fat(1,iat)+spring_const*dx
        atoms%fat(2,iat)=atoms%fat(2,iat)+spring_const*dy
        atoms%fat(3,iat)=atoms%fat(3,iat)+spring_const*dz
    enddo
    end associate
    end associate
    end associate
end subroutine cent2_force
!*****************************************************************************************
subroutine cal_potential_cent2(parini,ann_arr,atoms,cent)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    !local variables
    integer:: iat
    real(8):: epot_es, dx, dy, dz
    real(8):: hardness, spring_const
    if(trim(atoms%boundcond)/='bulk') then
        write(*,*) 'ERROR: CENT2 is ready only for bulk BC.'
        stop
    endif
    ann_arr%ener_ref=0.d0
    do iat=1,atoms%nat
        ann_arr%ener_ref=ann_arr%ener_ref+ann_arr%ann(atoms%itypat(iat))%ener_ref
    enddo
    atoms%epot=0.d0
    cent%qgrad=0.d0
    cent%rgrad=0.d0
    do iat=1,atoms%nat !summation over ions/electrons
        atoms%epot=atoms%epot+ann_arr%chi_o(iat)*(atoms%zat(iat)+atoms%qat(iat))
        hardness=ann_arr%ann(atoms%itypat(iat))%hardness
        atoms%epot=atoms%epot+0.5d0*hardness*(atoms%zat(iat)+atoms%qat(iat))**2
        cent%qgrad(iat)=cent%qgrad(iat)+ann_arr%chi_o(iat)+(atoms%zat(iat)+atoms%qat(iat))*hardness
        dx=cent%rel(1,iat)-atoms%rat(1,iat)
        dy=cent%rel(2,iat)-atoms%rat(2,iat)
        dz=cent%rel(3,iat)-atoms%rat(3,iat)
        spring_const=ann_arr%ann(atoms%itypat(iat))%spring_const
        atoms%epot=atoms%epot+0.5d0*spring_const*(dx**2+dy**2+dz**2)
        cent%rgrad(1,iat)=cent%rgrad(1,iat)+spring_const*dx
        cent%rgrad(2,iat)=cent%rgrad(2,iat)+spring_const*dy
        cent%rgrad(3,iat)=cent%rgrad(3,iat)+spring_const*dz
    enddo
    call cal_pot_with_bps(parini,ann_arr,atoms,cent,epot_es)
    atoms%epot=atoms%epot+epot_es
    atoms%epot=atoms%epot+ann_arr%ener_ref !CORRECT_IT
end subroutine cal_potential_cent2
!*****************************************************************************************
subroutine cal_pot_with_bps(parini,ann_arr,atoms,cent,epot_es)
    use mod_interface
    use mod_ann, only: typ_ann_arr
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
    !local variables
    real(8):: ehartree, pi
    real(8):: time1, time2, time3, time4, time5, time6, time7
    real(8):: sqrt_one_over_twopi
    integer:: iat
    !logical:: ewald
    pi=4.d0*atan(1.d0)
    sqrt_one_over_twopi=1.d0/sqrt(2.d0*pi)
    associate(nx=>cent%ewald_p3d%poisson_p3d%ngpx)
    associate(ny=>cent%ewald_p3d%poisson_p3d%ngpy)
    associate(nz=>cent%ewald_p3d%poisson_p3d%ngpz)
    !-------------------------------------------------------
    atoms%stress=0.d0
    !-------------------------------------------------------
    !open(unit=111,file="tinput",status='old')
    !read(111,*) ewald
    !close(111)
    !if(.not. ewald) then
    !    gw_ion_t=gw_ion
    !endif
    call cpu_time(time1)
    call put_gauss_to_grid(parini,atoms,cent)
    call cpu_time(time2)
    ehartree=0.d0
    call cal_hartree_pot_bps(cent%ewald_p3d,atoms,ehartree)
    epot_es=0.d0
    do iat=1,atoms%nat
        epot_es=epot_es-atoms%zat(iat)**2*sqrt_one_over_twopi/cent%gwit(iat)
    enddo
    epot_es=epot_es+ehartree
    call gauss_gradient(parini,'bulk',atoms%nat,cent%rel,atoms%cellvec,atoms%qat,cent%gwe, &
        cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%pot,cent%rgrad,cent%qgrad)

    !if(ewald) then
    call cal_shortrange_ewald(parini,ann_arr,atoms,cent,epot_es)
    !endif
    !gw_ion_t

    !write(61,'(f20.9)') epot_es
    !do iat=1,atoms%nat
    !    write(61,'(i4,4f10.5)') iat,rgrad(1,iat),rgrad(2,iat),rgrad(3,iat),qgrad(iat)
    !enddo
    !stop 'TESTING EWALD'

    end associate
    end associate
    end associate
end subroutine cal_pot_with_bps
!*****************************************************************************************
subroutine put_gauss_to_grid(parini,atoms,cent)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_cent
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_cent), intent(inout):: cent
    !local variables
    character(10):: bc
    integer:: nx, ny, nz
    nx=cent%ewald_p3d%poisson_p3d%ngpx
    ny=cent%ewald_p3d%poisson_p3d%ngpy
    nz=cent%ewald_p3d%poisson_p3d%ngpz
    bc=trim(atoms%boundcond)
    call gauss_grid(parini,bc,.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%zat, &
        cent%gwit,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
    call gauss_grid(parini,bc,.false.,atoms%nat,cent%rel,atoms%cellvec,atoms%qat, &
        cent%gwe,cent%ewald_p3d%rgcut,nx,ny,nz,cent%ewald_p3d%poisson_p3d%rho)
end subroutine put_gauss_to_grid
!*****************************************************************************************
subroutine gauss_grid(parini,bc,reset,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,rho)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    logical, intent(in):: reset
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: rho(ngx,ngy,ngz)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: cell(3) !dimensions of a smaller orthogonal cell for replication
    real(8):: vol
    real(8):: cvinv(3) !cell vectors of inverse coordinate, actual one at a time
    real(8):: htx, hty, htz
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hxx_g, hxy_g, hxz_g, hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: cvinv_norm
    real(8):: dmx, dmy, dmz, dmsq, gwsq_inv
    real(8):: xred, yred, zred
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: ncellx, ncelly, ncellz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igzysq
    integer:: iii, nlim, nlimsq, iix, iiy, iiz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz, nex, ney, nez
    integer:: finalx, igxs, igxf 
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8):: sqpi, tt
    integer:: istartx, istarty, istartz

    allocate(ratred(3,nat))
    call rxyz_cart2int_alborz(nat,cv,rxyz,ratred)
    do iat=1,nat
        xred=ratred(1,iat)
        yred=ratred(2,iat)
        zred=ratred(3,iat)
        if(xred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(yred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(zred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
        if(.not. (xred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(.not. (yred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(.not. (zred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
    enddo

    !reciprocal lattice to be used to determine the distance of corners of
    !the parallelepiped to its facets. Then those distances are used to
    !determine the number of grid points in each direction that are within
    !the cutoff of Gaussian function.
    call cell_vol(nat,cv,vol)
    vol=abs(vol)*nat
    call cross_product_alborz(cv(1,1),cv(1,2),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(3)=vol/cvinv_norm
    call cross_product_alborz(cv(1,2),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(1)=vol/cvinv_norm
    call cross_product_alborz(cv(1,1),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(2)=vol/cvinv_norm
    if(parini%iverbose>1) then
        write(*,*) 'cell  ', cell(1),cell(2),cell(3)
    endif
    htx=cell(1)/real(ngx,8)
    hty=cell(2)/real(ngy,8)
    htz=cell(3)/real(ngz,8)
    nbgx=int(rgcut/htx)+2
    nbgy=int(rgcut/hty)+2
    nbgz=int(rgcut/htz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    nlim = max(nagx,nagy,nagz)
    nlimsq = nlim**2
    !detemining the largest dimension for the pseudogrid.
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz

    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !if(trim(bc)=='bulk') then
    !    iii=0
    !elseif(trim(bc)=='slab') then
    !    iii=1
    !endif
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    sqpi=sqrt(pi)
    !-------------------------------------------------------
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqpi)**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1
        facqiat=fac*qat(iat)
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmx=ximg-rxyz(1,iat)
                    dmy=yimg-rxyz(2,iat)
                    dmz=zimg-rxyz(3,iat)
                    dmsq=dmx**2+dmy**2+dmz**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
                call vdexp(2*iix+1,exponentval(-iix),expval(-iix))
                wa(imgx-iix:imgx+iix,jgy,jgz)=wa(imgx-iix:imgx+iix,jgy,jgz)+facqiat*expval(-iix:iix)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    if(reset) then
        !if the input array of charge density does not contain any previous value
        !wanted to be preserved.
        rho = 0.d0
    endif
        !if the input array of charge density already some value that must be preserved.
    istartx = modulo((1-nagx-1),ngx)+1
    finalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-finalx+1
    istarty = modulo((-nagy),ngy)+1
    istartz = modulo((-nagz),ngz)+1
    iiz=istartz-1

    do igz=1-nagz,ngz+nagz
        iiy=istarty-1
        iiz=iiz+1
        if (iiz==ngz+1) iiz=1
        do igy=1-nagy,ngy+nagy
            iiy=iiy+1
            if (iiy==ngy+1) iiy=1
            iix=istartx-1
            do igx=1-nagx,igxs-1
                iix=iix+1
                rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(igx,igy,igz)
            enddo
            do igx=igxs,igxf-1,ngx
                do iix=1,ngx
                    jgx=igx+iix-1
                    rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(jgx,igy,igz)
                enddo
            enddo
            iix=0
            do igx=igxf,ngx+nagx
                iix=iix+1
                rho(iix,iiy,iiz)=rho(iix,iiy,iiz)+wa(igx,igy,igz)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    call f_free(wa)
end subroutine gauss_grid
!*****************************************************************************************
subroutine gauss_gradient(parini,bc,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,pot,rgrad,qgrad)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: pot(ngx,ngy,ngz)
    real(8), intent(out):: rgrad(3,nat), qgrad(nat)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: cell(3) !dimensions of a smaller orthogonal cell for replication
    real(8):: vol
    real(8):: cvinv(3) !cell vectors of inverse coordinate, actual one at a time
    real(8):: htx, hty, htz
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hxx_g, hxy_g, hxz_g, hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: cvinv_norm, vol_voxel, ttx, tty, ttz, ttq, tt1
    real(8):: dmsq, gwsq_inv
    real(8):: xred, yred, zred
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: ncellx, ncelly, ncellz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igyt, igzt, igzysq
    integer:: iii, nlim, nlimsq, iix, iiy, iiz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz, nex, ney, nez
    integer:: finalx, igxs, igxf 
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8), allocatable:: dmxarr(:), dmyarr(:), dmzarr(:)
    real(8):: sqpi, tt
    integer:: istartx, istarty, istartz

    allocate(ratred(3,nat))
    call rxyz_cart2int_alborz(nat,cv,rxyz,ratred)
    do iat=1,nat
        xred=ratred(1,iat)
        yred=ratred(2,iat)
        zred=ratred(3,iat)
        if(xred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(yred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(zred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
        if(.not. (xred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(.not. (yred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(.not. (zred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
    enddo

    !reciprocal lattice to be used to determine the distance of corners of
    !the parallelepiped to its facets. Then those distances are used to
    !determine the number of grid points in each direction that are within
    !the cutoff of Gaussian function.
    call cell_vol(nat,cv,vol)
    vol=abs(vol)*nat
    call cross_product_alborz(cv(1,1),cv(1,2),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(3)=vol/cvinv_norm
    call cross_product_alborz(cv(1,2),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(1)=vol/cvinv_norm
    call cross_product_alborz(cv(1,1),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(2)=vol/cvinv_norm
    if(parini%iverbose>1) then
        write(*,*) 'cell  ', cell(1),cell(2),cell(3)
    endif
    htx=cell(1)/real(ngx,8)
    hty=cell(2)/real(ngy,8)
    htz=cell(3)/real(ngz,8)
    nbgx=int(rgcut/htx)+2
    nbgy=int(rgcut/hty)+2
    nbgz=int(rgcut/htz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    nlim = max(nagx,nagy,nagz)
    nlimsq = nlim**2
    !detemining the largest dimension for the pseudogrid.
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz

    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
    allocate(dmxarr(-nbgx:nbgx),dmyarr(-nbgx:nbgx),dmzarr(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !if(trim(bc)=='bulk') then
    !    iii=0
    !elseif(trim(bc)=='slab') then
    !    iii=1
    !endif
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    !---------------------------------------------------------------------------
    istartx = modulo((1-nagx-1),ngx)+1
    finalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-finalx+1
    istarty = modulo((-nagy),ngy)+1
    istartz = modulo((-nagz),ngz)+1
    iiz=istartz-1

    do igz=1-nagz,ngz+nagz
        iiy=istarty-1
        iiz=iiz+1
        if (iiz==ngz+1) iiz=1
        do igy=1-nagy,ngy+nagy
            iiy=iiy+1
            if (iiy==ngy+1) iiy=1
            iix=istartx-1
            do igx=1-nagx,igxs-1
                iix=iix+1
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
            do igx=igxs,igxf-1,ngx
                do iix=1,ngx
                    jgx=igx+iix-1
                    wa(jgx,igy,igz)=pot(iix,iiy,iiz)
                enddo
            enddo
            iix=0
            do igx=igxf,ngx+nagx
                iix=iix+1
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
        enddo
    enddo
    !-------------------------------------------------------
    vol_voxel=vol/(ngx*ngy*ngz)
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqrt(pi))**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1
        facqiat=fac*qat(iat)
        ttq=0.d0
        ttx=0.d0
        tty=0.d0
        ttz=0.d0
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmxarr(igx)=ximg-rxyz(1,iat)
                    dmyarr(igx)=yimg-rxyz(2,iat)
                    dmzarr(igx)=zimg-rxyz(3,iat)
                    dmsq=dmxarr(igx)**2+dmyarr(igx)**2+dmzarr(igx)**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
                call vdexp( 2*iix+1, exponentval(-iix), expval(-iix) )
                do igx=-iix,iix
                    ttq=ttq+fac*expval(igx)*wa(igx+imgx,jgy,jgz)
                    tt1=facqiat*expval(igx)*wa(igx+imgx,jgy,jgz)*(2.d0*gwsq_inv)
                    ttx=ttx+tt1*dmxarr(igx)
                    tty=tty+tt1*dmyarr(igx)
                    ttz=ttz+tt1*dmzarr(igx)
                enddo
            enddo
        enddo
        qgrad(iat)=qgrad(iat)+ttq*vol_voxel
        rgrad(1,iat)=rgrad(1,iat)+ttx*vol_voxel
        rgrad(2,iat)=rgrad(2,iat)+tty*vol_voxel
        rgrad(3,iat)=rgrad(3,iat)+ttz*vol_voxel
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    deallocate(dmxarr,dmyarr,dmzarr)
    call f_free(wa)
end subroutine gauss_gradient
!*****************************************************************************************
subroutine cal_shortrange_ewald(parini,ann_arr,atoms,cent,epot_es)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(in):: atoms
    type(typ_cent), intent(inout):: cent
    real(8), intent(inout):: epot_es
    !local variables
    integer:: iat, jat, ib
    real(8):: alpha, epot_short, gama, dx, dy, dz, r, pi, vol, shift
    real(8):: sqrt_one_over_twopi, ee1, tt1, tt21, tt22, ttg, tt31, tt32
    real(8):: zat_tot
    !type(typ_atoms):: atoms_e
    pi=4.d0*atan(1.d0)
    sqrt_one_over_twopi=1.d0/sqrt(2.d0*pi)
    alpha=parini%alpha_ewald
    epot_short=0.d0
    do ib=1,cent%ewald_p3d%linked_lists%maxbound_rad
        iat=cent%ewald_p3d%linked_lists%bound_rad(1,ib)
        jat=cent%ewald_p3d%linked_lists%bound_rad(2,ib)
        !---------------------------------------------------
        dx=cent%ewald_p3d%pia_arr%pia(ib)%dr(1)
        dy=cent%ewald_p3d%pia_arr%pia(ib)%dr(2)
        dz=cent%ewald_p3d%pia_arr%pia(ib)%dr(3)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwi(jat)**2)
        epot_short=epot_short-atoms%zat(iat)*atoms%zat(jat)*erfc(gama*r)/r
        gama=1.d0/(sqrt(2.d0)*alpha)
        epot_short=epot_short+atoms%zat(iat)*atoms%zat(jat)*erfc(gama*r)/r
        !---------------------------------------------------
        dx=cent%ewald_p3d%pia_arr%pia(ib)%dr(1)+cent%rel(1,jat)-atoms%rat(1,jat)
        dy=cent%ewald_p3d%pia_arr%pia(ib)%dr(2)+cent%rel(2,jat)-atoms%rat(2,jat)
        dz=cent%ewald_p3d%pia_arr%pia(ib)%dr(3)+cent%rel(3,jat)-atoms%rat(3,jat)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jat)**2)
        epot_short=epot_short+atoms%zat(iat)*atoms%qat(jat)*erf(gama*r)/r
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt21=atoms%zat(iat)*atoms%qat(jat)*ee1
        tt31=atoms%zat(iat)*erf(gama*r)/r
        !-------------------------------------------
        gama=1.d0/sqrt(alpha**2+cent%gwe(jat)**2)
        epot_short=epot_short-atoms%zat(iat)*atoms%qat(jat)*erf(gama*r)/r
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt22=-atoms%zat(iat)*atoms%qat(jat)*ee1
        tt32=-atoms%zat(iat)*erf(gama*r)/r
        !-------------------------------------------
        cent%rgrad(1,jat)=cent%rgrad(1,jat)+(tt21+tt22)*dx
        cent%rgrad(2,jat)=cent%rgrad(2,jat)+(tt21+tt22)*dy
        cent%rgrad(3,jat)=cent%rgrad(3,jat)+(tt21+tt22)*dz
        cent%qgrad(jat)=cent%qgrad(jat)+tt31+tt32
        !---------------------------------------------------
        dx=-cent%ewald_p3d%pia_arr%pia(ib)%dr(1)+cent%rel(1,iat)-atoms%rat(1,iat)
        dy=-cent%ewald_p3d%pia_arr%pia(ib)%dr(2)+cent%rel(2,iat)-atoms%rat(2,iat)
        dz=-cent%ewald_p3d%pia_arr%pia(ib)%dr(3)+cent%rel(3,iat)-atoms%rat(3,iat)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gama=1.d0/sqrt(cent%gwi(jat)**2+cent%gwe(iat)**2)
        epot_short=epot_short+atoms%zat(jat)*atoms%qat(iat)*erf(gama*r)/r
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt21=atoms%zat(jat)*atoms%qat(iat)*ee1
        tt31=atoms%zat(jat)*erf(gama*r)/r
        !-------------------------------------------
        gama=1.d0/sqrt(alpha**2+cent%gwe(iat)**2)
        epot_short=epot_short-atoms%zat(jat)*atoms%qat(iat)*erf(gama*r)/r
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt22=-atoms%zat(jat)*atoms%qat(iat)*ee1
        tt32=-atoms%zat(jat)*erf(gama*r)/r
        !-------------------------------------------
        cent%rgrad(1,iat)=cent%rgrad(1,iat)+(tt21+tt22)*dx
        cent%rgrad(2,iat)=cent%rgrad(2,iat)+(tt21+tt22)*dy
        cent%rgrad(3,iat)=cent%rgrad(3,iat)+(tt21+tt22)*dz
        cent%qgrad(iat)=cent%qgrad(iat)+tt31+tt32
        !---------------------------------------------------
    enddo
    zat_tot=sum(atoms%zat(1:atoms%nat))
    call getvol_alborz(atoms%cellvec,vol)
    shift=0.d0
    do iat=1,atoms%nat
        shift=shift-(alpha**2-cent%gwi(iat)**2)*atoms%zat(iat)
    enddo
    shift=shift*pi/vol
    do iat=1,atoms%nat
        jat=iat
        dx=cent%rel(1,jat)-atoms%rat(1,iat)
        dy=cent%rel(2,jat)-atoms%rat(2,iat)
        dz=cent%rel(3,jat)-atoms%rat(3,iat)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        if(r>0.3d0) then
            write(*,'(a,es14.5,i6,1x,a)') 'ERROR: Center of electron far from atom: r= ', &
                r,iat,trim(atoms%sat(iat))
            stop
        endif
        if(r<0.1d0) then
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jat)**2)
            call erf_over_r_taylor(gama*r,tt1,ttg)
            epot_short=epot_short+atoms%zat(iat)*atoms%qat(jat)*(tt1*gama)
            tt21=atoms%zat(iat)*atoms%qat(jat)*gama**3*ttg
            tt31=atoms%zat(iat)*(tt1*gama)
            !-------------------------------------------
            gama=1.d0/sqrt(alpha**2+cent%gwe(jat)**2)
            call erf_over_r_taylor(gama*r,tt1,ttg)
            epot_short=epot_short-atoms%zat(iat)*atoms%qat(jat)*(tt1*gama)
            tt22=-atoms%zat(iat)*atoms%qat(jat)*gama**3*ttg
            tt32=-atoms%zat(iat)*(tt1*gama)
        else
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jat)**2)
            epot_short=epot_short+atoms%zat(iat)*atoms%qat(jat)*erf(gama*r)/r
            ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
            tt21=atoms%zat(iat)*atoms%qat(jat)*ee1
            tt31=atoms%zat(iat)*erf(gama*r)/r
            !-------------------------------------------
            gama=1.d0/sqrt(alpha**2+cent%gwe(jat)**2)
            epot_short=epot_short-atoms%zat(iat)*atoms%qat(jat)*erf(gama*r)/r
            ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
            tt22=-atoms%zat(iat)*atoms%qat(jat)*ee1
            tt32=-atoms%zat(iat)*erf(gama*r)/r
        endif
        !-------------------------------------------
        cent%rgrad(1,jat)=cent%rgrad(1,jat)+(tt21+tt22)*dx
        cent%rgrad(2,jat)=cent%rgrad(2,jat)+(tt21+tt22)*dy
        cent%rgrad(3,jat)=cent%rgrad(3,jat)+(tt21+tt22)*dz
        cent%qgrad(jat)=cent%qgrad(jat)+tt31+tt32
        cent%qgrad(jat)=cent%qgrad(jat)+shift
    enddo
    epot_es=epot_es+epot_short
end subroutine cal_shortrange_ewald
!*****************************************************************************************
subroutine gauss_force(parini,bc,nat,rxyz,cv,qat,gw,rgcut,ngx,ngy,ngz,pot,fat)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    character(*), intent(in):: bc
    integer, intent(in):: nat
    real(8), intent(in):: rxyz(3,nat)
    real(8), intent(in):: cv(3,3)
    real(8), intent(in):: qat(nat)
    real(8), intent(in):: gw(nat)
    real(8), intent(in):: rgcut
    integer, intent(in):: ngx, ngy, ngz
    real(8), intent(inout):: pot(ngx,ngy,ngz)
    real(8), intent(out):: fat(3,nat)
    !local variables
    !work arrays to save the values of one dimensional gaussian function.
    real(8):: pi
    real(8):: facqiat, fac
    real(8):: cell(3) !dimensions of a smaller orthogonal cell for replication
    real(8):: vol
    real(8):: cvinv(3) !cell vectors of inverse coordinate, actual one at a time
    real(8):: htx, hty, htz
    real(8):: hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz
    real(8):: hxx_g, hxy_g, hxz_g, hyx_g, hyy_g, hyz_g, hzx_g, hzy_g, hzz_g
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: cvinv_norm, vol_voxel, ttx, tty, ttz, ttq, tt1
    real(8):: dmsq, gwsq_inv
    real(8):: xred, yred, zred
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: ncellx, ncelly, ncellz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz, igyt, igzt, igzysq
    integer:: iii, nlim, nlimsq, iix, iiy, iiz
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz, nex, ney, nez
    integer:: finalx, igxs, igxf 
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: ratred(:,:)
    real(8), allocatable:: exponentval(:), expval(:)
    real(8), allocatable:: dmxarr(:), dmyarr(:), dmzarr(:)
    real(8):: sqpi, tt
    integer:: istartx, istarty, istartz

    allocate(ratred(3,nat))
    call rxyz_cart2int_alborz(nat,cv,rxyz,ratred)
    do iat=1,nat
        xred=ratred(1,iat)
        yred=ratred(2,iat)
        zred=ratred(3,iat)
        if(xred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(yred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(zred<0.d0) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
        if(.not. (xred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sx ',iat,xred
        if(.not. (yred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sy ',iat,yred
        if(.not. (zred<1.d0)) write(*,*) 'ATOM OUTSIDE: iat,sz ',iat,zred
    enddo

    !reciprocal lattice to be used to determine the distance of corners of
    !the parallelepiped to its facets. Then those distances are used to
    !determine the number of grid points in each direction that are within
    !the cutoff of Gaussian function.
    call cell_vol(nat,cv,vol)
    vol=abs(vol)*nat
    call cross_product_alborz(cv(1,1),cv(1,2),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(3)=vol/cvinv_norm
    call cross_product_alborz(cv(1,2),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(1)=vol/cvinv_norm
    call cross_product_alborz(cv(1,1),cv(1,3),cvinv)
    cvinv_norm=sqrt(cvinv(1)**2+cvinv(2)**2+cvinv(3)**2)
    cell(2)=vol/cvinv_norm
    if(parini%iverbose>1) then
        write(*,*) 'cell  ', cell(1),cell(2),cell(3)
    endif
    htx=cell(1)/real(ngx,8)
    hty=cell(2)/real(ngy,8)
    htz=cell(3)/real(ngz,8)
    nbgx=int(rgcut/htx)+2
    nbgy=int(rgcut/hty)+2
    nbgz=int(rgcut/htz)+2
    nagx=nbgx+1
    nagy=nbgy+1
    nagz=nbgz+1
    nlim = max(nagx,nagy,nagz)
    nlimsq = nlim**2
    !detemining the largest dimension for the pseudogrid.
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz

    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    allocate(exponentval(-nbgx:nbgx),expval(-nbgx:nbgx))
    allocate(dmxarr(-nbgx:nbgx),dmyarr(-nbgx:nbgx),dmzarr(-nbgx:nbgx))
   !  allocate(wa(1-nagx:ngx+nagx,1-nagy:ngy+nagy,1-nagz:ngz+nagz))

    hrxinv=real(ngx,8) !inverse of grid spacing in reduced coordinates
    hryinv=real(ngy,8) !inverse of grid spacing in reduced coordinates
    hrzinv=real(ngz,8) !inverse of grid spacing in reduced coordinates
    !if(trim(bc)=='bulk') then
    !    iii=0
    !elseif(trim(bc)=='slab') then
    !    iii=1
    !endif
    !-------------------------------------------------------
    pi=4.d0*atan(1.d0)
    !---------------------------------------------------------------------------
    istartx = modulo((1-nagx-1),ngx)+1
    finalx = modulo((ngx+nagx-1),ngx)+1
    igxs = 1-nagx+(ngx-istartx+1)
    igxf = ngx+nagx-finalx+1
    istarty = modulo((-nagy),ngy)+1
    istartz = modulo((-nagz),ngz)+1
    iiz=istartz-1

    do igz=1-nagz,ngz+nagz
        iiy=istarty-1
        iiz=iiz+1
        if (iiz==ngz+1) iiz=1
        do igy=1-nagy,ngy+nagy
            iiy=iiy+1
            if (iiy==ngy+1) iiy=1
            iix=istartx-1
            do igx=1-nagx,igxs-1
                iix=iix+1
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
            do igx=igxs,igxf-1,ngx
                do iix=1,ngx
                    jgx=igx+iix-1
                    wa(jgx,igy,igz)=pot(iix,iiy,iiz)
                enddo
            enddo
            iix=0
            do igx=igxf,ngx+nagx
                iix=iix+1
                wa(igx,igy,igz)=pot(iix,iiy,iiz)
            enddo
        enddo
    enddo
    !-------------------------------------------------------
    vol_voxel=vol/(ngx*ngy*ngz)
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqrt(pi))**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1
        facqiat=fac*qat(iat)
        ttx=0.d0
        tty=0.d0
        ttz=0.d0
        do igz=-nbgz,nbgz
            jgz=igz+imgz
            hzx_g=(jgz-1)*hzx
            hzy_g=(jgz-1)*hzy
            hzz_g=(jgz-1)*hzz
            do igy=-nbgy,nbgy
                igzysq=igz**2+igy**2
                if(igzysq>nlimsq) cycle
                jgy=igy+imgy
                hyx_g=(jgy-1)*hyx
                hyy_g=(jgy-1)*hyy
                hyz_g=(jgy-1)*hyz
                tt=nlimsq-igzysq
                iix=min(floor(sqrt(tt)),nbgx)
                do igx=-iix,iix
                    jgx=igx+imgx
                    ximg=(jgx-1)*hxx+hyx_g+hzx_g    
                    yimg=(jgx-1)*hxy+hyy_g+hzy_g
                    zimg=(jgx-1)*hxz+hyz_g+hzz_g
                    dmxarr(igx)=ximg-rxyz(1,iat)
                    dmyarr(igx)=yimg-rxyz(2,iat)
                    dmzarr(igx)=zimg-rxyz(3,iat)
                    dmsq=dmxarr(igx)**2+dmyarr(igx)**2+dmzarr(igx)**2
                    exponentval(igx)=-dmsq*gwsq_inv
                enddo
                call vdexp( 2*iix+1, exponentval(-iix), expval(-iix) )
                do igx=-iix,iix
                    tt1=facqiat*expval(igx)*wa(igx+imgx,jgy,jgz)*(2.d0*gwsq_inv)
                    ttx=ttx+tt1*dmxarr(igx)
                    tty=tty+tt1*dmyarr(igx)
                    ttz=ttz+tt1*dmzarr(igx)
                enddo
            enddo
        enddo
        fat(1,iat)=fat(1,iat)-ttx*vol_voxel
        fat(2,iat)=fat(2,iat)-tty*vol_voxel
        fat(3,iat)=fat(3,iat)-ttz*vol_voxel
    enddo
    !---------------------------------------------------------------------------
    deallocate(ratred)
    deallocate(exponentval,expval)
    deallocate(dmxarr,dmyarr,dmzarr)
    call f_free(wa)
end subroutine gauss_force
!*****************************************************************************************
subroutine cal_shortrange_ewald_force(parini,ann_arr,atoms,cent)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    !local variables
    integer:: iat, jat, ib
    real(8):: alpha, gama, dx, dy, dz, r, pi, vol, shift
    real(8):: sqrt_one_over_twopi, ee1, tt1, tt21, tt22, ttg, tt31, tt32
    real(8):: zat_tot
    !type(typ_atoms):: atoms_e
    pi=4.d0*atan(1.d0)
    sqrt_one_over_twopi=1.d0/sqrt(2.d0*pi)
    alpha=parini%alpha_ewald
    do ib=1,cent%ewald_p3d%linked_lists%maxbound_rad
        iat=cent%ewald_p3d%linked_lists%bound_rad(1,ib)
        jat=cent%ewald_p3d%linked_lists%bound_rad(2,ib)
        !---------------------------------------------------
        dx=cent%ewald_p3d%pia_arr%pia(ib)%dr(1)
        dy=cent%ewald_p3d%pia_arr%pia(ib)%dr(2)
        dz=cent%ewald_p3d%pia_arr%pia(ib)%dr(3)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwi(jat)**2)
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt21=atoms%zat(iat)*atoms%zat(jat)*ee1
        !---------------------------------------------------
        gama=1.d0/(sqrt(2.d0)*alpha)
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt22=-atoms%zat(iat)*atoms%zat(jat)*ee1
        !---------------------------------------------------
        atoms%fat(1,jat)=atoms%fat(1,jat)-(tt21+tt22)*dx
        atoms%fat(2,jat)=atoms%fat(2,jat)-(tt21+tt22)*dy
        atoms%fat(3,jat)=atoms%fat(3,jat)-(tt21+tt22)*dz
        atoms%fat(1,iat)=atoms%fat(1,iat)+(tt21+tt22)*dx
        atoms%fat(2,iat)=atoms%fat(2,iat)+(tt21+tt22)*dy
        atoms%fat(3,iat)=atoms%fat(3,iat)+(tt21+tt22)*dz
        !---------------------------------------------------
        dx=cent%ewald_p3d%pia_arr%pia(ib)%dr(1)+cent%rel(1,jat)-atoms%rat(1,jat)
        dy=cent%ewald_p3d%pia_arr%pia(ib)%dr(2)+cent%rel(2,jat)-atoms%rat(2,jat)
        dz=cent%ewald_p3d%pia_arr%pia(ib)%dr(3)+cent%rel(3,jat)-atoms%rat(3,jat)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jat)**2)
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt21=atoms%zat(iat)*atoms%qat(jat)*ee1
        !-------------------------------------------
        gama=1.d0/sqrt(alpha**2+cent%gwe(jat)**2)
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt22=-atoms%zat(iat)*atoms%qat(jat)*ee1
        !-------------------------------------------
        atoms%fat(1,iat)=atoms%fat(1,iat)+(tt21+tt22)*dx
        atoms%fat(2,iat)=atoms%fat(2,iat)+(tt21+tt22)*dy
        atoms%fat(3,iat)=atoms%fat(3,iat)+(tt21+tt22)*dz
        !---------------------------------------------------
        dx=-cent%ewald_p3d%pia_arr%pia(ib)%dr(1)+cent%rel(1,iat)-atoms%rat(1,iat)
        dy=-cent%ewald_p3d%pia_arr%pia(ib)%dr(2)+cent%rel(2,iat)-atoms%rat(2,iat)
        dz=-cent%ewald_p3d%pia_arr%pia(ib)%dr(3)+cent%rel(3,iat)-atoms%rat(3,iat)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gama=1.d0/sqrt(cent%gwi(jat)**2+cent%gwe(iat)**2)
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt21=atoms%zat(jat)*atoms%qat(iat)*ee1
        !-------------------------------------------
        gama=1.d0/sqrt(alpha**2+cent%gwe(iat)**2)
        ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
        tt22=-atoms%zat(jat)*atoms%qat(iat)*ee1
        !-------------------------------------------
        atoms%fat(1,jat)=atoms%fat(1,jat)+(tt21+tt22)*dx
        atoms%fat(2,jat)=atoms%fat(2,jat)+(tt21+tt22)*dy
        atoms%fat(3,jat)=atoms%fat(3,jat)+(tt21+tt22)*dz
        !---------------------------------------------------
    enddo
    do iat=1,atoms%nat
        jat=iat
        dx=cent%rel(1,jat)-atoms%rat(1,iat)
        dy=cent%rel(2,jat)-atoms%rat(2,iat)
        dz=cent%rel(3,jat)-atoms%rat(3,iat)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        if(r>0.3d0) then
            write(*,'(a,es14.5,i6,1x,a)') 'ERROR: Center of electron far from atom: r= ', &
                r,iat,trim(atoms%sat(iat))
            stop
        endif
        if(r<0.1d0) then
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jat)**2)
            call erf_over_r_taylor(gama*r,tt1,ttg)
            tt21=atoms%zat(iat)*atoms%qat(jat)*gama**3*ttg
            !-------------------------------------------
            gama=1.d0/sqrt(alpha**2+cent%gwe(jat)**2)
            call erf_over_r_taylor(gama*r,tt1,ttg)
            tt22=-atoms%zat(iat)*atoms%qat(jat)*gama**3*ttg
        else
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jat)**2)
            ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
            tt21=atoms%zat(iat)*atoms%qat(jat)*ee1
            !-------------------------------------------
            gama=1.d0/sqrt(alpha**2+cent%gwe(jat)**2)
            ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
            tt22=-atoms%zat(iat)*atoms%qat(jat)*ee1
        endif
        !-------------------------------------------
        atoms%fat(1,jat)=atoms%fat(1,jat)+(tt21+tt22)*dx
        atoms%fat(2,jat)=atoms%fat(2,jat)+(tt21+tt22)*dy
        atoms%fat(3,jat)=atoms%fat(3,jat)+(tt21+tt22)*dz
    enddo
end subroutine cal_shortrange_ewald_force
!*****************************************************************************************
subroutine erf_over_r_taylor(r,funcval,funcval_der)
    implicit none
    real(8), intent(in):: r
    real(8), intent(out):: funcval, funcval_der
    !local variables
    real(8):: pi, rsq, tt1
    real(8):: a0, a1, a2, a3, a4, a5, a6, a7, a8
    pi=4.d0*atan(1.d0)
    tt1=2.d0/sqrt(pi)
    a0=tt1*( 1.d0          )
    a1=tt1*(-1.d0/3.d0     )
    a2=tt1*( 1.d0/10.d0    )
    a3=tt1*(-1.d0/42.d0    )
    a4=tt1*( 1.d0/216.d0   )
    a5=tt1*(-1.d0/1320.d0  )
    a6=tt1*( 1.d0/9360.d0  )
    a7=tt1*(-1.d0/75600.d0 )
    a8=tt1*( 1.d0/685440.d0)
    rsq=r**2
    funcval=a0+rsq*(a1+rsq*(a2+rsq*(a3+rsq*(a4+rsq*(a5+rsq*(a6+rsq*(a7+rsq*a8)))))))
    funcval_der=(2.d0*a1+rsq*(4.d0*a2+rsq*(6.d0*a3+rsq*(8.d0*a4+rsq* &
                (1.d1*a5+rsq*(1.2d1*a6+rsq*(1.4d1*a7+rsq*a8*1.6d1)))))))
end subroutine erf_over_r_taylor
!*****************************************************************************************
