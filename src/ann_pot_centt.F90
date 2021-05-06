!*****************************************************************************************
subroutine cal_ann_centt(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr, typ_cent, convert_ann_epotd
    use mod_symfunc, only: typ_symfunc
    use mod_linked_lists, only: typ_pia_arr
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    !local variables
    type(typ_pia_arr):: pia_arr_tmp
    type(typ_cent):: cent
    integer:: iat, i, j, ng
    real(8):: epot_c, out_ann
    real(8):: time1, time2, time3, time4, time5, time6, time7, time8
    real(8):: tt1, tt2, tt3, fx_es, fy_es, fz_es, hinv(3,3), vol, fnet(3)
    real(8),allocatable :: gausswidth(:)
    call f_routine(id='cal_ann_centt')
    call update_ratp(atoms)
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
    allocate(gausswidth(atoms%nat))
    gausswidth(:)=1.d0 !TO_BE_CORRECTED
    if(trim(parini%psolver)/='pairsum') then
        cent%poisson%task_finit="alloc_rho:set_ngp"
        call init_hartree(parini,atoms,cent%poisson,gausswidth)
    endif
    deallocate(gausswidth)
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
            call convert_ann_epotd(ann_arr%ann(i),ann_arr%num(i),ann_arr%g_per_atom(1,iat))
            ann_arr%g_per_atom(1:ann_arr%num(1),iat)=ann_arr%g_per_atom(1:ann_arr%num(1),iat)*ann_arr%ann(i)%ampl_chi*ann_arr%ann(i)%prefactor_chi*(1.d0-tt1**2)
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
    call get_qat_from_chi_centt(parini,ann_arr,atoms,cent)
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
        !write(*,'(a,f8.3)') 'Timing:centt: initialize matrix          ',time2-time1
        !write(*,'(a,f8.3)') 'Timing:centt: calculation of symfunc     ',time3-time2
        !write(*,'(a,f8.3)') 'Timing:centt: neural network process     ',time4-time3
        !write(*,'(a,f8.3)') 'Timing:centt: linear equations solver    ',time5-time4
        !write(*,'(a,f8.3)') 'Timing:centt: force (SR term)            ',time6-time5
        !write(*,'(a,f8.3)') 'Timing:centt: energy (SR+LR), force (LR) ',time7-time6
        !write(*,'(a,f8.3)') 'Timing:centt: total time                 ',time7-time1
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
    !if(trim(atoms%boundcond)=='slab' .or. trim(atoms%boundcond)=='bulk') then
    if(trim(parini%psolver)/='pairsum') then
        call fini_hartree(parini,atoms,cent%poisson)
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
end subroutine cal_ann_centt
!*****************************************************************************************
subroutine get_qat_from_chi_centt(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    use yaml_output
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
    call init_centt(parini,ann_arr,atoms,cent)
    allocate(qgrad_old(atoms%nat),rgrad_old(3,atoms%nat))
    allocate(qat_old(atoms%nat),rel_old(3,atoms%nat))
    alpha_q=2.d-1*parini%alphax_q
    alpha_r=parini%alphax_r
    call yaml_sequence_open('Charge equilibration process')
    do istep=0,parini%nstep_cep
        if(parini%iverbose>=2) then
            call yaml_sequence(advance='no')
        endif
        call cal_potential_centt(parini,ann_arr,atoms,cent)
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
            !write(*,'(a,i5,es24.15,3es11.2,4f8.3)') 'cep: ', &
            !    istep,atoms%epot,de,gnrm,gnrm2,q1,qtot,alpha_r,alpha_q
            !write(51,'(i5,2f8.3)') istep,cent%rel(3,1)-atoms%rat(3,1),cent%rel(3,2)-atoms%rat(3,2)
            call yaml_mapping_open('cep',flow=.true.)
            call yaml_map('iter',istep,fmt='(i5)')
            call yaml_map('epot',atoms%epot,fmt='(es22.13)')
            call yaml_map('de',de,fmt='(es10.2)')
            call yaml_map('gnrm',gnrm,fmt='(es10.3)')
            call yaml_map('gnrm2',gnrm2,fmt='(es10.3)')
            call yaml_map('qtot',qtot,fmt='(es8.0)')
            call yaml_map('alpha_q',alpha_q,fmt='(es9.2)')
            call yaml_map('alpha_r',alpha_r,fmt='(es9.2)')
            !call yaml_map('stepsize',alpha/alphax,fmt='(es12.3)')
            !call yaml_map('qtot',qtot,fmt='(f10.3)')
            !call yaml_map('dpx',dipole(1),fmt='(f10.3)')
            !call yaml_map('dpy',dipole(2),fmt='(f10.3)')
            !call yaml_map('dpz',dipole(3),fmt='(f10.3)')
            call yaml_mapping_close()
        endif
        if(gnrm<parini%rgnrmtol .and. gnrm2<parini%qgnrmtol) then
            !write(*,'(a,i5,es24.15,2es11.2,2f8.3)') 'CEP converged: ', &
            !    istep,atoms%epot,gnrm,gnrm2,q1,qtot
            call yaml_sequence_close()
            call yaml_mapping_open('CEP',flow=.true.)
            call yaml_map('iter',istep,fmt='(i5)')
            call yaml_map('epot',atoms%epot,fmt='(es22.13)')
            call yaml_map('de',de,fmt='(es11.2)')
            call yaml_map('gnrm',gnrm,fmt='(es12.3)')
            call yaml_map('gnrm2',gnrm2,fmt='(es12.3)')
            call yaml_map('qtot',qtot,fmt='(es12.3)')
            call yaml_mapping_close()
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

    call centt_force(parini,ann_arr,atoms,cent)

    call charge_analysis(parini,atoms,ann_arr)
    call final_centt(cent)
    if(parini%iverbose>=2) then
        call yaml_map('charge on atoms',atoms%zat(1:atoms%nat)+atoms%qat(1:atoms%nat),fmt='(f10.5)')
    endif
end subroutine get_qat_from_chi_centt
!*****************************************************************************************
subroutine init_centt(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_cent
    use mod_atoms, only: typ_atoms
    use dynamic_memory
    use yaml_output
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
        !if(trim(parini%rng_type)=='only_for_tests') then
        !    call random_number_generator_simple(3,ttrand)
        !else
        !    call random_number(ttrand)
        !endif
        cent%rel(1:3,iat)=atoms%ratp(1:3,iat) !+(ttrand(1:3)-0.5d0)*2.d0*1.d-2
    enddo

    cent%poisson%linked_lists%rcut=parini%rcut_ewald
    !This linked list is used for the short range part of the Ewald.
    call call_linkedlist(parini,atoms,.false.,cent%poisson%linked_lists,cent%poisson%pia_arr)
    qtot=0.d0
    do iat=1,atoms%nat
        zion=ann_arr%ann(atoms%itypat(iat))%zion
        atoms%zat(iat)=zion
        atoms%qat(iat)=ann_arr%ann(atoms%itypat(iat))%qinit-zion
        qtot=qtot+atoms%qat(iat)+zion
        !write(*,'(i,3f8.2)') iat,atoms%zat(iat),atoms%qat(iat),atoms%zat(iat)+atoms%qat(iat)
    enddo
    do iat=1,atoms%nat
        atoms%qat(iat)=atoms%qat(iat)+(atoms%qtot-qtot)/real(atoms%nat,8)
    enddo
    qtot_ion=sum(atoms%zat(1:atoms%nat))
    qtot_ele=sum(atoms%qat(1:atoms%nat))
    qtot=qtot_ion+qtot_ele
    !write(*,'(a,3es14.5)') 'Initial total charges: ionic,electronic,net ', &
    !    qtot_ion,qtot_ele,qtot
    call yaml_mapping_open('Initial total charges',flow=.true.)
    call yaml_map('ionic',qtot_ion,fmt='(es9.1)')
    call yaml_map('electronic',qtot_ele,fmt='(es9.1)')
    call yaml_map('net',qtot,fmt='(es9.1)')
    call yaml_mapping_close()
    do iat=1,atoms%nat
        cent%gwi(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
        cent%gwe(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth
        cent%gwit(iat)=parini%alpha_ewald
    enddo
    alpha_elec=maxval(cent%gwe(1:atoms%nat))
    alpha_largest=sqrt(alpha_elec**2+parini%alpha_ewald**2)
    alpha_largest=max(alpha_largest,parini%alpha_ewald*sqrt(2.d0))
    rcut=cent%poisson%linked_lists%rcut
    error_erfc_over_r=erfc(rcut/alpha_largest)/rcut
    !write(*,*) 'short range at cut-off: ',error_erfc_over_r !CORRECT_IT
    call yaml_map('short range at cut-off',error_erfc_over_r) !CORRECT_IT
    !cent%poisson%rgcut=8.d0/0.529d0 !parini%rgcut_ewald*poisson%alpha !CORRECT_IT
    cent%poisson%rgcut=sqrt(-log(1.d-8))*alpha_largest
end subroutine init_centt
!*****************************************************************************************
subroutine final_centt(cent)
    use mod_ann, only: typ_cent
    use dynamic_memory
    implicit none
    type(typ_cent), intent(inout):: cent
    !local variables
    deallocate(cent%poisson%linked_lists%prime_bound)
    deallocate(cent%poisson%linked_lists%bound_rad)
    deallocate(cent%poisson%linked_lists%bound_ang)
    deallocate(cent%poisson%pia_arr%pia)
    deallocate(cent%rgrad)
    deallocate(cent%qgrad)
    deallocate(cent%rel)
    deallocate(cent%gwi)
    deallocate(cent%gwe)
    deallocate(cent%gwit)
end subroutine final_centt
!*****************************************************************************************
subroutine centt_force(parini,ann_arr,atoms,cent)
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
    if(trim(parini%psolver)/='pairsum') then
        cent%poisson%nat=atoms%nat
        cent%poisson%cv=atoms%cellvec
        cent%poisson%bc=atoms%boundcond
        cent%poisson%q(1:cent%poisson%nat)=atoms%zat(1:atoms%nat)
        cent%poisson%gw(1:cent%poisson%nat)=cent%gwit(1:atoms%nat)
        cent%poisson%rcart(1:3,1:cent%poisson%nat)=atoms%ratp(1:3,1:atoms%nat)
        call get_hartree_force(parini,cent%poisson,atoms)
        call cal_shortrange_ewald_force_centt(parini,ann_arr,atoms,cent)
    else
        call cal_centt_pairsum_force(parini,ann_arr,atoms,cent)
    endif
    do iat=1,atoms%nat !summation over ions/electrons
        dx=cent%rel(1,iat)-atoms%ratp(1,iat)
        dy=cent%rel(2,iat)-atoms%ratp(2,iat)
        dz=cent%rel(3,iat)-atoms%ratp(3,iat)
        spring_const=ann_arr%ann(atoms%itypat(iat))%spring_const
        atoms%fat(1,iat)=atoms%fat(1,iat)+spring_const*dx
        atoms%fat(2,iat)=atoms%fat(2,iat)+spring_const*dy
        atoms%fat(3,iat)=atoms%fat(3,iat)+spring_const*dz
    enddo
end subroutine centt_force
!*****************************************************************************************
subroutine cal_potential_centt(parini,ann_arr,atoms,cent)
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
    !if(trim(atoms%boundcond)/='bulk') then
    !    write(*,*) 'ERROR: CENTT is ready only for bulk BC.'
    !    stop
    !endif
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
        dx=cent%rel(1,iat)-atoms%ratp(1,iat)
        dy=cent%rel(2,iat)-atoms%ratp(2,iat)
        dz=cent%rel(3,iat)-atoms%ratp(3,iat)
        spring_const=ann_arr%ann(atoms%itypat(iat))%spring_const
        atoms%epot=atoms%epot+0.5d0*spring_const*(dx**2+dy**2+dz**2)
        cent%rgrad(1,iat)=cent%rgrad(1,iat)+spring_const*dx
        cent%rgrad(2,iat)=cent%rgrad(2,iat)+spring_const*dy
        cent%rgrad(3,iat)=cent%rgrad(3,iat)+spring_const*dz
    enddo
    if(trim(parini%psolver)=='pairsum') then
        call cal_centt_pot_pairsum(parini,ann_arr,atoms,cent,epot_es)
    else
        !write(*,*) 'AAAAAAAAAAAAA'
        !stop
        call cal_centt_pot_bps(parini,ann_arr,atoms,cent,epot_es)
    endif
    atoms%epot=atoms%epot+epot_es
    atoms%epot=atoms%epot+ann_arr%ener_ref !CORRECT_IT
end subroutine cal_potential_centt
!*****************************************************************************************
subroutine cal_centt_pot_pairsum(parini,ann_arr,atoms,cent,epot_es)
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
    integer:: iat, jat, iel, jel
    real(8):: dx, dy, dz, r !, tt1, tt2, tt3, tt4, tt5, tt6 ,tt7, ttf
    real(8):: pi, gama
    real(8):: sqrt_one_over_twopi, tt1, tt2, ttg, ee1, tt3
    pi=4.d0*atan(1.d0)
    sqrt_one_over_twopi=1.d0/sqrt(2.d0*pi)
    !fat=0.d0
    !write(*,*) cent%gwi(:)
    !write(*,*) cent%gwe(:)
    !stop 'BBB'
    epot_es=0.d0
    !-------------------------------------------------------
    !electronic self-interaction
    do iat=1,atoms%nat !summation over electrons
        epot_es=epot_es              +atoms%qat(iat)**2*sqrt_one_over_twopi/cent%gwe(iat)
        cent%qgrad(iat)=cent%qgrad(iat)+2.d0*atoms%qat(iat)*sqrt_one_over_twopi/cent%gwe(iat)
    enddo
    !-------------------------------------------------------
    !electron-electron interaction
    do iel=1,atoms%nat !summation over electrons
        do jel=iel+1,atoms%nat !summation over electrons
            dx=cent%rel(1,jel)-cent%rel(1,iel)
            dy=cent%rel(2,jel)-cent%rel(2,iel)
            dz=cent%rel(3,jel)-cent%rel(3,iel)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            gama=1.d0/sqrt(cent%gwe(iel)**2+cent%gwe(jel)**2)
            epot_es=epot_es              +atoms%qat(iel)*atoms%qat(jel)*erf(gama*r)/r
            ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
            tt2=atoms%qat(iel)*atoms%qat(jel)*ee1
            cent%rgrad(1,iel)=cent%rgrad(1,iel)-tt2*dx
            cent%rgrad(2,iel)=cent%rgrad(2,iel)-tt2*dy
            cent%rgrad(3,iel)=cent%rgrad(3,iel)-tt2*dz
            cent%rgrad(1,jel)=cent%rgrad(1,jel)+tt2*dx
            cent%rgrad(2,jel)=cent%rgrad(2,jel)+tt2*dy
            cent%rgrad(3,jel)=cent%rgrad(3,jel)+tt2*dz
            !write(*,'(a,2i5,4f10.5)') 'EE',iel,jel,tt2,dx,dy,dz
            !write(*,'(a,2i5,6f10.5)') 'EE',iel,jel,tt2,dx,dy,dz,atoms%qat(1),atoms%qat(2)
            !write(*,'(a,2i5,2f10.5)') 'EE',iel,jel,cent%rgrad(3,1),cent%rgrad(3,2)
            cent%qgrad(iel)=cent%qgrad(iel)+atoms%qat(jel)*erf(gama*r)/r
            cent%qgrad(jel)=cent%qgrad(jel)+atoms%qat(iel)*erf(gama*r)/r
        enddo
    enddo
    !-------------------------------------------------------
    !nuclei-nuclei interaction
    do iat=1,atoms%nat !summation over ions
        do jat=iat+1,atoms%nat !summation over ions
            dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
            dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
            dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwi(jat)**2)
            epot_es=epot_es              +atoms%zat(iat)*atoms%zat(jat)*erf(gama*r)/r
        enddo
    enddo
    !-------------------------------------------------------
    !electron-nuclei interaction
    do iat=1,atoms%nat !summation over ions
        do jel=1,atoms%nat !summation over electrons
            dx=cent%rel(1,jel)-atoms%ratp(1,iat)
            dy=cent%rel(2,jel)-atoms%ratp(2,iat)
            dz=cent%rel(3,jel)-atoms%ratp(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jel)**2)
            if(r<0.1d0) then
                call erf_over_r_taylor(gama*r,tt1,ttg)
                epot_es=epot_es              +atoms%zat(iat)*atoms%qat(jel)*(tt1*gama)
                tt2=atoms%zat(iat)*atoms%qat(jel)*gama**3*ttg
                tt3=atoms%zat(iat)*(tt1*gama)
            else
                epot_es=epot_es              +atoms%zat(iat)*atoms%qat(jel)*erf(gama*r)/r
                ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
                tt2=atoms%zat(iat)*atoms%qat(jel)*ee1
                tt3=atoms%zat(iat)*erf(gama*r)/r
            endif
            cent%qgrad(jel)=cent%qgrad(jel)+tt3
            cent%rgrad(1,jel)=cent%rgrad(1,jel)+tt2*dx
            cent%rgrad(2,jel)=cent%rgrad(2,jel)+tt2*dy
            cent%rgrad(3,jel)=cent%rgrad(3,jel)+tt2*dz
            !write(*,'(a,2i5,4f10.5)') 'EN',iat,jel,tt2,dx,dy,dz
            !write(*,'(a,2i5,2f10.5)') 'EN',iat,jel,cent%rgrad(3,1),cent%rgrad(3,2)
        enddo
    enddo
    !-------------------------------------------------------
    !epot=epot+ener_ref
end subroutine cal_centt_pot_pairsum
!*****************************************************************************************
subroutine cal_centt_pairsum_force(parini,ann_arr,atoms,cent)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_cent
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_cent), intent(inout):: cent
    !local variables
    integer:: iat, jat, iel, jel
    real(8):: dx, dy, dz, r !, tt1, tt2, tt3, tt4, tt5, tt6 ,tt7, ttf
    real(8):: pi, gama
    real(8):: sqrt_one_over_twopi, tt1, tt2, ttg, ee1, tt3
    pi=4.d0*atan(1.d0)
    sqrt_one_over_twopi=1.d0/sqrt(2.d0*pi)
    !fat=0.d0
    !write(*,*) cent%gwi(:)
    !write(*,*) cent%gwe(:)
    !stop 'BBB'
    !-------------------------------------------------------
    !nuclei-nuclei interaction
    do iat=1,atoms%nat !summation over ions
        do jat=iat+1,atoms%nat !summation over ions
            dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
            dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
            dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            !gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwi(jat)**2)
            !epot_es=epot_es              +atoms%zat(iat)*atoms%zat(jat)*erf(gama*r)/r
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwi(jat)**2)
            ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
            tt2=atoms%zat(iat)*atoms%zat(jat)*ee1
            !---------------------------------------------------
            atoms%fat(1,jat)=atoms%fat(1,jat)-tt2*dx
            atoms%fat(2,jat)=atoms%fat(2,jat)-tt2*dy
            atoms%fat(3,jat)=atoms%fat(3,jat)-tt2*dz
            atoms%fat(1,iat)=atoms%fat(1,iat)+tt2*dx
            atoms%fat(2,iat)=atoms%fat(2,iat)+tt2*dy
            atoms%fat(3,iat)=atoms%fat(3,iat)+tt2*dz
        enddo
    enddo
    !-------------------------------------------------------
    !electron-nuclei interaction
    do iat=1,atoms%nat !summation over ions
        do jel=1,atoms%nat !summation over electrons
            dx=cent%rel(1,jel)-atoms%ratp(1,iat)
            dy=cent%rel(2,jel)-atoms%ratp(2,iat)
            dz=cent%rel(3,jel)-atoms%ratp(3,iat)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            if((iat==jel) .and. r>0.3d0) then
                write(*,'(a,es14.5,i6,1x,a)') 'ERROR: Center of electron far from atom: r= ', &
                    r,iat,trim(atoms%sat(iat))
                stop
            endif
            gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwe(jel)**2)
            if(r<0.1d0) then
                call erf_over_r_taylor(gama*r,tt1,ttg)
                tt2=atoms%zat(iat)*atoms%qat(jel)*gama**3*ttg
            else
                ee1=(2.d0/sqrt(pi)*gama*exp(-gama**2*r**2)-erf(gama*r)/r)/r**2
                tt2=atoms%zat(iat)*atoms%qat(jel)*ee1
            endif
            !-------------------------------------------
            atoms%fat(1,iat)=atoms%fat(1,iat)+tt2*dx
            atoms%fat(2,iat)=atoms%fat(2,iat)+tt2*dy
            atoms%fat(3,iat)=atoms%fat(3,iat)+tt2*dz
        enddo
    enddo
    !-------------------------------------------------------
    !epot=epot+ener_ref
end subroutine cal_centt_pairsum_force
!*****************************************************************************************
subroutine cal_centt_pot_bps(parini,ann_arr,atoms,cent,epot_es)
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
    real(8):: ehartree, pi, ehartree_t
    real(8):: time1, time2, time3, time4, time5, time6, time7
    real(8):: sqrt_one_over_twopi
    integer:: iat
    real(8),allocatable :: gausswidth(:)
    !logical:: ewald
    pi=4.d0*atan(1.d0)
    sqrt_one_over_twopi=1.d0/sqrt(2.d0*pi)
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
    call put_centt_gauss_to_grid(parini,atoms,cent)
    call cpu_time(time2)
    ehartree=0.d0
    allocate(gausswidth(atoms%nat))
    gausswidth(:)=1.d0 !TO_BE_CORRECTED
    call get_hartree(parini,cent%poisson,atoms,gausswidth,ehartree)
    deallocate(gausswidth)
    epot_es=0.d0
    do iat=1,atoms%nat
        epot_es=epot_es-atoms%zat(iat)**2*sqrt_one_over_twopi/cent%gwit(iat)
    enddo
    epot_es=epot_es+ehartree
    if(.not. cent%poisson%initialized) then
        stop 'ERROR: calculate_forces_energy: poisson is not initialized!'
    endif
    cent%poisson%nat=atoms%nat
    cent%poisson%cv=atoms%cellvec
    cent%poisson%bc=atoms%boundcond
    cent%poisson%q(1:cent%poisson%nat)=atoms%qat(1:atoms%nat)
    cent%poisson%gw(1:cent%poisson%nat)=cent%gwe(1:atoms%nat)
    cent%poisson%rcart(1:3,1:cent%poisson%nat)=cent%rel(1:3,1:atoms%nat)
    cent%poisson%rgrad(1:3,1:cent%poisson%nat)=0.d0
    cent%poisson%qgrad(1:cent%poisson%nat)=0.d0
    call get_hartree_grad_rho(parini,cent%poisson,atoms,ehartree_t)
    cent%rgrad(1:3,1:cent%poisson%nat)=cent%rgrad(1:3,1:cent%poisson%nat)+cent%poisson%rgrad(1:3,1:cent%poisson%nat)
    cent%qgrad(1:cent%poisson%nat)=cent%qgrad(1:cent%poisson%nat)+cent%poisson%qgrad(1:cent%poisson%nat)

    !if(ewald) then
    call cal_centt_shortrange_ewald(parini,ann_arr,atoms,cent,epot_es)
    !endif
    !gw_ion_t

    !write(61,'(f20.9)') epot_es
    !do iat=1,atoms%nat
    !    write(61,'(i4,4f10.5)') iat,rgrad(1,iat),rgrad(2,iat),rgrad(3,iat),qgrad(iat)
    !enddo
    !stop 'TESTING EWALD'
end subroutine cal_centt_pot_bps
!*****************************************************************************************
subroutine put_centt_gauss_to_grid(parini,atoms,cent)
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
    nx=cent%poisson%ngpx
    ny=cent%poisson%ngpy
    nz=cent%poisson%ngpz
    if(.not. cent%poisson%initialized) then
        stop 'ERROR: put_centt_gauss_to_grid: poisson is not initialized!'
    endif
    cent%poisson%reset_rho=.true.
    cent%poisson%nat=atoms%nat
    cent%poisson%cv=atoms%cellvec
    cent%poisson%bc=atoms%boundcond
    cent%poisson%q(1:cent%poisson%nat)=atoms%zat(1:atoms%nat)
    cent%poisson%gw(1:cent%poisson%nat)=cent%gwit(1:atoms%nat)
    cent%poisson%rcart(1:3,1:cent%poisson%nat)=atoms%ratp(1:3,1:atoms%nat)
    call put_charge_density(parini,cent%poisson)
    if(.not. cent%poisson%initialized) then
        stop 'ERROR: put_centt_gauss_to_grid: poisson is not initialized!'
    endif
    cent%poisson%reset_rho=.false.
    cent%poisson%nat=atoms%nat
    cent%poisson%cv=atoms%cellvec
    cent%poisson%bc=atoms%boundcond
    cent%poisson%q(1:cent%poisson%nat)=atoms%qat(1:atoms%nat)
    cent%poisson%gw(1:cent%poisson%nat)=cent%gwe(1:atoms%nat)
    cent%poisson%rcart(1:3,1:cent%poisson%nat)=cent%rel(1:3,1:atoms%nat)
    call put_charge_density(parini,cent%poisson)
end subroutine put_centt_gauss_to_grid
!*****************************************************************************************
subroutine cal_centt_shortrange_ewald(parini,ann_arr,atoms,cent,epot_es)
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
    do ib=1,cent%poisson%linked_lists%maxbound_rad
        iat=cent%poisson%linked_lists%bound_rad(1,ib)
        jat=cent%poisson%linked_lists%bound_rad(2,ib)
        !---------------------------------------------------
        dx=cent%poisson%pia_arr%pia(ib)%dr(1)
        dy=cent%poisson%pia_arr%pia(ib)%dr(2)
        dz=cent%poisson%pia_arr%pia(ib)%dr(3)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gama=1.d0/sqrt(cent%gwi(iat)**2+cent%gwi(jat)**2)
        epot_short=epot_short-atoms%zat(iat)*atoms%zat(jat)*erfc(gama*r)/r
        gama=1.d0/(sqrt(2.d0)*alpha)
        epot_short=epot_short+atoms%zat(iat)*atoms%zat(jat)*erfc(gama*r)/r
        !---------------------------------------------------
        dx=cent%poisson%pia_arr%pia(ib)%dr(1)+cent%rel(1,jat)-atoms%ratp(1,jat)
        dy=cent%poisson%pia_arr%pia(ib)%dr(2)+cent%rel(2,jat)-atoms%ratp(2,jat)
        dz=cent%poisson%pia_arr%pia(ib)%dr(3)+cent%rel(3,jat)-atoms%ratp(3,jat)
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
        dx=-cent%poisson%pia_arr%pia(ib)%dr(1)+cent%rel(1,iat)-atoms%ratp(1,iat)
        dy=-cent%poisson%pia_arr%pia(ib)%dr(2)+cent%rel(2,iat)-atoms%ratp(2,iat)
        dz=-cent%poisson%pia_arr%pia(ib)%dr(3)+cent%rel(3,iat)-atoms%ratp(3,iat)
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
        dx=cent%rel(1,jat)-atoms%ratp(1,iat)
        dy=cent%rel(2,jat)-atoms%ratp(2,iat)
        dz=cent%rel(3,jat)-atoms%ratp(3,iat)
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
end subroutine cal_centt_shortrange_ewald
!*****************************************************************************************
subroutine cal_shortrange_ewald_force_centt(parini,ann_arr,atoms,cent)
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
    do ib=1,cent%poisson%linked_lists%maxbound_rad
        iat=cent%poisson%linked_lists%bound_rad(1,ib)
        jat=cent%poisson%linked_lists%bound_rad(2,ib)
        !---------------------------------------------------
        dx=cent%poisson%pia_arr%pia(ib)%dr(1)
        dy=cent%poisson%pia_arr%pia(ib)%dr(2)
        dz=cent%poisson%pia_arr%pia(ib)%dr(3)
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
        dx=cent%poisson%pia_arr%pia(ib)%dr(1)+cent%rel(1,jat)-atoms%ratp(1,jat)
        dy=cent%poisson%pia_arr%pia(ib)%dr(2)+cent%rel(2,jat)-atoms%ratp(2,jat)
        dz=cent%poisson%pia_arr%pia(ib)%dr(3)+cent%rel(3,jat)-atoms%ratp(3,jat)
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
        dx=-cent%poisson%pia_arr%pia(ib)%dr(1)+cent%rel(1,iat)-atoms%ratp(1,iat)
        dy=-cent%poisson%pia_arr%pia(ib)%dr(2)+cent%rel(2,iat)-atoms%ratp(2,iat)
        dz=-cent%poisson%pia_arr%pia(ib)%dr(3)+cent%rel(3,iat)-atoms%ratp(3,iat)
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
        dx=cent%rel(1,jat)-atoms%ratp(1,iat)
        dy=cent%rel(2,jat)-atoms%ratp(2,iat)
        dz=cent%rel(3,jat)-atoms%ratp(3,iat)
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
end subroutine cal_shortrange_ewald_force_centt
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
subroutine calc_multipoles_centt(parini,atoms,poisson,rel)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: rel(3,atoms%nat)
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat
    real(8):: cx, cy, cz, qtot, ztot
    real(8):: dpx, dpy, dpz
    real(8):: qpxx, qpyy, qpzz, qpxy, qpxz, qpyz
    real(8):: x, y, z, rsq
    cx=0.d0 ; cy=0.d0 ; cz=0.d0
    ztot=0.d0
    do iat=1,atoms%nat
        ztot=ztot+atoms%zat(iat)
        cx=cx+atoms%zat(iat)*atoms%ratp(1,iat)
        cy=cy+atoms%zat(iat)*atoms%ratp(2,iat)
        cz=cz+atoms%zat(iat)*atoms%ratp(3,iat)
    enddo
    cx=cx/ztot
    cy=cy/ztot
    cz=cz/ztot
    dpx=0.0 ; dpy=0.0 ; dpz=0.0
    qpxx=0.0 ; qpyy=0.0 ; qpzz=0.0
    qpxy=0.0 ; qpxz=0.0 ; qpyz=0.0
    do iat=1,atoms%nat
        !ionic part
        x=atoms%ratp(1,iat)-cx
        y=atoms%ratp(2,iat)-cy
        z=atoms%ratp(3,iat)-cz
        rsq=x**2+y**2+z**2
        qpxx=qpxx+(3.0*x*x-rsq)*atoms%zat(iat)
        qpyy=qpyy+(3.0*y*y-rsq)*atoms%zat(iat)
        qpzz=qpzz+(3.0*z*z-rsq)*atoms%zat(iat)
        qpxy=qpxy+(3.0*x*y)*atoms%zat(iat)
        qpxz=qpxz+(3.0*x*z)*atoms%zat(iat)
        qpyz=qpyz+(3.0*y*z)*atoms%zat(iat)
        !electronic part
        x=rel(1,iat)-cx
        y=rel(2,iat)-cy
        z=rel(3,iat)-cz
        rsq=x**2+y**2+z**2
        dpx=dpx+atoms%qat(iat)*x
        dpy=dpy+atoms%qat(iat)*y
        dpz=dpz+atoms%qat(iat)*z
        qpxx=qpxx+(3.0*x*x-rsq)*atoms%qat(iat)
        qpyy=qpyy+(3.0*y*y-rsq)*atoms%qat(iat)
        qpzz=qpzz+(3.0*z*z-rsq)*atoms%qat(iat)
        qpxy=qpxy+(3.0*x*y)*atoms%qat(iat)
        qpxz=qpxz+(3.0*x*z)*atoms%qat(iat)
        qpyz=qpyz+(3.0*y*z)*atoms%qat(iat)
    enddo
    poisson%dpm(1)=dpx ; poisson%dpm(2)=dpy ; poisson%dpm(3)=dpz
    poisson%qpm(1,1)=qpxx ; poisson%qpm(1,2)=qpxy ; poisson%qpm(1,3)=qpxz
    poisson%qpm(2,1)=qpxy ; poisson%qpm(2,2)=qpyy ; poisson%qpm(2,3)=qpyz
    poisson%qpm(3,1)=qpxz ; poisson%qpm(3,2)=qpyz ; poisson%qpm(3,3)=qpzz
    !write(*,'(a,3f8.4)') 'Dipole moments',dpx,dpy,dpz
    !write(*,'(a,3f8.4)') 'Quadrupole moments',qpxx,qpxy,qpxz
    !write(*,'(a,3f8.4)') 'Quadrupole moments',qpxy,qpyy,qpyz
    !write(*,'(a,3f8.4)') 'Quadrupole moments',qpxz,qpyz,qpzz
end subroutine calc_multipoles_centt
!*****************************************************************************************
subroutine calc_multipoles_grid_centt(parini,atoms,poisson)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_poisson
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    type(typ_poisson), intent(inout):: poisson
    !local variables
    integer:: iat, ix, iy, iz
    real(8):: cx, cy, cz, qtot, ztot
    real(8):: dpx, dpy, dpz
    real(8):: qpxx, qpyy, qpzz, qpxy, qpxz, qpyz
    real(8):: x, y, z, rsq, vol_voxel, rho
    cx=0.d0 ; cy=0.d0 ; cz=0.d0
    ztot=0.d0
    do iat=1,atoms%nat
        ztot=ztot+atoms%zat(iat)
        cx=cx+atoms%zat(iat)*atoms%ratp(1,iat)
        cy=cy+atoms%zat(iat)*atoms%ratp(2,iat)
        cz=cz+atoms%zat(iat)*atoms%ratp(3,iat)
    enddo
    cx=cx/ztot
    cy=cy/ztot
    cz=cz/ztot
    dpx=0.0 ; dpy=0.0 ; dpz=0.0
    qpxx=0.0 ; qpyy=0.0 ; qpzz=0.0
    qpxy=0.0 ; qpxz=0.0 ; qpyz=0.0
    do iz=1,poisson%ngpz
        do iy=1,poisson%ngpy
            do ix=1,poisson%ngpx
                x=(ix-1)*poisson%hgrid(1,1)-cx
                y=(iy-1)*poisson%hgrid(2,2)-cy
                z=(iz-1)*poisson%hgrid(3,3)-cz
                rsq=x**2+y**2+z**2
                rho=poisson%rho(ix,iy,iz)
                dpx=dpx+x*rho
                dpy=dpy+y*rho
                dpz=dpz+z*rho
                qpxx=qpxx+(3.0*x*x-rsq)*rho
                qpyy=qpyy+(3.0*y*y-rsq)*rho
                qpzz=qpzz+(3.0*z*z-rsq)*rho
                qpxy=qpxy+(3.0*x*y)*rho
                qpxz=qpxz+(3.0*x*z)*rho
                qpyz=qpyz+(3.0*y*z)*rho
            enddo
        enddo
    enddo
    vol_voxel=poisson%hgrid(1,1)*poisson%hgrid(2,2)*poisson%hgrid(3,3)
    dpx=dpx*vol_voxel
    dpy=dpy*vol_voxel
    dpz=dpz*vol_voxel
    qpxx=qpxx*vol_voxel
    qpyy=qpyy*vol_voxel
    qpzz=qpzz*vol_voxel
    qpxy=qpxy*vol_voxel
    qpxz=qpxz*vol_voxel
    qpyz=qpyz*vol_voxel
    do iat=1,atoms%nat
        !ionic part
        x=atoms%ratp(1,iat)-cx
        y=atoms%ratp(2,iat)-cy
        z=atoms%ratp(3,iat)-cz
        rsq=x**2+y**2+z**2
        qpxx=qpxx+(3.0*x*x-rsq)*atoms%zat(iat)
        qpyy=qpyy+(3.0*y*y-rsq)*atoms%zat(iat)
        qpzz=qpzz+(3.0*z*z-rsq)*atoms%zat(iat)
        qpxy=qpxy+(3.0*x*y)*atoms%zat(iat)
        qpxz=qpxz+(3.0*x*z)*atoms%zat(iat)
        qpyz=qpyz+(3.0*y*z)*atoms%zat(iat)
    enddo
    poisson%dpm(1)=dpx ; poisson%dpm(2)=dpy ; poisson%dpm(3)=dpz
    poisson%qpm(1,1)=qpxx ; poisson%qpm(1,2)=qpxy ; poisson%qpm(1,3)=qpxz
    poisson%qpm(2,1)=qpxy ; poisson%qpm(2,2)=qpyy ; poisson%qpm(2,3)=qpyz
    poisson%qpm(3,1)=qpxz ; poisson%qpm(3,2)=qpyz ; poisson%qpm(3,3)=qpzz
    !write(*,'(a,3f8.4)') 'Dipole moments',dpx,dpy,dpz
    !write(*,'(a,3f8.4)') 'Quadrupole moments',qpxx,qpxy,qpxz
    !write(*,'(a,3f8.4)') 'Quadrupole moments',qpxy,qpyy,qpyz
    !write(*,'(a,3f8.4)') 'Quadrupole moments',qpxz,qpyz,qpzz
end subroutine calc_multipoles_grid_centt
!*****************************************************************************************
