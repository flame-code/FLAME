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
    real(8), allocatable:: grad1_p1(:,:)
    real(8), allocatable:: grad1_p2(:,:)
    real(8), allocatable:: grad2_p1(:)
    real(8), allocatable:: grad2_p2(:)
    allocate(grad1_p1(3,atoms%nat),source=0.d0)
    allocate(grad1_p2(3,atoms%nat),source=0.d0)
    allocate(grad2_p1(atoms%nat),source=0.d0)
    allocate(grad2_p2(atoms%nat),source=0.d0)
    atoms%fat=0.d0
    atoms%epot=0.d0
    do iat=1,atoms%nat !summation over ions/electrons
        atoms%epot=atoms%epot+ann_arr%chi_o(iat)*(atoms%zat(iat)+atoms%qat(iat))
        hardness=ann_arr%ann(atoms%itypat(iat))%hardness
        atoms%epot=atoms%epot+0.5d0*hardness*(atoms%zat(iat)+atoms%qat(iat))**2
        grad2_p1(iat)=grad2_p1(iat)+ann_arr%chi_o(iat)+(atoms%zat(iat)+atoms%qat(iat))*hardness
        dx=rel(1,iat)-atoms%rat(1,iat)
        dy=rel(2,iat)-atoms%rat(2,iat)
        dz=rel(3,iat)-atoms%rat(3,iat)
        spring_const=ann_arr%ann(atoms%itypat(iat))%spring_const
        atoms%epot=atoms%epot+0.5d0*spring_const*(dx**2+dy**2+dz**2)
        grad1_p1(1,iat)=grad1_p1(1,iat)+spring_const*dx
        grad1_p1(2,iat)=grad1_p1(2,iat)+spring_const*dy
        grad1_p1(3,iat)=grad1_p1(3,iat)+spring_const*dz
    enddo
    !call cal_electrostatic_cent2(nat,atoms%rat,rel,atoms%qat,epot_es, &
    !    atoms%fat,grad1_p1,grad2_p1,epot_atom)
    call cal_pot_with_bps(ann_arr,atoms,rel,epot_es_bps,grad1_p2,grad2_p2)
    epot_es=epot_es_bps
    do iat=1,atoms%nat
        grad1(1,iat)=grad1_p1(1,iat)+grad1_p2(1,iat)
        grad1(2,iat)=grad1_p1(2,iat)+grad1_p2(2,iat)
        grad1(3,iat)=grad1_p1(3,iat)+grad1_p2(3,iat)
        grad2(iat)=grad2_p1(iat)+grad2_p2(iat)
    enddo
    atoms%epot=atoms%epot+epot_es
    !epot=epot+ener_ref
    deallocate(grad1_p1)
    deallocate(grad1_p2)
    deallocate(grad2_p1)
    deallocate(grad2_p2)
end subroutine cal_potential_cent2
!*****************************************************************************************
subroutine cal_pot_with_bps(ann_arr,atoms,rel,epot_es,grad1,grad2)
    use mod_interface
    use mod_ann, only: typ_ann_arr
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    use dynamic_memory
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    real(8), intent(in):: rel(3,atoms%nat)
    real(8), intent(inout):: epot_es, grad1(3,atoms%nat), grad2(atoms%nat)
    !local variables
    real(8):: cell(3)
    type(typ_parini):: parini
    !type(typ_ewald_p3d):: ewald_p3d_rough
    type(typ_ewald_p3d):: ewald_p3d
    real(8):: ehartree, error, pi
    real(8):: time1, time2, time3, time4, time5, time6, time7
    real(8), allocatable:: potref(:,:,:)
    real(8), allocatable:: rel_t(:,:)
    real(8), allocatable:: ratred(:,:), fat(:,:), fat_m(:,:)
    real(8), allocatable:: gw_ion(:), gw(:), eqd(:), qat_tot(:)
    real(8):: ehartree_kwald, stress(3,3), celldv(3,3), stress_m(3,3)
    real(8):: x, y, z, v1, v2
    integer:: ix, iy, iz, iat
    pi=4.d0*atan(1.d0)
    associate(nx=>ewald_p3d%poisson_p3d%ngpx)
    associate(ny=>ewald_p3d%poisson_p3d%ngpy)
    associate(nz=>ewald_p3d%poisson_p3d%ngpz)
    rel_t=f_malloc([1.to.3,1.to.atoms%nat],id='rel_t')
    ewald_p3d%poisson_p3d%ngpx=50
    ewald_p3d%poisson_p3d%ngpy=50
    ewald_p3d%poisson_p3d%ngpz=50
    ewald_p3d%poisson_p3d%rho=f_malloc([1.to.nx,1.to.ny,1.to.nz],id='rho')
    ewald_p3d%poisson_p3d%pot=f_malloc([1.to.nx,1.to.ny,1.to.nz],id='pot')
    potref=f_malloc([1.to.nx,1.to.ny,1.to.nz],id='potref')
    rel_t(1:3,1:atoms%nat)=rel(1:3,1:atoms%nat)
    !call put_in_cell(atoms,rel_t,cell)
    cell(1)=atoms%cellvec(1,1)
    cell(2)=atoms%cellvec(2,2)
    cell(3)=atoms%cellvec(3,3)
    ewald_p3d%hgx=cell(1)/nx
    ewald_p3d%hgy=cell(2)/ny
    ewald_p3d%hgz=cell(3)/nz
    ewald_p3d%rgcut=6.d0/0.529d0 !parini%rgcut_ewald*ewald_p3d%alpha
    ewald_p3d%nbgpx=int(ewald_p3d%rgcut/ewald_p3d%hgx)+2
    ewald_p3d%nbgpy=int(ewald_p3d%rgcut/ewald_p3d%hgy)+2
    ewald_p3d%nbgpz=int(ewald_p3d%rgcut/ewald_p3d%hgz)+2
    ewald_p3d%nagpx=ewald_p3d%nbgpx+1
    ewald_p3d%nagpy=ewald_p3d%nbgpy+1
    ewald_p3d%nagpz=ewald_p3d%nbgpz+1

    !-------------------------------------------------------
    atoms%stress=0.d0
    stress=0.d0
    stress_m=0.d0
    allocate(gw(atoms%nat),gw_ion(atoms%nat))
    allocate(ratred(3,atoms%nat))
    allocate(fat(3,atoms%nat),eqd(atoms%nat),qat_tot(atoms%nat))
    allocate(fat_m(3,atoms%nat))
    do iat=1,atoms%nat
        gw_ion(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth_ion
        gw(iat)=ann_arr%ann(atoms%itypat(iat))%gausswidth
    enddo
    qat_tot(1:atoms%nat)=atoms%zat(1:atoms%nat)+atoms%qat(1:atoms%nat)
    !-------------------------------------------------------
    call cpu_time(time1)
    call put_gauss_to_grid(parini,atoms,rel_t,gw_ion,gw,ewald_p3d,potref)

    call cpu_time(time2)
    call construct_ewald_bps(parini,atoms,ewald_p3d)
    call cpu_time(time3)
    call cal_hartree_pot_bps(ewald_p3d,atoms,ehartree)


    !call longerange_forces(parini,atoms%boundcond,.true. ,atoms%nat,atoms%rat,qat_tot,gw_ion,ewald_p3d,fat_m,stress_m) !CORRECT_IT
    !call longerange_forces(parini,atoms%boundcond,.true. ,atoms%nat,atoms%rat,zat,gw_ion,ewald_p3d,fat_m,stress_m)
    !call longerange_forces(parini,atoms%boundcond,.false.,atoms%nat,rel      ,qat,gw    ,ewald_p3d,fat_m,stress_m)
    !write(*,*) 'cell(1)*cell(2)*cell(3) ',cell(1)*cell(2)*cell(3)
    write(*,*) 'VOLUME ',(cell(1)*cell(2)*cell(3))
    stress_m(1:3,1:3)=stress_m(1:3,1:3)/(cell(1)*cell(2)*cell(3))

    write(*,'(a,3es14.5)') 'stress (BPS) ',atoms%stress(1,1),atoms%stress(1,2),atoms%stress(1,3)
    write(*,'(a,3es14.5)') 'stress (BPS) ',atoms%stress(2,1),atoms%stress(2,2),atoms%stress(2,3)
    write(*,'(a,3es14.5)') 'stress (BPS) ',atoms%stress(3,1),atoms%stress(3,2),atoms%stress(3,3)
    write(*,*)
    write(*,'(a,es14.5)') 'stress (BPS) ',atoms%stress(1,1)*(cell(1)*cell(2)*cell(3))

    write(*,'(a,f20.10)') 'ehartree ',ehartree
    write(*,'(a,2f20.10,es14.5)') 'forces ',fat(1,1),fat_m(1,1),fat(1,1)-fat_m(1,1)
    write(*,'(a,2f20.10,es14.5)') 'stress ',atoms%stress(1,1),stress(1,1),atoms%stress(1,1)-stress(1,1)
    stop 'STOPPED TO COMPARE STRESS'
    !!!! call cpu_time(time4)
    !!!! call cal_grad_long(atoms,qat,rel_t,ewald_p3d,grad1,grad2)
    !!!! call cpu_time(time5)
    !!!! call destruct_ewald_bps(ewald_p3d)
    !!!! call cpu_time(time6)
    !!!! epot_es=0.d0
    !!!! call cal_shortrange_ewald(atoms,qat,rel_t,epot_es,grad1,grad2)
    !!!! call cpu_time(time7)
    !!!! write(*,'(a,7f6.2)') 'TIME ',time2-time1,time3-time2,time4-time3,time5-time4, &
    !!!!                              time6-time5,time7-time6,time7-time1

    !do iat=1,atoms%nat
    !    epot_es=epot_es-zat(iat)**2/(gw_ion(iat)*sqrt(2.d0*pi))
    !enddo
    write(*,*) 'ehartree ',ehartree
    epot_es=epot_es+ehartree
    !write(*,*) 'RMSE ',ewald_p3d%hgx,error
    call f_free(potref)
    call f_free(ewald_p3d%poisson_p3d%pot)
    call f_free(ewald_p3d%poisson_p3d%rho)
    call f_free(atoms%rat)
    call f_free(rel_t)
    end associate
    end associate
    end associate
end subroutine cal_pot_with_bps
!*****************************************************************************************
subroutine put_gauss_to_grid(parini,atoms,rel,gw_ion,gw,ewald_p3d,potref)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_electrostatics, only: typ_ewald_p3d
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    real(8), intent(in):: rel(3,atoms%nat)
    real(8), intent(in):: gw_ion(atoms%nat)
    real(8), intent(in):: gw(atoms%nat)
    type(typ_ewald_p3d), intent(inout):: ewald_p3d
    real(8), intent(inout):: potref(ewald_p3d%poisson_p3d%ngpx,ewald_p3d%poisson_p3d%ngpy,ewald_p3d%poisson_p3d%ngpz)
    !local variables
    integer:: ix, iy, iz, iat
    real(8):: x, y, z, r, rsq, pi, ttg, factor, tt1, qtot, alpha
    associate(nx=>ewald_p3d%poisson_p3d%ngpx)
    associate(ny=>ewald_p3d%poisson_p3d%ngpy)
    associate(nz=>ewald_p3d%poisson_p3d%ngpz)
    !character(4):: bc
    pi=4.d0*atan(1.d0)
    alpha=2.d0
    ewald_p3d%poisson_p3d%rho=0.d0
    !bc='bulk'
    !call putgaussgrid(parini,atoms%boundcond,.true. ,atoms%nat,atoms%rat,atoms%zat,gw_ion,ewald_p3d)
    !call putgaussgrid(parini,atoms%boundcond,.false.,atoms%nat,rel      ,atoms%qat,gw    ,ewald_p3d)
    call gauss_grid(parini,'bulk',.true.,atoms%nat,atoms%rat,atoms%cellvec,atoms%zat,gw_ion, &
        ewald_p3d%rgcut,nx,ny,nz,ewald_p3d%poisson_p3d%rho)
    call gauss_grid(parini,'bulk',.false.,atoms%nat,atoms%rat,atoms%cellvec,atoms%qat,gw    , &
        ewald_p3d%rgcut,nx,ny,nz,ewald_p3d%poisson_p3d%rho)
    end associate
    end associate
    end associate
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
    real(8):: hrxinv, hryinv, hrzinv
    real(8):: cvinv_norm
    real(8):: dmx, dmy, dmz, dmsq, gwsq_inv
    real(8):: xred, yred, zred
    real(8):: ximg, yimg, zimg
    integer:: imgx, imgy, imgz
    integer:: ncellx, ncelly, ncellz
    integer:: iat, igx, igy, igz, jgx, jgy, jgz
    integer:: iii
    integer:: nbgx, nbgy, nbgz, nagx, nagy, nagz, nex, ney, nez
    integer:: ilgx, ilgy, ilgz, irgx, irgy, irgz
    real(8), allocatable:: wa(:,:,:)
    real(8), allocatable:: wm(:,:,:)
    real(8), allocatable:: ratred(:,:)

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
    !detemining the largest dimension for the pseudogrid.
    hxx=cv(1,1)/ngx ; hxy=cv(2,1)/ngx ; hxz=cv(3,1)/ngx
    hyx=cv(1,2)/ngy ; hyy=cv(2,2)/ngy ; hyz=cv(3,2)/ngy
    hzx=cv(1,3)/ngz ; hzy=cv(2,3)/ngz ; hzz=cv(3,3)/ngz

    wa=f_malloc0([1-nagx.to.ngx+nagx,1-nagy.to.ngy+nagy,1-nagz.to.ngz+nagz],id='wa')
    wm=f_malloc0([1.to.ngx,1.to.ngy,1.to.ngz],id='wm')

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
    !-------------------------------------------------------
    do iat=1,nat
        gwsq_inv=1.d0/gw(iat)**2
        fac=1.d0/(gw(iat)*sqrt(pi))**3
        imgx=nint(ratred(1,iat)*hrxinv)+1
        imgy=nint(ratred(2,iat)*hryinv)+1
        imgz=nint(ratred(3,iat)*hrzinv)+1
        facqiat=fac*qat(iat)
        do igz=-nbgz,nbgz
            jgz=imgz+igz
            do igy=-nbgy,nbgy
                jgy=imgy+igy
                do igx=-nbgx,nbgx
                    jgx=imgx+igx
                    ximg=(jgx-1)*hxx+(jgy-1)*hyx+(jgz-1)*hzx
                    yimg=(jgx-1)*hxy+(jgy-1)*hyy+(jgz-1)*hzy
                    zimg=(jgx-1)*hxz+(jgy-1)*hyz+(jgz-1)*hzz
                    dmx=ximg-rxyz(1,iat)
                    dmy=yimg-rxyz(2,iat)
                    dmz=zimg-rxyz(3,iat)
                    dmsq=dmx**2+dmy**2+dmz**2
                    wa(jgx,jgy,jgz)=wa(jgx,jgy,jgz)+facqiat*exp(-dmsq*gwsq_inv)
                enddo
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    ncellx=nagx/ngx
    ncelly=nagy/ngy
    ncellz=nagz/ngz
    ilgx=-ncellx*ngx+1 ; irgx=(ncellx+1)*ngx
    ilgy=-ncelly*ngy+1 ; irgy=(ncelly+1)*ngy
    ilgz=-ncellz*ngz+1 ; irgz=(ncellz+1)*ngz
    !---------------------------------------------------------------------------
    !extended box constructed by integer number of cells
    nex=max(ngx,irgx-ilgx+1)
    ney=max(ngy,irgy-ilgy+1)
    nez=max(ngz,irgz-ilgz+1)
    !wrap around grid points that are outside the extended box in into the extended box,
    !these grid points do not form a complete cell.
    do igz=1-nagz,ngz+nagz
        do igy=1-nagy,ilgy-1
            do igx=1-nagx,ilgx-1
                wa(igx+nex,igy+ney,igz)=wa(igx+nex,igy+ney,igz)+wa(igx,igy,igz)
            enddo
            do igx=ilgx,irgx
                wa(igx,igy+ney,igz)=wa(igx,igy+ney,igz)+wa(igx,igy,igz)
            enddo
            do igx=irgx+1,ngx+nagx
                wa(igx-nex,igy+ney,igz)=wa(igx-nex,igy+ney,igz)+wa(igx,igy,igz)
            enddo
        enddo
        do igy=ilgy,irgy
            do igx=1-nagx,ilgx-1
                wa(igx+nex,igy,igz)=wa(igx+nex,igy,igz)+wa(igx,igy,igz)
            enddo
            do igx=irgx+1,ngx+nagx
                wa(igx-nex,igy,igz)=wa(igx-nex,igy,igz)+wa(igx,igy,igz)
            enddo
        enddo
        do igy=irgy+1,ngy+nagy
            do igx=1-nagx,ilgx-1
                wa(igx+nex,igy-ney,igz)=wa(igx+nex,igy-ney,igz)+wa(igx,igy,igz)
            enddo
            do igx=ilgx,irgx
                wa(igx,igy-ney,igz)=wa(igx,igy-ney,igz)+wa(igx,igy,igz)
            enddo
            do igx=irgx+1,ngx+nagx
                wa(igx-nex,igy-ney,igz)=wa(igx-nex,igy-ney,igz)+wa(igx,igy,igz)
            enddo
        enddo
    enddo
    do igz=1-nagz,ilgz-1
        do igy=ilgy,irgy
            do igx=ilgx,irgx
                wa(igx,igy,igz+nez)=wa(igx,igy,igz+nez)+wa(igx,igy,igz)
            enddo
        enddo
    enddo
    do igz=irgz+1,ngz+nagz
        do igy=ilgy,irgy
            do igx=ilgx,irgx
                wa(igx,igy,igz-nez)=wa(igx,igy,igz-nez)+wa(igx,igy,igz)
            enddo
        enddo
    enddo
    !---------------------------------------------------------------------------
    if(ncellx==0 .and. ncelly==0 .and. ncellz==0) then
        do igz=1,ngz
            do igy=1,ngy
                do igx=1,ngx
                    wm(igx,igy,igz)=wa(igx,igy,igz)
                enddo
            enddo
        enddo
    else
        !wrap around grid points which form a complete cell and are outside the main cell.
        do igz=ilgz,irgz
            jgz=modulo(igz-1,ngz)+1
            do igy=ilgy,irgy
                jgy=modulo(igy-1,ngy)+1
                do igx=ilgx,irgx
                    jgx=modulo(igx-1,ngx)+1
                    wm(jgx,jgy,jgz)=wm(jgx,jgy,jgz)+wa(igx,igy,igz)
                enddo
            enddo
        enddo
    endif
    !---------------------------------------------------------------------------
    if(reset) then
        !if the input array of charge density does not contain any previous value
        !wanted to be preserved.
        do igz=1,ngz
            do igy=1,ngy
                do igx=1,ngx
                    rho(igx,igy,igz)=wm(igx,igy,igz)
                enddo
            enddo
        enddo
    else
        !if the input array of charge density already some value that must be preserved.
        do igz=1,ngz
            do igy=1,ngy
                do igx=1,ngx
                    rho(igx,igy,igz)=rho(igx,igy,igz)+wm(igx,igy,igz)
                enddo
            enddo
        enddo
    endif
    !---------------------------------------------------------------------------
    deallocate(ratred)
    call f_free(wa)
    call f_free(wm)
end subroutine gauss_grid
!*****************************************************************************************
