!*****************************************************************************************
subroutine cal_ann_tb(parini,partb,atoms,ann_arr,symfunc,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_tightbinding, only: typ_partb
    use mod_potl, only: potl_typ
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, typ_symfunc, typ_ekf
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_ekf), intent(inout):: ekf
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc), intent(inout):: symfunc
    type(typ_partb), intent(inout):: partb
    !local variables
    type(potl_typ):: pplocal
    real(8), allocatable:: hgen(:,:), dhgen(:,:)
    integer:: iat, jat, ng, i, j, k, isat, ib, nb, ixyz
    real(8):: hgen_der(4,1:atoms%nat,1:atoms%nat)  , ttxyz !derivative of 
    real(8):: epotn, tt, epotdh, c, dx, dy, dz, r, rsq, hbar, fc, dfc, tt1
    !call f_routine(id='cal_ann_tb')
    atoms%fat=0.d0
    partb%paircut=ann_arr%rcut
    !partb%dedh=f_malloc0([1.to.4,1.to.atoms%nat,1.to.atoms%nat],id='partb%dedh')
    !partb%hgenall0=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall0')
    !partb%hgenall1=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall1')
    !partb%hgenall2=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall2')
    !partb%hgenall3=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall3')
    !partb%dhgenall0=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall0')
    !partb%dhgenall1=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall1')
    !partb%dhgenall2=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall2')
    !partb%dhgenall3=f_malloc0([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall3')
    allocate(partb%dedh(4,atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%hgenall0(atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%hgenall1(atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%hgenall2(atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%hgenall3(atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%dhgenall0(atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%dhgenall1(atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%dhgenall2(atoms%nat,atoms%nat),source=0.d0)
    allocate(partb%dhgenall3(atoms%nat,atoms%nat),source=0.d0)
    if(trim(ann_arr%event)=='train') then
        !The following is allocated with ekf%num(1), this means number of
        !nodes in the input layer is the same for all atom types.
        !Therefore, it must be fixed later.
        !g_per_atom=f_malloc([1.to.ekf%num(1),1.to.atoms%nat],id='g_per_atom') !HERE
        do i=1,ann_arr%n
            call convert_x_ann(ekf%num(i),ekf%x(ekf%loc(i)),ann_arr%ann(i))
        enddo
    endif
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    endif
    !if(symfunc%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
    nb=symfunc%linked_lists%maxbound_rad!/2
    !hgen=f_malloc([1.to.4,1.to.nb],id='hgen')
    !dhgen=f_malloc([1.to.4,1.to.nb],id='dhgen')
    allocate(hgen(4,nb))
    allocate(dhgen(4,nb))
    if(trim(ann_arr%event)=='train') then
        !ann_arr%g_per_bond=f_malloc([1.to.ekf%num(1),1.to.4,1.to.nb],id='ann_arr%g_per_bond') !HERE
        allocate(ann_arr%g_per_bond(ekf%num(1),4,nb))
    endif
    over_i: do i=1,4
        over_ib: do ib=1,nb
            ng=ann_arr%ann(i)%nn(0)
            ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,ib)
            if(trim(ann_arr%event)=='train') then
            !write(*,*) "symfunc", ng, nb, symfunc%y(1,ib)
                call cal_architecture_der(ann_arr%ann(i),hgen(i,ib))
                !write(*,*) "hopping", hgen(i,ib)
                call convert_ann_epotd(ann_arr%ann(i),ekf%num(i),ann_arr%g_per_bond(1,i,ib))
                !write(*,*) "dhda", ann_arr%g_per_bond(1,i,ib)
            elseif(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
                call cal_architecture(ann_arr%ann(i),hgen(i,ib))          
                iat=symfunc%linked_lists%bound_rad(1,ib)
                jat=symfunc%linked_lists%bound_rad(2,ib)
                dhgen(i,ib)=0.d0
                do j=1,ann_arr%ann(i)%nn(0)
                    dhgen(i,ib)=dhgen(i,ib) + ann_arr%ann(i)%d(j)*symfunc%y0d_bond(j,ib)
                enddo
            else
                stop 'ERROR: in cal_ann_tb undefined content for ann_arr%event'
            endif
        enddo over_ib
    enddo over_i
    do ib=1,nb
        iat=symfunc%linked_lists%bound_rad(1,ib)
        jat=symfunc%linked_lists%bound_rad(2,ib)
        dx=atoms%rat(1,iat)-atoms%rat(1,jat)
        dy=atoms%rat(2,iat)-atoms%rat(2,jat)
        dz=atoms%rat(3,iat)-atoms%rat(3,jat)
        rsq=dx**2+dy**2+dz**2
        r=sqrt(rsq)
        fc=(1.d0-rsq/partb%paircut**2)**3
        dfc=-6.d0*r*(1.d0-rsq/partb%paircut**2)**2/partb%paircut**2
        do i=1,4
            hbar=hgen(i,ib)
            hgen(i,ib)=hbar*fc
            dhgen(i,ib)=dhgen(i,ib)*fc+hbar*dfc
            do j=1, ekf%num(1)
                ann_arr%g_per_bond(j,i,ib)=fc*ann_arr%g_per_bond(j,i,ib)
            enddo
        enddo
        
        partb%hgenall0(iat,jat)=hgen(1,ib)
        partb%hgenall1(iat,jat)=hgen(2,ib)
        partb%hgenall2(iat,jat)=hgen(3,ib)
        partb%hgenall3(iat,jat)=hgen(4,ib)
        partb%hgenall0(jat,iat)=partb%hgenall0(iat,jat)
        partb%hgenall1(jat,iat)=partb%hgenall1(iat,jat)
        partb%hgenall2(jat,iat)=partb%hgenall2(iat,jat)
        partb%hgenall3(jat,iat)=partb%hgenall3(iat,jat)
        partb%dhgenall0(iat,jat)=dhgen(1,ib)
        partb%dhgenall1(iat,jat)=dhgen(2,ib)
        partb%dhgenall2(iat,jat)=dhgen(3,ib)
        partb%dhgenall3(iat,jat)=dhgen(4,ib)
        partb%dhgenall0(jat,iat)=partb%dhgenall0(iat,jat)
        partb%dhgenall1(jat,iat)=partb%dhgenall1(iat,jat)
        partb%dhgenall2(jat,iat)=partb%dhgenall2(iat,jat)
        partb%dhgenall3(jat,iat)=partb%dhgenall3(iat,jat)
        !if(trim(ann_arr%event)=='evalu' .and. icall>0) then
        !    write(100+icall2,'(a,5e26.17)') 'hgen-N',hgen(1,ib),hgen(2,ib),hgen(3,ib),hgen(4,ib),r
        !endif
    enddo
        partb%event=ann_arr%event
        call lenoskytb_ann(partb,atoms,atoms%nat,c)
        atoms%epot=atoms%epot+atoms%nat*ann_arr%ann(atoms%itypat(1))%ener_ref
        if(trim(ann_arr%event)=='train') then
            ekf%g=0.d0
            do i=1,4
                do j=1,ekf%num(1)
                    tt1=0.d0
                    do ib=1,nb
                    iat=symfunc%linked_lists%bound_rad(1,ib)
                    jat=symfunc%linked_lists%bound_rad(2,ib)
                    !write(*,'(a,5i4,2es14.5,i4)') 'DEDH',i,j,ib,iat,jat,partb%dedh(i,iat,jat),ann_arr%g_per_bond(j,i,ib),ekf%loc(i)+j-1
                    tt1=tt1+partb%dedh(i,iat,jat)*ann_arr%g_per_bond(j,i,ib)
                    enddo
                    ekf%g(ekf%loc(i)+j-1)=tt1
                enddo
            enddo
        endif
    tt=(ann_arr%ann(1)%ebounds(2)-ann_arr%ann(1)%ebounds(1))/2.d0
    atoms%epot=((atoms%epot+1.d0)*tt+ann_arr%ann(1)%ebounds(1)) !*atoms%nat
    !call f_free(hgen)
    !call f_free(dhgen)
    deallocate(hgen)
    deallocate(dhgen)
    if(trim(ann_arr%event)=='train') then
        !call f_free(ann_arr%g_per_bond)
        !call f_free(partb%dedh)
        deallocate(ann_arr%g_per_bond)
        deallocate(partb%dedh)
    endif
    !call f_release_routine()
end subroutine cal_ann_tb
!*****************************************************************************************
subroutine lenoskytb_ann(partb,atoms,natsi,count_md)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_frame, only: CLSMAXATOM, CLSMAXNEIGHB
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    use mod_const, only: ha2ev, bohr2ang
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
    !local variables
    integer, save:: icall=0, firstcall=1
    type(potl_typ), save:: pplocal
    type(typ_partb), save:: partb_pair
    real(8), allocatable:: fat(:,:)
    integer:: iat
    if(firstcall==1) then
        write(*,'(a)') 'Reading spline potential coeff.cls'
        call prmst38c(partb_pair,pplocal) !Reads potential 
        firstcall=0
    endif
    partb%usepairpot=partb_pair%usepairpot
    partb%paircut=partb_pair%paircut 

    icall=icall+1
    !write(*,'(a,f)') 'paircut= ',partb%paircut
    call lenoskytb_init(partb,atoms,natsi)
    count_md=count_md+1.d0
    !PRINT SOME WARNINGS
    if(atoms%nat>CLSMAXATOM) then
        write(*,'(a,i,a,i)') 'WARNING Number of atoms = ',atoms%nat,' is greater than CLSMAXATOM = ',CLSMAXATOM
        write(*,'(a)') 'THIS MAY CAUSE FAILURE OR UNPREDICTABLE BEHAVIOR'
        write(*,'(a)') 'INCREASE CLSMAXATOM in defines.h and RECOMPILE'
    endif
    if(atoms%nat>CLSMAXNEIGHB) then
        write(*,'(a,i,a,i)') 'WARNING Number of atoms = ',atoms%nat,' is greater than CLSMAXNEIGHB = ',CLSMAXNEIGHB
        write(*,'(a)') 'THIS MAY CAUSE FAILURE OR UNPREDICTABLE BEHAVIOR'
        write(*,'(a)') 'INCREASE CLSMAXNEIGHB in defines.h and RECOMPILE'
    endif
    if(natsi>atoms%nat) write(*,'(a)') 'WARNING natsi = ',natsi,' is greater than number of atoms = ',atoms%nat
    if(atoms%nat == 0) write(*,'(a)') 'WARNING lenoskytb called with zero atoms'
    if(atoms%nat < 0) write(*,'(a)') 'WARNING lenoskytb called with negative number of atoms'
    !write(*,*) 'natsi ',natsi
    call totalenergy(partb,atoms,natsi,pplocal) 

    allocate(fat(3,atoms%nat))
    do iat=1,atoms%nat
        fat(1,iat)=atoms%fat(1,iat)
        fat(2,iat)=atoms%fat(2,iat)
        fat(3,iat)=atoms%fat(3,iat)
        atoms%fat(1,iat)=0.d0
        atoms%fat(2,iat)=0.d0
        atoms%fat(3,iat)=0.d0
    enddo

    call pairenergy(partb,atoms,pplocal,natsi)
    atoms%epot=atoms%epot+partb%pairen
    !write(61,*) atoms%epot
    !stop 'BBBBBBBBB'

    do iat=1,atoms%nat
        atoms%fat(1,iat)=fat(1,iat)+atoms%fat(1,iat)
        atoms%fat(2,iat)=fat(2,iat)+atoms%fat(2,iat)
        atoms%fat(3,iat)=fat(3,iat)+atoms%fat(3,iat)
    enddo
    atoms%epot=atoms%epot+partb%pairen
    deallocate(fat)

    call lenoskytb_final(partb)
end subroutine lenoskytb_ann
!*****************************************************************************************
