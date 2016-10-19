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
    integer:: iat, jat, ng, i, j, k, isat, ib, nb
    real(8):: hgen_der(4,1:atoms%nat,1:atoms%nat)   !derivative of 
    real(8):: epotn, tt, epotdh, c
    call f_routine(id='cal_ann_tb')
    partb%paircut=ann_arr%rcut
    partb%hgenall0=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall0')
    partb%hgenall1=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall1')
    partb%hgenall2=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall2')
    partb%hgenall3=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%hgenall3')
    partb%dhgenall0=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall0')
    partb%dhgenall1=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall1')
    partb%dhgenall2=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall2')
    partb%dhgenall3=f_malloc([1.to.atoms%nat,1.to.atoms%nat],id='partb%dhgenall3')
    if(trim(ann_arr%event)=='train') then
        ekf%gc=f_malloc([1.to.ekf%num(1),1.to.4],id='ekf%gc') !HERE
    endif
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    endif
    if(symfunc%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
    nb=symfunc%linked_lists%maxbound_rad/2
    hgen=f_malloc([1.to.4,1.to.nb],id='hgen')
    dhgen=f_malloc([1.to.4,1.to.nb],id='dhgen')
    over_i: do i=1,4
        over_ib: do ib=1,nb
            ng=ann_arr%ann(i)%nn(0)
            !if(ann_arr%compute_symfunc) then
            !    ann_arr%ann(i)%y(1:ng,0)=ann_arr%yall(1:ng,ib)
            !else
                ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,ib)
            !endif
            if(trim(ann_arr%event)=='train') then
                call cal_architecture_der(ann_arr%ann(i),hgen(i,ib))
                call convert_ann_epotd(ann_arr%ann(i),ekf%num(i),ekf%gc(1,i))
            elseif(trim(ann_arr%event)=='evalu') then
                call cal_architecture(ann_arr%ann(i),hgen(i,ib))
                !call cal_architecture_force(ann_arr%ann(i),atoms%nat,hgen(i,iat,jat),atoms%fat(1,iat))
                !r=sqrt((atoms%rat(1,iat))**2+(atoms%rat(2,iat))**2+(atoms%rat(3,iat))**2)
                !dhgen(i,iat,jat)=(atoms%rat(1,iat)/r)*atoms%fat(1,iat)
            else
                stop 'ERROR: undefined content for ann_arr%event'
            endif
        enddo over_ib
    enddo over_i
    do ib=1,nb
        iat=symfunc%linked_lists%bound_rad(1,ib)
        jat=symfunc%linked_lists%bound_rad(2,ib)
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
        !write(*,'(a,4es19.10)') 'hgen-A ',hgen(1,1),hgen(2,1),hgen(3,1),hgen(4,1)
    enddo
        if(trim(ann_arr%event)=='train') then
            partb%event=ann_arr%event
        endif
        call lenoskytb_ann(partb,atoms,atoms%nat,c)
        atoms%epot=atoms%epot-0.23d0
        !atoms%epot=atoms%epot-0.2208033067776594d0
        write(*,*) 'energy ',atoms.epot
        if(trim(ann_arr%event)=='train') then
            ekf%g=0.d0
            do i=1,4
                !do ib=1,nb
                    do j=1,ekf%num(1)
                        ekf%g(ekf%loc(i)+j-1)=partb%dedh(i)*ekf%gc(j,i)
                        !write(*,'(a,2i5,3es14.5)') 'GGG ',i,j,ekf%g(ekf%loc(i)+j-1),partb%dedh(i),ekf%gc(j,i)
                    enddo
                !enddo
            enddo
        endif
    tt=(ann_arr%ann(1)%ebounds(2)-ann_arr%ann(1)%ebounds(1))/2.d0
    atoms%epot=((atoms%epot+1.d0)*tt+ann_arr%ann(1)%ebounds(1)) !*atoms%nat
    call f_free(hgen)
    call f_free(dhgen)
    if(trim(ann_arr%event)=='train') then
        call f_free(ekf%gc)
    endif
    call f_release_routine()
end subroutine cal_ann_tb
!*****************************************************************************************
subroutine lenoskytb_ann(partb,atoms,natsi,count_md)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_frame, only: CLSMAXATOM, CLSMAXNEIGHB
    use mod_potl, only: potl_typ
    use mod_tightbinding, only: typ_partb, lenosky
    implicit none
    type(typ_partb), intent(inout):: partb
    type(typ_atoms), intent(inout):: atoms
    integer, intent(in):: natsi
    real(8), intent(inout):: count_md
    !local variables
    integer, save:: icall=0
    type(potl_typ), save:: pplocal
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
    call lenoskytb_final(partb)
end subroutine lenoskytb_ann
!*****************************************************************************************
