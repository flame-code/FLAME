!*****************************************************************************************
subroutine cal_ann_atombased(parini,atoms,symfunc,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
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
    integer:: iat, i, j, ng, jat, ib
    real(8):: epoti, tt,vol
    real(8):: ttx, tty, ttz
    real(8):: sxx, sxy, sxz, syx, syy, syz, szx, szy, szz
    call f_routine(id='cal_ann_atombased')
    if(ann_arr%compute_symfunc) then
        call symmetry_functions(parini,ann_arr,atoms,symfunc,.true.)
    else
        symfunc%linked_lists%rcut=ann_arr%rcut
        symfunc%linked_lists%triplex=.true.
        call call_linkedlist(parini,atoms,.true.,symfunc%linked_lists,pia_arr_tmp)
    endif
    if(parini%save_symfunc_behnam) then
        ann_arr%cal_force=.false.
    endif
    i=1
    atoms%epot=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    over_iat: do iat=1,atoms%nat
        ng=ann_arr%ann(i)%nn(0)
        !if(ann_arr%compute_symfunc) then
        !    ann_arr%ann(i)%y(1:ng,0)=ann_arr%yall(1:ng,iat)
        !else
            ann_arr%ann(i)%y(1:ng,0)=symfunc%y(1:ng,iat)
        !endif
        if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
            call cal_architecture(ann_arr%ann(i),epoti)
            if(ann_arr%cal_force) then
            do ib=symfunc%linked_lists%prime_bound(iat),symfunc%linked_lists%prime_bound(iat+1)-1
                jat=symfunc%linked_lists%bound_rad(2,ib)
                ttx=0.d0 ; tty=0.d0 ; ttz=0.d0
                sxx=0.d0 ; sxy=0.d0 ; sxz=0.d0
                syx=0.d0 ; syy=0.d0 ; syz=0.d0
                szx=0.d0 ; szy=0.d0 ; szz=0.d0
                do j=1,ann_arr%ann(i)%nn(0)
                    ttx=ttx+ann_arr%ann(i)%d(j)*symfunc%y0d(j,1,ib)
                    tty=tty+ann_arr%ann(i)%d(j)*symfunc%y0d(j,2,ib)
                    ttz=ttz+ann_arr%ann(i)%d(j)*symfunc%y0d(j,3,ib)
                    sxx=sxx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,1,ib)
                    sxy=sxy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,2,ib)
                    sxz=sxz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,3,ib)
                    syx=syx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,4,ib)
                    syy=syy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,5,ib)
                    syz=syz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,6,ib)
                    szx=szx+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,7,ib)
                    szy=szy+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,8,ib)
                    szz=szz+ann_arr%ann(i)%d(j)*symfunc%y0dr(j,9,ib)
                enddo !over j
                atoms%fat(1,jat)=atoms%fat(1,jat)-ttx
                atoms%fat(2,jat)=atoms%fat(2,jat)-tty
                atoms%fat(3,jat)=atoms%fat(3,jat)-ttz
                atoms%fat(1,iat)=atoms%fat(1,iat)+ttx
                atoms%fat(2,iat)=atoms%fat(2,iat)+tty
                atoms%fat(3,iat)=atoms%fat(3,iat)+ttz
                atoms%stress(1,1)=atoms%stress(1,1)-sxx 
                atoms%stress(2,1)=atoms%stress(2,1)-syx
                atoms%stress(3,1)=atoms%stress(3,1)-szx
                atoms%stress(1,2)=atoms%stress(1,2)-sxy
                atoms%stress(2,2)=atoms%stress(2,2)-syy
                atoms%stress(3,2)=atoms%stress(3,2)-szy
                atoms%stress(1,3)=atoms%stress(1,3)-sxz
                atoms%stress(2,3)=atoms%stress(2,3)-syz
                atoms%stress(3,3)=atoms%stress(3,3)-szz
            enddo !over jat
            endif
        elseif(trim(ann_arr%event)=='train') then
            call cal_architecture_der(ann_arr%ann(i),epoti)
            call convert_ann_epotd(ann_arr%ann(i),ann_arr%num(i),ann_arr%g_per_atom(1,iat))
        else
            stop 'ERROR: undefined content for ann_arr%event'
        endif
        atoms%epot=atoms%epot+epoti
    enddo over_iat
    call cell_vol(atoms%nat,atoms%cellvec,vol)
    vol=vol*real(atoms%nat,8)
    atoms%stress(1:3,1:3)=atoms%stress(1:3,1:3)/vol
!    tt=(ann_arr%ann(1)%ebounds(2)-ann_arr%ann(1)%ebounds(1))/2.d0
!    !opt_ann%epot=(opt_ann%epot*atoms%nat+1.d0)*tt+ann_arr%ann(1)%ebounds(1)
!    opt_ann%epot=((opt_ann%epot+1.d0)*tt+ann_arr%ann(1)%ebounds(1))
!    tt=abs(opt_ann%epot-atoms%epot)/atoms%nat
    deallocate(symfunc%linked_lists%prime_bound)
    deallocate(symfunc%linked_lists%bound_rad)
    deallocate(symfunc%linked_lists%bound_ang)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        if(ann_arr%cal_force) then
        call f_free(symfunc%y)
        call f_free(symfunc%y0d)
        call f_free(symfunc%y0dr)
        endif
    endif
    call f_release_routine()
end subroutine cal_ann_atombased
!*****************************************************************************************
