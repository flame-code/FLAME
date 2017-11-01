!*****************************************************************************************
subroutine symmetry_functions(parini,ann_arr,atoms,symfunc,apply_gbounds)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc
    use mod_atoms, only: typ_atoms
    use time_profiling
    use mod_timing , only: TCAT_SYMFUNC_COMPUT
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms), intent(inout):: atoms
    type(typ_symfunc), intent(inout):: symfunc
    logical, intent(in):: apply_gbounds
    !local variables
    integer:: i0, iat, jat, isat, i, ib, istat, ng, ig
    real(8):: gleft
    call f_routine(id='symmetry_functions')
    call f_timing(TCAT_SYMFUNC_COMPUT,'ON')
    !-----------------------------------------------------------------
    !first index of "y0d" is for number of symmetry function
    !second index of "y0d" is for x,y,z
    !3th index of "y0d" is for number of atoms
    !-----------------------------------------------------------------
    !write(*,*) ann_arr%ann(1)%nn(0),atoms%nat,allocated(ann_arr%y0d)
    !stop
bondbased: if(parini%bondbased_ann) then
    call symmetry_functions_driver_bond(parini,ann_arr,atoms,symfunc)
        do ib=1,symfunc%linked_lists%maxbound_rad
        !do iat=1,atoms%nat
        !i=atoms%itypat(iat)
        !do jat=1,atoms%nat
            if(apply_gbounds) then
            do i0=1,ann_arr%ann(1)%nn(0)
                gleft=ann_arr%ann(1)%gbounds(1,i0)
                !normalization of y
                !ann_arr%yall_bond(i0,iat,jat)=(ann_arr%yall_bond(i0,iat,jat)-gleft)*ann_arr%ann(1)%two_over_gdiff(i0)-1.d0
                symfunc%y(i0,ib)=(symfunc%y(i0,ib)-gleft)*ann_arr%ann(1)%two_over_gdiff(i0)-1.d0
                !normalization of y0d
                !ann_arr%y0d_bond(i0,1:3,1:atoms%nat,1:atoms%nat)=ann_arr%y0d_bond(i0,1:3,1:atoms%nat,1:atoms%nat)*ann_arr%ann(1)%two_over_gdiff(i0)
                symfunc%y0d_bond(i0,ib)=symfunc%y0d_bond(i0,ib)*ann_arr%ann(1)%two_over_gdiff(i0)
            enddo
            endif
            if(ann_arr%ann(1)%nn(0)/=ann_arr%ann(1)%ng1+ann_arr%ann(1)%ng2+ann_arr%ann(1)%ng3+ann_arr%ann(1)%ng4) then
                write(*,'(a,2i)') 'ERROR: inconsistency between # of input nodes and # of symmetry functions:', &
                ann_arr%ann(1)%nn(0),ann_arr%ann(1)%ng1+ann_arr%ann(1)%ng2+ann_arr%ann(1)%ng3+ann_arr%ann(1)%ng4
                stop
            endif
        !enddo
        enddo
else bondbased
    !allocate(ann_arr%yall(ann_arr%ann(1)%nn(0),atoms%nat),stat=istat,source=0.d0)
    !if(istat/=0) stop 'ERROR: unable to allocate array ann_arr%yall.'
    !do iat=1,atoms%nat
    !    isat=atoms%itypat(iat)
    !    ann_arr%yall(1:ann_arr%ann(isat)%nn(0),iat)=0.d0
    !enddo
    !do ib=1,atoms%maxbound_rad
    !    ann_arr%y0dr(1:ann_arr%ann(isat)%nn(0),1:9,ib)=0.d0
    !    ann_arr%y0d(1:ann_arr%ann(isat)%nn(0),1,ib)=0.d0
    !    ann_arr%y0d(1:ann_arr%ann(isat)%nn(0),2,ib)=0.d0
    !    ann_arr%y0d(1:ann_arr%ann(isat)%nn(0),3,ib)=0.d0
    !enddo
    if (parini%symfunc_type_ann=='behler') then 
        call symmetry_functions_driver(parini,ann_arr,atoms,symfunc)
    else
        call symmetry_functions_driver_stefan(parini,ann_arr,atoms,symfunc)
    endif
    !i0=0
    !call symmetry_functions_g02(ann,iat,atoms,i0)
    !call symmetry_functions_g04(ann,iat,atoms,i0)
    !call symmetry_functions_g05(ann,iat,atoms,i0)
    !call symmetry_functions_g06(ann,iat,atoms,i0)
    if(apply_gbounds) then
    do iat=1,atoms%nat
    isat=atoms%itypat(iat)
        do i0=1,ann_arr%ann(isat)%nn(0)
            gleft=ann_arr%ann(isat)%gbounds(1,i0)
            !normalization of y
            symfunc%y(i0,iat)=(symfunc%y(i0,iat)-gleft)*ann_arr%ann(isat)%two_over_gdiff(i0)-1.d0
        enddo
    enddo
    do ib=1,symfunc%linked_lists%maxbound_rad
    iat=symfunc%linked_lists%bound_rad(1,ib)
    isat=atoms%itypat(iat)
        do i0=1,ann_arr%ann(isat)%nn(0)
            gleft=ann_arr%ann(isat)%gbounds(1,i0)
            !normalization of y0d
            symfunc%y0d(i0,1:3,ib)=symfunc%y0d(i0,1:3,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
            !normalization of y0dr
            symfunc%y0dr(i0,1:9,ib)=symfunc%y0dr(i0,1:9,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
        enddo
    enddo
    endif
    !---------------------------------------------
    do iat=1,atoms%nat
    isat=atoms%itypat(iat)
    ng=ann_arr%ann(isat)%nn(0)
    if(ng/=ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ann_arr%ann(isat)%ng3+ann_arr%ann(isat)%ng4+ann_arr%ann(isat)%ng5+ann_arr%ann(isat)%ng6) then
        write(*,'(a,2i)') &
            'ERROR: inconsistency between # of input nodes and # of symmetry functions:', &
            ann_arr%ann(isat)%nn(0),ann_arr%ann(isat)%ng1+ann_arr%ann(isat)%ng2+ann_arr%ann(isat)%ng3+ann_arr%ann(isat)%ng4+ann_arr%ann(isat)%ng5+ann_arr%ann(isat)%ng6
        stop
    endif
    !symfunc%y(1:ng,iat)=ann_arr%yall(1:ng,iat)
    enddo
    !associate(nat=>atoms%nat)
    !associate(nb=>symfunc%linked_lists%maxbound_rad)
    !do iat=1,atoms%nat
    !    do ig=1,ng
    !        symfunc%y(ig,iat)=ann_arr%yall(ig,iat)
    !    enddo
    !enddo
    !do ib=1,nb
    !    do i=1,3
    !        do ig=1,ng
    !            symfunc%y0d(ig,i,ib)=ann_arr%y0d(ig,i,ib)
    !        enddo
    !    enddo
    !    do i=1,9
    !        do ig=1,ng
    !            symfunc%y0dr(ig,i,ib)=ann_arr%y0dr(ig,i,ib)
    !        enddo
    !    enddo
    !enddo
    !end associate
    !end associate
endif bondbased
    call f_timing(TCAT_SYMFUNC_COMPUT,'OF')
    call f_release_routine()
end subroutine symmetry_functions
!*****************************************************************************************
