!*****************************************************************************************
module mod_train
    use mod_opt_ann, only: typ_opt_ann, typ_cost_object
    use mod_atoms, only: typ_atoms_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_ann, only: typ_ann_arr
    use mod_refdata, only: typ_refdata
    implicit none
    private
    public:: ann_train
contains
!*****************************************************************************************
subroutine ann_train(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann, ann_lm
    use mod_ann_io_yaml, only: write_ann_all_yaml
    use mod_processors, only: iproc
    use mod_train_chi, only: ekf_rivals_fitchi
    use mod_train_cent1, only: ekf_rivals_cent1
    use mod_train_cent2, only: ekf_rivals_cent2
    use mod_train_atombased, only: ekf_rivals_atombased
    use mod_train_dvec, only: ekf_rivals_dvec
    use mod_refdata, only: typ_refdata
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr):: ann_arr
    type(typ_opt_ann):: opt_ann
    type(typ_refdata):: refdata
    call f_routine(id='ann_train')
    !-------------------------------------------------------
    call init_ann_train(parini,ann_arr)
    if(trim(parini%optimizer_ann)=='rivals') then
        if(trim(parini%approach_ann)=='cent1') then
            call ekf_rivals_cent1(parini,ann_arr,opt_ann,refdata)
        elseif(trim(parini%approach_ann)=='cent2') then
            call ekf_rivals_cent2(parini,ann_arr,opt_ann,refdata)
        elseif(trim(parini%approach_ann)=='dvec') then
            call ekf_rivals_dvec(parini,ann_arr,opt_ann,refdata)
        elseif(trim(parini%approach_ann)=='atombased') then
            call ekf_rivals_atombased(parini,ann_arr,opt_ann,refdata)
        else
            stop 'ERROR: unknown approach_ann is used with ekf_rivals!'
        endif
    else if(trim(parini%optimizer_ann)=='rivals_fitchi') then
        call ekf_rivals_fitchi(parini,ann_arr,opt_ann,refdata)
    else if(trim(parini%optimizer_ann)=='lm') then
        stop 'ERROR: call to ann_lm needs some change due to the changes!'
        call ann_lm(parini,ann_arr,refdata,opt_ann)
    else
        write(*,*) 'ERROR: unknown optimzer in ANN training'
    endif
    if(iproc==0) then
        if( ann_arr%exists_yaml_file) then
            call write_ann_all_yaml(parini,ann_arr,-1)
        else
            call write_ann_all(parini,ann_arr,-1)
        endif
    endif
    call fini_ann_train(parini,ann_arr,opt_ann,refdata)
    call f_release_routine()
end subroutine ann_train
!*****************************************************************************************
subroutine init_ann_train(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_processors, only: iproc
    use yaml_output
    use futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    character(30):: fnout
    integer:: ierr
    !---------------------------------------------
    if(iproc==0) then
        !write(fnout,'(a12,i3.3)') 'err_train',iproc
        fnout="train_output.yaml"
        ann_arr%iunit=f_get_free_unit(10**5)
        call yaml_set_stream(unit=ann_arr%iunit,filename=trim(fnout),&
             record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
        if (ierr/=0) then
           call yaml_warning('Failed to create'//trim(fnout)//', error code='//trim(yaml_toa(ierr)))
        end if
        call yaml_release_document(ann_arr%iunit)
        call yaml_sequence_open('training iterations',unit=ann_arr%iunit)
    endif
end subroutine init_ann_train
!*****************************************************************************************
subroutine fini_ann_train(parini,ann_arr,opt_ann,refdata)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann
    use mod_atoms, only: typ_atoms_arr, atom_deallocate_old
    use mod_refdata, only: typ_refdata
    use mod_processors, only: iproc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_refdata), intent(inout):: refdata
    !local variables
    integer:: iconf
    if(iproc==0) then
        call yaml_close_stream(unit=ann_arr%iunit)
    endif

    call ann_arr%fini_ann_arr()
    call opt_ann%fini_opt_ann()

    do iconf=1,refdata%atoms_train%nconf
        call f_free(refdata%symfunc_train%symfunc(iconf)%y)
        call atom_deallocate_old(refdata%atoms_train%atoms(iconf))
    enddo

    do iconf=1,refdata%atoms_valid%nconf
        call f_free(refdata%symfunc_valid%symfunc(iconf)%y)
        call atom_deallocate_old(refdata%atoms_valid%atoms(iconf))
    enddo

    do iconf=1,refdata%atoms_train%nconf
        call refdata%symfunc_train%symfunc(iconf)%fini_symfunc()
    enddo
    deallocate(refdata%symfunc_train%symfunc)

    do iconf=1,refdata%atoms_valid%nconf
        call refdata%symfunc_valid%symfunc(iconf)%fini_symfunc()
    enddo
    deallocate(refdata%symfunc_valid%symfunc)

    deallocate(refdata%atoms_train%conf_inc)
    deallocate(refdata%atoms_valid%conf_inc)
end subroutine fini_ann_train
!*****************************************************************************************
subroutine apply_gbounds_bond(parini,ann_arr,atoms_arr,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    integer:: iconf, iat, jat, ig, ib, i0, isat
    real(8):: tt
    do iconf=1,atoms_arr%nconf
        if(atoms_arr%atoms(iconf)%ntypat>1) stop 'ERROR: this part not ready for ntypat>1'
        !if(symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
        do ib=1,symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad
            iat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(1,ib)
            jat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(2,ib)
            if(iat>jat) cycle
            do ig=1,symfunc_arr%symfunc(iconf)%ng
                tt=symfunc_arr%symfunc(iconf)%y(ig,ib)
                tt=(tt-ann_arr%ann(1)%gbounds(1,ig))*ann_arr%ann(1)%two_over_gdiff(ig)-1.d0
                symfunc_arr%symfunc(iconf)%y(ig,ib)=tt
            enddo
        enddo
        if(parini%save_symfunc_force_ann) then
            do ib=1,symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad
                iat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(1,ib)
                jat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(2,ib)
                if(iat>jat) cycle
                isat=atoms_arr%atoms(iconf)%itypat(iat)
                do i0=1,ann_arr%ann(isat)%nn(0)
                    !normalization of y0d
                    tt=ann_arr%ann(isat)%two_over_gdiff(i0)
                    symfunc_arr%symfunc(iconf)%y0d_bond(i0,ib)=symfunc_arr%symfunc(iconf)%y0d_bond(i0,ib)*tt
                    !symfunc_arr%symfunc(iconf)%y0d(i0,1,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,1,ib)*tt
                    !symfunc_arr%symfunc(iconf)%y0d(i0,2,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,2,ib)*tt
                    !symfunc_arr%symfunc(iconf)%y0d(i0,3,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,3,ib)*tt
                    !normalization of y0dr
                    !symfunc%y0dr(i0,1:9,ib)=symfunc%y0dr(i0,1:9,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                enddo
            enddo
        endif
    enddo
end subroutine apply_gbounds_bond
!*****************************************************************************************
end module mod_train
!*****************************************************************************************
