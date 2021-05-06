!*****************************************************************************************
module mod_callback_ann
    use mod_atoms, only: typ_atoms_arr
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_opt_ann, only: typ_opt_ann
    implicit none
    type(typ_parini), pointer:: parini_t
    type(typ_atoms_arr), pointer:: atoms_train_t
    type(typ_atoms_arr), pointer:: atoms_valid_t
    type(typ_atoms_arr), pointer:: atoms_smplx_t
    type(typ_ann_arr), pointer:: ann_arr_t
    type(typ_opt_ann), pointer:: opt_ann_t
    type(typ_symfunc_arr), pointer:: symfunc_train_t
    type(typ_symfunc_arr), pointer:: symfunc_valid_t
end module mod_callback_ann
!*****************************************************************************************
module mod_train
    implicit none
    private
    public:: ann_train
contains
!*****************************************************************************************
subroutine ann_train(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_opt_ann, only: typ_opt_ann, ekf_rivals, ekf_behler, ann_lm
    use mod_opt_ann, only: set_annweights
    use mod_atoms, only: typ_atoms_arr, typ_atoms
    use mod_processors, only: iproc
    use mod_callback_ann
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr):: ann_arr
    type(typ_opt_ann):: opt_ann
    type(typ_atoms_arr), target:: atoms_train
    type(typ_atoms_arr), target:: atoms_valid
    type(typ_atoms_arr):: atoms_smplx
    type(typ_symfunc_arr), target:: symfunc_train
    type(typ_symfunc_arr), target:: symfunc_valid
    real(8):: time1, time2, time3
    logical:: file_exists
    call f_routine(id='ann_train')
    !-------------------------------------------------------
    !Reading configurations and their energies and forces
    inquire(file="list_posinp_train.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_train.yaml',atoms_train)
    else
        call read_data_old(parini,'list_posinp_train',atoms_train)
    endif
    inquire(file="list_posinp_valid.yaml",exist=file_exists)
    if(file_exists) then
        call read_data_yaml(parini,'list_posinp_valid.yaml',atoms_valid)
    else
        call read_data_old(parini,'list_posinp_valid',atoms_valid)
    endif
    if(trim(parini%approach_ann)=='centt') then
        call read_data_yaml(parini,'list_posinp_smplx.yaml',atoms_smplx)
    endif
    call init_ann_train(parini,ann_arr,opt_ann,atoms_train,atoms_valid)
    !-------------------------------------------------------
    if(iproc==0) then
        call yaml_map('number of ANN wights',opt_ann%n)
        call yaml_map('number of training data points',atoms_train%nconf)
        call yaml_map('number of validating data points',atoms_valid%nconf)
        !write(*,'(a,i)') 'number of ANN wights:             ',opt_ann%n
        !write(*,'(a,i)') 'number of training data points:   ',atoms_train%nconf
        !write(*,'(a,i)') 'number of validating data points: ',atoms_valid%nconf
    endif
    call set_conf_inc_random(parini,atoms_train)
    call set_conf_inc_random(parini,atoms_valid)
    call prepare_atoms_arr(parini,ann_arr,atoms_train)
    call prepare_atoms_arr(parini,ann_arr,atoms_valid)
    call prepare_atoms_arr(parini,ann_arr,atoms_smplx)
    !allocate(atoms_train%inclusion(atoms_train%nconf),source=0)
    !-------------------------------------------------------
    call cpu_time(time1)
    call set_gbounds(parini,ann_arr,atoms_valid,'bounds_valid',symfunc_valid)
    call cpu_time(time2)
    call set_gbounds(parini,ann_arr,atoms_train,'bounds_train',symfunc_train)
    call cpu_time(time3)
    !write(*,'(a,2f10.1)') 'TIMING: evaluation symmetry functions: ',time2-time1,time3-time2
    call yaml_map('TIMING evaluation symmetry functions for valid',time2-time1)
    call yaml_map('TIMING evaluation symmetry functions for train',time3-time2)
    !-------------------------------------------------------------------------------------
    !IMPORTANT: The following must be done after set_gbounds is called for training set.
    !if(trim(parini%symfunc)/='do_not_save') then
    if(parini%bondbased_ann) then
        call apply_gbounds_bond(parini,ann_arr,atoms_valid,symfunc_valid)
        call apply_gbounds_bond(parini,ann_arr,atoms_train,symfunc_train)
    else
        call apply_gbounds_atom(parini,ann_arr,atoms_valid,symfunc_valid)
        call apply_gbounds_atom(parini,ann_arr,atoms_train,symfunc_train)
    endif
    !endif
    !-------------------------------------------------------------------------------------
    call set_annweights(parini,opt_ann,ann_arr)
    if(trim(parini%approach_ann)=='centt') then
        call set_single_atom_energy(parini,ann_arr,opt_ann)
    endif

    ann_arr%compute_symfunc=.false.
    !if(parini%prefit_ann .and. trim(parini%approach_ann)=='centt') then
    if(parini%prefit_ann ) then
        !call prefit_cent_ener_ref(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
        call prefit_cent(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    endif

    if(trim(parini%approach_ann)=='centt' .and. parini%prefit_centt_ann) then
        call centt_simplex(parini,ann_arr,atoms_smplx,opt_ann)
    endif
    atoms_train_t=>atoms_train
    atoms_valid_t=>atoms_valid
    symfunc_train_t=>symfunc_train
    symfunc_valid_t=>symfunc_valid
    if(trim(parini%optimizer_ann)=='behler') then
        call ekf_behler(parini,ann_arr,opt_ann)
    else if(trim(parini%optimizer_ann)=='rivals') then
        call ekf_rivals(parini,ann_arr,opt_ann)
    !else if(trim(parini%optimizer_ann)=='rivals_tmp') then
    !    call ekf_rivals_tmp(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    else if(trim(parini%optimizer_ann)=='lm') then
        call ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,opt_ann)
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
    call fini_ann_train(parini,ann_arr,opt_ann,atoms_train,atoms_valid,symfunc_train,symfunc_valid)

    call f_release_routine()
end subroutine ann_train
!*****************************************************************************************
subroutine set_single_atom_energy(parini,ann_arr,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann, convert_opt_x_ann_arr
    use mod_atoms, only: typ_atoms, atom_allocate_old, atom_deallocate_old, set_rat_iat
    use mod_symfunc, only: typ_symfunc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    integer:: ityp
    real(8):: t_ener_ref
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
    !call atom_copy_old(atoms_train%atoms(1),atoms,'atoms_arr%atoms(iconf)->atoms')
    call atom_allocate_old(atoms,1,0,0)
    !write(*,*) atoms%nat
    !write(*,*) atoms%sat(:)
    !write(*,*) atoms%itypat(:)
    call set_rat_iat(atoms,1,(/1.d0,1.d0,1.d0/))
    atoms%cellvec(1:3,1:3)=0.d0
    atoms%cellvec(1,1)=10.d0
    atoms%cellvec(2,2)=10.d0
    atoms%cellvec(3,3)=10.d0
    atoms%boundcond='free'
    ann_arr%event='potential'
    ann_arr%compute_symfunc=.true.
    !call yaml_sequence_open('total charge of single atom in adjusting ener_ref')
    do ityp=1,parini%ntypat
        !call yaml_sequence(advance='no')
        atoms%sat(1)=parini%stypat(ityp)
        atoms%itypat(1)=parini%ltypat(ityp)
        t_ener_ref=ann_arr%ann(atoms%itypat(1))%ener_ref
        ann_arr%ann(atoms%itypat(1))%ener_ref=0.d0
        call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
        ann_arr%ann(atoms%itypat(1))%ener_ref=t_ener_ref-atoms%epot
        call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
        !call yaml_map('type',trim(atoms%sat(1)))
        !call yaml_map('charge',atoms%zat(1)+atoms%qat(1))
        !write(*,'(a,f)') 'Adjusting ener_ref: total charge= ',atoms%zat(1)+atoms%qat(1)
    enddo
    !call yaml_sequence_close()
    call atom_deallocate_old(atoms)
end subroutine set_single_atom_energy
!*****************************************************************************************
subroutine centt_simplex(parini,ann_arr,atoms_smplx,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann
    use mod_atoms, only: typ_atoms_arr
    use mod_callback_ann
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in), target:: parini
    type(typ_ann_arr), intent(inout), target:: ann_arr
    type(typ_opt_ann), intent(inout), target:: opt_ann
    type(typ_atoms_arr), intent(inout), target:: atoms_smplx
    !local variables
    real(8):: vertices(10,11), fval(11)
    real(8):: step, ftol
    integer:: ndim, iter, i
    !external:: cal_rmse_force_centt
    !external:: cal_rmse_energy_centt
    ndim=10
    ftol=parini%ftol_ann
    step=0.d0
    !vertices(1,1)=ann_arr%ann(2)%chi0-ann_arr%ann(1)%chi0
    vertices(1,1)=ann_arr%ann(1)%chi0
    vertices(2,1)=ann_arr%ann(2)%chi0

    !vertices(2,1)=ann_arr%ann(1)%zion
    !vertices(3,1)=ann_arr%ann(2)%zion

    vertices(3,1)=ann_arr%ann(1)%hardness
    vertices(4,1)=ann_arr%ann(2)%hardness

    vertices(5,1)=ann_arr%ann(1)%spring_const
    vertices(6,1)=ann_arr%ann(2)%spring_const

    vertices(7,1)=ann_arr%ann(1)%gausswidth_ion
    vertices(8,1)=ann_arr%ann(2)%gausswidth_ion

    vertices(9,1)=ann_arr%ann(1)%ener_ref
    vertices(10,1)=ann_arr%ann(2)%ener_ref

    do i=2,ndim-1
        vertices(1:ndim,i)=vertices(1:ndim,1)
        vertices(i-1,i)=vertices(i-1,i)+1.d-2
    enddo
    vertices(1:ndim,ndim)=vertices(1:ndim,1)
    vertices(ndim-1,ndim)=vertices(ndim-1,ndim)+1.d-3
    vertices(1:ndim,ndim+1)=vertices(1:ndim,1)
    vertices(ndim,ndim+1)=vertices(ndim,ndim+1)+1.d-3

    atoms_smplx_t=>atoms_smplx
    parini_t=>parini
    ann_arr_t=>ann_arr
    opt_ann_t=>opt_ann
    !call simplex(vertices,fval,step,ndim,ftol,cal_rmse_force_centt,iter)
    call simplex(vertices,fval,step,ndim,ftol,cal_rmse_energy_centt,iter)
end subroutine centt_simplex
!*****************************************************************************************
subroutine cal_rmse_force_centt(ndim,vertex,rmse_force_centt)
    use mod_callback_ann, only: atoms_smplx=>atoms_smplx_t, parini=>parini_t
    use mod_callback_ann, only: ann_arr=>ann_arr_t, opt_ann=>opt_ann_t
    use mod_atoms, only: typ_atoms, atom_copy_old 
    use mod_symfunc, only: typ_symfunc
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: vertex(ndim)
    real(8), intent(out) :: rmse_force_centt
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    real(8):: rmse
    integer:: iat, iconf, nat_tot
    ann_arr%ann(1)%chi0=-vertex(1)/2.d0
    ann_arr%ann(2)%chi0= vertex(1)/2.d0

    ann_arr%ann(1)%zion=vertex(2)
    ann_arr%ann(2)%zion=vertex(3)

    ann_arr%ann(1)%hardness=vertex(4)
    ann_arr%ann(2)%hardness=vertex(5)

    ann_arr%ann(1)%spring_const=vertex(6)
    ann_arr%ann(2)%spring_const=vertex(7)

    ann_arr%ann(1)%gausswidth_ion=vertex(8)
    ann_arr%ann(2)%gausswidth_ion=vertex(9)

    write(*,'(5(a,2f7.4,4x))') &
        'CHI= ',ann_arr%ann(1)%chi0,ann_arr%ann(2)%chi0, &
        'ZION= ',ann_arr%ann(1)%zion,ann_arr%ann(2)%zion, &
        'HARD= ',ann_arr%ann(1)%hardness,ann_arr%ann(2)%hardness, &
        'K= ',ann_arr%ann(1)%spring_const,ann_arr%ann(2)%spring_const, &
        'BETA= ',ann_arr%ann(1)%gausswidth_ion,ann_arr%ann(2)%gausswidth_ion

    nat_tot=0
    rmse=0.d0
    ann_arr%event='potential'
    ann_arr%compute_symfunc=.true.
    do iconf=1,atoms_smplx%nconf
        call atom_copy_old(atoms_smplx%atoms(iconf),atoms,'atoms_smplx%atoms(iconf)->atoms')
        call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
        do iat=1,atoms%nat
            rmse=rmse+(atoms%fat(1,iat)-atoms_smplx%atoms(iconf)%fat(1,iat))**2 &
                     +(atoms%fat(2,iat)-atoms_smplx%atoms(iconf)%fat(2,iat))**2 &
                     +(atoms%fat(3,iat)-atoms_smplx%atoms(iconf)%fat(3,iat))**2
        enddo
        nat_tot=nat_tot+atoms%nat
    enddo
    rmse_force_centt=sqrt(rmse/real(3*nat_tot,8))
    write(*,*) 'rmse_force_centt ',rmse_force_centt
end subroutine cal_rmse_force_centt
!*****************************************************************************************
subroutine cal_rmse_energy_centt(ndim,vertex,rmse_energy_centt)
    use mod_callback_ann, only: atoms_smplx=>atoms_smplx_t, parini=>parini_t
    use mod_callback_ann, only: ann_arr=>ann_arr_t, opt_ann=>opt_ann_t
    use mod_atoms, only: typ_atoms, atom_copy_old
    use mod_symfunc, only: typ_symfunc
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: vertex(ndim)
    real(8), intent(out) :: rmse_energy_centt
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    real(8):: rmse
    integer:: iat, iconf
    !ann_arr%ann(1)%chi0=-vertex(1)/2.d0
    !ann_arr%ann(2)%chi0= vertex(1)/2.d0
    ann_arr%ann(1)%chi0=vertex(1)
    ann_arr%ann(2)%chi0=vertex(2)

    !ann_arr%ann(1)%zion=vertex(2)
    !ann_arr%ann(2)%zion=vertex(3)

    ann_arr%ann(1)%hardness=vertex(3)
    ann_arr%ann(2)%hardness=vertex(4)

    ann_arr%ann(1)%spring_const=vertex(5)
    ann_arr%ann(2)%spring_const=vertex(6)

    ann_arr%ann(1)%gausswidth_ion=vertex(7)
    ann_arr%ann(2)%gausswidth_ion=vertex(8)

    ann_arr%ann(1)%ener_ref=vertex(9)
    ann_arr%ann(2)%ener_ref=vertex(10)

    write(*,'(5(a,2f7.4,4x),a,2f10.4)') &
        'CHI= ',ann_arr%ann(1)%chi0,ann_arr%ann(2)%chi0, &
        'ZION= ',ann_arr%ann(1)%zion,ann_arr%ann(2)%zion, &
        'HARD= ',ann_arr%ann(1)%hardness,ann_arr%ann(2)%hardness, &
        'K= ',ann_arr%ann(1)%spring_const,ann_arr%ann(2)%spring_const, &
        'BETA= ',ann_arr%ann(1)%gausswidth_ion,ann_arr%ann(2)%gausswidth_ion, &
        'E0= ',ann_arr%ann(1)%ener_ref,ann_arr%ann(2)%ener_ref

    !call set_single_atom_energy(parini,ann_arr,opt_ann)
    rmse=0.d0
    ann_arr%event='potential'
    ann_arr%compute_symfunc=.true.
    do iconf=1,atoms_smplx%nconf
        call atom_copy_old(atoms_smplx%atoms(iconf),atoms,'atoms_smplx%atoms(iconf)->atoms')
        call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
        do iat=1,atoms%nat
            rmse=rmse+(atoms%epot-atoms_smplx%atoms(iconf)%epot)**2
        enddo
    enddo
    rmse_energy_centt=sqrt(rmse/real(atoms_smplx%nconf,8))
    write(*,*) 'rmse_energy_centt ',rmse_energy_centt
end subroutine cal_rmse_energy_centt
!*****************************************************************************************
subroutine init_ann_train(parini,ann_arr,opt_ann,atoms_train,atoms_valid)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, set_number_of_ann, init_ann_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_opt_ann, only: typ_opt_ann, init_opt_ann
    use mod_processors, only: iproc
    use yaml_output
    use futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_atoms_arr), intent(inout):: atoms_train, atoms_valid
    !local variables
    character(30):: fnout
    character (50)::fname
    integer:: ierr
    call set_number_of_ann(parini,ann_arr)
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in ann_train'
    call yaml_map('number of ann',ann_arr%nann)
    !write(*,*) 'Here', ann_arr%nann
    allocate(ann_arr%ann(ann_arr%nann))
    ann_arr%approach=trim(parini%approach_ann)
    fname = trim(parini%stypat(1))//'.ann.input.yaml'
    inquire(file=trim(fname),exist=ann_arr%exists_yaml_file)
    if( ann_arr%exists_yaml_file) then
        if (parini%restart_param) then
            call read_ann_yaml(parini,ann_arr)
        else
            call read_input_ann_yaml(parini,iproc,ann_arr)
        endif
    else
        call read_input_ann(parini,iproc,ann_arr)
    endif
    !---------------------------------------------
    call init_ann_arr(ann_arr)
    if(trim(ann_arr%approach)=='cent3') then
        call init_opt_ann(3*atoms_train%nconf,opt_ann,ann_arr)
    else
        call init_opt_ann(atoms_train%nconf,opt_ann,ann_arr)
    endif
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
subroutine fini_ann_train(parini,ann_arr,opt_ann,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, fini_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_opt_ann, only: typ_opt_ann, fini_opt_ann
    use mod_atoms, only: typ_atoms_arr
    use mod_processors, only: iproc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_atoms_arr), intent(inout):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    !local variables
    integer:: iconf
    !call f_free(opt_ann%x)

    if(iproc==0) then
        !close(11)
        !close(12)
        call yaml_close_stream(unit=ann_arr%iunit)
    endif

    call fini_ann_arr(ann_arr)
    call fini_opt_ann(opt_ann)

    do iconf=1,atoms_train%nconf
        call f_free(symfunc_train%symfunc(iconf)%y)
    enddo

    do iconf=1,atoms_valid%nconf
        call f_free(symfunc_valid%symfunc(iconf)%y)
    enddo

    deallocate(atoms_train%conf_inc)
    deallocate(atoms_valid%conf_inc)
    !deallocate(atoms_train%inclusion)
end subroutine fini_ann_train
!*****************************************************************************************
subroutine set_conf_inc_random(parini,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_utils
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: iconf, irand
    real(8):: tt
    if(parini%nconf_rmse==0) then
        write(*,*) 'ERROR: parini%nconf_rmse=0'
        stop
    endif
    if(parini%nconf_rmse>=atoms_arr%nconf) then
        allocate(atoms_arr%conf_inc(atoms_arr%nconf),source=.true.)
        atoms_arr%nconf_inc=atoms_arr%nconf
        return
    endif
    allocate(atoms_arr%conf_inc(atoms_arr%nconf),source=.false.)
    atoms_arr%nconf_inc=parini%nconf_rmse
    irand=0
    do
        if(irand==atoms_arr%nconf_inc) exit
        if(trim(parini%rng_type)=='only_for_tests') then
            call random_number_generator_simple(tt)
        else
            call random_number(tt)
        endif
        tt=tt*real(atoms_arr%nconf)
        iconf=int(tt)+1
        if(atoms_arr%conf_inc(iconf)) cycle
        atoms_arr%conf_inc(iconf)=.true.
        irand=irand+1
    enddo
end subroutine set_conf_inc_random
!*****************************************************************************************
subroutine apply_gbounds_atom(parini,ann_arr,atoms_arr,symfunc_arr)
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
    integer:: iconf, iat, ig, ib, i, i0, isat
    real(8):: tt
    do iconf=1,atoms_arr%nconf
        do iat=1,atoms_arr%atoms(iconf)%nat
            i=atoms_arr%atoms(iconf)%itypat(iat)
            do ig=1,symfunc_arr%symfunc(iconf)%ng
                tt=symfunc_arr%symfunc(iconf)%y(ig,iat)
                tt=(tt-ann_arr%ann(i)%gbounds(1,ig))*ann_arr%ann(i)%two_over_gdiff(ig)-1.d0
                symfunc_arr%symfunc(iconf)%y(ig,iat)=tt
            enddo
        enddo
        if(parini%save_symfunc_force_ann) then
            do ib=1,symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad
                iat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(1,ib)
                isat=atoms_arr%atoms(iconf)%itypat(iat)
                do i0=1,ann_arr%ann(isat)%nn(0)
                    !normalization of y0d
                    tt=ann_arr%ann(isat)%two_over_gdiff(i0)
                    symfunc_arr%symfunc(iconf)%y0d(i0,1,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,1,ib)*tt
                    symfunc_arr%symfunc(iconf)%y0d(i0,2,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,2,ib)*tt
                    symfunc_arr%symfunc(iconf)%y0d(i0,3,ib)=symfunc_arr%symfunc(iconf)%y0d(i0,3,ib)*tt
                    !normalization of y0dr
                    !symfunc%y0dr(i0,1:9,ib)=symfunc%y0dr(i0,1:9,ib)*ann_arr%ann(isat)%two_over_gdiff(i0)
                enddo
            enddo
        endif
    enddo
end subroutine apply_gbounds_atom
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
subroutine prepare_atoms_arr(parini,ann_arr,atoms_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms_arr, update_ratp, update_rat
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables
    integer:: iconf, iat, i
    real(8), allocatable:: ratred(:,:)
    do iconf=1,atoms_arr%nconf
        if(trim(atoms_arr%atoms(iconf)%boundcond)=='bulk') then
        allocate(ratred(1:3,1:atoms_arr%atoms(iconf)%nat))
        call update_ratp(atoms_arr%atoms(iconf))
        call rxyz_cart2int_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,atoms_arr%atoms(iconf)%ratp,ratred)
        call backtocell_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,ratred)
        call rxyz_int2cart_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,ratred,atoms_arr%atoms(iconf)%ratp)
        call update_rat(atoms_arr%atoms(iconf),upall=.true.)
        deallocate(ratred)
        endif
        do iat=1,atoms_arr%atoms(iconf)%nat
            do i=1,ann_arr%nann
                if(trim(atoms_arr%atoms(iconf)%sat(iat))==trim(parini%stypat(i))) then
                    atoms_arr%atoms(iconf)%itypat(iat)=parini%ltypat(i)
                    exit
                endif
            enddo
        enddo
    enddo
end subroutine prepare_atoms_arr
!*****************************************************************************************
subroutine set_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables 
    integer:: iconf
    if(.not. allocated(symfunc_arr%symfunc)) then
        symfunc_arr%nconf=atoms_arr%nconf
        allocate(symfunc_arr%symfunc(symfunc_arr%nconf))
    endif
    !write(*,'(a,i3,i6)') 'iproc,nconf ',iproc,atoms_arr%nconf
    do iconf=1,atoms_arr%nconf
        symfunc_arr%symfunc(iconf)%ng=ann_arr%ann(1)%nn(0) !HERE
        symfunc_arr%symfunc(iconf)%nat=atoms_arr%atoms(iconf)%nat
    enddo
    configuration: do iconf=1+iproc,atoms_arr%nconf,nproc
        if(trim(parini%symfunc)/='read') then
            call symmetry_functions(parini,ann_arr,atoms_arr%atoms(iconf),symfunc_arr%symfunc(iconf),.false.)
            if(.not. parini%save_symfunc_force_ann) then
                call f_free(symfunc_arr%symfunc(iconf)%y0d)
            endif
            call f_free(symfunc_arr%symfunc(iconf)%y0dr)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%prime_bound)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_rad)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_ang)
        !elseif(trim(parini%symfunc)/='read') then
        !    stop 'ERROR: arini%symfunc contains none of the three acceptable possibilies'
        endif
        if(trim(parini%symfunc)=='write') then
            call write_symfunc(parini,iconf,atoms_arr,strmess,symfunc_arr)
        elseif(trim(parini%symfunc)=='read') then
            call read_symfunc(parini,iconf,ann_arr,atoms_arr,strmess,symfunc_arr)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%prime_bound)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_rad)
            deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_ang)
        !elseif(trim(parini%symfunc)=='do_not_save') then
        !    call f_free(symfunc_arr%symfunc(iconf)%y0d)
        !    call f_free(symfunc_arr%symfunc(iconf)%y0dr)
        !    deallocate(symfunc_arr%symfunc(iconf)%linked_lists%prime_bound)
        !    deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_ang)
        endif
    enddo configuration
    call save_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
end subroutine set_gbounds
!*****************************************************************************************
subroutine write_symfunc(parini,iconf,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr, update_ratp
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iconf
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    real(8), allocatable:: wa(:)
    integer:: nwa
    character(30):: filename
    integer:: i, ig, iat, ios, n, ib
    !Symmetry functions are written into files to be used for
    !subsequent training runs.
    if(trim(strmess)=='bounds_train') then
        write(filename,'(a25,i5.5)') '../symfunc/train.symfunc.',iconf
    elseif(trim(strmess)=='bounds_valid') then
        write(filename,'(a25,i5.5)') '../symfunc/valid.symfunc.',iconf
    else
        stop 'ERROR: invalid content in strmess in gset_bounds '
    endif
    open(unit=311,file=trim(filename),status='replace',form='unformatted', &
        access='stream',iostat=ios)
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning file ',trim(filename)
        stop
    endif
    associate(nb=>symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad)
    associate(ng=>symfunc_arr%symfunc(iconf)%ng)
    associate(nat=>atoms_arr%atoms(iconf)%nat)
    if(parini%save_symfunc_force_ann) then
        nwa=3+nat*(3+ng)+ng*3*nb
    else
        nwa=3+nat*(3+ng)
    endif
    allocate(wa(1:nwa))
    wa(1)=real(nat,8)
    wa(2)=real(ng,8)
    wa(3)=real(nb,8)
    call update_ratp(atoms_arr%atoms(iconf))
    do iat=1,atoms_arr%atoms(iconf)%nat
        wa(3+iat*3-2)=atoms_arr%atoms(iconf)%ratp(1,iat)
        wa(3+iat*3-1)=atoms_arr%atoms(iconf)%ratp(2,iat)
        wa(3+iat*3-0)=atoms_arr%atoms(iconf)%ratp(3,iat)
    enddo
    n=3+3*atoms_arr%atoms(iconf)%nat
    bondbased: if(parini%bondbased_ann) then
        !-------------------------- bond symmetry functions ---------------------------
        write(*,*) 'ERROR: writing/reading symmetry function values from files not'
        write(*,*) 'working yet, for several reasons for example allocation of wa'
        stop
        !----------------------------------------------------------------------------
    else
        do iat=1,atoms_arr%atoms(iconf)%nat
            do ig=1,symfunc_arr%symfunc(iconf)%ng
                n=n+1
                wa(n)=symfunc_arr%symfunc(iconf)%y(ig,iat)
            enddo
        enddo
        if(parini%save_symfunc_force_ann) then
            do ib=1,nb
                do i=1,3
                    do ig=1,symfunc_arr%symfunc(iconf)%ng
                        n=n+1
                        wa(n)=symfunc_arr%symfunc(iconf)%y0d(ig,i,ib)
                    enddo
                enddo
            enddo
        !else
        !    call f_free(symfunc_arr%symfunc(iconf)%y0d)
        endif
    endif bondbased
    end associate
    end associate
    end associate
    
    write(311,iostat=ios) wa
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure writing to file ',trim(filename)
        stop
    endif
    !write(311,'(2i6)') atoms_arr%atoms(iconf)%nat,symfunc_arr%symfunc(iconf)%ng
    !do iat=1,atoms_arr%atoms(iconf)%nat
    !    do ig=1,symfunc_arr%symfunc(iconf)%ng
    !        write(311,'(es25.16)') symfunc_arr%symfunc(iconf)%y(ig,iat)
    !    enddo
    !enddo
    close(311)
    deallocate(wa)
end subroutine write_symfunc
!*****************************************************************************************
subroutine read_symfunc(parini,iconf,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr, update_ratp
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iconf
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    real(8), allocatable:: wa(:)
    integer:: nwa
    character(30):: filename
    integer:: i, ig, iat, ios, nat_t, ng_t, nb_t, n, ib
    type(typ_pia_arr):: pia_arr_tmp
    real(8):: ttx, tty, ttz
    real(8):: eps=epsilon(1.d0)
    character(100):: smsg
    !Symmetry functions which are previously calculated and written by
    !some other run is going to be read from files
    if(trim(strmess)=='bounds_train') then
        write(filename,'(a25,i5.5)') '../symfunc/train.symfunc.',iconf
    elseif(trim(strmess)=='bounds_valid') then
        write(filename,'(a25,i5.5)') '../symfunc/valid.symfunc.',iconf
    else
        stop 'ERROR: invalid content in strmess in gset_bounds '
    endif
    !open(unit=311,file=trim(filename),status='old',iostat=ios)
    open(unit=311,file=trim(filename),status='old',form='unformatted', &
        access='stream',iostat=ios)
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning file ',trim(filename)
        stop
    endif
    !read(311,*) nat_t,ng_t
    associate(nb=>symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad)
    associate(ng=>symfunc_arr%symfunc(iconf)%ng)
    associate(nat=>atoms_arr%atoms(iconf)%nat)
    symfunc_arr%symfunc(iconf)%linked_lists%rcut=ann_arr%rcut
    symfunc_arr%symfunc(iconf)%linked_lists%triplex=.true.
    call call_linkedlist(parini,atoms_arr%atoms(iconf),.true.,symfunc_arr%symfunc(iconf)%linked_lists,pia_arr_tmp)
    deallocate(pia_arr_tmp%pia)
    symfunc_arr%symfunc(iconf)%y=f_malloc0((/1.to.ng,1.to.nat/),id='symfunc%y')
    if(parini%save_symfunc_force_ann) then
        symfunc_arr%symfunc(iconf)%y0d=f_malloc0((/1.to.ng,1.to.3,1.to.nb/),id='symfunc%y0d')
        nwa=3+nat*(3+ng)+ng*3*nb
    else
        nwa=3+nat*(3+ng)
    endif
    wa=f_malloc([1.to.nwa],id='wa')
    read(311,iostat=ios) wa
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure reading from file ',trim(filename)
        stop
    endif
    nat_t=nint(wa(1))
    ng_t=nint(wa(2))
    nb_t=nint(wa(3))
    if(parini%save_symfunc_behnam) then
    eps=eps*1.d4
    if(nat_t/=nat .or. ng_t/=ng) then
        write(*,'(a,7i6)') 'ERROR: inconsistent nat or ng ',iconf,nat_t,nat,ng_t,ng,nb_t,nb
        stop
    endif
    else
    if(nat_t/=nat .or. ng_t/=ng .or. nb_t/=nb) then
        write(*,'(a,7i6)') 'ERROR: inconsistent nat or ng ',iconf,nat_t,nat,ng_t,ng,nb_t,nb
        stop
    endif
    endif
    call update_ratp(atoms_arr%atoms(iconf))
    do iat=1,nat
        ttx=abs(wa(3+iat*3-2)-atoms_arr%atoms(iconf)%ratp(1,iat))
        tty=abs(wa(3+iat*3-1)-atoms_arr%atoms(iconf)%ratp(2,iat))
        ttz=abs(wa(3+iat*3-0)-atoms_arr%atoms(iconf)%ratp(3,iat))
        if(ttx>eps .or. tty>eps .or. ttz>eps) then
            smsg='ERROR: inconsistency of configuration in symmetry functions file. '
            !write(*,*) wa(3+iat*3-2),wa(3+iat*3-1),wa(3+iat*3-0)
            !write(*,*) atoms_arr%atoms(iconf)%rat(1,iat),atoms_arr%atoms(iconf)%rat(2,iat),atoms_arr%atoms(iconf)%rat(3,iat)
            !write(*,*) 'IAT',iat
            write(*,'(a,3es14.5,i7,a,i5)') trim(smsg),ttx,tty,ttz,iconf,trim(atoms_arr%fn(iconf)),atoms_arr%lconf(iconf)
            stop
        endif
    enddo
    n=3+3*atoms_arr%atoms(iconf)%nat
    bondbased: if(parini%bondbased_ann) then
        !--------------------- bond symmetry functions -------------------------------
        write(*,*) 'ERROR: writing/reading symmetry function values from files not'
        write(*,*) 'working yet, for several reasons for example allocation of wa'
        stop
        !do iat=1,atoms_arr%atoms(iconf)%nat
        !do jat=1,atoms_arr%atoms(iconf)%nat
        !    do ig=1,symfunc_arr%symfunc(iconf)%ng
        !        n=n+1
        !        symfunc_arr%symfunc(iconf)%y_bond(ig,iat,jat)=wa(n)
        !    enddo
        !enddo
        !enddo
        !-----------------------------------------------------------------------------
    else
        do iat=1,nat
            do ig=1,ng
                n=n+1
                symfunc_arr%symfunc(iconf)%y(ig,iat)=wa(n)
            enddo
        enddo
        if(parini%save_symfunc_force_ann) then
            do ib=1,nb
                do i=1,3
                    do ig=1,ng
                        n=n+1
                        symfunc_arr%symfunc(iconf)%y0d(ig,i,ib)=wa(n)
                    enddo
                enddo
            enddo
        endif
    endif bondbased            
    end associate
    end associate
    end associate
    !do iat=1,atoms_arr%atoms(iconf)%nat
    !    do ig=1,symfunc_arr%symfunc(iconf)%ng
    !        read(311,*) symfunc_arr%symfunc(iconf)%y(ig,iat)
    !    enddo
    !enddo
    close(311)
    call f_free(wa)
    write(*,*) "Reading symmetry functions done."
end subroutine read_symfunc
!*****************************************************************************************
subroutine save_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables
    integer:: iconf, ib, ig, i, iat, jat, i0
    real(8), allocatable:: gminarr(:,:), gmaxarr(:,:) !, poll_period
    integer, allocatable:: iatmin(:,:), iatmax(:,:), iconfmin(:,:), iconfmax(:,:)
    integer:: ibmin(350), ibmax(350)
    integer:: ngmax
    ngmax=350
    allocate(gminarr(1:ngmax,1:parini%ntypat))
    gminarr=huge(1.d20)
    allocate(gmaxarr(1:ngmax,1:parini%ntypat))
    gmaxarr=-huge(1.d20)
    allocate(iatmin(1:ngmax,1:parini%ntypat))
    iatmin=0.d0
    allocate(iatmax(1:ngmax,1:parini%ntypat))
    iatmax=0.d0
    allocate(iconfmin(1:ngmax,1:parini%ntypat))
    iconfmin=0.d0
    allocate(iconfmax(1:ngmax,1:parini%ntypat))
    iconfmax=0.d0
    ibmin(1:350)=0 ; ibmax(1:350)=0
    do iconf=1,atoms_arr%nconf
        !if(mod(iconf-1,nproc)==iproc) cycle
        !write(41,'(i6,i3)',advance='no') mod(iconf-1,nproc),iproc
        if(parini%bondbased_ann) then
            !-------------------------------- bond symmetry function --------------------------------
            !if(symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad/=2) stop 'ERROR: correct next line'
            do ib=1,symfunc_arr%symfunc(iconf)%linked_lists%maxbound_rad
                iat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(1,ib)
                jat=symfunc_arr%symfunc(iconf)%linked_lists%bound_rad(2,ib)
                if(iat>jat) cycle
                do ig=1,symfunc_arr%symfunc(iconf)%ng
                    if(symfunc_arr%symfunc(iconf)%y(ig,ib)<gminarr(ig,1)) then
                        gminarr(ig,1)=symfunc_arr%symfunc(iconf)%y(ig,ib)
                        ibmin(ig)=ib
                        iconfmin(ig,1)=iconf
                    endif
                    if(symfunc_arr%symfunc(iconf)%y(ig,ib)>gmaxarr(ig,1)) then
                        gmaxarr(ig,1)=symfunc_arr%symfunc(iconf)%y(ig,ib)
                        ibmax(ig)=ib
                        iconfmax(ig,1)=iconf
                    endif
                enddo
            enddo
            !-----------------------------------------------------------------------------------------
        else
            do iat=1,atoms_arr%atoms(iconf)%nat
                i=atoms_arr%atoms(iconf)%itypat(iat)
                do ig=1,symfunc_arr%symfunc(iconf)%ng
                    if(symfunc_arr%symfunc(iconf)%y(ig,iat)<gminarr(ig,i)) then
                        gminarr(ig,i)=symfunc_arr%symfunc(iconf)%y(ig,iat)
                        iatmin(ig,i)=iat
                        iconfmin(ig,i)=iconf
                    endif
                    if(symfunc_arr%symfunc(iconf)%y(ig,iat)>gmaxarr(ig,i)) then
                        gmaxarr(ig,i)=symfunc_arr%symfunc(iconf)%y(ig,iat)
                        iatmax(ig,i)=iat
                        iconfmax(ig,i)=iconf
                    endif
                enddo
            enddo
        endif
    enddo
    call yaml_mapping_open('symfunc bounds')
    do i=1,ann_arr%nann
    do ig=1,ann_arr%ann(1)%nn(0) !HERE
        if(iproc==0) then
        call yaml_sequence(advance='no')
        if(parini%bondbased_ann) then
            call yaml_mapping_open(trim(strmess),flow=.true.)
            !write(*,'(2(i7,2i4,es20.10),1x,a)') &
            !    iconfmin(ig,1),atoms_arr%atoms(iconfmin(ig,1))%nat,ibmin(ig),gminarr(ig,1), &
            !    iconfmax(ig,1),atoms_arr%atoms(iconfmax(ig,1))%nat,ibmax(ig),gmaxarr(ig,1),trim(strmess)

            call yaml_map('iconfmin',iconfmin(ig,1),fmt='(i7)')
            call yaml_map('nat',atoms_arr%atoms(iconfmin(ig,1))%nat,fmt='(i4)')
            call yaml_map('ibmin',ibmin(ig),fmt='(i4)')
            call yaml_map('gminarr',gminarr(ig,1),fmt='(es20.10)')
            call yaml_map('iconfmax',iconfmax(ig,1),fmt='(i7)')
            call yaml_map('nat',atoms_arr%atoms(iconfmax(ig,1))%nat,fmt='(i4)')
            call yaml_map('ibmax',ibmax(ig),fmt='(i4)')
            call yaml_map('gmaxarr',gmaxarr(ig,1),fmt='(es20.10)')
            call yaml_map('strmess',trim(strmess))

            !write(*,'(2(a50,i6,1x))') trim(atoms_arr%fn(iconfmin(ig,1))),atoms_arr%lconf(iconfmin(ig,1)), &
            !    trim(atoms_arr%fn(iconfmax(ig,1))),atoms_arr%lconf(iconfmax(ig,1))

            call yaml_map('fn_min',trim(atoms_arr%fn(iconfmin(ig,1))))
            call yaml_map('lconf_min',atoms_arr%lconf(iconfmin(ig,1)),fmt='(i6)')
            call yaml_map('fn_max',trim(atoms_arr%fn(iconfmax(ig,1))))
            call yaml_map('lconf_max',atoms_arr%lconf(iconfmax(ig,1)),fmt='(i6)')

            call yaml_mapping_close()
        else
            call yaml_mapping_open(trim(strmess),flow=.true.)
            !write(*,'(2(i7,2i4,es20.10),1x,a)') &
            !    iconfmin(ig,i),atoms_arr%atoms(iconfmin(ig,i))%nat,ibmin(ig),gminarr(ig,i), &
            !    iconfmax(ig,i),atoms_arr%atoms(iconfmax(ig,i))%nat,ibmax(ig),gmaxarr(ig,i),trim(strmess)

            call yaml_map('iconfmin',iconfmin(ig,i),fmt='(i7)')
            call yaml_map('nat',atoms_arr%atoms(iconfmin(ig,i))%nat,fmt='(i4)')
            call yaml_map('ibmin',ibmin(ig),fmt='(i4)')
            call yaml_map('gminarr',gminarr(ig,i),fmt='(es20.10)')
            call yaml_map('iconfmax',iconfmax(ig,i),fmt='(i7)')
            call yaml_map('nat',atoms_arr%atoms(iconfmax(ig,i))%nat,fmt='(i4)')
            call yaml_map('ibmax',ibmax(ig),fmt='(i4)')
            call yaml_map('gmaxarr',gmaxarr(ig,i),fmt='(es20.10)')
            call yaml_map('strmess',trim(strmess))

            !write(*,'(2(a50,i6,1x))') trim(atoms_arr%fn(iconfmin(ig,i))),atoms_arr%lconf(iconfmin(ig,i)), &
            !    trim(atoms_arr%fn(iconfmax(ig,i))),atoms_arr%lconf(iconfmax(ig,i))

            call yaml_map('fn_min',trim(atoms_arr%fn(iconfmin(ig,i))))
            call yaml_map('lconf_min',atoms_arr%lconf(iconfmin(ig,i)),fmt='(i6)')
            call yaml_map('fn_max',trim(atoms_arr%fn(iconfmax(ig,i))))
            call yaml_map('lconf_max',atoms_arr%lconf(iconfmax(ig,i)),fmt='(i6)')

            call yaml_mapping_close()
        endif
        endif
    enddo
    enddo
    call yaml_mapping_close()

    if(parini%bondbased_ann) then
        !------------------------------ bond symmetry function ----------------------------------------
        !IMPORTANT: the upper bound of the loop over i needs to be corrected for multicomponent systems.
        if(atoms_arr%atoms(1)%ntypat>1) stop 'ERROR: this part not ready for ntypat>1'
        do i=1,ann_arr%nann
        do i0=1,ann_arr%ann(i)%nn(0)
            if(gminarr(i0,1)==0.d0) then
                ann_arr%ann(i)%gbounds(1,i0)=-epsilon(1.d0)
            else
                ann_arr%ann(i)%gbounds(1,i0)=gminarr(i0,1)
            endif
            if(gmaxarr(i0,1)==0.d0) then
                ann_arr%ann(i)%gbounds(2,i0)=epsilon(1.d0)
            else
                ann_arr%ann(i)%gbounds(2,i0)=gmaxarr(i0,1)
            endif
            write(*,*) 'gbounds', ann_arr%ann(i)%gbounds(2,i0),ann_arr%ann(i)%gbounds(1,i0)
            ann_arr%ann(i)%two_over_gdiff(i0)=2.d0/(ann_arr%ann(i)%gbounds(2,i0)-ann_arr%ann(i)%gbounds(1,i0))
        enddo
        enddo
        ! ----------------------------------------------------------------------------------------------
    else
        do i=1,ann_arr%nann
        do i0=1,ann_arr%ann(i)%nn(0)
            !if(abs(gminarr(i0))<epsilon(1.d0)
            if(gminarr(i0,i)==0.d0) then
                ann_arr%ann(i)%gbounds(1,i0)=-epsilon(1.d0)
            else
                ann_arr%ann(i)%gbounds(1,i0)=gminarr(i0,i)
            endif
            if(gmaxarr(i0,i)==0.d0) then
                ann_arr%ann(i)%gbounds(2,i0)=epsilon(1.d0)
            else
                ann_arr%ann(i)%gbounds(2,i0)=gmaxarr(i0,i)
            endif
            ann_arr%ann(i)%two_over_gdiff(i0)=2.d0/(ann_arr%ann(i)%gbounds(2,i0)-ann_arr%ann(i)%gbounds(1,i0))
        enddo
        enddo
    endif
    !if(trim(parini%symfunc)=='do_not_save') then
    !    do iconf=1,atoms_arr%nconf
    !        call f_free(symfunc_arr%symfunc(iconf)%y)
    !        deallocate(symfunc_arr%symfunc(iconf)%linked_lists%bound_rad)
    !    enddo
    !endif
    deallocate(gminarr)
    deallocate(gmaxarr)
    deallocate(iatmin)
    deallocate(iatmax)
    deallocate(iconfmin)
    deallocate(iconfmax)
end subroutine save_gbounds
!*****************************************************************************************
subroutine randomize_data_order(atoms_arr)
    use mod_atoms, only: typ_atoms_arr, typ_atoms, atom_copy_old, atom_deallocate_old
    implicit none
    type(typ_atoms_arr), intent(inout):: atoms_arr
    !local variables 
    integer:: i1, i2, iconf !, nat_t
    real(8):: t1, t2 !, rat_t(3,1000), epot_t
    type(typ_atoms):: atoms_t
    !if(natmax>1000) stop 'ERROR: natmax>1000'
    iconf=0
    do
        call random_number(t1)
        call random_number(t2)
        i1=int(t1*real(atoms_arr%nconf,8))+1
        i2=int(t2*real(atoms_arr%nconf,8))+1
        if(i1<1 .or. i1>atoms_arr%nconf) stop 'ERROR: something wrong with i1'
        if(i2<1 .or. i2>atoms_arr%nconf) stop 'ERROR: something wrong with i2'
        if(i1==i2) cycle
        iconf=iconf+1
        call atom_copy_old(atoms_arr%atoms(i1),atoms_t,'atoms_arr%atoms(i1)->atoms_t')
        call atom_copy_old(atoms_arr%atoms(i2),atoms_arr%atoms(i1),'atoms_arr%atoms(i2)->atoms_arr%atoms(i1)')
        call atom_copy_old(atoms_t,atoms_arr%atoms(i2),'atoms_t->atoms_arr%atoms(i2)')

        !nat_t=natarr(i1)
        !rat_t(1:3,1:nat_t)=ratall(1:3,1:nat_t,i1)
        !epot_t=epotall(i1)

        !natarr(i1)=natarr(i2)
        !ratall(1:3,1:natarr(i2),i1)=ratall(1:3,1:natarr(i2),i2)
        !epotall(i1)=epotall(i2)

        !natarr(i2)=nat_t
        !ratall(1:3,1:nat_t,i2)=rat_t(1:3,1:nat_t)
        !epotall(i2)=epot_t

        if(iconf>1*atoms_arr%nconf) exit
    enddo
    call atom_deallocate_old(atoms_t)
end subroutine randomize_data_order
!*****************************************************************************************
subroutine set_ref_energy(parini,atoms_train,atoms_ref,ind)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(in):: atoms_train
    !type(typ_atoms_arr), intent(in):: atoms_valid
    type(typ_atoms_arr), intent(inout):: atoms_ref
    integer, intent(out):: ind(200)
    !local variables
    !real(8), allocatable:: 
    integer:: iconf, mat
    atoms_ref%nconf=200
    allocate(atoms_ref%atoms(atoms_ref%nconf))
    do mat=1,200
        atoms_ref%atoms(mat)%epot=huge(1.d0)
    enddo
    write(*,*) 'atoms_train%nconf ',atoms_train%nconf
    ind(1:200)=0
    do iconf=1,atoms_train%nconf
        write(*,*) 'nat ',atoms_train%atoms(iconf)%nat
        mat=atoms_train%atoms(iconf)%nat
        if(atoms_train%atoms(iconf)%epot<atoms_ref%atoms(mat)%epot) then
            ind(mat)=iconf
            atoms_ref%atoms(mat)%epot=atoms_train%atoms(iconf)%epot
        endif
    enddo
    do mat=1,200
        if(ind(mat)==0) cycle
        iconf=ind(mat)
        call atom_copy_old(atoms_train%atoms(iconf),atoms_ref%atoms(mat),'copy to atoms_ref')
    enddo
end subroutine set_ref_energy
!*****************************************************************************************
end module mod_train
!*****************************************************************************************
subroutine ann_evaluate_all(parini,iter,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    use mod_callback_ann, only: atoms_train=>atoms_train_t
    use mod_callback_ann, only: atoms_valid=>atoms_valid_t
    use mod_callback_ann, only: symfunc_train=>symfunc_train_t
    use mod_callback_ann, only: symfunc_valid=>symfunc_valid_t
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
    call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
end subroutine ann_evaluate_all
!*****************************************************************************************
subroutine ann_evaluate(parini,iter,ann_arr,symfunc_arr,atoms_arr,data_set)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    use mod_opt_ann, only: typ_opt_ann
    use mod_processors, only: iproc
    use futile
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: data_set
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    type(typ_opt_ann):: opt_ann
    real(8):: rmse, errmax, tt, pi
    real(8):: frmse, ttx, tty, ttz, ppx, ppy, ppz, tt1, tt2, tt3, ttn, tta
    integer:: iconf, ierrmax, iat, nat_tot, nconf_force
    real(8):: time1=0.d0
    real(8):: time2=0.d0
    real(8), save:: time_p=0.d0
    real(8):: dtime1, dtime2
    integer:: ilarge1, ilarge2, ilarge3, iunit, ios
    !character(28):: frmt1='(i6,5f10.3,i7,i5,3i6,a40,i6)'
    !character(28):: frmt2='(i6,5e10.1,i7,i5,3i6,a40,i6)'
    !character(28):: frmt
    character(28):: fmt_main
    character(28):: fmt1='(f10.3)'
    character(28):: fmt2='(e10.1)'
    character(15):: filename
    character(len=8):: str_key
    call cpu_time(time1)
    pi=4.d0*atan(1.d0)
    rmse=0.d0
    frmse=0.d0
    ttn=0.d0
    tta=0.d0
    nat_tot=0
    nconf_force=0
    errmax=0.d0
    ierrmax=0
    ilarge1=0
    ilarge2=0
    ilarge3=0
    ann_arr%event='evalu'
    if(parini%save_symfunc_behnam) then
        ann_arr%compute_symfunc=.false.
    else
        ann_arr%compute_symfunc=.true.
    endif
    if(parini%print_energy) then
        write(filename,'(a12,i3.3)') 'detailed_err',iter
        iunit=f_get_free_unit(10**5)
        if(trim(data_set)=="train") then
            open(unit=iunit,file=trim(filename),status='unknown',iostat=ios)
        elseif(trim(data_set)=="valid") then
            open(unit=iunit,file=trim(filename),status='old',access='append',iostat=ios)
        endif
        if(ios/=0) then
            write(*,'(a,a)') 'ERROR: failure openning file: ',trim(filename)
            stop
        endif
        !write(iunit,'(a2,a44,4a23)') "#", " ","E_dft","E_ann","E_dft-E_ann/atom (Ha)","E_dft-E_ann (eV)" 
    endif
    configuration: do iconf=1,atoms_arr%nconf
        if(.not. atoms_arr%conf_inc(iconf)) cycle
        call atom_copy_old(atoms_arr%atoms(iconf),atoms,'atoms_arr%atoms(iconf)->atoms')
        if(parini%save_symfunc_behnam) then
            call cal_ann_main(parini,atoms,symfunc_arr%symfunc(iconf),ann_arr,opt_ann)
        else
            call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
        endif
        !if(iter==parini%nstep_opt_ann) then
        !    write(40+ifile,'(2i6,2es24.15,es14.5)') iconf,atoms_arr%atoms(iconf)%nat, &
        !        atoms_arr%atoms(iconf)%epot/atoms_arr%atoms(iconf)%nat,atoms%epot/atoms_arr%atoms(iconf)%nat, &
        !        (atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        !endif
        tt=abs(atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        !HERE
        if(parini%print_energy) then
            write(iunit,'(i7,es14.5,a40,i6,a)') iconf,tt,trim(atoms_arr%fn(iconf)),atoms_arr%lconf(iconf),trim(data_set)
        endif
        if(tt>1.d-2) ilarge1=ilarge1+1
        if(tt>1.d-3) ilarge2=ilarge2+1
        if(tt>1.d-4) ilarge3=ilarge3+1
        if(tt>errmax) then
            errmax=tt
            ierrmax=iconf
        endif
        rmse=rmse+tt**2
        !if(tt>1.d-2) then
        !    atoms_arr%inclusion(iconf)=0
        !else
        !    atoms_arr%inclusion(iconf)=1
        !endif
        !write(22,'(a,i5.5)') 'configuration ',iconf
        !if(atoms%nat<=parini%nat_force) then
        nat_tot=nat_tot+atoms%nat
        nconf_force=nconf_force+1
        do iat=1,atoms%nat
            ttx=atoms_arr%atoms(iconf)%fat(1,iat)
            tty=atoms_arr%atoms(iconf)%fat(2,iat)
            ttz=atoms_arr%atoms(iconf)%fat(3,iat)
            ppx=atoms%fat(1,iat)
            ppy=atoms%fat(2,iat)
            ppz=atoms%fat(3,iat)
            !write(41,'(2i6,6f7.3)') iconf,iat,ttx,tty,ttz,ppx,ppy,ppz
            tt1=sqrt(ttx**2+tty**2+ttz**2)
            tt2=sqrt(ppx**2+ppy**2+ppz**2)
            tt3=(ppx-ttx)**2+(ppy-tty)**2+(ppz-ttz)**2
            frmse=frmse+tt3
            tt3=sqrt(tt3)
        enddo
        ttn=ttn+ann_arr%fchi_norm
        tta=tta+ann_arr%fchi_angle
        !write(44,'(2i7,4es14.5)') iter,iconf,ann_arr%fchi_norm,ann_arr%fchi_angle,ttn/nconf_force,tta/nconf_force
        !endif
    enddo configuration
    rmse=sqrt(rmse/real(atoms_arr%nconf_inc,8))
    if(nconf_force==0) nconf_force=1
    ttn=ttn/real(nconf_force,8)
    tta=tta/real(nconf_force,8)
    if(nat_tot==0) nat_tot=1
    frmse=sqrt(frmse/real(3*nat_tot,8))
    if(iproc==0) then
        rmse=rmse*1.d3
        errmax=errmax*1.d3
        if(rmse>99999.d0) then
            !frmt=frmt2
            fmt_main=fmt2
        else
            !frmt=frmt1
            fmt_main=fmt1
        endif

        !write(ifile,frmt) iter,rmse,ttn,tta,frmse,errmax,ierrmax,atoms_arr%atoms(ierrmax)%nat, &
        !    ilarge1,ilarge2,ilarge3,trim(atoms_arr%fn(ierrmax)),atoms_arr%lconf(ierrmax)
        write(str_key,'(a)') trim(data_set)
        call yaml_mapping_open(trim(str_key),flow=.true.,unit=ann_arr%iunit)
        call yaml_map('iter',iter,unit=ann_arr%iunit)
        call yaml_map('rmse',rmse,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('ttn',ttn,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('tta',tta,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('frmse',frmse,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('errmax',errmax,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('ierrmax',ierrmax,unit=ann_arr%iunit)
        call yaml_map('nat',atoms_arr%atoms(ierrmax)%nat,unit=ann_arr%iunit)
        call yaml_map('ilarge1',ilarge1,unit=ann_arr%iunit)
        call yaml_map('ilarge2',ilarge2,unit=ann_arr%iunit)
        call yaml_map('ilarge3',ilarge3,unit=ann_arr%iunit)
        call yaml_map('fn_ierrmax',trim(atoms_arr%fn(ierrmax)),unit=ann_arr%iunit)
        call yaml_map('lconf_ierrmax',atoms_arr%lconf(ierrmax),unit=ann_arr%iunit)
        call yaml_mapping_close(unit=ann_arr%iunit)
    endif
    call cpu_time(time2)
    dtime1=time1-time_p
    dtime2=time2-time1
    !write(*,'(a,2f20.2)') 'TIME ',dtime1,dtime2
    time_p=time2
    if(parini%print_energy) then
        close(iunit)
    endif
    ann_arr%compute_symfunc=.false.
end subroutine ann_evaluate
!*****************************************************************************************
