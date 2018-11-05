!*****************************************************************************************
module mod_callback
    use mod_atoms, only: typ_atoms_arr
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_ekf, only: typ_ekf
    type(typ_parini), pointer:: parini_t
    type(typ_atoms_arr), pointer:: atoms_smplx_t
    type(typ_ann_arr), pointer:: ann_arr_t
    type(typ_ekf), pointer:: ekf_t
end module mod_callback
!*****************************************************************************************
module mod_train
    implicit none
    private
    public:: ann_train !, convert_x_ann, convert_ann_epotd
contains
!*****************************************************************************************
subroutine ann_train(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_ekf, only: typ_ekf, ekf_rivals, ekf_rivals_tmp, ekf_behler
    use mod_atoms, only: typ_atoms_arr, typ_atoms
    use mod_processors, only: iproc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr):: ann_arr
    type(typ_ekf):: ekf
    type(typ_atoms_arr):: atoms_train
    type(typ_atoms_arr):: atoms_valid
    type(typ_atoms_arr):: atoms_smplx
    type(typ_symfunc_arr):: symfunc_train
    type(typ_symfunc_arr):: symfunc_valid
    integer:: iconf, ia, ityp
    real(8):: time1, time2, time3, t_ener_ref
    logical:: file_exists
    call f_routine(id='ann_train')
    call init_ann_train(parini,ann_arr,ekf)
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
    if(trim(parini%approach_ann)=='cent2') then
        call read_data_yaml(parini,'list_posinp_smplx.yaml',atoms_smplx)
    endif
    !-------------------------------------------------------
    if(iproc==0) then
        call yaml_map('number of ANN wights',ekf%n)
        call yaml_map('number of training data points',atoms_train%nconf)
        call yaml_map('number of validating data points',atoms_valid%nconf)
        !write(*,'(a,i)') 'number of ANN wights:             ',ekf%n
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
    call set_ebounds(ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    !-------------------------------------------------------
    ekf%x=f_malloc([1.to.ekf%n],id='ekf%x')
    if (.not. parini%restart_param) then
        call set_annweights(parini,ekf)
        if(trim(parini%approach_ann)=='cent2') then
            do ia=1,ann_arr%n
                ekf%x(ekf%loc(ia)+ekf%num(1)-1)=0.d0
                !write(*,*) 'XXX ',ia,ekf%loc(ia)+ekf%num(1)-1
            enddo
        endif
    else
        do ia=1,ann_arr%n
            call convert_ann_x(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
        enddo
    endif
    if(trim(parini%approach_ann)=='cent2') then
        call set_single_atom_energy(parini,ann_arr,ekf)
    endif

    ann_arr%compute_symfunc=.false.
    !if(parini%prefit_ann .and. trim(parini%approach_ann)=='cent2') then
    if(parini%prefit_ann ) then
        !call prefit_cent_ener_ref(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
        call prefit_cent(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    endif

    if(trim(parini%approach_ann)=='cent2' .and. parini%prefit_cent2_ann) then
        call cent2_simplex(parini,ann_arr,atoms_smplx,ekf)
    endif
    if(trim(parini%optimizer_ann)=='behler') then
        call ekf_behler(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    else if(trim(parini%optimizer_ann)=='rivals') then
        call ekf_rivals(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    else if(trim(parini%optimizer_ann)=='rivals_tmp') then
        call ekf_rivals_tmp(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,ekf)
    else if(trim(parini%optimizer_ann)=='lm') then
        call ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    else
        write(*,*) 'ERROR: unknown optimzer in ANN training'
    endif

    !call convert_x_ann(ekf%n,ekf%x,ann_arr) !HERE
    if(iproc==0) then
        if( ann_arr%exists_yaml_file) then
            call write_ann_all_yaml(parini,ann_arr,-1)
        else
            call write_ann_all(parini,ann_arr,-1)
        endif
    endif
    call final_ann_train(parini,ann_arr,ekf,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    !stop 'AAAAAAAAAAAAAAAAAAA'

    call f_release_routine()
end subroutine ann_train
!*****************************************************************************************
subroutine set_single_atom_energy(parini,ann_arr,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr !, typ_symfunc_arr
    use mod_ekf, only: typ_ekf
    use mod_atoms, only: typ_atoms, atom_allocate_old, atom_deallocate_old
    use mod_ann, only: typ_symfunc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_ekf), intent(inout):: ekf
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    integer:: ia, ityp
    real(8):: t_ener_ref
    do ia=1,ann_arr%n
        call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
    enddo
    !call atom_copy_old(atoms_train%atoms(1),atoms,'atoms_arr%atoms(iconf)->atoms')
    call atom_allocate_old(atoms,1,0,0)
    !write(*,*) atoms%nat
    !write(*,*) atoms%sat(:)
    !write(*,*) atoms%itypat(:)
    atoms%rat(1,1)=1.d0 ; atoms%rat(2,1)=1.d0 ; atoms%rat(3,1)=1.d0
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
        call cal_ann_main(parini,atoms,symfunc,ann_arr,ekf)
        ann_arr%ann(atoms%itypat(1))%ener_ref=t_ener_ref-atoms%epot
        call cal_ann_main(parini,atoms,symfunc,ann_arr,ekf)
        !call yaml_map('type',trim(atoms%sat(1)))
        !call yaml_map('charge',atoms%zat(1)+atoms%qat(1))
        !write(*,'(a,f)') 'Adjusting ener_ref: total charge= ',atoms%zat(1)+atoms%qat(1)
    enddo
    !call yaml_sequence_close()
    call atom_deallocate_old(atoms)
end subroutine set_single_atom_energy
!*****************************************************************************************
subroutine cent2_simplex(parini,ann_arr,atoms_smplx,ekf)
    !use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr !, typ_symfunc_arr
    use mod_ekf, only: typ_ekf
    use mod_atoms, only: typ_atoms_arr
    use mod_callback
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in), target:: parini
    type(typ_ann_arr), intent(inout), target:: ann_arr
    type(typ_ekf), intent(inout), target:: ekf
    type(typ_atoms_arr), intent(inout), target:: atoms_smplx
    !local variables
    real(8):: vertices(10,11), fval(11)
    real(8):: step, ftol
    integer:: ndim, iter, i
    !external:: cal_rmse_force_cent2
    !external:: cal_rmse_energy_cent2
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
    ekf_t=>ekf
    !call simplex(vertices,fval,step,ndim,ftol,cal_rmse_force_cent2,iter)
    call simplex(vertices,fval,step,ndim,ftol,cal_rmse_energy_cent2,iter)
end subroutine cent2_simplex
!*****************************************************************************************
subroutine cal_rmse_force_cent2(ndim,vertex,rmse_force_cent2)
    use mod_interface
    use mod_callback, only: atoms_smplx=>atoms_smplx_t, parini=>parini_t
    use mod_callback, only: ann_arr=>ann_arr_t, ekf=>ekf_t
    use mod_atoms, only: typ_atoms, atom_copy_old 
    use mod_ann, only: typ_symfunc
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: vertex(ndim)
    real(8), intent(out) :: rmse_force_cent2
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
        call cal_ann_main(parini,atoms,symfunc,ann_arr,ekf)
        do iat=1,atoms%nat
            rmse=rmse+(atoms%fat(1,iat)-atoms_smplx%atoms(iconf)%fat(1,iat))**2 &
                     +(atoms%fat(2,iat)-atoms_smplx%atoms(iconf)%fat(2,iat))**2 &
                     +(atoms%fat(3,iat)-atoms_smplx%atoms(iconf)%fat(3,iat))**2
        enddo
        nat_tot=nat_tot+atoms%nat
    enddo
    rmse_force_cent2=sqrt(rmse/real(3*nat_tot,8))
    write(*,*) 'rmse_force_cent2 ',rmse_force_cent2
end subroutine cal_rmse_force_cent2
!*****************************************************************************************
subroutine cal_rmse_energy_cent2(ndim,vertex,rmse_energy_cent2)
    use mod_interface
    use mod_callback, only: atoms_smplx=>atoms_smplx_t, parini=>parini_t
    use mod_callback, only: ann_arr=>ann_arr_t, ekf=>ekf_t
    use mod_atoms, only: typ_atoms, atom_copy_old
    use mod_ann, only: typ_symfunc
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: vertex(ndim)
    real(8), intent(out) :: rmse_energy_cent2
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

    !call set_single_atom_energy(parini,ann_arr,ekf)
    rmse=0.d0
    ann_arr%event='potential'
    ann_arr%compute_symfunc=.true.
    do iconf=1,atoms_smplx%nconf
        call atom_copy_old(atoms_smplx%atoms(iconf),atoms,'atoms_smplx%atoms(iconf)->atoms')
        call cal_ann_main(parini,atoms,symfunc,ann_arr,ekf)
        do iat=1,atoms%nat
            rmse=rmse+(atoms%epot-atoms_smplx%atoms(iconf)%epot)**2
        enddo
    enddo
    rmse_energy_cent2=sqrt(rmse/real(atoms_smplx%nconf,8))
    write(*,*) 'rmse_energy_cent2 ',rmse_energy_cent2
end subroutine cal_rmse_energy_cent2
!*****************************************************************************************
subroutine init_ann_train(parini,ann_arr,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_ekf, only: typ_ekf
    use mod_processors, only: iproc
    use yaml_output
    use futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_ekf), intent(inout):: ekf
    !local variables
    integer:: ialpha, i, ios, ia
    character(30):: fnout
    character (50)::fname
    integer:: ierr
    ann_arr%n=parini%ntypat
    if(parini%bondbased_ann) then
        ann_arr%n=4
    endif
    if(ann_arr%n==0) stop 'ERROR: number of type of atoms zero in ann_train'
    call yaml_map('number of ann',ann_arr%n)
    !write(*,*) 'Here', ann_arr%n
    allocate(ann_arr%ann(ann_arr%n))
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
    ekf%num(1:10)=0
    ekf%n=0
    call yaml_sequence_open('EKF') !,flow=.true.)
    do i=1,ann_arr%n
        do ialpha=1,ann_arr%ann(i)%nl
            ekf%num(i)=ekf%num(i)+(ann_arr%ann(i)%nn(ialpha-1)+1)*ann_arr%ann(i)%nn(ialpha)
        enddo
        ekf%loc(i)=ekf%n+1
        ekf%n=ekf%n+ekf%num(i)
        call yaml_sequence(advance='no')
        call yaml_map('iann',i)
        call yaml_map('loc',ekf%loc(i))
        call yaml_map('num',ekf%num(i))
        call yaml_map('n',ekf%n)
        !write(*,'(a,3i5)') 'EKF: ',ekf%loc(i),ekf%num(i),ekf%n
    enddo
    call yaml_sequence_close()
    call ann_allocate(ekf,ann_arr)
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
        !fnout='err_train'
        !open(unit=11,file=fnout,status='replace',iostat=ios)
        !if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
        !!write(fnout,'(a12,i3.3)') 'err_valid',iproc
        !fnout='err_valid'
        !open(unit=12,file=fnout,status='replace',iostat=ios)
        !if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
    endif
end subroutine init_ann_train
!*****************************************************************************************
subroutine final_ann_train(parini,ann_arr,ekf,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_ekf, only: typ_ekf
    use mod_atoms, only: typ_atoms_arr
    use mod_processors, only: iproc
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_ekf), intent(inout):: ekf
    type(typ_atoms_arr), intent(inout):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    !local variables
    integer:: iconf
    call f_free(ekf%x)

    if(iproc==0) then
        !close(11)
        !close(12)
        call yaml_close_stream(unit=ann_arr%iunit)
    endif

    call ann_deallocate(ann_arr)

    do iconf=1,atoms_train%nconf
        call f_free(symfunc_train%symfunc(iconf)%y)
    enddo

    do iconf=1,atoms_valid%nconf
        call f_free(symfunc_valid%symfunc(iconf)%y)
    enddo

    deallocate(atoms_train%conf_inc)
    deallocate(atoms_valid%conf_inc)
    !deallocate(atoms_train%inclusion)
end subroutine final_ann_train
!*****************************************************************************************
subroutine set_conf_inc_random(parini,atoms_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
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
    do irand=1,atoms_arr%nconf_inc
        call random_number(tt)
        tt=tt*real(atoms_arr%nconf)
        iconf=int(tt)+1
        atoms_arr%conf_inc(iconf)=.true.
    enddo
end subroutine set_conf_inc_random
!*****************************************************************************************
subroutine apply_gbounds_atom(parini,ann_arr,atoms_arr,symfunc_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
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
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
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
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_atoms, only: typ_atoms_arr
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
        call rxyz_cart2int_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,atoms_arr%atoms(iconf)%rat,ratred)
        call backtocell_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,ratred)
        call rxyz_int2cart_alborz(atoms_arr%atoms(iconf)%nat,atoms_arr%atoms(iconf)%cellvec,ratred,atoms_arr%atoms(iconf)%rat)
        deallocate(ratred)
        endif
        do iat=1,atoms_arr%atoms(iconf)%nat
            do i=1,ann_arr%n
                if(trim(atoms_arr%atoms(iconf)%sat(iat))==trim(parini%stypat(i))) then
                    atoms_arr%atoms(iconf)%itypat(iat)=parini%ltypat(i)
                    exit
                endif
            enddo
        enddo
    enddo
end subroutine prepare_atoms_arr
!*****************************************************************************************
subroutine set_ebounds(ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    use mod_interface
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(in):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    !local variables
    integer:: iconf, nat
    real(8):: epot, tt
    ann_arr%ann(1)%ebounds(1)= 1.d20
    ann_arr%ann(1)%ebounds(2)=-1.d20
    do iconf=1,atoms_train%nconf
        tt=atoms_train%atoms(iconf)%epot/atoms_train%atoms(iconf)%nat
        ann_arr%ann(1)%ebounds(1)=min(ann_arr%ann(1)%ebounds(1),tt)
        ann_arr%ann(1)%ebounds(2)=max(ann_arr%ann(1)%ebounds(2),tt)
    enddo
    if(atoms_train%nconf==1) then
        tt=atoms_train%atoms(1)%epot !/atoms_train%atoms(1)%nat
        ann_arr%ann(1)%ebounds(1)=(tt-1.d0)/atoms_train%atoms(1)%nat
        ann_arr%ann(1)%ebounds(2)=(tt+1.d0)/atoms_train%atoms(1)%nat
    endif
    ann_arr%ann(1)%ebounds(1)=-1.d0
    ann_arr%ann(1)%ebounds(2)= 1.d0
    !write(*,'(a,2es14.5)') 'ebounds: ',ann_arr%ann(1)%ebounds(1),ann_arr%ann(1)%ebounds(2)
    tt=2.d0/(ann_arr%ann(1)%ebounds(2)-ann_arr%ann(1)%ebounds(1))
    do iconf=1,atoms_train%nconf
        nat=1 !atoms_train%atoms(iconf)%nat
        epot=atoms_train%atoms(iconf)%epot
        symfunc_train%symfunc(iconf)%epot=(epot/nat-ann_arr%ann(1)%ebounds(1))*tt-1.d0
        !write(*,'(a,i6,f8.3)') 'train: ',iconf,symfunc_train%symfunc(iconf)%epot
    enddo
    do iconf=1,atoms_valid%nconf
        nat=1 !atoms_valid%atoms(iconf)%nat
        epot=atoms_valid%atoms(iconf)%epot
        symfunc_valid%symfunc(iconf)%epot=(epot/nat-ann_arr%ann(1)%ebounds(1))*tt-1.d0
        !write(*,'(a,i6,f8.3)') 'valid: ',iconf,symfunc_valid%symfunc(iconf)%epot
    enddo
end subroutine set_ebounds
!*****************************************************************************************
subroutine set_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc, mpi_comm_abz
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: strmess
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    !local variables 
    integer:: iconf
    integer:: jproc, ierr, ireq, ireq_tmp, nreq, mreq, ii !, itry
    integer, allocatable:: ireqarr(:)
    logical:: flag
#if defined(MPI)
    include 'mpif.h'
    integer:: status_mpi(MPI_STATUS_SIZE)
#endif
#if defined(MPI)
    associate(MPI_DP=>MPI_DOUBLE_PRECISION)
    if(nproc>1) then
        nreq=atoms_arr%nconf/nproc
        if(mod(atoms_arr%nconf,nproc)>iproc) nreq=nreq+1
        nreq=atoms_arr%nconf-nreq
        ireqarr=f_malloc([1.to.nreq],id='ireqarr')
    endif
#endif
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
#if defined(MPI)
        if(nproc>1) then
        do jproc=0,nproc-1
            if(jproc==iproc) cycle
            associate(ntot=>symfunc_arr%symfunc(iconf)%ng*symfunc_arr%symfunc(iconf)%nat)
            call MPI_ISEND(symfunc_arr%symfunc(iconf)%y(1,1),ntot,MPI_DP,jproc,iconf,mpi_comm_abz,ireq_tmp,ierr)
            end associate
        enddo
        endif
#endif
    enddo configuration
#if defined(MPI)
    if(nproc>1) then
    ireq=0
    do iconf=1,atoms_arr%nconf
        jproc=mod(iconf-1,nproc)
        if(jproc==iproc) cycle
        ireq=ireq+1
        associate(ntot=>symfunc_arr%symfunc(iconf)%ng*symfunc_arr%symfunc(iconf)%nat)
        call MPI_IRECV(symfunc_arr%symfunc(iconf)%y(1,1),ntot,MPI_DP,jproc,iconf,mpi_comm_abz,ireqarr(ireq),ierr)
        end associate
    enddo
    call MPI_BARRIER(mpi_comm_abz,ierr)

    !itry=0
    ireq=0
    mreq=nreq
    do
        !itry=itry+1
        ireq=ireq+1
        flag=.false.
        call MPI_TEST(ireqarr(ireq),flag,status_mpi,ierr)
        if(flag) then
            !write(41,'(i3,i6,l2,4i6)') status_mpi(MPI_SOURCE),status_mpi(MPI_TAG),flag,mreq,nreq,ireq,itry
            ii=ireqarr(ireq)
            ireqarr(ireq)=ireqarr(mreq)
            ireqarr(mreq)=ii
            mreq=mreq-1
        endif
        if(mreq==0) exit
        if(ireq>=mreq) ireq=0
    enddo
    call MPI_BARRIER(mpi_comm_abz,ierr)
    endif !end of if for nproc>1
#endif
    call save_gbounds(parini,ann_arr,atoms_arr,strmess,symfunc_arr)
#if defined(MPI)
    if(nproc>1) then
        call f_free(ireqarr)
    endif
    end associate
#endif
end subroutine set_gbounds
!*****************************************************************************************
subroutine write_symfunc(parini,iconf,atoms_arr,strmess,symfunc_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc, mpi_comm_abz
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
    do iat=1,atoms_arr%atoms(iconf)%nat
        wa(3+iat*3-2)=atoms_arr%atoms(iconf)%rat(1,iat)
        wa(3+iat*3-1)=atoms_arr%atoms(iconf)%rat(2,iat)
        wa(3+iat*3-0)=atoms_arr%atoms(iconf)%rat(3,iat)
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
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc, mpi_comm_abz
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
    do iat=1,nat
        ttx=abs(wa(3+iat*3-2)-atoms_arr%atoms(iconf)%rat(1,iat))
        tty=abs(wa(3+iat*3-1)-atoms_arr%atoms(iconf)%rat(2,iat))
        ttz=abs(wa(3+iat*3-0)-atoms_arr%atoms(iconf)%rat(3,iat))
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
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_linked_lists, only: typ_pia_arr
    use mod_processors, only: iproc, nproc, mpi_comm_abz
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
        !associate(n=>symfunc_arr%symfunc(iconf)%ng*symfunc_arr%symfunc(iconf)%nat)
        !call MPI_RECV(symfunc_arr%symfunc(iconf)%y(1,1),n,MPI_DP,mod(iconf-1,nproc),iconf,mpi_comm_abz,status_mpi,ierr)
        !end associate
        !write(41,'(i6,i3)') status_mpi(MPI_TAG),status_mpi(MPI_SOURCE)
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
    do i=1,ann_arr%n
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
        do i=1,ann_arr%n
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
        do i=1,ann_arr%n
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
    use mod_interface
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
end module mod_train
!*****************************************************************************************
subroutine convert_x_ann(n,x,ann)
    use mod_interface
    use mod_ann, only: typ_ann
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: x(n)
    type(typ_ann), intent(inout):: ann
    !local variables
    integer:: i, j, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                l=l+1
                ann%a(i,j,ialpha)=x(l)
            enddo
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            ann%b(i,ialpha)=x(l)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_x_ann
!*****************************************************************************************
subroutine convert_ann_x(n,x,ann)
    use mod_interface
    use mod_ann, only: typ_ann
    implicit none
    integer, intent(in):: n
    real(8), intent(inout):: x(n)
    type(typ_ann), intent(inout):: ann
    !local variables
    integer:: i, j, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do j=1,ann%nn(ialpha)
            do i=1,ann%nn(ialpha-1)
                l=l+1
                x(l)=ann%a(i,j,ialpha)
            enddo
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            x(l)=ann%b(i,ialpha)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_ann_x
!*****************************************************************************************
subroutine convert_ann_epotd(ann,n,epotd)
    use mod_interface
    use mod_ann, only: typ_ann
    implicit none
    type(typ_ann), intent(in):: ann
    integer, intent(in):: n
    real(8), intent(inout):: epotd(n)
    !local variables
    integer:: i, ij, l, ialpha
    l=0
    do ialpha=1,ann%nl
        do ij=1,ann%nn(ialpha)*ann%nn(ialpha-1)
            l=l+1
            epotd(l)=ann%ad(ij,ialpha)
        enddo
        do i=1,ann%nn(ialpha)
            l=l+1
            epotd(l)=ann%bd(i,ialpha)
        enddo
    enddo
    if(l/=n) stop 'ERROR: l/=n'
end subroutine convert_ann_epotd
!*****************************************************************************************
