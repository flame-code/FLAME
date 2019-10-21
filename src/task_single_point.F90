!*****************************************************************************************
subroutine single_point_task(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, typ_file_info, set_ndof, atom_deallocate
    use mod_potential, only: fcalls, perfstatus, potential
    use mod_acf, only: acf_read_new, acf_write
    use mod_yaml_conf, only: write_yaml_conf, read_yaml_conf
    use mod_processors, only: iproc
    use mod_const, only: ev2ha, ang2bohr
    use yaml_output
    implicit none
    type(typ_parini), intent(inout):: parini !poscar_getsystem must be called from parser
    !so parini can be intent(in) in future.
    !local variables
    type(typ_atoms_arr):: atoms_arr
    type(typ_file_info):: file_info
    real(8):: tt1, tt2, fxyz(3)
    integer:: iconf, iat
    logical:: acf_exists, yaml_exists
    inquire(file='posinp.yaml',exist=yaml_exists)
    inquire(file='posinp.acf',exist=acf_exists)
    if(yaml_exists) then
        call read_yaml_conf(parini,'posinp.yaml',10000,atoms_arr)
    elseif(acf_exists) then
        call acf_read_new(parini,'posinp.acf',10000,atoms_arr)
    else
        atoms_arr%nconf=1
        allocate(atoms_arr%atoms(atoms_arr%nconf))
        call read_poscar_for_single_point(parini,atoms_arr%atoms(1))
    endif
    do iconf=1,atoms_arr%nconf
        call set_ndof(atoms_arr%atoms(iconf))
    enddo
    potential=trim(parini%potential_potential)
    file_info%filename_positions='posout.acf'
    file_info%print_force=parini%print_force_single_point
    file_info%file_position='new'
    if(trim(parini%frmt_single_point)/='unknown') then
        file_info%frmt=trim(parini%frmt_single_point)
    endif
   
    if (parini%usesocket) then
        call netsock_task(parini)
    else
        do iconf=1,atoms_arr%nconf
            if(trim(potential)/='netsock' .or. iconf==1) then 
                call init_potential_forces(parini,atoms_arr%atoms(iconf))
            endif
            if(iconf==1) then 
                call yaml_sequence_open('epot of all configurations')
            endif
            call yaml_sequence(advance='no')
            call cal_potential_forces(parini,atoms_arr%atoms(iconf))
            !call atoms_all%fatall(1:3,1:atoms_all%atoms%nat,iconf)=atoms_all%atoms%fat(1:3,1:atoms_all%atoms%nat)
            if(iconf==1) then
                tt1=0.d0
                tt2=0.d0
            else
                tt1=atoms_arr%atoms(iconf)%epot-atoms_arr%atoms(1)%epot
                tt2=atoms_arr%atoms(iconf)%epot-atoms_arr%atoms(iconf-1)%epot
            endif
            call yaml_map('conf. number',iconf,fmt='(i6.6)')
            call yaml_map('epot',atoms_arr%atoms(iconf)%epot,fmt='(e24.15)')
            call yaml_map('diff w.r.t. first conf',tt1,fmt='(f12.7)')
            call yaml_map('diff w.r.t. prev. conf',tt2,fmt='(f12.7)')
            !write(*,'(a,i7,e24.15,2f10.5)') 'EPOT',iconf,atoms_arr%atoms(iconf)%epot,tt1,tt2
            !if(parini%print_force_single_point) then
            !    do iat=1,atoms_all%atoms%nat
            !        fxyz(1)=atoms_all%atoms%fat(1,iat)
            !        fxyz(2)=atoms_all%atoms%fat(2,iat)
            !        fxyz(3)=atoms_all%atoms%fat(3,iat)
            !        write(*,'(3f12.8)') fxyz(1),fxyz(2),fxyz(3)
            !    enddo
            !endif
            if(trim(potential)/='netsock' .or. iconf==atoms_arr%nconf) then 
                call final_potential_forces(parini,atoms_arr%atoms(iconf))
            endif
            if (iconf==2)  file_info%file_position='append'
            if(yaml_exists) then
                file_info%filename_positions='posout.yaml'
                call write_yaml_conf(file_info,atoms=atoms_arr%atoms(iconf),strkey='posout')
            elseif(acf_exists) then
                file_info%filename_positions='posout.acf'
                call acf_write(file_info,atoms=atoms_arr%atoms(iconf),strkey='posout')
            endif
        enddo
        do iconf=1,atoms_arr%nconf
            call atom_deallocate(atoms_arr%atoms(iconf))
        enddo
    endif
    call yaml_sequence_close()
    !call atom_all_deallocate(atoms_all,ratall=.true.,fatall=.true.,epotall=.true.)
end subroutine single_point_task
!*****************************************************************************************
subroutine read_poscar_for_single_point(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, atom_allocate_old, update_rat
    use global, only: units
    implicit none
    type(typ_parini), intent(inout):: parini !poscar_getsystem must be called from parser
    !so parini can be intent(in) in future.
    type(typ_atoms):: atoms
    !local variables
    real(8), allocatable:: xred(:,:)
    real(8), allocatable:: fcart(:,:)
    logical, allocatable:: fixat(:)
    integer, allocatable:: fragarr(:)
    real(8):: strten(6), printval1, printval2
    logical:: fixlat(7), readfix, readfrag
    integer:: iat
    character(40):: filename
    logical:: file_exists
    filename='posinp.vasp'
    inquire(file=trim(filename),exist=file_exists)
    if (.not. file_exists) then 
        filename='POSCAR'
        inquire(file=trim(filename),exist=file_exists)
        if (.not. file_exists) stop "VASP file not found"
    endif
    write (*,*) "Reading structure from ",trim(filename)

    call poscar_getsystem(parini,trim(filename))
    allocate(xred(3,parini%nat),source=0.d0)
    allocate(fcart(3,parini%nat),source=0.d0)
    if(.not.allocated(fixat)) allocate(fixat(parini%nat),source=.false.)
    if(.not.allocated(fragarr)) allocate(fragarr(parini%nat),source=0)
    atoms%nat=parini%nat
    atoms%boundcond='bulk'
    call atom_allocate_old(atoms,parini%nat,0,0)
    call read_atomic_file_poscar(filename,atoms%nat,units,xred,atoms%cellvec,fcart,strten, &
        fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
    call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,xred,atoms%ratp)
    call update_rat(atoms,upall=.true.)
    do iat=1,parini%nat
        atoms%sat(iat)=trim(parini%char_type(parini%typat_global(iat)))
    enddo
    deallocate(fixat)
    deallocate(xred)
    deallocate(fcart)
    deallocate(fragarr)
    if(allocated(parini%znucl)) deallocate(parini%znucl)
    if(allocated(parini%char_type)) deallocate(parini%char_type)
    if(allocated(parini%amu)) deallocate(parini%amu)
    if(allocated(parini%rcov)) deallocate(parini%rcov)
    if(allocated(parini%typat_global)) deallocate(parini%typat_global)
end subroutine read_poscar_for_single_point
!*****************************************************************************************
