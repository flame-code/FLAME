!*****************************************************************************************
! IMPORTANT NOTE: 
! The structure written by the subroutine "lammps_write" needs the rotated structure. 
! At the moment, it supposes that the structure sent to this subroutine is
! already rotated. In case of changing the format of posinp.acf, the
! lammps_write should be changed in a way it does the rotation iteself.
!-----------------------------------------------------------------------------------------
subroutine lammps_task(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr, atom_copy, atom_deallocate, set_typat
    use mod_atoms, only: get_rat, set_rat
    use mod_yaml_conf, only: write_yaml_conf, read_yaml_conf
    use mod_potential, only: potcode
    use mod_potential, only: init_potential_forces
    use mod_potential, only: fini_potential_forces
#if defined(HAVE_LAMMPS)
    use mpi
    use LAMMPS
    use mod_acf, only: acf_read_new
    use callback !, only: atoms
    use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int, C_FUNPTR
#endif
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
#if defined(HAVE_LAMMPS)
    type (C_ptr) :: lmp
    !character(100):: str_run
    integer:: iconf
    type(typ_atoms_arr):: atoms_arr
    logical:: acf_exists, yaml_exists
    real(8), allocatable:: rat(:,:)
    real(8):: dproj(6), rotmat(3,3)
    !-------------------------------------------------------
    potcode=trim(parini%potential_potential)
    !-------------------------------------------------------
    call copy_parini_for_lammps(parini)
    inquire(file='posinp.yaml',exist=yaml_exists)
    inquire(file='posinp.acf',exist=acf_exists)
    if(yaml_exists) then
        call read_yaml_conf(parini,'posinp.yaml',1000,atoms_arr)
        allocate(rat(3,atoms_arr%atoms(1)%nat))
        call get_rat(atoms_arr%atoms(1),rat)
        call latvec2dproj(dproj,atoms_arr%atoms(1)%cellvec,rotmat,rat,atoms_arr%atoms(1)%nat)
        call set_rat(atoms_arr%atoms(1),rat,setall=.true.)
        deallocate(rat)
    elseif(acf_exists) then
        call acf_read_new(parini,'posinp.acf',1000,atoms_arr)
    else
        stop 'ERROR: this task reads only yaml and acf.'
    endif
    call atom_copy(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(1)->atoms')
    call set_typat(atoms)
    call init_potential_forces(parini,atoms)
    call lammps_write(parini,atoms)
    call lammps_open_no_mpi ('lmp -log log.lammps', lmp)
    call lammps_file (lmp, 'in.lammps.1')
    call lammps_set_callback(lmp)
    call lammps_file (lmp, 'in.lammps.2')
    call lammps_set_external_vector_length(lmp,2)
    ! nrun in flame_in.yaml cannot have a smaller value than 1
    !if(parini%nrun_lammps<1) then
    !    stop 'ERROR: parini%nrun_lammps<1 in lammps_task'
    !endif
    !write(str_run,'(1a,1x,1i8)') 'run',parini%nrun_lammps
    !lammps_command is the LAMMPS function that does the
    !task we have asked, e.g. molecular dynamics,
    !therefore, lammps_command will not be left until
    !the required task is completed.
    !call lammps_command (lmp,str_run)
    call lammps_close (lmp)
    call fini_potential_forces(parini,atoms)
    do iconf=1,atoms_arr%nconf
        call atom_deallocate(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
    call atom_deallocate(atoms)
#else
    stop 'ERROR: FLAME is compiled without linking with LAMMPS.'
#endif
end subroutine lammps_task
!*****************************************************************************************
subroutine lammps_write(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    implicit none
    integer:: iat,j
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms
    integer,allocatable:: id_type(:)
    allocate(id_type(1:atoms%nat))
    do iat=1,atoms%nat
        do j=1,atoms%ntypat
            if(trim(atoms%sat(iat))==trim(atoms%stypat(j))) id_type(iat)=j
        enddo
    enddo

    open(unit=10,file="lammps_posinp",status='replace')
    write(10,'(a)') "# Position data file"
    write(10,*)
    write(10,'(i8,a)') atoms%nat," atoms"
    write(10,'(i3,a)') atoms%ntypat, " atom types"
    write(10,*)
    write(10,'(a,es24.15,a)') "0.000000",atoms%cellvec(1,1),"   xlo xhi"
    write(10,'(a,es24.15,a)') "0.000000",atoms%cellvec(2,2),"   ylo yhi"
    write(10,'(a,es24.15,a)') "0.000000",atoms%cellvec(3,3),"   zlo zhi"
    write(10,*)
    write(10,'(3es24.15,3x,a)') atoms%cellvec(1,2),atoms%cellvec(1,3),atoms%cellvec(2,3),"xy xz yz"
    write(10,*)
    write(10,'(a)')"Atoms"
    write(10,*)

    call update_ratp(atoms)
    do iat= 1,atoms%nat
        write(10,'(i8,i3,3es24.15)') iat, id_type(iat),atoms%ratp(1,iat),atoms%ratp(2,iat),atoms%ratp(3,iat)
    enddo
    close(10)
    deallocate(id_type)
end subroutine lammps_write
!*****************************************************************************************
