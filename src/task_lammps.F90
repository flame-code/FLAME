!*****************************************************************************************
subroutine lammps_task(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
#if defined(HAVE_LAMMPS)
    use mpi
    use LAMMPS
    use callback !, only: atoms
    use, intrinsic :: ISO_C_binding, only : C_double, C_ptr, C_int, C_FUNPTR
#endif
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
#if defined(HAVE_LAMMPS)
    type (C_ptr) :: lmp
    character(100):: str_run
    integer:: iconf
    type(typ_atoms_arr):: atoms_arr
    call copy_parini_for_lammps(parini)
    call acf_read_new(parini,'posinp.acf',10000,atoms_arr)
    call atom_copy(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(1)->atoms')
    call lammps_open_no_mpi ('lmp -log log.lammps', lmp)
    call lammps_file (lmp, 'in.lammps')
    call lammps_set_callback(lmp)
    call lammps_set_external_vector_length(lmp,2)
    if(parini%nrun_lammps<1) then
        stop 'ERROR: parini%nrun_lammps<1 in lammps_task'
    endif
    write(str_run,'(a,1x,i)') 'run',parini%nrun_lammps
    call lammps_command (lmp,str_run)
    call lammps_close (lmp)
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
