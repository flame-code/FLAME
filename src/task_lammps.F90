!*****************************************************************************************
subroutine lammps_task(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_potential, only: potential
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
    !-------------------------------------------------------
    potential=trim(parini%potential_potential)
    !-------------------------------------------------------
    call copy_parini_for_lammps(parini)
    call acf_read_new(parini,'posinp.acf',10000,atoms_arr)
    call atom_copy(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(1)->atoms')
    call lammps_write(parini,atoms)
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
subroutine lammps_write(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    integer:: iat,j
    type(typ_parini), intent(in):: parini
    type(typ_atoms):: atoms
    integer,allocatable:: id_type(:)
    call set_typat(atoms)
    allocate(id_type(1:atoms%nat))
    do iat=1,atoms%nat
        do j=1,atoms%ntypat
            if(trim(atoms%sat(iat))==trim(atoms%stypat(j))) id_type(iat)=j
        enddo
    enddo

    open(unit=10,file="lammps_posinp",status='replace')
    write(10,'(a)') "# Position data file"
    write(10,*)
    write(10,'(i4,a)') atoms%nat," atoms"
    write(10,'(i1,a)') parini%ntypat, " atom types"
    write(10,*)
    write(10,'(a,f14.7,a)') "0.000000",atoms%cellvec(1,1),"   xlo xhi"
    write(10,'(a,f14.7,a)') "0.000000",atoms%cellvec(2,2),"   ylo yhi"
    write(10,'(a,f14.7,a)') "0.000000",atoms%cellvec(3,3),"   zlo zhi"
    write(10,*)
    write(10,'(f8.6,f11.6,f11.6,3x,a)') atoms%cellvec(1,2),atoms%cellvec(1,3),atoms%cellvec(2,3),"xy xz yz"
    write(10,*)
    write(10,'(a)')"Atoms"
    write(10,*)

    do iat= 1,atoms%nat
        write(10,'(i4,i2,3f15.10)') iat, id_type(iat),atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    enddo
    deallocate(id_type)
end subroutine lammps_write
!*****************************************************************************************
