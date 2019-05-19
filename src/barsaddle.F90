!*****************************************************************************************
subroutine bar_saddle(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use dictionaries
    use dynamic_memory
    use yaml_output
    use mod_processors, only: nproc, iproc
    use mod_atoms, only: typ_atoms, atom_deallocate, typ_atoms_arr
    use mod_atoms, only: atom_copy, get_rat
    use mod_potential, only: potential
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer :: iconf, infocode
    real(8) :: etot,fnoise
    type(typ_atoms) :: atoms
    type(typ_atoms_arr) :: atoms_arr
    real(8), dimension(:,:), allocatable :: fxyz
    real(8), dimension(:,:), pointer :: rxyz1, rxyz2
    potential=trim(parini%potential_potential)
    !The following reads a maximum of two configurations but
    !I am sending you a posinp.yaml that includes one configuration.
    call read_yaml_conf(parini,'posinp.yaml',2,atoms_arr)
    call atom_copy(atoms_arr%atoms(1),atoms,'atoms_arr%atoms(1)=>atoms')
    do iconf=1,atoms_arr%nconf
        call atom_deallocate(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
    allocate(rxyz1(3,atoms%nat),rxyz2(3,atoms%nat))
    allocate(fxyz(3,atoms%nat))
    call get_rat(atoms,rxyz1)
    call init_potential_forces(parini,atoms)
    !-------------------------------------------------------
    !MAX, call your routine here.
    !The following line is an example to get energy and forces
    call call_bigdft(nproc,iproc,atoms,rxyz1,etot,fxyz,fnoise,infocode,parini)
    call yaml_map('epot',etot,fmt='(e24.15)')
    !-------------------------------------------------------
    call final_potential_forces(parini,atoms)
    deallocate(rxyz1,rxyz2)
    deallocate(fxyz)
    call atom_deallocate(atoms)
end subroutine bar_saddle
!*****************************************************************************************
