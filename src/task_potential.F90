!*****************************************************************************************
subroutine alborz_as_potential_init(nat,sat)
    use mod_alborz_as_potential, only: parini, parres, file_ini, atoms
    use mod_potential, only: potential
    use mod_atoms, only: atom_allocate_old
    implicit none
    integer, intent(in):: nat
    character(5), intent(in):: sat(nat)
    !local variables
    integer:: iat
    if(trim(parini%task)/='minhocao') then
        call alborz_init(parini,parres,file_ini)
    endif
    potential=trim(parini%potential_potential)
    call atom_allocate_old(atoms,nat,0,0)
    do iat=1,atoms%nat
        atoms%sat(iat)=sat(iat)
    enddo

    call init_potential_forces(parini,atoms)
end subroutine alborz_as_potential_init
!*****************************************************************************************
subroutine alborz_as_potential_get(boundcond,nat,cellvec,rat,sat,fat,epot,stress)
    use mod_alborz_as_potential, only: parini, atoms
    use mod_atoms, only: set_rat
    implicit none
    character(*), intent(in):: boundcond
    integer, intent(in):: nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3)
    character(5), intent(in):: sat(nat)
    real(8), intent(out):: fat(3,nat), epot, stress(3,3)
    !local variables
    integer:: iat
    atoms%boundcond=trim(boundcond)
    if(atoms%nat/=nat) then
        write(*,'(a,2i6)') 'ERROR: atoms%nat=/nat in alborz_as_potential_get',atoms%nat,nat
    endif
    call set_rat(atoms,rat,setall=.true.)
    atoms%cellvec(1:3,1:3)=cellvec(1:3,1:3)
    call cal_potential_forces(parini,atoms)
    epot=atoms%epot
    do iat=1,nat
        fat(1,iat)=atoms%fat(1,iat)
        fat(2,iat)=atoms%fat(2,iat)
        fat(3,iat)=atoms%fat(3,iat)
    enddo
    stress(1:3,1:3)=atoms%stress(1:3,1:3)
end subroutine alborz_as_potential_get
!*****************************************************************************************
subroutine alborz_as_potential_final
    use mod_alborz_as_potential, only: parini, file_ini, atoms
    use mod_atoms, only: atom_deallocate_old
    implicit none
    call atom_deallocate_old(atoms)
    if(trim(parini%task)/='minhocao') then
        call alborz_final(parini,file_ini)
    endif
end subroutine alborz_as_potential_final
!*****************************************************************************************
