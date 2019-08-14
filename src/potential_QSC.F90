!*****************************************************************************************
subroutine init_potential_forces_qsc(atoms_t)
    use mod_potential, only: l_atomic_num
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat
    !character(5):: sat(nat)
    !local variables
    integer:: iat
    if(atoms_t%nat<1) stop 'ERROR: nat<1 in init_potential_forces_qsc'
    if(atoms_t%nat>1000) stop 'ERROR: nat>1000 in init_potential_forces_qsc'
    !if(.not. allocated(l_atomic_num)) allocate(l_atomic_num(nat))
    do iat=1,atoms_t%nat
        if(trim(atoms_t%sat(iat))=='Ni') then
            l_atomic_num(iat)=28
        elseif(trim(atoms_t%sat(iat))=='Cu') then 
            l_atomic_num(iat)=29
        elseif(trim(atoms_t%sat(iat))=='Rh') then 
            l_atomic_num(iat)=45
        elseif(trim(atoms_t%sat(iat))=='Pd') then 
            l_atomic_num(iat)=46
        elseif(trim(atoms_t%sat(iat))=='Ag') then 
            l_atomic_num(iat)=47
        elseif(trim(atoms_t%sat(iat))=='Ir') then 
            l_atomic_num(iat)=77
        elseif(trim(atoms_t%sat(iat))=='Pt') then 
            l_atomic_num(iat)=78
        elseif(trim(atoms_t%sat(iat))=='Au') then 
            l_atomic_num(iat)=79
        else
            write(*,'(2a)') 'ERROR: QSC parameters not available for element ',trim(atoms_t%sat(iat))
            stop
        endif
    enddo
    if(abs(atoms_t%cellvec(2,1))>1.d-12 .or. abs(atoms_t%cellvec(3,1))>1.d-12 .or. & 
       abs(atoms_t%cellvec(1,2))>1.d-12 .or. abs(atoms_t%cellvec(3,2))>1.d-12 .or. & 
       abs(atoms_t%cellvec(1,3))>1.d-12 .or. abs(atoms_t%cellvec(2,3))>1.d-12) then
        stop 'ERROR: QSC is only for Orthorhombic cells'
    endif
    !cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
end subroutine init_potential_forces_qsc
!*****************************************************************************************
subroutine cal_potential_forces_qsc(atoms)
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    use mod_potential, only: l_atomic_num
    use mod_const, only: bohr2ang, ev2ha
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat
    real(8):: cv(3,3), tt
    real(8), allocatable:: rat(:,:)
    call update_ratp(atoms)
    do iat=1,atoms%nat
        atoms%ratp(1,iat)=modulo(modulo(atoms%ratp(1,iat),atoms%cellvec(1,1)),atoms%cellvec(1,1))
        atoms%ratp(2,iat)=modulo(modulo(atoms%ratp(2,iat),atoms%cellvec(2,2)),atoms%cellvec(2,2))
        atoms%ratp(3,iat)=modulo(modulo(atoms%ratp(3,iat),atoms%cellvec(3,3)),atoms%cellvec(3,3))
    enddo
    call update_rat(atoms,upall=.true.)
    allocate(rat(3,atoms%nat))
    cv(1:3,1:3)=atoms%cellvec(1:3,1:3)*bohr2ang
    do iat=1,atoms%nat
        rat(1,iat)=atoms%ratp(1,iat)*bohr2ang
        rat(2,iat)=atoms%ratp(2,iat)*bohr2ang
        rat(3,iat)=atoms%ratp(3,iat)*bohr2ang
    enddo
    call qsc_energy_forces(atoms%nat,rat,l_atomic_num,atoms%fat,atoms%epot,cv);
    atoms%epot=atoms%epot*ev2ha
    tt=bohr2ang*ev2ha
    do iat=1,atoms%nat
        atoms%fat(1,iat)=atoms%fat(1,iat)*tt
        atoms%fat(2,iat)=atoms%fat(2,iat)*tt
        atoms%fat(3,iat)=atoms%fat(3,iat)*tt
    enddo
    deallocate(rat)
end subroutine cal_potential_forces_qsc
!*****************************************************************************************
subroutine final_potential_forces_qsc
    implicit none
    !if(allocated(l_atomic_num)) deallocate(l_atomic_num)
    call qsc_energy_forces_final();
end subroutine final_potential_forces_qsc
!*****************************************************************************************
