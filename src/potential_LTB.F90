!*****************************************************************************************
subroutine init_lenosky_tb(atoms_t)
    use mod_potential, only: cell, natsi
    use mod_atoms, only: typ_atoms
    use mod_const, only: bohr2ang
    implicit none
    type(typ_atoms), intent(in):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat
    !character(5):: sat(nat)
    !local variables
    integer:: iat
    if(atoms_t%nat==0) stop 'ERROR: nat=0 in init_lenosky_tb' 
    natsi=0
    do iat=1,atoms_t%nat
        if(trim(atoms_t%sat(iat))=='Si') then
            natsi=natsi+1
        elseif(trim(atoms_t%sat(iat))=='H') then
        elseif(trim(atoms_t%sat(iat))=='C') then
        else
            write(*,'(a)') 'ERROR: this is Lensoky TB for Si and H.'
            write(*,'(2a)') 'ERROR: the unacceptable element is ',trim(atoms_t%sat(iat))
            stop
        endif
    enddo
    if(abs(atoms_t%cellvec(2,1))>1.d-12 .or. abs(atoms_t%cellvec(3,1))>1.d-12 .or. & 
       abs(atoms_t%cellvec(1,2))>1.d-12 .or. abs(atoms_t%cellvec(3,2))>1.d-12 .or. & 
       abs(atoms_t%cellvec(1,3))>1.d-12 .or. abs(atoms_t%cellvec(2,3))>1.d-12) then
        stop 'ERROR: LenoskyTB is only for Orthorhombic cells'
   endif
   cell(1)=atoms_t%cellvec(1,1)*bohr2ang
   cell(2)=atoms_t%cellvec(2,2)*bohr2ang
   cell(3)=atoms_t%cellvec(3,3)*bohr2ang
end subroutine init_lenosky_tb
!*****************************************************************************************
subroutine lenosky_tb(atoms)
    use mod_interface
    use mod_atoms, only: typ_atoms
    use mod_potential, only: cell, natsi
    use mod_const, only: bohr2ang, ev2ha
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat
    real(8):: count_t, tt, cellvec(3,3)
    real(8), allocatable:: rat(:,:)
    !if(iproc==0) call sleep(1)
    count_t=0.d0
    do iat=1,atoms%nat
        atoms%rat(1,iat)=modulo(modulo(atoms%rat(1,iat),atoms%cellvec(1,1)),atoms%cellvec(1,1))
        atoms%rat(2,iat)=modulo(modulo(atoms%rat(2,iat),atoms%cellvec(2,2)),atoms%cellvec(2,2))
        atoms%rat(3,iat)=modulo(modulo(atoms%rat(3,iat),atoms%cellvec(3,3)),atoms%cellvec(3,3))
    enddo
    allocate(rat(3,atoms%nat))
    do iat=1,atoms%nat
        rat(1,iat)=atoms%rat(1,iat)
        rat(2,iat)=atoms%rat(2,iat)
        rat(3,iat)=atoms%rat(3,iat)
        atoms%rat(1,iat)=atoms%rat(1,iat)*bohr2ang
        atoms%rat(2,iat)=atoms%rat(2,iat)*bohr2ang
        atoms%rat(3,iat)=atoms%rat(3,iat)*bohr2ang
    enddo
    cellvec(1:3,1:3)=atoms%cellvec(1:3,1:3)
    atoms%cellvec(1:3,1:3)=atoms%cellvec(1:3,1:3)*bohr2ang
    call lenoskytb_alborz(atoms,natsi,count_t)
    write(11,*) natsi, count_t
    atoms%cellvec(1:3,1:3)=cellvec(1:3,1:3)
    atoms%epot=atoms%epot*ev2ha
    tt=bohr2ang*ev2ha
    do iat=1,atoms%nat
        atoms%rat(1,iat)=rat(1,iat)
        atoms%rat(2,iat)=rat(2,iat)
        atoms%rat(3,iat)=rat(3,iat)
        atoms%fat(1,iat)=atoms%fat(1,iat)*tt
        atoms%fat(2,iat)=atoms%fat(2,iat)*tt
        atoms%fat(3,iat)=atoms%fat(3,iat)*tt
    enddo
    deallocate(rat)
end subroutine lenosky_tb
!*****************************************************************************************
