!*****************************************************************************************
subroutine init_lenosky_tb(atoms_t)
    use mod_potential, only: cell, natsi
    use mod_atoms, only: typ_atoms
    use mod_tightbinding, only: lenosky
    implicit none
    type(typ_atoms), intent(in):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat
    !character(5):: sat(nat)
    !local variables
    integer:: iat
    lenosky=.true.
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
   cell(1)=atoms_t%cellvec(1,1)
   cell(2)=atoms_t%cellvec(2,2)
   cell(3)=atoms_t%cellvec(3,3)
end subroutine init_lenosky_tb
!*****************************************************************************************
subroutine lenosky_tb(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    use mod_potential, only: cell, natsi
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat
    real(8):: count_t
    !if(iproc==0) call sleep(1)
    count_t=0.d0
    call update_ratp(atoms)
    do iat=1,atoms%nat
        atoms%ratp(1,iat)=modulo(modulo(atoms%ratp(1,iat),atoms%cellvec(1,1)),atoms%cellvec(1,1))
        atoms%ratp(2,iat)=modulo(modulo(atoms%ratp(2,iat),atoms%cellvec(2,2)),atoms%cellvec(2,2))
        atoms%ratp(3,iat)=modulo(modulo(atoms%ratp(3,iat),atoms%cellvec(3,3)),atoms%cellvec(3,3))
    enddo
    call update_rat(atoms,upall=.true.)
    call lenoskytb_alborz(parini,atoms,natsi,count_t)
end subroutine lenosky_tb
!*****************************************************************************************
