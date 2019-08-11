!*****************************************************************************************
module mod_velocity
    implicit none
    private
    public:: set_velocities
contains
!*****************************************************************************************
subroutine set_velocities(atoms, ekin_arg)
    use mod_atoms, only: typ_atoms, update_ratp
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: i, iat
    real(8):: ekin, sclvel, t1,ekin_target
    real(8), optional::ekin_arg
    call random_number(atoms%vat)
    atoms%vat(1:3,1:atoms%nat)=atoms%vat(1:3,1:atoms%nat)-0.5d0

    !write(*,*) atoms%vat(1:3,1)
    !write(*,*) atoms%vat(1:3,2)
    !write(*,*) atoms%vat(1:3,3)
    !vcm(1:3)=0.d0
    !do iat=1,atoms%nat
    !    t1=atoms%amass(iat)
    !    vcm(1:3)=vcm(1:3)+t1*atoms%vat(1:3,iat)
    !enddo
    !vcm(1:3)=vcm(1:3)/totmass
    !do iat=1,atoms%nat
    !    atoms%vat(1:3,iat)=atoms%vat(1:3,iat)-vcm(1:3)
    !enddo
    !nmd=100
    if( present(ekin_arg)) then
        ekin_target=ekin_arg
    else
       ekin_target=0.5d0*8.d-1*atoms%nat
    endif
    if( present(ekin_arg)) then
        call elim_moment_mass(atoms%nat,atoms%vat,atoms%amass)
    else
        call elim_moment_alborz(atoms%nat,atoms%vat)
        call update_ratp(atoms)
        call elim_torque_reza_alborz(atoms%nat,atoms%ratp,atoms%vat)
    endif
    do i=1,10
        ekin=0.d0
        do iat=1,atoms%nat
            t1=atoms%amass(iat)
            ekin=ekin+t1*(atoms%vat(1,iat)**2+atoms%vat(2,iat)**2+atoms%vat(3,iat)**2)
        enddo
        ekin=0.5d0*ekin
        sclvel=dsqrt(ekin_target/ekin)
        atoms%vat(1:3,1:atoms%nat)=atoms%vat(1:3,1:atoms%nat)*sclvel
        if( present(ekin_arg)) then
            call elim_moment_mass(atoms%nat,atoms%vat,atoms%amass)
        else
            call elim_moment_alborz(atoms%nat,atoms%vat)
            call update_ratp(atoms)
            call elim_torque_reza_alborz(atoms%nat,atoms%ratp,atoms%vat)
        endif
        write(*,*) 'i,ekin ',i,ekin
    enddo
end subroutine set_velocities
!*****************************************************************************************
end module mod_velocity
!*****************************************************************************************
