module interface_alborz
    implicit none
    integer:: icount_alborz=0
    character(5):: sat(1000)
end module interface_alborz
subroutine call_to_alborz_init(parini,nat)
    use mod_parini, only: typ_parini
    use interface_alborz
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: nat
    integer:: iat
    do iat=1,nat
        sat(iat)=trim(adjustl(parini%char_type(parini%typat_global(iat))))
    enddo
    call alborz_as_potential_init(nat,sat)
end subroutine call_to_alborz_init
subroutine call_to_alborz_get(boundcond,nat,latvec,xred,fcart,energy,strten)
    use interface_alborz
    implicit none
    character(*), intent(in):: boundcond
    integer, intent(in):: nat
    real(8), intent(in):: latvec(3,3), xred(3,nat)
    real(8), intent(inout):: fcart(3,nat), energy, strten(6)
    integer:: iat
    real(8), allocatable:: rat(:,:)
    real(8):: stress(3,3), vol
    allocate(rat(3,nat))
    call rxyz_int2cart_alborz(nat,latvec,xred,rat)
    call alborz_as_potential_get(boundcond,nat,latvec,rat,sat,fcart,energy,stress)
    strten(1) = -stress(1,1)
    strten(2) = -stress(2,2)
    strten(3) = -stress(3,3)
    strten(6) = -stress(2,1)
    strten(5) = -stress(3,1)
    strten(4) = -stress(3,2)
    call getvol_alborz(latvec,vol)
    strten=strten/vol
    deallocate(rat)
end subroutine call_to_alborz_get
subroutine call_to_alborz_final
    implicit none
    call alborz_as_potential_final
end subroutine call_to_alborz_final
