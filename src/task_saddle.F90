!*****************************************************************************************
subroutine task_saddle(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    if(trim(parini%method_saddle)=='dimer') then
        call dimer_method(parini)
    else if(trim(parini%method_saddle)=='splined_saddle') then
        call splined_saddle(parini)
    else if(trim(parini%method_saddle)=='bar_saddle') then
        call bar_saddle(parini)
    else
        write(*,'(2a)') 'ERROR: unknown parini%method_saddle ',trim(parini%method_saddle)
        stop
    endif
end subroutine task_saddle
!*****************************************************************************************
