!*****************************************************************************************
subroutine miscellaneous_task(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    if(trim(parini%subtask_misc)=='solve_poisson') then
        call solve_poisson(parini)
    elseif(trim(parini%subtask_misc)=='linkedlist_test') then
        call linkedlist_test(parini)
    else
        write(*,'(2a)') 'ERROR: unknown parini%subtask_misc ',trim(parini%subtask_misc)
        stop
    endif
end subroutine miscellaneous_task
!*****************************************************************************************
