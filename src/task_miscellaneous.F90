!*****************************************************************************************
subroutine miscellaneous_task(parini)
    use mod_parini, only: typ_parini
    use mod_qat_target, only: get_qat_target
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    if(trim(parini%subtask_misc)=='solve_poisson') then
        call solve_poisson(parini)
    elseif(trim(parini%subtask_misc)=='linkedlist_test') then
        call linkedlist_test(parini)
    elseif(trim(parini%subtask_misc)=='fit_elecpot') then
        call subtask_fit_elecpot(parini)
    elseif(trim(parini%subtask_misc)=='test_free_bps') then
        call test_free_BPS(parini)
    elseif(trim(parini%subtask_misc)=='get_qat_target') then
        call get_qat_target(parini)
    else
        write(*,'(2a)') 'ERROR: unknown parini%subtask_misc ',trim(parini%subtask_misc)
        stop
    endif
end subroutine miscellaneous_task
!*****************************************************************************************
