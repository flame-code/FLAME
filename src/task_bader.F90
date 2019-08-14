!*****************************************************************************************
subroutine task_bader(parini)
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    if(trim(parini%approach_bader)=='ongrid') then
        call bader_ongrid(parini)
    else if(trim(parini%approach_bader)=='neargrid') then
        call bader_neargrid(parini)
    else if(trim(parini%approach_bader)=='weight') then
        call bader_weight(parini)
    else
        write(*,'(2a)') 'ERROR: unknown parini%approach_bader ',trim(parini%approach_bader)
        stop
    endif
end subroutine task_bader
!*****************************************************************************************
