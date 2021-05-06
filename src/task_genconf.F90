!*****************************************************************************************
subroutine task_genconf(parini)
    use mod_genconf, only: typ_genconf
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    type(typ_genconf):: genconf

    genconf%cal_pot=parini%cal_pot_genconf
    genconf%nat_add=parini%nat_add_genconf
    if(trim(parini%subtask_genconf)=='diatomic') then
        call genconf_diatomic(parini,genconf)
    else if(trim(parini%subtask_genconf)=='trimer') then
        call genconf_trimer(parini,genconf)
    else if(trim(parini%subtask_genconf)=='genrandom') then
        call genrandom(parini,genconf)
    else if(trim(parini%subtask_genconf)=='rangrow') then
        call rangrow(parini,genconf)
    else
        write(*,'(2a)') 'ERROR: unknown parini%subtask_genconf, parini%subtask_genconf=', &
            trim(parini%subtask_genconf)
        stop
    endif
end subroutine task_genconf
!*****************************************************************************************
