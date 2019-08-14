!*****************************************************************************************
program alborz
    use mod_task, only: typ_file_ini
    use mod_parini, only: typ_parini
    use mod_alborz_as_potential, only: parini_of_potential=>parini
    use yaml_output
    implicit none
    type(typ_file_ini):: file_ini
    type(typ_parini):: parini
    type(typ_parini):: parres
    call alborz_init(parini,parres,file_ini)
    call yaml_flush_document()
    !-----------------------------------------------------------------
    if(trim(parini%task)=='minhopp') then
        call minimahopping(parini)
    elseif(trim(parini%task)=='minhocao') then
        parini_of_potential=parini
        call task_minhocao(parini,parres)
    elseif(trim(parini%task)=='geopt') then
        call geopt(parini)
    elseif(trim(parini%task)=='saddle') then
        call task_saddle(parini)
    elseif(trim(parini%task)=='dynamics') then
        call dynamics(parini)
    elseif(trim(parini%task)=='conf_comp') then
        call conf_comp(parini)
    elseif(trim(parini%task)=='ann') then
        call task_ann(parini)
    elseif(trim(parini%task)=='genconf') then
        call task_genconf(parini)
    elseif(trim(parini%task)=='testforces') then
        call task_testforces(parini)
    elseif(trim(parini%task)=='single_point') then
        call single_point_task(parini)
    elseif(trim(parini%task)=='netsock') then
        call netsock_task(parini)
    elseif(trim(parini%task)=='bader') then
        call task_bader(parini)
    elseif(trim(parini%task)=='lammps') then
        call lammps_task(parini)
    elseif(trim(parini%task)=='phonon') then
        call cal_hessian_4p(parini)
    elseif(trim(parini%task)=='misc') then
        call miscellaneous_task(parini)
    elseif(trim(parini%task)=='potential') then
        write(*,'(a)') 'ERROR: task potential is not supposed to be called here'
        stop
    else
        write(*,'(2a)') 'ERROR: unknown parini%task name ',trim(parini%task)
    endif
    !-----------------------------------------------------------------
    call alborz_final(parini,file_ini)
end program alborz
!*****************************************************************************************
