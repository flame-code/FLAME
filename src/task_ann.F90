!*****************************************************************************************
subroutine task_ann(parini)
    use mod_parini, only: typ_parini
    use mod_train, only: ann_train
    use mod_fit_bf_cent2, only: get_basis_functions_cent2
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    if(trim(parini%subtask_ann)=='train') then
        call ann_train(parini)
    else if(trim(parini%subtask_ann)=='evaluate') then
        call ann_evaluate_subtask(parini)
    else if(trim(parini%subtask_ann)=='best_symfunc') then
        call ann_best_symfunc(parini)
    else if(trim(parini%subtask_ann)=='gen_symmetry_function') then
        call ann_gen_symmetry_function(parini)
    else if(trim(parini%subtask_ann)=='check_symmetry_function') then
        call ann_check_symmetry_function(parini)
    else if(trim(parini%subtask_ann)=='fit_bf_cent2') then
        call get_basis_functions_cent2(parini)
    else
        write(*,'(2a)') 'ERROR: unknown parini%subtask_ann ',trim(parini%subtask_ann)
        stop
    endif
end subroutine task_ann
!*****************************************************************************************
