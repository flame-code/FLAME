!*****************************************************************************************
subroutine ann_evaluate_subtask(parini)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann
    use mod_symfunc, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann):: ann
    type(typ_symfunc_arr):: symfunc_arr
    integer:: iat, i
    real(8):: epot, epoti
    type(typ_atoms_arr):: atoms_arr
    !real(8), allocatable:: ratall(:,:,:), epotall(:)
    stop 'ERROR: this subtask is NOT ready yet.'
    !-------------------------------------------------------
    !call read_data_old(parini,'../list_posinp_test',atoms_arr)
    !call read_ann(parini,'ann.param',ann)
!HERE    call ann_evaluate(0,ann,symfunc_arr,atoms_arr,6,.true.)
    !-------------------------------------------------------
end subroutine ann_evaluate_subtask
!*****************************************************************************************
