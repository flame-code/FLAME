!*****************************************************************************************
module mod_trial_energy
    implicit none
    private
    public:: trial_energy_allocate_old, trial_energy_deallocate_old
    public:: trial_energy_copy_old
    public:: trial_energy_allocate, trial_energy_deallocate
    type, public:: typ_trial_energy
        private
        integer, public:: ntrial=0
        real(8), public:: ehartree_scn_excl
        real(8), allocatable, public:: energy(:)
        real(8), allocatable, public:: disp(:,:)
        integer, allocatable, public:: iat_list(:)
    end type typ_trial_energy
contains
!*****************************************************************************************
subroutine trial_energy_allocate_old(ntrial,trial_energy)
    implicit none
    integer, intent(in):: ntrial
    type(typ_trial_energy), pointer, intent(out):: trial_energy
    !local variables
    integer:: istat
    if(associated(trial_energy)) then
    if(allocated(trial_energy%energy)) stop 'ERROR: trial_energy%energy is already allocated'
    if(allocated(trial_energy%disp)) stop 'ERROR: trial_energy%disp is already allocated'
    if(allocated(trial_energy%iat_list)) stop 'ERROR: trial_energy%iat_list is already allocated'
    endif
    allocate(trial_energy,stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating trial_energy'
    trial_energy%ntrial=ntrial
    allocate(trial_energy%energy(ntrial))
    allocate(trial_energy%iat_list(ntrial))
    allocate(trial_energy%disp(3,ntrial))
end subroutine trial_energy_allocate_old
!*****************************************************************************************
subroutine trial_energy_deallocate_old(trial_energy)
    implicit none
    type(typ_trial_energy), pointer, intent(inout):: trial_energy
    !local variables
    integer:: istat
    if(allocated(trial_energy%energy)) deallocate(trial_energy%energy)
    if(allocated(trial_energy%disp)) deallocate(trial_energy%disp)
    if(allocated(trial_energy%iat_list)) deallocate(trial_energy%iat_list)
    deallocate(trial_energy,stat=istat)
    if(istat/=0) stop 'ERROR: failure deallocating trial_energy'
    nullify(trial_energy)
end subroutine trial_energy_deallocate_old
!*****************************************************************************************
subroutine trial_energy_copy_old(trial_energy_in,trial_energy_out)
    implicit none
    type(typ_trial_energy), pointer, intent(in):: trial_energy_in
    type(typ_trial_energy), pointer, intent(out):: trial_energy_out
    !local variables
    integer:: istat
    if(.not. associated(trial_energy_in)) then
        stop 'ERROR: trial_energy%energy is already allocated'
    endif
    if(associated(trial_energy_out)) then
        call trial_energy_deallocate_old(trial_energy_out)
    endif
    call trial_energy_allocate_old(trial_energy_in%ntrial,trial_energy_out)
    trial_energy_out%energy=trial_energy_in%energy
    trial_energy_out%disp=trial_energy_in%disp
    trial_energy_out%iat_list=trial_energy_in%iat_list
end subroutine trial_energy_copy_old
!*****************************************************************************************
subroutine trial_energy_allocate(ntrial,trial_energy)
    use dynamic_memory
    implicit none
    integer, intent(in):: ntrial
    type(typ_trial_energy), pointer, intent(out):: trial_energy
    !local variables
    integer:: istat
    if(associated(trial_energy)) then
    if(allocated(trial_energy%energy)) stop 'ERROR: trial_energy%energy is already allocated'
    if(allocated(trial_energy%disp)) stop 'ERROR: trial_energy%disp is already allocated'
    if(allocated(trial_energy%iat_list)) stop 'ERROR: trial_energy%iat_list is already allocated'
    endif
    allocate(trial_energy,stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating trial_energy'
    trial_energy%ntrial=ntrial
    trial_energy%energy=f_malloc0([1.to.ntrial],id='trial_energy%energy')
    trial_energy%iat_list=f_malloc0([1.to.ntrial],id='trial_energy%iat_list')
    trial_energy%disp=f_malloc0([1.to.3,1.to.ntrial],id='trial_energy%disp')
end subroutine trial_energy_allocate
!*****************************************************************************************
subroutine trial_energy_deallocate(trial_energy)
    use dynamic_memory
    implicit none
    type(typ_trial_energy), pointer, intent(inout):: trial_energy
    !local variables
    integer:: istat
    if(allocated(trial_energy%energy)) call f_free(trial_energy%energy)
    if(allocated(trial_energy%disp)) call f_free(trial_energy%disp)
    if(allocated(trial_energy%iat_list)) call f_free(trial_energy%iat_list)
    deallocate(trial_energy,stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating trial_energy'
    nullify(trial_energy)
end subroutine trial_energy_deallocate
!*****************************************************************************************
end module mod_trial_energy
!*****************************************************************************************
