!*****************************************************************************************
module mod_processors
    implicit none
    private
    public:: get_proc_stake
    integer, public:: mpi_group_world
    integer, public:: mpi_group_abz
    !integer:: mpi_comm_abz
    integer, public:: nproc_world=1 !number of processors, default is set for serial.
    integer, public:: iproc_world=0 !index of the current processor, default is set for serial.
    integer, public:: nproc_pot=-1
    integer, public:: iproc_pot=-1
    integer, public:: nproc=1 
    integer, public:: iproc=0 
    !integer, allocatable:: iproc_type_all(:)
    !integer, allocatable:: iproc_list_abz(:)
    integer, public:: iproc_type_all(0:999)
    integer, public:: iproc_list_abz(0:999)
    integer, public:: iproc_list_pot(0:999)
    integer, public:: imaster=0
    logical, public:: parallel=.true.
contains
!*****************************************************************************************
subroutine get_proc_stake(mpi_env,n,is,ie)
    use mod_flm_futile
    implicit none
    type(mpi_environment), intent(in):: mpi_env
    integer, intent(in):: n
    integer, intent(out):: is, ie
    !local variables
    integer:: m, mproc
    if(mpi_env%nproc>1) then
        m=n/mpi_env%nproc
        is=mpi_env%iproc*m+1
        mproc=mod(n,mpi_env%nproc)
        is=is+max(0,mpi_env%iproc-mpi_env%nproc+mproc)
        if(mpi_env%iproc>mpi_env%nproc-mproc-1) m=m+1
        ie=is+m-1
    else
        is=1
        ie=n
    endif
end subroutine get_proc_stake
!*****************************************************************************************
end module mod_processors
!*****************************************************************************************
