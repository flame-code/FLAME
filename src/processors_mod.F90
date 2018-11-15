!*****************************************************************************************
module mod_processors
    implicit none
    integer:: mpi_group_world
    integer:: mpi_group_abz
    integer:: mpi_comm_abz
    integer:: nproc_world=1 !number of processors, default is set for serial.
    integer:: iproc_world=0 !index of the current processor, default is set for serial.
    integer:: nproc_pot=-1
    integer:: iproc_pot=-1
    integer:: nproc=1 
    integer:: iproc=0 
    !integer, allocatable:: iproc_type_all(:)
    !integer, allocatable:: iproc_list_abz(:)
    integer:: iproc_type_all(0:999)
    integer:: iproc_list_abz(0:999)
    integer:: iproc_list_pot(0:999)
    integer:: imaster=0
    logical:: parallel=.true.
end module mod_processors
!*****************************************************************************************
