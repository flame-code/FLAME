!*****************************************************************************************
subroutine initprocessors
    use mod_processors, only: parallel, nproc, iproc, iproc_world, nproc_world, imaster, &
        iproc_type_all, iproc_list_abz, mpi_group_world, mpi_group_abz, mpi_comm_abz
    use yaml_output
    implicit none
    !local variables
    integer:: ierr, ios, nproc_t, jproc, kproc
    integer:: iproc_type=1
#if defined(MPI)
    include 'mpif.h'
    call MPI_INIT(ierr)
    CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD,iproc_world,ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc_world,ierr)
    call yaml_mapping_open('mpi started',flow=.true.)
    call yaml_map('iproc_world',iproc_world)
    call yaml_map('nproc_world',nproc_world)
    !write(*,'(a,2i4)') 'mpi started: iproc_world,nproc_world ',iproc_world,nproc_world
    if(iproc_world==imaster) then
        open(unit=21,file='master.dat',status='replace',iostat=ios)
        if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning master.dat';stop;endif
        write(21,*) imaster
        close(21)
    endif
    !allocate(iproc_type_all(0:nproc_world-1))
    iproc_type=1
    call MPI_ALLGATHER(iproc_type,1,MPI_INTEGER,iproc_type_all,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    !do jproc=0,nproc_world-1
    !    write(*,'(2(a,i3.3),5x,i)') 'proc',iproc_world,'PROC',jproc,iproc_type_all(jproc)
    !enddo
    nproc_t=0
    do jproc=0,nproc_world-1
        if(iproc_type_all(jproc)==1) nproc_t=nproc_t+1
    enddo
    !write(*,*) 'nproc_t ',nproc_t
    !allocate(iproc_list_abz(0:nproc_t-1))
    kproc=0
    do jproc=0,nproc_world-1
        if(iproc_type_all(jproc)==1) then
            iproc_list_abz(kproc)=jproc
            kproc=kproc+1
        endif
    enddo
    call MPI_COMM_GROUP(MPI_COMM_WORLD,mpi_group_world,ierr)
    call MPI_GROUP_INCL(mpi_group_world,nproc_t,iproc_list_abz,mpi_group_abz,ierr)
    call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_group_abz,mpi_comm_abz,ierr)
    call MPI_COMM_RANK(mpi_comm_abz,iproc,ierr)
    call MPI_COMM_SIZE(mpi_comm_abz,nproc,ierr)
    if(nproc/=nproc_t .or. nproc>nproc_world) then
        write(*,'(a)') 'ERROR: nproc/=nproc_t .or. nproc>nproc_world'
        call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
    endif
    !iproc=iproc_world
    !nproc=nproc_world
    call yaml_map('iproc',iproc)
    call yaml_map('nproc',nproc)
    call yaml_mapping_close()
    !write(*,'(a,2i4)') 'mpi started: iproc,nproc ',iproc,nproc
    parallel=.true. !parallel MPI set parallel=.true., for serial parallel=.false.
#else
    iproc_world=0
    nproc_world=1
    iproc=0
    nproc=1
    imaster=0
    return
#endif
end subroutine initprocessors
!*****************************************************************************************
subroutine finalizeprocessors
    use mod_processors, only: parallel, mpi_comm_abz, mpi_group_abz, mpi_group_world, &
        iproc_world, nproc_world
    use yaml_output
    implicit none
    !local variables
    integer:: ierr
#if defined(MPI)
    call MPI_COMM_FREE(mpi_comm_abz,ierr)
    call MPI_GROUP_FREE(mpi_group_abz,ierr)
    call MPI_GROUP_FREE(mpi_group_world,ierr)
    !deallocate(iproc_type_all,iproc_list_abz)
    call MPI_FINALIZE(ierr)
    call yaml_mapping_open('mpi finalized',flow=.true.)
    call yaml_map('iproc_world',iproc_world)
    call yaml_map('nproc_world',nproc_world)
    call yaml_mapping_close()
    !write(*,'(a,2i4)') 'mpi finalized: iproc_world,nproc_world ',iproc_world,nproc_world
#endif
end subroutine finalizeprocessors
!*****************************************************************************************
