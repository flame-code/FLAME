!*****************************************************************************************
subroutine init_potential_forces_mpmd(atoms_t)
    use mod_potential, only: ntypat, natarr, stypat, cellvec, l_atomic_num
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat
    !character(5):: sat(nat)
    !local variables
    integer:: iat
    logical:: typeisnew
    integer:: itypat, ntot
#if !defined(MPI)
    stop 'ERROR: this potential works only with MPI while it is compiled without MPI'
#else
    if(atoms_t%nat==0) stop 'ERROR: nat=0 in init_potential_forces_mpmd'
    if(atoms_t%nat>1000) stop 'ERROR: nat>1000 in init_potential_forces_mpmd'
    cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
    ntypat=1
    stypat(1)=trim(atoms_t%sat(1))
    natarr(1)=1
    do iat=2,atoms_t%nat
        !write(*,*) 'TYPE ',trim(atoms_t%sat(iat)),trim(stypat(ntypat))
        typeisnew=.true.
        do itypat=1,ntypat
            if(trim(atoms_t%sat(iat))==trim(stypat(itypat))) then
                typeisnew=.false.
                exit
            endif
        enddo
        if(typeisnew) then
            ntypat=ntypat+1
            natarr(ntypat)=1
            stypat(ntypat)=trim(atoms_t%sat(iat))
            !comment=trim(comment)//" "//trim(atoms_t%sat(iat))
        else
            natarr(itypat)=natarr(itypat)+1
        endif
    enddo
    itypat=0
    ntot=0
    do iat=1,atoms_t%nat
        if(iat>ntot) then
            itypat=itypat+1
            ntot=ntot+natarr(itypat)
        endif
        !write(*,*) 'REZA ',iat,itypat
        l_atomic_num(iat)=itypat
    enddo
    call mpmd_init
#endif
end subroutine init_potential_forces_mpmd
!*****************************************************************************************
subroutine cal_potential_forces_mpmd(atoms)
    use mod_const
    use mod_atoms, only: typ_atoms, update_ratp, update_rat
    use mod_processors, only: iproc_pot
    use mod_potential, only: cellvec, l_atomic_num !, fcalls
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat, i, ierr
    integer:: ipbc=1
    integer:: failed
    character(1024):: cwd
    integer:: icwd(1024)
    integer(8):: nat_long
    logical:: flag
    real(8):: poll_period
#if defined(MPI)
    include 'mpif.h'
    integer:: status_mpi(MPI_STATUS_SIZE)
#endif
#if !defined(MPI)
    stop 'ERROR: this potential works only with MPI while it is compiled without MPI'
#else
    !if(iproc==0) call sleep(1)
    nat_long=int(atoms%nat,8)
    call getcwd(cwd)
    do i=1,1024
        icwd(i)=ichar(cwd(i:i))
    enddo
    cellvec = cellvec * bohr2ang 
    call update_ratp(atoms)
    atoms%ratp = atoms%ratp * bohr2ang
    call MPI_SEND(atoms%nat   ,1          ,MPI_INTEGER         ,iproc_pot,0,MPI_COMM_WORLD,ierr)
    call MPI_SEND(l_atomic_num,atoms%nat  ,MPI_INTEGER         ,iproc_pot,0,MPI_COMM_WORLD,ierr)
    call MPI_SEND(atoms%ratp   ,3*atoms%nat,MPI_DOUBLE_PRECISION,iproc_pot,0,MPI_COMM_WORLD,ierr)
    call MPI_SEND(cellvec     ,9          ,MPI_DOUBLE_PRECISION,iproc_pot,0,MPI_COMM_WORLD,ierr)
    call MPI_SEND(ipbc        ,1          ,MPI_INTEGER         ,iproc_pot,0,MPI_COMM_WORLD,ierr)
    call MPI_SEND(icwd        ,1024       ,MPI_INTEGER         ,iproc_pot,0,MPI_COMM_WORLD,ierr)

    !if (poll_period > 0.0) {
    !    while (MPI::COMM_WORLD.Iprobe(potentialRank, 0) == false) {
    !        usleep((useconds_t)(poll_period/1000000.0));
    !    }
    !}
    poll_period=1.d-2
    if(poll_period>0.d0) then
        do 
            call MPI_IPROBE(iproc_pot,0,MPI_COMM_WORLD,flag,status_mpi,ierr)
            if(flag) exit
            call wrapper_usleep(poll_period)
        enddo
    endif


    call MPI_RECV(failed     ,1          ,MPI_INTEGER         ,iproc_pot,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(atoms%epot ,1          ,MPI_DOUBLE_PRECISION,iproc_pot,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
    call MPI_RECV(atoms%fat  ,3*atoms%nat,MPI_DOUBLE_PRECISION,iproc_pot,0,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierr)
     atoms%fat = atoms%fat / (ha2ev/bohr2ang)
     atoms%epot = atoms%epot / ha2ev
     cellvec = cellvec / bohr2ang
     atoms%ratp = atoms%ratp/ bohr2ang
     call update_rat(atoms)

    !do iat=1,atoms%nat
    !    write(101+int(fcalls),'(3f14.6)') fat(1,iat),fat(2,iat),fat(3,iat)
    !enddo
    !do iat=1,atoms%nat
    !    write(201+int(fcalls),'(3f14.6)') atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
    !call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
#endif
end subroutine cal_potential_forces_mpmd
!*****************************************************************************************
subroutine mpmd_init
    use mod_processors, only: iproc_list_pot, iproc_type_all, nproc_world, mpi_group_world, &
        nproc_pot, iproc_pot, nproc, iproc
    implicit none
    integer:: jproc, ierr, kproc
    integer:: mpi_group_pot
    integer:: mpi_comm_pot
    integer:: ioffset, nsize_pot_group
#if defined(MPI)
    include 'mpif.h'
#endif
#if !defined(MPI)
    stop 'ERROR: this potential works only with MPI while it is compiled without MPI'
#else
    nproc_pot=0
    do jproc=0,nproc_world-1
        if(iproc_type_all(jproc)==2) nproc_pot=nproc_pot+1
    enddo
    if(nproc_pot/=nproc_world-nproc) then
        write(*,'(a)') 'ERROR: nproc_pot/=nproc_world-nproc'
        call MPI_ABORT(MPI_COMM_WORLD,0,ierr)
    endif
    kproc=0
    do jproc=0,nproc_world-1
        if(iproc_type_all(jproc)==2) then
            iproc_list_pot(kproc)=jproc
            kproc=kproc+1
        endif
    enddo
    if(mod(nproc_pot,nproc)/=0) then
        write(*,'(a)') 'WARNING: mod(nproc_pot,nproc)/=0'
    endif
    nsize_pot_group=nproc_pot/nproc
    write(*,'(a,2i4)') 'POTENTIAL: npotentials,potential_group_size ',nproc_pot,nsize_pot_group
    do jproc=0,nproc-1
        ioffset=jproc*nsize_pot_group
        call MPI_GROUP_INCL(mpi_group_world,nsize_pot_group,iproc_list_pot(ioffset),mpi_group_pot,ierr)
        call MPI_COMM_CREATE(MPI_COMM_WORLD,mpi_group_pot,mpi_comm_pot,ierr)
    enddo
    iproc_pot=iproc_list_pot(iproc*nsize_pot_group)
#endif
end subroutine mpmd_init
!*****************************************************************************************
