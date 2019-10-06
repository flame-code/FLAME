!*****************************************************************************************
subroutine init_potential_forces_bigdft(atoms_t)
    use mod_potential, only: cellvec, sat
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_atoms), intent(in):: atoms_t
    !real(8):: cv_t(3,3)
    !integer:: nat
    !character(5):: sat_t(nat)
    !local variables
    integer:: iat
    if(atoms_t%nat==0) stop 'ERROR: nat=0 in init_potential_forces_bigdft'
    cellvec(1:3,1:3)=atoms_t%cellvec(1:3,1:3)
    if(allocated(sat)) then
        stop 'ERROR: sat already allocated in init_potential_forces_bigdft'
    else
        allocate(sat(atoms_t%nat))
    endif
    do iat=1,atoms_t%nat
        sat(iat)=trim(atoms_t%sat(iat))
    enddo
end subroutine init_potential_forces_bigdft
!*****************************************************************************************
subroutine final_potential_forces_bigdft
    use mod_potential, only: sat
    implicit none
    !local variables
    if(.not. allocated(sat)) then
        stop 'ERROR: sat not allocated in final_potential_forces_bigdft'
    else
        deallocate(sat)
    endif
end subroutine final_potential_forces_bigdft
!*****************************************************************************************
subroutine cal_potential_forces_bigdft(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_potential, only: perfstatus, perfstatus_old, fcalls, bigdft_restart
    use mod_processors, only: iproc
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    character(6):: dir
    character(50):: filename
    character(100):: command
    !real(8):: stress(3,3)
    integer:: iat, itry
    logical:: success
    write(dir,'(a3,i3.3)') 'tmp',iproc
    filename=dir//'/posinp.xyz'
    call update_ratp(atoms)
    call writexyz_bigdft(filename,atoms%nat,atoms%ratp,'angstroemd0')
    if(fcalls<1.d0 .or. trim(perfstatus)/=trim(perfstatus_old)) then
        bigdft_restart=.false.
    else
        bigdft_restart=.true.
    endif
    perfstatus_old=trim(perfstatus)
    do itry=1,2
    if(trim(perfstatus)=='fast') then
        command='cp -f input.dft.fast '//dir//'/input.dft'
        call system(command)
    elseif(trim(perfstatus)=='normal') then
        command='cp -f input.dft.normal '//dir//'/input.dft'
        call system(command)
    elseif(trim(perfstatus)=='accurate') then
        command='cp -f input.dft.accurate '//dir//'/input.dft'
        call system(command)
    else
        stop 'ERROR: perfstatus is not set properly.'
    endif
    if(bigdft_restart) then
        command='sed -i "s/template1/2/g" '//dir//'/input.dft'
        call system(command)
    else
        command='sed -i "s/template1/0/g" '//dir//'/input.dft'
        call system(command)
    endif
    !command='cp KPOINTS '//dir//'/KPOINTS'
    !call system(command)
    call system("sleep 1")
    command='cd '//dir//'; ./brun ; cd ..'
    call system(command)
    call system("sleep 1")
    filename=dir//'/forces_posinp.xyz'
    !write(*,*) filename
    call get_output_bigdft(iproc,filename,atoms%nat,atoms%fat,atoms%epot,success)
    if(success) exit
    bigdft_restart=.false.
    enddo !end of loop over itry
    if(.not. success) then
        write(*,'(a)') 'ERROR: failed to find all requested variables from '
        write(*,'(a,i3)') '       forces_posinp.xyz by MHM iproc= ',iproc
        stop
    endif
    !do iat=1,atoms%nat
    !    write(*,*) atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
end subroutine cal_potential_forces_bigdft
!*****************************************************************************************
subroutine writexyz_bigdft(filename,nat,rat,comment)
    use mod_potential, only: sat, cellvec
    implicit none
    integer:: nat,iat
    real(8)::rat(3,nat),x,y,z,cellx,celly,cellz
    character(*)::filename,comment
    character(5)::an
    open(unit=1358, file=filename,status='replace')
    write(1358,'(i6,2x,a)')   nat,trim(comment)
    !write(1358,'(a)') '   1.66728000E+01   1.66728000E+01   1.66728000E+01'
    cellx=cellvec(1,1) !*0.529177d0
    celly=cellvec(2,2) !*0.529177d0
    cellz=cellvec(3,3) !*0.529177d0
    write(1358,'(a,3es20.10)') 'free ',cellx,celly,cellz
    do iat=1,nat
        !if(nzx(imass(iat))==1 ) an=' H'
        !if(nzx(imass(iat))==14) an='Si'
        !an=' C '
        x=rat(1,iat) !*0.529177d0
        y=rat(2,iat) !*0.529177d0
        z=rat(3,iat) !*0.529177d0
        write(1358,'(a5,2x,3es25.15)') trim(sat(iat)),x,y,z
    enddo
    close(1358)
end subroutine writexyz_bigdft
!*****************************************************************************************
subroutine get_output_bigdft(iproc,filename,nat,fat,epot,success)
    !Since its a single call, we only have forces from one configuration
    implicit none
    integer, intent(in):: iproc, nat
    character(*):: filename
    real(8):: fat(3,nat),epot
    logical:: success
    !local variables
    integer:: ioserr, i, iat, nat_t
    character(5):: ttsat
    character(20):: units_t
    character(250):: strline
    epot=1.d10
    fat=1.d10
    open(unit=32,file=filename,status='old',iostat=ioserr)
    if(ioserr/=0) stop 'ERROR: failure openning forces_posinp.xyz'
    iat=0
    do i=1,2*nat+3
        read(32,'(a)',end=99) strline
        !write(*,'(a,i3,2x,a)') 'IPROC ',iproc,trim(strline)
        if(i==1) then
            read(strline,*) nat_t,units_t,epot
            if(nat_t/=nat) stop 'ERROR: Has nat in forces_posinp.xyz wrong value?'
        elseif(i>2+nat+1) then
            iat=iat+1
            read(strline,*) ttsat,fat(1,iat),fat(2,iat),fat(3,iat)
        endif

    enddo
    99 continue
    close(32)
    !write(*,'(a,i3,2x,2es26.17)') 'IPROC ',iproc,epot,fat(3,nat)
    if(epot==1.d10 .or.fat(3,nat)==1.d10) then
        success=.false.
    else
        success=.true.
        epot=epot*27.2114d0
        do iat=1,nat
            fat(1,iat)=fat(1,iat)*(27.2114d0/0.529177d0)
            fat(2,iat)=fat(2,iat)*(27.2114d0/0.529177d0)
            fat(3,iat)=fat(3,iat)*(27.2114d0/0.529177d0)
        enddo
    endif
end subroutine get_output_bigdft
!*****************************************************************************************
