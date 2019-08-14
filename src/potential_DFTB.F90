!*****************************************************************************************
subroutine init_potential_forces_dftb(atoms)
    !use mod_potential, only: cellvec, sat
    use mod_atoms, only: typ_atoms
    use mod_processors, only: iproc
    implicit none
    type(typ_atoms), intent(in):: atoms
    !real(8):: cv_t(3,3)
    !integer:: nat
    !character(5):: sat_t(nat)
    !local variables
    integer:: iat
    !character(7):: dir
    character(2):: dir
    character(50):: filename
    logical:: dftb_input_exists
    if(atoms%nat==0) stop 'ERROR: nat=0 in init_potential_forces_dftb'
    !write(dir,'(a4,i3.3)') 'proc',iproc
    dir='./'
    filename=dir//'/dftb_in.hsd'
    inquire(file=trim(filename),exist=dftb_input_exists)
    if(dftb_input_exists) then
        write(*,'(a)') 'potential dftb initialized.'
    else
        write(*,*) 'ERROR: DFTB input file ',trim(filename),' does not exist.'
        stop
    endif
    !cellvec(1:3,1:3)=atoms%cellvec(1:3,1:3)
    !if(allocated(sat)) then
    !    stop 'ERROR: sat already allocated in init_potential_forces_dftb'
    !else
    !    allocate(sat(atoms%nat))
    !endif
    !do iat=1,atoms%nat
    !    sat(iat)=trim(atoms%sat(iat))
    !enddo
end subroutine init_potential_forces_dftb
!*****************************************************************************************
subroutine final_potential_forces_dftb
    use mod_potential, only: sat
    implicit none
    !local variables
    write(*,'(a)') 'potential dftb finalized.'
    !if(.not. allocated(sat)) then
    !    stop 'ERROR: sat not allocated in final_potential_forces_dftb'
    !else
    !    deallocate(sat)
    !endif
end subroutine final_potential_forces_dftb
!*****************************************************************************************
subroutine cal_potential_forces_dftb(atoms)
    use mod_atoms, only: typ_atoms, set_typat
    use mod_const, only: bohr2ang
    !use mod_potential, only:
    use mod_processors, only: iproc
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    !character(7):: dir
    character(2):: dir
    character(50):: filename
    character(100):: command
    !real(8):: stress(3,3)
    integer:: iat, itry, i
    logical:: success
    call set_typat(atoms)
    !write(*,*) atoms%ntypat
    !write(*,*) atoms%stypat(1:atoms%ntypat)
    !write(*,*) atoms%ltypat(1:atoms%ntypat)
    !stop
    !write(dir,'(a4,i3.3)') 'proc',iproc
    dir='./'
    filename=dir//'/posinp.gen'
    !call writexyz_dftb(filename,atoms)
    !-------------------------------------------------------
    open(unit=87,file="input_geometry.gen")
    if(trim(atoms%boundcond)=='free') then
        write(87,'(i5,a)') atoms%nat, " C"
    elseif(trim(atoms%boundcond)=='wire' .or. trim(atoms%boundcond)=='slab' .or. &
        trim(atoms%boundcond)=='bulk') then
        write(87,'(i5,a)') atoms%nat, " S"
    else
        write(*,*) 'ERROR: unknown BC ',trim(atoms%boundcond)
        stop
    endif
    write(87,*) (trim(atoms%stypat(i))//" ", i=1,atoms%ntypat)
    if(trim(atoms%boundcond)=='free') then
        do iat=1,atoms%nat
            write(87,'(i5,1x,i5,3(1x,es25.15))') iat,atoms%itypat(iat),atoms%ratp(:,iat)*bohr2ang
        end do
    elseif(trim(atoms%boundcond)=='wire' .or. trim(atoms%boundcond)=='slab' .or. &
        trim(atoms%boundcond)=='bulk') then
        do iat=1,atoms%nat
            write(87,'(i5,1x,i5,3(1x,es25.15))') iat,atoms%itypat(iat),atoms%ratp(:,iat)*bohr2ang
        end do
        write(87,'(3(1x,es25.15))') 0.d0,0.d0,0.d0
        write(87,'(3(1x,es25.15))') atoms%cellvec(:,1)*bohr2ang
        write(87,'(3(1x,es25.15))') atoms%cellvec(:,2)*bohr2ang
        write(87,'(3(1x,es25.15))') atoms%cellvec(:,3)*bohr2ang
    else
        write(*,*) 'ERROR: unknown BC ',trim(atoms%boundcond)
        stop
    endif
    close(87)
    !-------------------------------------------------------
    !call system("sleep 1")
    !command='cd '//dir//'; ~/dftb+ ; cd ..'
    command="./runjob.sh"
    call system(trim(command))
    filename=dir//'/detailed.out'
    !write(*,*) filename
    call get_output_dftb_alborz(filename,atoms,success)
    !if(success) exit
    if(.not. success) then
        write(*,'(3a,i3.3)') 'ERROR: failed to find requested variables from ', &
            trim(filename),' by iproc= ',iproc
        stop
    endif
    !write(*,'(f20.10)') atoms%epot
    !do iat=1,atoms%nat
    !    write(81,'(3es24.15)') atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
    !enddo
    !stop
    !do iat=1,atoms%nat
    !    write(*,*) atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
end subroutine cal_potential_forces_dftb
!*****************************************************************************************
subroutine writexyz_dftb(filename,atoms)
    use mod_atoms, only: typ_atoms, get_rat_iat
    use mod_const, only: bohr2ang
    use mod_processors, only: iproc
    !use mod_potential, only: sat, cellvec
    implicit none
    character(*):: filename !,comment
    type(typ_atoms), intent(in):: atoms
    !local variables
    integer:: iat, i
    real(8):: x, y, z, cellx, celly, cellz, xyz(3)
    character(1):: bc
    !character(5):: an
    open(unit=1358,file=trim(filename),status='replace')
    if(trim(atoms%boundcond)=='free') then
        bc='C'
    else
        write(*,'(a,i3.3,a)') 'ERROR: proc= ',iproc,' since lattice vectors are not '
        write(*,'(a)') 'written in dftb gen file, currently only free BCs can be used '
        stop
        bc='F' !This is for bulk,slab,wire BCs
    endif
    write(1358,'(i6,2x,a1)')   atoms%nat,bc
    !write(1358,'(a)') '   1.66728000E+01   1.66728000E+01   1.66728000E+01'
    !cellx=cellvec(1,1) !*0.529177d0
    !celly=cellvec(2,2) !*0.529177d0
    !cellz=cellvec(3,3) !*0.529177d0
    !write(1358,'(a,3es20.10)') 'free ',cellx,celly,cellz
    do i=1,atoms%ntypat
        write(1358,'(1x,a,1x)',advance='no') trim(atoms%stypat(i))
    enddo
    write(1358,*)
    do iat=1,atoms%nat
        !if(nzx(imass(iat))==1 ) an=' H'
        !if(nzx(imass(iat))==14) an='Si'
        !an=' C '
        call get_rat_iat(atoms,iat,xyz)
        x=xyz(1)*bohr2ang !*0.529177d0
        y=xyz(2)*bohr2ang !*0.529177d0
        z=xyz(3)*bohr2ang !*0.529177d0
        write(1358,'(i5,i3,3x,3es23.14)') iat,atoms%itypat(iat),x,y,z
    enddo
    close(1358)
end subroutine writexyz_dftb
!*****************************************************************************************
subroutine get_output_dftb_alborz(filename,atoms,success)
    use mod_atoms, only: typ_atoms
    implicit none
    character(*), intent(in):: filename
    type(typ_atoms), intent(inout):: atoms
    logical, intent(out):: success
    !local variables
    integer:: iat, n, k, l, m, ios, nat_cell
    !real(8):: strten(6),latvec(3,3),xred(3,nat),str_matrix(3,3),vol,a(3,3),scaling
    character(250)::all_line
    logical:: got_energy, got_force
    got_energy=.false.
    got_force=.false.
    !atoms%fat(1:3,1:atoms%nat)=1.d10
    !strten=1.d10
    !nat_cell=0
    open(unit=32,file="detailed.out")
    open(unit=32,file=trim(filename),status='old',iostat=ios)
    if(ios/=0) then
        write(*,'(2a)') 'ERROR: failure openning ',trim(filename)
        stop
    endif
    do while(.true.)
        read(32,'(a250)',end=99) all_line
        !!write(*,*) all_line
        n=len_trim(all_line)
        k=index(all_line(1:n),"Total Mermin free energy")
        if(k.ne.0) then
            m=len_trim(all_line)
            l=scan(all_line(1:m),":",.true.)
            read(all_line(l+1:m),*) atoms%epot
            got_energy=.true.
            cycle
        endif
        k=index(all_line(1:n),"Total Forces")
        if(k.ne.0) then
            do iat=1,atoms%nat
                read(32,*) atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
            enddo
            got_force=.true.
            cycle
        endif
        !k=index(all_line(1:n),"Total stress tensor")
        !if(k.ne.0) then
        !    do iat=1,3
        !        read(32,*) str_matrix(:,iat)
        !    enddo
        !    strten(1)=-str_matrix(1,1)
        !    strten(2)=-str_matrix(2,2)
        !    strten(3)=-str_matrix(3,3)
        !    strten(6)=-str_matrix(1,2)
        !    strten(5)=-str_matrix(1,3)
        !    strten(4)=-str_matrix(2,3)
        !    cycle
        !endif
    enddo
    99 continue 
    if(got_energy .and. got_force) then
        success=.true.
    else
        success=.false.
    endif
    close(32)
    !if(bc==2) strten=0.d0
    !if(energy==1.d10 .or. strten(1)==1.d10.or.fcart(1,1)==1.d10) stop "Could not find all requested variables"
    !Transform all to bohr
    !energy=energy!/real(nat_cell,8)*real(nat,8)/Ha_eV
    !strten=strten!/HaBohr3_GPa
    !strten=0.d0
    !fcart=fcart!/Ha_eV*Bohr_Ang

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!    !Since its a single call, we only have forces from one configuration
!    implicit none
!    integer, intent(in):: iproc, nat
!    character(*):: filename
!    real(8):: fat(3,nat),epot
!    logical:: success
!    !local variables
!    integer:: ioserr, i, iat, nat_t
!    character(5):: ttsat
!    character(20):: units_t
!    character(250):: strline
!    epot=1.d10
!    fat=1.d10
!    iat=0
!    do i=1,2*nat+3
!        read(32,'(a)',end=99) strline
!        !write(*,'(a,i3,2x,a)') 'IPROC ',iproc,trim(strline)
!        if(i==1) then
!            read(strline,*) nat_t,units_t,epot
!            if(nat_t/=nat) stop 'ERROR: Has nat in forces_posinp.xyz wrong value?'
!        elseif(i>2+nat+1) then
!            iat=iat+1
!            read(strline,*) ttsat,fat(1,iat),fat(2,iat),fat(3,iat)
!        endif
!
!    enddo
!    99 continue
!    close(32)
!    !write(*,'(a,i3,2x,2es26.17)') 'IPROC ',iproc,epot,fat(3,nat)
!    if(epot==1.d10 .or.fat(3,nat)==1.d10) then
!        success=.false.
!    else
!        success=.true.
!        epot=epot*27.2114d0
!        do iat=1,nat
!            fat(1,iat)=fat(1,iat)*(27.2114d0/0.529177d0)
!            fat(2,iat)=fat(2,iat)*(27.2114d0/0.529177d0)
!            fat(3,iat)=fat(3,iat)*(27.2114d0/0.529177d0)
!        enddo
!    endif
end subroutine get_output_dftb_alborz
!*****************************************************************************************
