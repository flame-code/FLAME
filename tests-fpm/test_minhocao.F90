!*****************************************************************************************
subroutine test_minhocao
    use iso_fortran_env, only: error_unit, output_unit
    use mod_parini, only: typ_parini
    use mod_flm_futile
    use mod_colors, only: green_passed, red_failed
    implicit none
    type(typ_parini):: parini, parres
    integer:: nat, iconf
    real(kind=8):: strten(6), printval1, printval2, errmax
    real(kind=8):: cellvec(3,3), cellvec_ref(3,3)
    logical:: fixlat(7), readfix, readfrag, file_exists
    character(len=256):: filename, command
    integer, allocatable:: fragarr(:)
    real(kind=8), allocatable:: xred(:,:), fcart(:,:)
    real(kind=8), allocatable:: rat(:,:), rat_ref(:,:)
    logical, allocatable:: fixat(:)
    call set_dict_parini_default(parini)
    parini%dict_user=>dict_new()
    call dict_set(parini%dict_user//'main'//'task','minhocao')
    call dict_set(parini%dict_user//'main'//'types','C Si')
    call dict_set(parini%dict_user//'main'//'nat',6)
    call dict_set(parini%dict_user//'main'//'typat','1*1 5*2')
    call dict_set(parini%dict_user//'main'//'pressure',5.d0)
    call dict_set(parini%dict_user//'main'//'verbose',0)
    !call dict_set(parini%dict_user//'main'//'rng_type','only_for_tests')
    call dict_set(parini%dict_user//'geopt'//'nit',2000)
    call dict_set(parini%dict_user//'geopt'//'method','fire')
    call dict_set(parini%dict_user//'geopt'//'dt_start',10.d0)
    call dict_set(parini%dict_user//'geopt'//'fmaxtol',1.d-5)
    call dict_set(parini%dict_user//'potential'//'potential','tersoff')
    call dict_update(parini%dict,parini%dict_user)
    call yaml_dict_dump(parini%dict)
    call yaml_get_all_parameters(parini)
    call dict_free(parini%dict)
    nullify(parini%dict)
    call dict_free(parini%dict_user)
    nullify(parini%dict_user)
    !-------------------------------------------------------
    parini%mpi_env%nproc=1
    parini%mpi_env%iproc=0
    parini%nmd_dynamics=1000
    inquire(file='global.mon',exist=file_exists)
    if(file_exists) then
        call execute_command_line('rm -f global.mon')
    endif
    call minhocao_write_input_files_for_test()
    call task_minhocao(parini,parres)
    nat=6
    allocate(fragarr(nat))
    allocate(xred(3,nat))
    allocate(rat(3,nat),rat_ref(3,nat))
    allocate(fcart(3,nat))
    allocate(fixat(nat))

    call minhocao_check_outputs_for_test()
    !-----------------------------------------------------------------
    filename='tests-fpm/data/minhocao_C1Si5_poscur.ascii'
    call read_atomic_file_ascii(trim(filename),nat,'bohr',xred,cellvec_ref,fcart, \
        strten,fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
    call rxyz_int2cart_alborz(nat,cellvec_ref,xred,rat_ref)
    filename='poscur.ascii'
    call read_atomic_file_ascii(trim(filename),nat,'bohr',xred,cellvec,fcart, \
        strten,fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
    call rxyz_int2cart_alborz(nat,cellvec,xred,rat)
    errmax=maxval(abs(rat_ref(1:3,1:nat)-rat(1:3,1:nat)))
    if(errmax<1.d-6) then
        write(output_unit,'(2a)') green_passed,' in test_minhocao: poscur.ascii'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_minhocao for poscur.ascii,   errmax=',errmax
        call exit(1)
    end if
    !-----------------------------------------------------------------
    do iconf=1,4
    write(filename,'(a,i1,a)') 'tests-fpm/data/minhocao_C1Si5_poslow0000',iconf,'.ascii'
    call read_atomic_file_ascii(trim(filename),nat,'bohr',xred,cellvec_ref,fcart, \
        strten,fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
    call rxyz_int2cart_alborz(nat,cellvec_ref,xred,rat_ref)
    write(filename,'(a,i1,a)') 'poslow0000',iconf,'.ascii'
    call read_atomic_file_ascii(trim(filename),nat,'bohr',xred,cellvec,fcart, \
        strten,fixat,fixlat,readfix,fragarr,readfrag,printval1,printval2)
    call rxyz_int2cart_alborz(nat,cellvec,xred,rat)
    !write(*,*) 'AAA',sqrt(sum((rat_all_ref(1:3,1:nat,iconf)-rat_all(1:3,1:nat,iconf))**2))
    !write(*,*) 'AAA',maxval(abs(rat_all_ref(1:3,1:nat,iconf)-rat_all(1:3,1:nat,iconf)))
    errmax=maxval(abs(rat_ref(1:3,1:nat)-rat(1:3,1:nat)))
    if(errmax<1.d-6) then
        write(output_unit,'(2a,i1)') green_passed,' in test_minhocao: poslow0000',iconf
    else
        write(error_unit,'(2a,i1,a,es14.5)') red_failed,' in test_minhocao for poslow0000',iconf,'     errmax=',errmax
        call exit(1)
    end if
    enddo
    !-----------------------------------------------------------------
    do iconf=1,4
        write(filename,'(a,i1,a)') 'poslow0000',iconf,'.ascii'
        command='rm -f '//trim(filename)
        call execute_command_line(trim(command))
    enddo
    call execute_command_line('rm -f strten.bin')
    call execute_command_line('rm -f poslow.bin')
    call execute_command_line('rm -f poscur.ascii')
    call execute_command_line('rm -f latvec.bin')
    call execute_command_line('rm -f ioput')
    call execute_command_line('rm -f global.mon')
    call execute_command_line('rm -f fp.bin')
    call execute_command_line('rm -f fcart.bin')
    call execute_command_line('rm -f earr.dat')
end subroutine test_minhocao
!*****************************************************************************************
subroutine minhocao_write_input_files_for_test()
    use mod_flm_futile
    implicit none
    !local variables
    integer:: iunit
    iunit=f_get_free_unit(10**5)
    open(unit=iunit,file='poscur.ascii',status='replace')
    write(iunit,'(a)') ' 6'
    write(iunit,'(a)') ' 3.241 -1.622  5.735'
    write(iunit,'(a)') '-1.617 -0.956  6.097'
    write(iunit,'(a)') ' 0.655  2.657  0.181  C'
    write(iunit,'(a)') '-0.967  5.344  0.180 Si'
    write(iunit,'(a)') '-0.966  3.154  0.955 Si'
    write(iunit,'(a)') ' 0.658 -0.706  4.595 Si'
    write(iunit,'(a)') '-0.964  3.068  3.264 Si'
    write(iunit,'(a)') '-0.962  1.102  4.548 Si'
    close(iunit)
    open(unit=iunit,file='ioput',status='replace')
    write(iunit,'(a)') '7.E-03  3000.0 5.E7  ediff, temperature, maximal temperature'
    close(iunit)
    open(unit=iunit,file='earr.dat',status='replace')
    write(iunit,'(a)') ' 0  4 #No. of minima already found, no. of minima to be found'
    write(iunit,'(a)') ' 1.E-04  0.15  #delta_enthalpy, delta_fingerprint'
    close(iunit)
end subroutine minhocao_write_input_files_for_test
!*****************************************************************************************
subroutine minhocao_check_outputs_for_test()
    use iso_fortran_env, only: error_unit, output_unit
    use mod_colors, only: green_passed, red_failed
    use mod_flm_futile
    implicit none
    !local variables
    integer:: iunit, iline, ii, iline_errmax, lstr, iline_difference
    real(kind=8):: rr, ener, ener_ref, errmax
    character(len=256):: filename, str_global(6), str_global_ref(6)
    character(len=256):: str_line, str_g, str_g_ref
    logical:: global_mon
    iunit=f_get_free_unit(10**5)
    filename='tests-fpm/data/minhocao_C1Si5_global.mon'
    open(unit=iunit,file=trim(filename),status='old')
    do iline=1,6
        read(iunit,'(a)') str_global_ref(iline)
    enddo
    close(iunit)
    filename='global.mon'
    open(unit=iunit,file=trim(filename),status='old')
    do iline=1,6
        read(iunit,'(a)') str_global(iline)
    enddo
    close(iunit)
    global_mon=.true.
    iline_difference=0
    do iline=1,6
        str_line=str_global(iline)
        lstr=len_trim(str_line)
        str_g=str_line(46:lstr)
        str_line=str_global_ref(iline)
        lstr=len_trim(str_line)
        str_g_ref=str_line(46:lstr)
        if(trim(str_g)/=trim(str_g_ref)) then
            global_mon=.false.
            iline_difference=iline
            exit
        endif
    enddo
    if(global_mon) then
        write(output_unit,'(2a)') green_passed,' in test_minhocao for global.mon right side'
    else
        write(error_unit,'(2a,i1)') red_failed,' in test_minhocao for global.mon right side iline= ',iline_difference
        call exit(1)
    end if
    errmax=0.d0
    iline_errmax=0
    do iline=1,6
        read(str_global(iline),*) ii,rr,ener
        read(str_global_ref(iline),*) ii,rr,ener_ref
        if(errmax<abs(ener-ener_ref)) then
            errmax=abs(ener-ener_ref)
            iline_errmax=iline
        endif
    enddo
    if(errmax<1.d-7) then
        write(output_unit,'(2a)') green_passed,' in test_minhocao for global.mon'
    else
        write(error_unit,'(2a,i1,a,es14.5)') red_failed,' in test_minhocao for global.mon iline',iline_errmax,'     errmax=',errmax
        call exit(1)
    end if
end subroutine minhocao_check_outputs_for_test
!*****************************************************************************************
