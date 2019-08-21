!*****************************************************************************************
subroutine alborz_init(parini,parres,file_ini)
    use mod_processors, only: iproc, mpi_comm_abz, imaster
    use mod_task, only: typ_file_ini, time_start
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: read_file_input
    use futile
#ifndef __GFORTRAN__ 
    use ifport
#endif
    use time_profiling
    !use mod_timing , only: TCAT_ALBORZ_INIT_FINAL
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    type(typ_parini), intent(inout):: parres
    !local variables
    integer:: istat, ierr
    character(len=*), parameter:: filename='flame_log.yaml'
    logical:: flib_profiling
    call f_lib_initialize()
    inquire(file="NO_FLIB_PROFILING",exist=flib_profiling)
    if(flib_profiling) then
        call f_malloc_set_status()!profiling_depth=0)
    endif
    call f_routine(id='alborz_init')
    !-----------------------------------------------------------------
    !call alborz_initialize_timing_categories()
    !call f_timing(TCAT_ALBORZ_INIT_FINAL,'ON')
    !-----------------------------------------------------------------
    parini%iunit=f_get_free_unit(10**5)
    call yaml_set_stream(unit=parini%iunit,filename=trim(filename),&
         record_length=92,istat=ierr,setdefault=.false.,tabbing=0,position='rewind')
    if (ierr/=0) then
       call yaml_warning('Failed to create'//trim(filename)//', error code='//trim(yaml_toa(ierr)))
    end if
    call yaml_release_document(parini%iunit)
    call yaml_set_default_stream(parini%iunit,ierr)
    !call yaml_get_default_stream(unit_log)
    call yaml_new_document()
    !call yaml_invoice_example()
    call flm_print_logo(parini)
    !-----------------------------------------------------------------
    istat=getcwd(parini%cwd)
    if(istat/=0) stop 'ERROR: could not get CWD'
    call cpu_time(time_start)
    !parsing all blocks in input.ini
    inquire(file="flame_in.yaml",exist=parini%exists_yaml_file)
    if(parini%exists_yaml_file) then
        call yaml_get_parameters(parini)
    else
        call yaml_comment('flame_in.yaml does not exists, trying input.ini')
        allocate(file_ini%file_lines(file_ini%nline_max)) !,comment_line(nline_max))
        file_ini%file_lines(1:file_ini%nline_max)=' '
        !reading file input.ini into array file_lines which later it will be parsed.
        call read_file_input(file_ini)
        allocate(file_ini%stat_line_is_read(file_ini%nline),source=.false.)
        call get_main_parameters(file_ini,parini)
        call set_atomc_types_info(parini)
        call get_minhopp_parameters(file_ini,parini)
        call get_geopt_parameters(file_ini,parini)
        call get_geopt_prec_parameters(file_ini,parini)
        call get_saddle_opt_parameters(file_ini,parini)
        call get_saddle_parameters(file_ini,parini)
        call get_potential_parameters(file_ini,parini)
        call get_ann_parameters(file_ini,parini)
        call get_dynamics_parameters(file_ini,parini)
        call get_bader_parameters(file_ini,parini)
        call get_genconf_parameters(file_ini,parini)
        call get_conf_comp_parameters(file_ini,parini)
        call get_testforces_parameters(file_ini,parini)
        call get_single_point_parameters(file_ini,parini)
        call get_ewald_parameters(file_ini,parini)
        call get_misc_parameters(file_ini,parini)
    endif
    if(trim(parini%task)/='minhocao') then
        call initprocessors !start MPI in parallel version.
    endif
    if(trim(parini%task)/='potential') then
        call init_random_seed(parini)
    endif
    parres=parini
    !-----------------------------------------------------------------
    !call f_timing(TCAT_ALBORZ_INIT_FINAL,'OF')
    call f_release_routine()
end subroutine alborz_init
!*****************************************************************************************
!subroutine alborz_initialize_timing_categories
!    use yaml_output
!    use time_profiling
!    use mod_timing, only: dict_timing_info, TCAT_PSOLVER, TCAT_SYMFUNC_COMPUT
!    use mod_timing, only: TCAT_ALBORZ_INIT_FINAL
!    use wrapper_mpi, only: mpi_initialize_timing_categories, wmpi_init_thread, mpifinalize
!    use wrapper_mpi, only: mpisize, mpirank
!    use mod_processors, only: iproc
!    implicit none
!    integer:: iverbose=3
!    character(len=*), parameter :: pscpt1='Alborz init-fin'
!    character(len=*), parameter :: pscpt2='potential'
!    call f_timing_category_group(pscpt1,'Alborz computations total time')
!    call f_timing_category_group(pscpt2,'potential computations total time')
!    !-----------------------------------------------------------------
!    call f_timing_category('Alborz initialize-finalize',pscpt1,'computation of init-fin',TCAT_ALBORZ_INIT_FINAL)
!    !-----------------------------------------------------------------
!    call f_timing_category('psolver',pscpt2,'computation of Psolver',TCAT_PSOLVER)
!    call f_timing_category('symfunc compute',pscpt2,'computation of symmetry functions',TCAT_SYMFUNC_COMPUT)
!    !-----------------------------------------------------------------
!    call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=iverbose>2)
!end subroutine alborz_initialize_timing_categories
!*****************************************************************************************
subroutine alborz_final(parini,file_ini)
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini, time_start, time_end
    use mod_processors, only: iproc
    use yaml_output
    use time_profiling
    !use mod_timing , only: dict_timing_info, TCAT_ALBORZ_INIT_FINAL
    use dynamic_memory
    implicit none
    type(typ_parini), intent(inout):: parini
    type(typ_file_ini), intent(inout):: file_ini
    !local variables
    call f_routine(id='alborz_final')
    call f_free(parini%qt)
    call f_free(parini%at)
    !call f_timing(TCAT_ALBORZ_INIT_FINAL,'ON')
    if(.not. parini%exists_yaml_file) then
        deallocate(file_ini%file_lines,file_ini%stat_line_is_read) !,comment_line)
    endif
    call cpu_time(time_end)
    call yaml_mapping_open('CPU time',flow=.true.)
    call yaml_map('iproc',iproc,fmt='(i4)')
    call yaml_map('hrs',(time_end-time_start)/3600.d0,fmt='(e15.3)')
    call yaml_map('min',(time_end-time_start)/60.d0,fmt='(e15.3)')
    call yaml_map('sec',(time_end-time_start),fmt='(e15.3)')
    call yaml_mapping_close()
    !write(*,'(a,1x,i4,e15.3)') 'CPU time: iproc,time(hrs)',iproc,(time_end-time_start)/3600.d0
    !write(*,'(a,1x,i4,e15.3)') 'CPU time: iproc,time(min)',iproc,(time_end-time_start)/60.d0
    !write(*,'(a,1x,i4,e15.3)') 'CPU time: iproc,time(sec)',iproc,(time_end-time_start)
    if(trim(parini%task)/='minhocao') then
        call finalizeprocessors
    endif
    !-----------------------------------------------------------------
    !call f_timing(TCAT_ALBORZ_INIT_FINAL,'OF')
    !call f_timing_stop(dict_info=dict_timing_info)
    !-----------------------------------------------------------------
    call yaml_flush_document(unit=parini%iunit)
    !ASK LUIGI OR DAMIEN ABOUT THE PUZZLING ISSUE OF WHICH MUST COME FIRST:
    !yaml_close_stream or f_lib_finalize
    !call yaml_close_stream(unit=parini%iunit)
    call f_release_routine()
    call f_lib_finalize()
end subroutine alborz_final
!*****************************************************************************************
subroutine init_random_seed(parini)
    use mod_processors, only: iproc, mpi_comm_abz, imaster
    use mod_parini, only: typ_parini
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer:: istat, nseed, ierr, i
    integer, allocatable:: iseed(:)
#if defined(MPI)
    include 'mpif.h'
#endif
    call random_seed(size=nseed)
    allocate(iseed(nseed))
    if(iproc==0) then
        if(parini%iseed>=0) then
            !seed is set in input.ini
            iseed(1)=parini%iseed
            if(nseed==2) iseed(2)=iseed(1)+7
            if(nseed==3) iseed(3)=iseed(2)+17
        else if(parini%iseed==-1) then
            !seed must be set using system time as requested in input.ini
            call random_seed
            call random_seed(get=iseed)
        else
            !using default values
            call random_seed(get=iseed)
        endif
    endif
#if defined(MPI)
    if(trim(parini%task)/='minhocao') then 
         call MPI_BARRIER(mpi_comm_abz,ierr)
         call MPI_BCAST(iseed,nseed,MPI_INTEGER,imaster,mpi_comm_abz,ierr)
    endif
#endif
    iseed(1:nseed)=iseed(1:nseed)+iproc*11
    call random_seed(put=iseed)
    call yaml_mapping_open('random number generator',flow=.true.)
    call yaml_map('iproc',iproc,fmt='(i4)')
    if(trim(parini%rng_type)/='only_for_tests') then
        call yaml_map('seed',iseed) !,fmt='(i12)')
    endif
    call yaml_mapping_close()
    deallocate(iseed)
end subroutine init_random_seed
!*****************************************************************************************
subroutine set_atomc_types_info(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: sat_to_iatom
    use yaml_output
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: ntypat, iatomnum(20), i, j, ind, ind_tmp(1)
    character(5):: stypat(20)='unknown'
    if(trim(parini%types_main)=='unknown') then
        write(*,'(2a)') 'ERROR: keyword types is missing in block [main]: ', &
            trim(parini%types_main)
        stop
    endif
    iatomnum(1:20)=1000
    call count_words(parini%types_main,ntypat)
    read(parini%types_main,*) stypat(1:ntypat)
    do i=1,ntypat
        call sat_to_iatom(stypat(i),iatomnum(i))
    enddo
    !-------------------------------------------------------
    !The following commented lines were used in the past
    !when type of atoms were sorted according to atom number.
    !do i=1,ntypat
    !    ind_tmp(1:1)=minloc(iatomnum)
    !    ind=ind_tmp(1)
    !    do j=1,i-1
    !        if(stypat(ind)==parini%stypat(j)) then
    !            write(*,'(a,2a5)') 'ERROR: repeated types in types of block [main]:', &
    !                stypat(ind),parini%stypat(j)
    !            stop
    !        endif
    !    enddo
    !    parini%ltypat(i)=i
    !    parini%iatomnum(i)=iatomnum(ind)
    !    parini%stypat(i)=stypat(ind)
    !    iatomnum(ind)=1000
    !enddo
    !-------------------------------------------------------
    do i=1,ntypat
        parini%ltypat(i)=i
        parini%iatomnum(i)=iatomnum(i)
        parini%stypat(i)=stypat(i)
    enddo
    do i=1,ntypat
        if(parini%iatomnum(i)==-1) stop 'ERROR: atomic number not set properly.'
    enddo
    parini%ntypat=ntypat
    call yaml_mapping_open('system info',flow=.true.)
    call yaml_map('ntypat',parini%ntypat,fmt='(i2)')
    do i=1,ntypat
        call yaml_map(trim(parini%stypat(i)),parini%ltypat(i),fmt='(i2)')
    enddo
    call yaml_mapping_close()
    !write(*,'(a,i2)') 'ntypat ',parini%ntypat
    !do i=1,ntypat
    !    write(*,'(a5,i2,a5)') 'type ',parini%ltypat(i),trim(parini%stypat(i))
    !enddo
end subroutine set_atomc_types_info
!*****************************************************************************************
subroutine flm_print_logo(parini)
    use mod_parini, only: typ_parini
    use futile
    implicit none
    type(typ_parini), intent(inout):: parini
    call yaml_mapping_open('Code logo')
    call yaml_scalar('"__________________ Fully-Loaded Atomistic Modeling Environment')
    call yaml_scalar('        .                                                     ')     
    call yaml_scalar('       .M                                                     ')
    call yaml_scalar('      ,MM                                                     ')
    call yaml_scalar('      MM:                                                     ')
    call yaml_scalar('  .   YMM,                                                    ')
    call yaml_scalar('  M   `MMM,     .                                             ')
    call yaml_scalar('  M.   `MMM    .M                                             ')
    call yaml_scalar('  MM,  ,MMM   ,MM                                             ')
    call yaml_scalar('  `MM, MMM`  ,MM` .                                           ')
    call yaml_scalar('  ,MMM./MMMM.MMM, M                                           ')
    call yaml_scalar('  MMMMMM MMMMMMMMMMI  FFFFFF LL        AA    M        M EEEEEE')
    call yaml_scalar('  MMMMMM   MMMMMMMMM  F      LL       A  A   MM      MM E     ')
    call yaml_scalar('  `MMMM     MMMMMMM`  FFFFFF LL      A    A  M M    M M EEEEEE')
    call yaml_scalar('   /MMMMM   MMMMMM`   F      LL      AAAAAA  M  M  M  M E     ')
    call yaml_scalar('    MMMMMM  MMMMM`    F      LLLLLL A      A M   MM   M EEEEEE')     
    call yaml_scalar('________________________________________ www.flame-code.org   "')
    call yaml_mapping_close()
    call yaml_map('Reference Paper','To Be Added Later.')
    !call yaml_map('Version Number',package_version)
    call yaml_map('Timestamp of this run',yaml_date_and_time_toa())
end subroutine flm_print_logo
!*****************************************************************************************
!subroutine yaml_invoice_example()
!  use yaml_output
!  use yaml_strings, only: yaml_date_toa
!  implicit none
!
!  !call yaml_set_stream(tabbing=0)
!  call yaml_comment('Yaml Invoice Example',hfill='-')
!  call yaml_map('invoice',34843)
!  call yaml_map('date',trim(yaml_date_toa()))
!  call yaml_mapping_open('bill-to',label='id001')
!   call yaml_map('given','Chris')
!   call yaml_mapping_open('address')
!      call yaml_mapping_open('lines')
!      call yaml_scalar('458 Walkman Dr.')
!      call yaml_scalar('Suite #292')
!      call yaml_mapping_close()
!   call yaml_mapping_close()
!  call yaml_mapping_close()
!  call yaml_map('ship_to','*id001')
!  
!  !next step: sequence elements
!  call yaml_sequence_open('product')
!  !call yaml_sequence_open()
!    call yaml_sequence(advance='no')
!!    call yaml_mapping_open()
!      call yaml_map('sku','BL394D')
!      call yaml_map('quantity',4)
!      call yaml_map('description','Basketball')
!      call yaml_map('price',450.,fmt='(f6.2)')
!!    call yaml_mapping_close()
!    !call yaml_newline() !new line in a flow 
!     call yaml_sequence(advance='no')
!!     call yaml_mapping_open()
!     call yaml_map('sku','BL4438H')
!     call yaml_map('quantity',1)
!     call yaml_map('description','Super Hoop')
!     call yaml_map('price',2392.,fmt='(f8.2)')
!!     call yaml_mapping_close()
!    call yaml_sequence_close()
!    !final part
!    call yaml_map('tax',251.42,fmt='(f6.2)')
!    call yaml_map('total',4443.52d0,fmt='(f6.2)') !wrong format on purpose
!    call yaml_map('comments','Late afternoon is best. Backup contact is Nancy Billsmer @ 338-4338.')
!
!      !call yaml_mapping_close()
!
!end subroutine yaml_invoice_example
