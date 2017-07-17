!*****************************************************************************************
subroutine alborz_init(parini,file_ini)
    use mod_interface
    use mod_processors, only: iproc, mpi_comm_abz, imaster
    use mod_task, only: typ_file_ini, time_start
    use mod_parini, only: typ_parini
    use ifport
    use time_profiling
    use mod_timing , only: TCAT_ALBORZ_INIT_FINAL
    use dynamic_memory
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: istat
    call f_lib_initialize()
    call f_routine(id='alborz_init')
    !-----------------------------------------------------------------
    call alborz_initialize_timing_categories()
    call f_timing(TCAT_ALBORZ_INIT_FINAL,'ON')
    !-----------------------------------------------------------------
    istat=getcwd(parini%cwd)
    if(istat/=0) stop 'ERROR: could not get CWD'
    call cpu_time(time_start)
    !parsing all blocks in input.ini
    inquire(file="flame_in.yaml",exist=parini%exists_yaml_file)
    if(parini%exists_yaml_file) then
        call yaml_get_parameters(parini)
    else
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
        call get_saddle_1s_opt_parameters(file_ini,parini)
        call get_saddle_1s_parameters(file_ini,parini)
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
    !-----------------------------------------------------------------
    call f_timing(TCAT_ALBORZ_INIT_FINAL,'OF')
    call f_release_routine()
end subroutine alborz_init
!*****************************************************************************************
subroutine alborz_initialize_timing_categories
    use yaml_output
    use time_profiling
    use mod_timing, only: dict_timing_info, TCAT_PSOLVER, TCAT_SYMFUNC_COMPUT
    use mod_timing, only: TCAT_ALBORZ_INIT_FINAL
    use wrapper_mpi, only: mpi_initialize_timing_categories, wmpi_init_thread, mpifinalize
    use wrapper_mpi, only: mpisize, mpirank
    use mod_processors, only: iproc
    implicit none
    integer:: iverbose=3
    character(len=*), parameter :: pscpt1='Alborz init-fin'
    character(len=*), parameter :: pscpt2='potential'
    call f_timing_category_group(pscpt1,'Alborz computations total time')
    call f_timing_category_group(pscpt2,'potential computations total time')
    !-----------------------------------------------------------------
    call f_timing_category('Alborz initialize-finalize',pscpt1,'computation of init-fin',TCAT_ALBORZ_INIT_FINAL)
    !-----------------------------------------------------------------
    call f_timing_category('psolver',pscpt2,'computation of Psolver',TCAT_PSOLVER)
    call f_timing_category('symfunc compute',pscpt2,'computation of symmetry functions',TCAT_SYMFUNC_COMPUT)
    !-----------------------------------------------------------------
    call f_timing_reset(filename='time.yaml',master=iproc==0,verbose_mode=iverbose>2)
end subroutine alborz_initialize_timing_categories
!*****************************************************************************************
subroutine alborz_final(parini,file_ini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_task, only: typ_file_ini, time_start, time_end
    use mod_processors, only: iproc
    use yaml_output
    use time_profiling
    use mod_timing , only: dict_timing_info, TCAT_ALBORZ_INIT_FINAL
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_file_ini), intent(inout):: file_ini
    !local variables
    call f_routine(id='alborz_final')
    call f_timing(TCAT_ALBORZ_INIT_FINAL,'ON')
    if(.not. parini%exists_yaml_file) then
        deallocate(file_ini%file_lines,file_ini%stat_line_is_read) !,comment_line)
    endif
    call cpu_time(time_end)
    write(*,'(a,1x,i4,e15.3)') 'CPU time: iproc,time(hrs)',iproc,(time_end-time_start)/3600.d0
    write(*,'(a,1x,i4,e15.3)') 'CPU time: iproc,time(min)',iproc,(time_end-time_start)/60.d0
    write(*,'(a,1x,i4,e15.3)') 'CPU time: iproc,time(sec)',iproc,(time_end-time_start)
    if(trim(parini%task)/='minhocao') then
        call finalizeprocessors
    endif
    !-----------------------------------------------------------------
    call f_timing(TCAT_ALBORZ_INIT_FINAL,'OF')
    call f_timing_stop(dict_info=dict_timing_info)
    !-----------------------------------------------------------------
    call f_release_routine()
    call f_lib_finalize()
end subroutine alborz_final
!*****************************************************************************************
subroutine init_random_seed(parini)
    use mod_interface
    use mod_processors, only: iproc, mpi_comm_abz, imaster
    use mod_parini, only: typ_parini
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    integer:: istat, nseed, ierr
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
    if(nseed==1) then
        write(*,'(a,i4,i12)') 'iproc,seed ',iproc,iseed(1)
    else if(nseed==2) then
        write(*,'(a,i4,2i12)') 'iproc,seed ',iproc,iseed(1),iseed(2)
    else if(nseed==3) then
        write(*,'(a,i4,3i12)') 'iproc,seed ',iproc,iseed(1),iseed(2),iseed(3)
    else
        stop 'ERROR: seed size greater than 3 not considered in Alborz.'
    endif
    deallocate(iseed)
end subroutine init_random_seed
!*****************************************************************************************
subroutine set_atomc_types_info(parini)
    use mod_interface
    use mod_parini, only: typ_parini
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
    do i=1,ntypat
        ind_tmp(1:1)=minloc(iatomnum)
        ind=ind_tmp(1)
        do j=1,i-1
            if(stypat(ind)==parini%stypat(j)) then
                write(*,'(a,2a5)') 'ERROR: repeated types in types of block [main]:', &
                    stypat(ind),parini%stypat(j)
                stop
            endif
        enddo
        parini%ltypat(i)=i
        parini%iatomnum(i)=iatomnum(ind)
        parini%stypat(i)=stypat(ind)
        iatomnum(ind)=1000
    enddo
    do i=1,ntypat
        if(parini%iatomnum(i)==-1) stop 'ERROR: atomic number not set properly.'
    enddo
    parini%ntypat=ntypat
    write(*,'(a,i2)') 'ntypat ',parini%ntypat
    do i=1,ntypat
        write(*,'(a5,i2,a5)') 'type ',parini%ltypat(i),trim(parini%stypat(i))
    enddo
end subroutine set_atomc_types_info
!*****************************************************************************************
