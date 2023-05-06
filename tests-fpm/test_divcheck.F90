!*****************************************************************************************
subroutine test_divcheck()
    use iso_fortran_env, only: error_unit, output_unit
    use mod_parini, only: typ_parini
    use mod_colors, only: green_passed, red_failed
    use mod_flm_futile
    implicit none
    !local variables
    real(8):: errmax
    type(typ_parini):: parini
    integer:: nd, ios, ii, jj, id, i
    character(len=256):: tts1, tts2, str
    real(8), allocatable:: d(:), dref(:)
    integer, allocatable:: iconf(:), jconf(:)
    parini%iverbose=1
    parini%symfunc_type_ann='behler'
    parini%types_main='C N'
    call set_atomc_types_info(parini)
    parini%mpi_env%nproc=1
    parini%mpi_env%iproc=0
    open(unit=761,file='list_posinp_check.yaml',status='replace',iostat=ios)
    write(761,'(a)') 'files:'
    write(761,'(a)') ' - tests-fpm/data/poslow.yaml'
    close(761)
    call ann_check_symmetry_function(parini,'tests-fpm/data')
    nd=3
    allocate(d(nd),dref(nd),iconf(nd),jconf(nd))
    open(unit=761,file='tests-fpm/data/distall',status='old',iostat=ios)
    do id=1,nd
        read(761,*) tts1,iconf(id),tts2,jconf(id),dref(id)
        !write(*,*) trim(tts1),iconf(id),trim(tts2),jconf(id),dref(id)
    enddo
    close(761)
    open(unit=761,file='distall',status='old',iostat=ios)
    do id=1,nd
        read(761,'(a)') str
        do i=1,len_trim(str)
            if(str(i:i)=='/') str(i:i)="-"
        enddo
        !write(*,'(a)') trim(str)
        read(str,*) tts1,ii,tts2,jj,d(id)
        !write(*,*) trim(tts1),ii,trim(tts2),jj,d(id)
        if(ii/=iconf(id) .or. jj/=jconf(id)) then
            write(error_unit,'(2a)') red_failed,' in test_divcheck: '
            write(error_unit,'(a,2i6)') 'iconf,jconf should be: ',iconf(id),jconf(id)
            write(error_unit,'(a,2i6)') '       current values: ',ii,jj
            call exit(1)
        endif
    enddo
    close(761)
    errmax=maxval(abs(d-dref))
    if(errmax<1.d-13) then
        write(output_unit,'(2a)') green_passed,' in test_divcheck: errmax'
    else
        write(error_unit,'(2a,es14.5)') red_failed,' in test_divcheck: errmax=  ',errmax
        call exit(1)
    end if
    call execute_command_line("rm -f list_posinp_check.yaml")
    call execute_command_line("rm -f distall")
    call execute_command_line("rm -f incompatible")
    deallocate(d,dref,iconf,jconf)
end subroutine test_divcheck
!*****************************************************************************************
