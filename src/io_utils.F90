!*****************************************************************************************
subroutine read_list_files_yaml(fname,nfiles_max,fn_list,nfiles)
    use dictionaries
    use yaml_parse
    use dynamic_memory
    implicit none
    character(len=*), intent(in):: fname
    integer, intent(in):: nfiles_max
    character(len=256), intent(out):: fn_list(nfiles_max)
    integer, intent(out):: nfiles
    !local variales
    type(dictionary), pointer :: dict=>null()
    type(dictionary), pointer :: subdict=>null()
    type(dictionary), pointer :: subsubdict=>null()
    character, dimension(:), allocatable :: fbuf
    integer(kind = 8) :: cbuf, cbuf_len
    integer:: ifile
    character(len=256):: path
    call getFileContent(cbuf,cbuf_len,fname,len_trim(fname))
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')
    call copyCBuffer(fbuf,cbuf,cbuf_len)
    call freeCBuffer(cbuf)
    !then parse the user's input file
    call yaml_parse_from_char_array(dict,fbuf)
    call f_free_str(1,fbuf)
    !call dict_copy(subdictdict_user,dict//0)
    subdict=>dict//0
    if(.not. has_key(subdict,"files")) then
        write(*,'(2a)') 'ERROR: no files key in yaml file ',trim(fname)
        stop
    endif
    !subsubdict=>subdict//"files"
    !nfiles=dict_len(subsubdict)
    !write(*,*) 'LENGTH ',nfiles
    !do ifile=0,nfiles-1
    !    fn_list(ifile+1)=subsubdict//ifile
    !enddo
    subsubdict=>dict_iter(subdict//"files")
    ifile=0
    do while(associated(subsubdict))
        ifile=ifile+1
        fn_list(ifile)=dict_value(subsubdict)
        !write(*,*) 'LENGTH ',trim(fn_list(ifile))
        subsubdict=>dict_next(subsubdict)
    end do
    nfiles=ifile
    if(has_key(subdict,"path")) then
        path=dict_value(subdict//"path")
        !write(*,*) trim(path)
        do ifile=1,nfiles
            fn_list(ifile)=trim(path)//trim(fn_list(ifile))
        enddo
    endif
    nullify(subsubdict)
    nullify(subdict)
    call dict_free(dict)
    nullify(dict)
end subroutine read_list_files_yaml
!*****************************************************************************************
