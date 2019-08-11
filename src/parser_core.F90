!*****************************************************************************************
module mod_parser_ini
    implicit none
    private
    public:: read_file_input, get_header_location, split_line, get_one_param
contains
!*****************************************************************************************
subroutine read_file_input(file_ini)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    !local variables
    integer:: ios_err, iline
    integer:: istart1, istart2, i
    character(256):: str_tmp
    open(unit=1,file='input.ini',status='old',iostat=ios_err)
    if(ios_err/=0) then
        write(*,'(a)') 'ERROR: failure openning input.ini'
        stop
    endif
    do iline=1,file_ini%nline_max
        read(1,'(a)',iostat=ios_err) str_tmp
        if(ios_err>0) then
            write(*,'(a)') 'WARNING: input.ini seems to be at least partially binary'
        elseif(ios_err<0) then
            exit
        endif
        do i=1,256
            if(str_tmp(i:i)==char(9)) then
                !write(*,*) 'YES ',i
                str_tmp(i:i)=' '
            endif
        enddo
        if(len_trim(str_tmp)>254) then
            write(*,*) 'ERROR: line ',iline,' in input.ini is too long'
            stop
        endif
        str_tmp=adjustl(trim(str_tmp))
        !call convertlower(str_tmp)
        !if(tr_tmp(1:1)=='#') then
        !    comment_line(iline)=.true.
        !else
        !    comment_line(iline)=.false.
        !endif
        istart1=scan(str_tmp,'=')
        istart2=scan(str_tmp,'=',.true.)
        if(istart2>0 .and. istart2>istart1) then
            write(*,'(a47,i8)') 'ERROR: input.ini: too many equal sign in line #',iline
            stop
        endif
        if(istart1>0) then
            if((istart1>1 .and. istart1<254)) then
                if(str_tmp(istart1-1:istart1-1)/=' ') then
                    str_tmp=str_tmp(1:istart1-1)//' '//str_tmp(istart1:254)
                    istart1=istart1+1
                endif
                if(str_tmp(istart1+1:istart1+1)/=' ') then
                    str_tmp=str_tmp(1:istart1)//' '//str_tmp(istart1+1:255)
                endif
            else
                write(*,'(a59,i8)') 'ERROR: input.ini: improper position of equal sign in line #',iline
            endif
        endif
        file_ini%file_lines(iline)=str_tmp
    enddo
    file_ini%nline=iline-1
    close(1)
end subroutine read_file_input
!*****************************************************************************************
subroutine get_header_location(file_ini,str_header)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    character(*), intent(in):: str_header
    !local variables
    integer:: iline
    !determining the line number of the header
    file_ini%iline_header=0
    do iline=1,file_ini%nline
        file_ini%str_tmp=trim(file_ini%file_lines(iline))
        if(trim(file_ini%str_tmp)==trim(str_header)) then
            file_ini%iline_header=iline
            !write(*,*) file_ini%iline_header
            exit
        endif
    enddo
    if(file_ini%iline_header==0) return
    !determining the line number of the next header
    file_ini%iline_next_header=file_ini%nline+1
    do iline=file_ini%iline_header+1,file_ini%nline
        file_ini%str_tmp=adjustl(file_ini%file_lines(iline))
        if(file_ini%str_tmp(1:1)=='[') then
            file_ini%iline_next_header=iline
            exit
        endif
    enddo
end subroutine get_header_location
!*****************************************************************************************
subroutine split_line(file_ini)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    !local variables
    character(1):: ch
    character(256):: str_line_tmp
    integer:: ios_err, i_eqsign
    read(file_ini%file_lines(file_ini%iline),*,iostat=ios_err) file_ini%str_par,ch
    if(ios_err/=0 .or. ch/='=') then
        file_ini%str_par=''
        file_ini%str_val=''
    else
        i_eqsign=scan(file_ini%file_lines(file_ini%iline),'=')
        if(i_eqsign>1 .and. i_eqsign<254) then
            str_line_tmp=file_ini%file_lines(file_ini%iline)
            str_line_tmp(1:i_eqsign)=' '
            read(str_line_tmp,'(a)',iostat=ios_err) file_ini%str_val
            if(ios_err/=0) file_ini%str_val=''
        else
            file_ini%str_val=''
        endif
    endif
    !write(*,'(3(a,x))') trim(file_ini%str_par),trim(ch),trim(file_ini%str_val)
end subroutine split_line
!*****************************************************************************************
subroutine get_one_param(file_ini,var_name,int_var,real_var,char_var,char_line_var,log_var)
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    character(*), intent(in):: var_name
    integer, optional, intent(out):: int_var
    real(8), optional, intent(out):: real_var
    character(*), optional, intent(out):: char_var
    character(*), optional, intent(out):: char_line_var
    logical, optional, intent(out):: log_var
    !local variables
    integer:: ios_err
    logical:: something_present
    character(256):: str_m
    if(file_ini%stat_line_is_read(file_ini%iline)) return
    something_present=.false.
    if(trim(file_ini%str_par)==trim(var_name)) then
        if(.not. something_present .and. present(int_var)) then
            read(file_ini%str_val,*,iostat=ios_err) int_var
            something_present=.true.
        endif
        if(.not. something_present .and. present(real_var)) then
            read(file_ini%str_val,*,iostat=ios_err) real_var
            something_present=.true.
        endif
        if(.not. something_present .and. present(char_var)) then
            read(file_ini%str_val,*,iostat=ios_err) char_var
            something_present=.true.
        endif
        if(.not. something_present .and. present(char_line_var)) then
            read(file_ini%str_val,'(a)',iostat=ios_err) char_line_var
            something_present=.true.
        endif
        if(.not. something_present .and. present(log_var)) then
            read(file_ini%str_val,*,iostat=ios_err) log_var
            something_present=.true.
        endif
        if(ios_err/=0) then
            str_m='ERROR: input.ini: improper content for '//trim(var_name)//' at line #'
            write(*,'(a100,i8)') trim(str_m),file_ini%iline
            stop
        endif
        file_ini%stat_line_is_read(file_ini%iline)=.true.
    endif
end subroutine get_one_param
!*****************************************************************************************
end module mod_parser_ini
!*****************************************************************************************
