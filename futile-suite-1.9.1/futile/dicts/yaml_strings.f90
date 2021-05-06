!> @file
!! Define the modules (yaml_strings and yaml_output) and the methods to write yaml output
!! yaml: Yet Another Markeup Language (ML for Human)
!! @author
!!    Copyright (C) 2011-2015 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Module defining all yaml strings for output.
!! This module must be only used by the module yaml_output
module yaml_strings
  use f_precisions
  implicit none

  private

  integer, parameter :: max_value_length=95 !< Not a parameter in order to be used by C bindings but constant

  character(len=*), parameter :: yaml_int_fmt  = '(i0)'                      !< Default format for integer
  character(len=*), parameter :: yaml_real_fmt = '(1pe18.9)'                 !< Default format for single
  character(len=*), parameter :: yaml_dble_fmt = '(1pg26.16e3)'!'(1pe25.17)' !< Default format for double
  character(len=*), parameter :: yaml_char_fmt = '(a)'                       !< Default format for strings

  !>escape sequences for pretty printing
  character(len=*), parameter :: escape_normal=char(27)//"[m"
  character(len=*), parameter :: escape_bold=char(27)//"[0;1m"
  character(len=*), parameter :: escape_blink=char(27)//"[0;5m"

  character(len=*), parameter, public :: YAML_INFINITY='.inf'

  !> structure containing the string and its length
  !! for the moment implement it basically, we might then
  !! identifty a strategy to allocate the string according to the needs
  type, public :: f_string
     character(len=4*max_value_length) :: msg
  end type f_string

  interface yaml_toa
     module procedure yaml_itoa,yaml_litoa,yaml_ftoa,yaml_dtoa,yaml_ltoa,yaml_ctoa
     module procedure yaml_dvtoa,yaml_ivtoa,yaml_cvtoa,yaml_ztoa,yaml_zvtoa,yaml_lvtoa,yaml_rvtoa
     module procedure yaml_livtoa
  end interface

  interface cnv_fmt  !< Give the default format corresponding to the nature of the data
     module procedure fmt_i,fmt_r,fmt_d,fmt_a,fmt_li
  end interface

  interface f_strcpy
     module procedure f_strcpy,f_strcpy_str
  end interface

  interface operator(.eqv.)
     module procedure string_equivalence
  end interface operator(.eqv.)

  interface operator(.neqv.)
     module procedure string_inequivalence
  end interface operator(.neqv.)

  interface operator(//)
     module procedure string_and_integer,string_and_double,string_and_simple
     module procedure string_and_long,string_and_msg,msg_and_string
     module procedure integer_and_string,integer_and_msg,msg_and_msg
  end interface

  interface operator(+)
     module procedure combine_strings,attach_ci,attach_cli,attach_cd,combine_msg
     module procedure attach_c_msg,attach_msg_c,attach_lic
  end interface

  interface assignment(=)
     module procedure msg_to_string,string_to_msg
  end interface assignment(=)

  interface operator(**)
     module procedure yaml_itoa_fmt,yaml_litoa_fmt,yaml_dtoa_fmt,yaml_ctoa_fmt
  end interface operator(**)

  !Public routines
  public :: yaml_toa, buffer_string, align_message, shiftstr,yaml_date_toa
  public :: yaml_date_and_time_toa,yaml_time_toa,is_atoi,is_atof,is_atol,is_atoli
  public :: read_fraction_string,f_strcpy
  public:: yaml_bold,yaml_blink,rstrip,f_char_ptr,convert_f_char_ptr
  public :: operator(.eqv.),operator(.neqv.),operator(+),operator(//),operator(**),assignment(=)

contains

  pure function fmt_li(data)
    implicit none
    integer(f_long), intent(in) :: data
    character(len=len(yaml_int_fmt)) :: fmt_li
    !local variables
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    fmt_li=yaml_int_fmt
  end function fmt_li

  pure function fmt_i(data)
    implicit none
    integer(f_integer), intent(in) :: data
    character(len=len(yaml_int_fmt)) :: fmt_i
    !local variables
    integer :: kindt
    kindt=kind(data) !to remove compilation warning
    fmt_i=yaml_int_fmt
  end function fmt_i

  pure function fmt_r(data)
    implicit none
    real(f_simple), intent(in) :: data
    character(len=len(yaml_real_fmt)) :: fmt_r
    !local variables
    integer :: kindt
    kindt=kind(data) !to remove compilation warning

    fmt_r=yaml_real_fmt
  end function fmt_r

  pure function fmt_d(data)
    implicit none
    real(f_double), intent(in) :: data
    character(len=len(yaml_dble_fmt)) :: fmt_d
    !local variables
    integer :: kindt
    kindt=kind(data) !to remove compilation warning

    fmt_d=yaml_dble_fmt
  end function fmt_d

  pure function fmt_a(data)
    implicit none
    character(len=*), intent(in) :: data
    character(len=len(yaml_char_fmt)) :: fmt_a
    !local variables
    integer :: kindt
    kindt=kind(data) !to remove compilation warning

    fmt_a=yaml_char_fmt
  end function fmt_a

  !> Write the strings as if they were written by write
  pure subroutine f_strcpy(dest,src)
    implicit none
    character(len=*), intent(out) :: dest
    character(len=*), intent(in) :: src
    !local variables
    integer :: i,n
    !dest=repeat(' ',len(dest))
    n=min(len(src),len(dest))
    do i=1,n
       dest(i:i)=src(i:i)
    end do
    do i=n+1,len(dest)
       dest(i:i)=' '
    end do

  end subroutine f_strcpy

  pure subroutine f_strcpy_str(dest,src)
    implicit none
    character(len=*), intent(out) :: dest
    type(f_string), intent(in) :: src

    call f_strcpy(dest=dest,src=src%msg)

  end subroutine f_strcpy_str

  pure function yaml_escape(str,escape_sequence) result(bstr)
    implicit none
    character(len=*), intent(in) :: str,escape_sequence
    character(len=max_value_length) :: bstr
    !local variables
    integer :: ipos
    !open the bold field
    ipos=0
    call buffer_string(bstr,len(bstr),escape_sequence,ipos)
    !copy the input string
    call buffer_string(bstr,len(bstr),trim(str),ipos)
    !then copy the normal escape sequence
    call buffer_string(bstr,len(bstr),escape_normal,ipos,&
         back=ipos+len(escape_normal) > len(bstr))
    bstr(ipos+1:len(bstr))=' '
  end function yaml_escape

  !> boldify a string (to be used only) if the terminal is of tty type
  pure function yaml_bold(str) result(bstr)
    implicit none
    character(len=*), intent(in) :: str
    character(len=max_value_length) :: bstr
    bstr=yaml_escape(str,escape_bold)
  end function yaml_bold

  !> boldify a string (to be used only) if the terminal is of tty type
  pure function yaml_blink(str) result(bstr)
    implicit none
    character(len=*), intent(in) :: str
    character(len=max_value_length) :: bstr
    bstr=yaml_escape(str,escape_blink)
  end function yaml_blink

  !> Add a buffer to a string and increase its length
  pure subroutine buffer_string(string,string_lgt,buffer,string_pos,back,istat)
    implicit none
    integer, intent(in) :: string_lgt                   !< Length of the string towrite
    character(len=string_lgt), intent(inout) :: string  !< String towrite
    integer, intent(inout) :: string_pos                !< Position to add buffer into string and for the next.
    character(len=*), intent(in) :: buffer              !< Buffer to add
    logical, optional, intent(in) :: back               !< Add string from the end
    integer, optional, intent(out) :: istat             !< Error status (if present otherwise truncate output if error)
    !local variables
    integer :: lgt_add

    if (present(istat)) istat=0 !no errors

    lgt_add=len(buffer)
    !do not copy strings which are too long if istat is present
    if (lgt_add+string_pos > string_lgt) then
       !try to eliminate trailing spaces
       lgt_add=len_trim(buffer)
    end if

    if (lgt_add+string_pos > string_lgt) then
       !ic corrections August 23rd 2013
       !If the string is too long, truncate the length to the length of String: string_lgt-string_pos
       !write in stderr that a problem is produced (not compatible with pure procedures
       !write(0,*)'WARNING: Length of Buffer is too long: ',lgt_add
       !write(0,*)'WARNING: Missing string: ',buffer(string_lgt-string_pos:lgt_add)
       if (present(istat)) then
          istat=-1
          return
       else if (lgt_add+string_pos > string_lgt) then
!          write(*,*)'#ERROR (buffer string): string too long'
!          write(*,*)'#Initial String: ',string(1:string_pos)
!          write(*,*)'#Buffer: ',trim(buffer)
!          write(*,*)'#String position: ',string_pos
!          write(*,*)'#Length of Buffer: ',lgt_add
!          write(*,*)'#String limit: ',string_lgt
          lgt_add=string_lgt-string_pos-1
!          write(*,*)'#Buffer shortened into: ',buffer(1:lgt_add)
          !stop
       end if
    end if

    if (lgt_add==0) return
    if (present(back)) then
       if (back) then
          call shiftstr(string,lgt_add)
          string(1:lgt_add)=buffer(1:lgt_add)
       else
          string(string_pos+1:string_pos+lgt_add)=buffer(1:lgt_add)
       end if
    else
       string(string_pos+1:string_pos+lgt_add)=buffer(1:lgt_add)
    end if

    string_pos=string_pos+lgt_add

  end subroutine buffer_string


  !> Add the spaces necessary to align the first occurrence of a given anchor
  !! into a tabular value. Can be done either by moving rigidly the message or
  !! by adding spaces between the anchor and the rest of the message
  subroutine align_message(rigid,maxlen,tabval,anchor,message)
    implicit none
    logical, intent(in) :: rigid
    integer, intent(in) :: maxlen
    integer, intent(in) :: tabval
    character(len=*), intent(in) :: anchor
    character(len=maxlen), intent(inout) :: message
    !local variables
    integer :: iscpos,ishift

    !cannot align, tabular too far
    if (tabval>maxlen) return

    iscpos=index(message,anchor)
    ishift=tabval-iscpos
    if (rigid) then
       call shiftstr(message,ishift)
    else
       message=message(1:iscpos-1)//repeat(' ',ishift)//anchor//&
            message(iscpos+1:maxlen-ishift)  ! shift right
    end if

  end subroutine align_message


  !> Convert integer to character
  pure function yaml_itoa(data,fmt) result(str)
    implicit none
    integer(f_integer), intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_itoa

  !> Convert longinteger to character
  pure function yaml_litoa(data,fmt) result(str)
    implicit none
    integer(f_long), intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_litoa

  !> Convert float to character
  pure function yaml_ftoa(data,fmt) result(str)
    implicit none
    real(f_simple), intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_ftoa

  !> Convert double to character
  pure function yaml_dtoa(data,fmt) result(str)
    implicit none
    real(f_double), intent(in) :: data
    include 'yaml_toa-inc.f90'
  end function yaml_dtoa

  !> character to character, only for genericity
  pure function yaml_ctoa(d,fmt)
    implicit none
    character(len=*), intent(in) :: d
    character(len=len(d)) :: yaml_ctoa
    character(len=*), optional, intent(in) :: fmt

    if (present(fmt)) then
       write(yaml_ctoa,fmt) trim(d)
    else
       call f_strcpy(src=d,dest=yaml_ctoa)
       !yaml_ctoa(1:max_value_length)=trim(d)
    end if
  end function yaml_ctoa


  !> Convert double complex to character
  !! use python notation for yaml complex
  pure function yaml_ztoa(z,fmt)
    implicit none
    double complex, intent(in) :: z
    character(len=max_value_length) :: yaml_ztoa
    character(len=*), optional, intent(in) :: fmt
    !local variables
    integer :: ipos,rpos
    character(len=max_value_length) :: zr,zi
    double complex :: ztmp
    double precision, dimension(2) :: zeta

    yaml_ztoa=repeat(' ',max_value_length)
    ztmp=z
    zeta=transfer(ztmp,zeta)
    zr=yaml_ztoa
    zi=yaml_ztoa

    if (present(fmt)) then
       write(zr,fmt) zeta(1)
       write(zi,fmt) zeta(2)
    else
       write(zr,yaml_dble_fmt) zeta(1)
       write(zi,yaml_dble_fmt) zeta(2)
    end if

    zr=yaml_adjust(zr,clean=.not. present(fmt))
    zi=yaml_adjust(zi,clean=.not. present(fmt))
    rpos=len(trim(zr))
    ipos=min(len(trim(zi)),max_value_length-rpos-2)

    yaml_ztoa(1:rpos)=zr(1:rpos)
    if (zeta(2) >= 0.d0) then
       yaml_ztoa(rpos+1:rpos+2)='+'
       rpos=rpos+1
    end if
    yaml_ztoa(rpos+1:rpos+ipos)=zi(1:ipos)
    yaml_ztoa(rpos+ipos+1:rpos+ipos+2)='j'
  end function yaml_ztoa

  !> Convert logical to character
  pure function yaml_ltoa(l,fmt)
    implicit none
    logical, intent(in) :: l
    character(len=max_value_length) :: yaml_ltoa
    character(len=*), optional, intent(in) :: fmt

    yaml_ltoa=repeat(' ',max_value_length)

    if (present(fmt)) then
       write(yaml_ltoa,fmt) l
    else
       if (l) then
          write(yaml_ltoa,'(a3)')'Yes'
       else
          write(yaml_ltoa,'(a3)')'No'
       end if
    end if

    yaml_ltoa=yaml_adjust(yaml_ltoa)
  end function yaml_ltoa

  !> Convert vector of double to character
  pure function yaml_dvtoa(vec,fmt) result(vec_toa)
    implicit none
    real(kind=8), dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_dvtoa

  !> Convert vector of double complex to character
  pure function yaml_zvtoa(vec,fmt) result(vec_toa)
    implicit none
    double complex, dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_zvtoa

  !> Convert vector of integer to character
  !! @warning Truncate if too long
  pure function yaml_ivtoa(vec,fmt) result(vec_toa)
    implicit none
    integer(kind=4), dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_ivtoa

  !> Convert vector of characters to a chain of characters
  pure function yaml_cvtoa(vec,fmt) result(vec_toa)
    implicit none
    character(len=*), dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_cvtoa

  pure function yaml_lvtoa(vec,fmt) result(vec_toa)
    implicit none
    logical, dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_lvtoa

  pure function yaml_rvtoa(vec,fmt) result(vec_toa)
    implicit none
    real, dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_rvtoa

  pure function yaml_livtoa(vec,fmt) result(vec_toa)
    implicit none
    integer(kind=8), dimension(:), intent(in) :: vec
    include 'yaml_toa-arr-inc.f90'
  end function yaml_livtoa

  !> Yaml Spaced format for Date and Time
  function yaml_date_and_time_toa(values,zone)
    implicit none
    logical, optional, intent(in) :: zone
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_date_and_time_toa
    !local variables
    character(len=*), parameter :: &
         deffmt='i4.4,"-",i2.2,"-",i2.2," ",i2.2,":",i2.2,":",i2.2,".",i3.3'
    logical :: zon
    integer :: zonhrs,zonmin
    integer, dimension(8) :: vals
    character(len=4) :: sgn

    zon=.false.
    if (present(zone)) zon=zone

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    if (zon) then
       zonmin=abs(mod(vals(4),60))
       zonhrs=abs(vals(4)/60)
       if (vals(4) < 0) then
          sgn='" -"'
       else
          sgn='" +"'
       end if
       write(yaml_date_and_time_toa,'('//deffmt//','//sgn//',i2.2,":",i2.2)')&
            vals(1:3),vals(5:8),zonhrs,zonmin

    else
       write(yaml_date_and_time_toa,'('//deffmt//')')vals(1:3),vals(5:8)
    end if

    !There is no - sign so we skip this step (TD)
    !yaml_date_and_time_toa=yaml_adjust(yaml_date_and_time_toa)

  end function yaml_date_and_time_toa

  function yaml_date_toa(values)
    ! Yaml Spaced format for Date
    implicit none
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_date_toa
    !local variables
    integer, dimension(8) :: vals

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    write(yaml_date_toa,'(i4.4,"-",i2.2,"-",i2.2)')vals(1:3)

    yaml_date_toa=yaml_adjust(yaml_date_toa,clean=.false.)

  end function yaml_date_toa

  function yaml_time_toa(values)
    implicit none
    integer, dimension(8), optional, intent(in) :: values
    character(len=max_value_length) :: yaml_time_toa
    !local variables
    integer, dimension(8) :: vals

    if (present(values)) then
       vals=values
    else
       call date_and_time(values=vals)
    end if

    write(yaml_time_toa,'(i2.2,":",i2.2,":",i2.2,".",i3.3)')vals(5:8)

    yaml_time_toa=yaml_adjust(yaml_time_toa,clean=.false.)

  end function yaml_time_toa

  pure function yaml_adjust(str,clean)
    implicit none
    character(len=*), intent(in) :: str
    logical, intent(in), optional :: clean
    character(len=max_value_length) :: yaml_adjust
    !local variables
    logical :: clean0

    clean0=.true.
    if (present(clean)) clean0=clean

    yaml_adjust=adjustl(str)

    if (clean0) yaml_adjust=clean_zero(yaml_adjust)

    !put a space if there is no sign
    if (yaml_adjust(1:1)/='-') then
       call shiftstr(yaml_adjust,1)
    else
       call shiftstr(yaml_adjust,0)
    end if

  end function yaml_adjust

  pure function clean_zero(str)
    implicit none
    character(len=*), intent(in) :: str
    character(len=max_value_length) :: clean_zero
    !local variables
    integer :: idot,iexpo,i

    !first fill with all the values up to the dot if it exist
    idot=scan(str,'.')
    if (idot==0) then
       !no dot, nothing to clean
       clean_zero(1:max_value_length)=str
    else
       !first find the position of the end of the string
!       iend=len_trim(str)
       !then search for the position of the exponent or of the space if present
       iexpo=scan(str(idot+2:),'eE ')+idot+1
       !print *,'there',iexpo,'str',str(idot+2:)
       if (iexpo==idot+1) iexpo=len(str)+1
       i=iexpo
       find_last_zero: do while(i > idot+1) !first digit after comma always stays
          i=i-1
          if (str(i:i) /= '0') exit find_last_zero
       end do find_last_zero
       clean_zero(1:i)=str(1:i)
       !print *,'here',i,clean_zero(1:i),'iexpo',iexpo,str(iexpo:)
       !then copy the exponent
       if (str(iexpo:) /= 'E+00' .and. str(iexpo:) /= 'e+00' .and. str(iexpo:) /= 'E+000' .and. &
            str(iexpo:) /= 'e+000') then
          clean_zero(i+1:max_value_length)=str(iexpo:)
       else
          clean_zero(i+1:max_value_length)=' '
       end if
       !try to put at the old position a termination character
!       clean_zero(iend:iend)=char(0)
    end if
  end function clean_zero


  !>find if a string is an integer
  !! use the portable mode described in
  !! http://flibs.sourceforge.net/fortran_aspects.html#check_integers
  pure function is_atoi(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    integer :: ierr,ival
    character(len=20) :: form

    !fill the string describing the format to be used for reading
    !use the trimmed string and the yaml_toa function as i0 can add extra zeros in the specifications
    write(form,'(a20)')'(i'//adjustl(trim(yaml_toa(len_trim(str),fmt='(i17)')))//')'
    read(str,trim(form),iostat=ierr)ival
    yes=ierr==0
  end function is_atoi

  !>find if a string is a long integer
  !! use the portable mode described in
  !! http://flibs.sourceforge.net/fortran_aspects.html#check_integers
  !! note that this function also gives positive answer if the character fits with ddefault integer type
  !! therefore care should be taken in the usage (use only when long is needed)
  pure function is_atoli(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    integer :: ierr
    integer(kind=8) :: ival
    character(len=20) :: form

    !fill the string describing the format to be used for reading
    !use the trimmed string and the yaml_toa function as i0 can add extra zeros in the specifications
    write(form,'(a20)')'(i'//adjustl(trim(yaml_toa(len_trim(str),fmt='(i17)')))//')'
    read(str,trim(form),iostat=ierr)ival
    yes=ierr==0
  end function is_atoli


  !>check if str contains a floating point number.
  !!note that in principle this function gives positive answer also
  !!if the number in str is an integer. Therefore is_atoi should be used to check before
  pure function is_atof(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    integer :: ierr,is,ie
    double precision :: rval

    ie=len_trim(str)
    is=scan(trim(str),' ')+1
    yes=scan(str(is:ie),' ') ==0 !there is no other space in the string
    if (yes) then
       read(str(is:ie),*,iostat=ierr)rval
       yes=ierr==0 .and. str(max(ie,1):ie)/='/' !the slash terminator is not allowed
    end if

  end function is_atof

  !> check if str contains a logical in yaml specification (Yes=.true. and No=.false.)
  pure function is_atol(str) result(yes)
    implicit none
    character(len=*), intent(in) :: str
    logical :: yes
    !local variables
    integer :: is,ie
    !fill the string describing the format to be used for reading
    !use the trimmed string and the yaml_toa function as i0 can add extra zeros in the specifications
    ie=len(trim(str))
    is=max(scan(trim(str),' '),1)
    yes=scan(str(is:ie),' ') ==0 !there is no other space in the string
    if (yes) yes= (ie-is+1==3 .and. any(str(is:ie) == ['Yes', 'yes', 'YES'])) &
       &     .or. (ie-is+1==2 .and. any(str(is:ie) == ['No', 'no', 'NO'])) &
       &     .or. (ie-is+1==4 .and. any(str(is:ie) == ['True', 'true', 'TRUE'])) &
       &     .or. (ie-is+1==5 .and. any(str(is:ie) == ['False', 'false', 'FALSE']))
  end function is_atol

  !> Read a real or real/real, real:real
  !! Here the fraction is indicated by the ':' or '/'
  !! The problem is that / is a separator for Fortran
  pure subroutine read_fraction_string(string,var,ierror)
    implicit none
    !Arguments
    character(len=*), intent(in) :: string
    double precision, intent(out) :: var
    integer, intent(out) :: ierror
    !Local variables
    character(len=256) :: tmp
    integer :: num,den,pfr,psp

    !First look at the first blank after trim
    tmp(1:len(tmp))=trim(adjustl(string))
    psp = scan(tmp,' ')
    !see whether there is a fraction in the string
    if(psp==0) psp=len(tmp)
    pfr = scan(tmp(1:psp),':')
    if (pfr == 0) pfr = scan(tmp(1:psp),'/')
    !It is not a fraction
    if (pfr == 0) then
       read(tmp(1:psp),*,iostat=ierror) var
    else
       read(tmp(1:pfr-1),*,iostat=ierror) num
       read(tmp(pfr+1:psp),*,iostat=ierror) den
       if (ierror == 0) var=dble(num)/dble(den)
    end if
    !Value by defaut
    if (ierror /= 0) var = huge(1.d0)
  END SUBROUTINE read_fraction_string

  pure function string_inequivalence(a,b) result(notok)
    implicit none
    character(len=*), intent(in) :: a,b
    logical :: notok
    notok = .not. string_equivalence(a,b)
  end function string_inequivalence

  pure function string_equivalence(a,b) result(ok)
    implicit none
    character(len=*), intent(in) :: a,b
    logical :: ok
    !local variables
    integer :: ln
    !to be equivalent the two strings must have already the same length
    ln= len_trim(adjustl(a))
    ok= ln == len_trim(adjustl(b))
    if (.not. ok .or. ln ==0) return
    ok=case_insensitive_equiv(trim(adjustl(a)),trim(adjustl(b)))
  end function string_equivalence

  !> Compare two strings (case-insensitive). Blanks are relevant!
  pure function case_insensitive_equiv(stra,strb)
    implicit none
    character(len=*), intent(in) :: stra,strb
    logical :: case_insensitive_equiv
    !Local variables
    integer :: i,ica,icb,ila,ilb,ilength
    ila=len(stra)
    ilb=len(strb)
    ilength=min(ila,ilb)
    ica=ichar(stra(1:1))
    icb=ichar(strb(1:1))
    case_insensitive_equiv=(modulo(ica-icb,32) == 0) .and. (ila==ilb)
    do i=2,ilength
       ica=ichar(stra(i:i))
       icb=ichar(strb(i:i))
       case_insensitive_equiv=case_insensitive_equiv .and. &
            &   (modulo(ica-icb,32) == 0)
       if (.not. case_insensitive_equiv) exit
    end do

  end function case_insensitive_equiv

  !> modifies string to remove substring from the right
  pure subroutine rstrip(string,substring)
    implicit none
    character(len=*), intent(in) :: substring
    character(len=*), intent(inout) :: string
    !local variables
    integer :: ipos
    
    ipos=index(string,substring,back=.true.)
    if (ipos > 0) string=string(1:ipos-1)

  end subroutine rstrip

  !define the strings which combine them, without the need of using trim or adjustl specifications
  pure function combine_strings(a,b) result(c)
    implicit none
    character(len=*), intent(in) :: a
    character(len=*), intent(in) :: b
    character(len=len_trim(a)+len_trim(b)) :: c
    c=trim(a)//trim(adjustl(b))
  end function combine_strings

  !define the strings which combine them, without the need of using trim or adjustl specifications
  pure function combine_msg(a,b) result(c)
    implicit none
    type(f_string), intent(in) :: a
    type(f_string), intent(in) :: b
    character(len=len_trim(a%msg)+len_trim(b%msg)) :: c
    c=trim(a%msg)//trim(adjustl(b%msg))
  end function combine_msg

  !define the strings which combine them, without the need of using trim or adjustl specifications
  pure function attach_msg_c(a,b) result(c)
    implicit none
    type(f_string), intent(in) :: a
    character(len=*), intent(in) :: b
    character(len=len_trim(a%msg)+len_trim(b)) :: c
    c=trim(a%msg)//trim(adjustl(b))
  end function attach_msg_c

  !define the strings which combine them, without the need of using trim or adjustl specifications
  pure function attach_c_msg(a,b) result(c)
    implicit none
    character(len=*), intent(in) :: a
    type(f_string), intent(in) :: b
    character(len=len_trim(a)+len_trim(b%msg)) :: c
    c=trim(a)//trim(adjustl(b%msg))
  end function attach_c_msg

  pure function string_and_integer(a,num) result(c)
    implicit none
    integer(f_integer), intent(in) :: num
    character(len=*), intent(in) :: a
    type(f_string) :: c
    call f_strcpy(c%msg,a//trim(yaml_toa(num)))
  end function string_and_integer

  pure function integer_and_string(a,num) result(c)
    implicit none
    integer(f_integer), intent(in) :: a
    character(len=*), intent(in) :: num
    type(f_string) :: c
    call f_strcpy(c%msg,trim(adjustl(yaml_toa(a)))//trim(num))
  end function integer_and_string

  pure function integer_and_msg(a,num) result(c)
    implicit none
    integer(f_integer), intent(in) :: a
    type(f_string), intent(in) :: num
    type(f_string) :: c
    call f_strcpy(c%msg,trim(adjustl(yaml_toa(a)))//trim(num%msg))
  end function integer_and_msg

  pure function string_and_long(a,num) result(c)
    implicit none
    integer(f_long), intent(in) :: num
    character(len=*), intent(in) :: a
    type(f_string) :: c
    call f_strcpy(c%msg,trim(adjustl(a))//trim(yaml_toa(num)))
  end function string_and_long

  pure function string_and_double(a,num) result(c)
    implicit none
    real(f_double), intent(in) :: num
    character(len=*), intent(in) :: a
    type(f_string) :: c
    call f_strcpy(c%msg,trim(adjustl(a))//trim(yaml_toa(num)))
  end function string_and_double

  pure function string_and_simple(a,num) result(c)
    implicit none
    real(f_simple), intent(in) :: num
    character(len=*), intent(in) :: a
    type(f_string) :: c
    call f_strcpy(c%msg,trim(adjustl(a))//trim(yaml_toa(num)))
  end function string_and_simple


  pure function string_and_msg(a,num) result(c)
    implicit none
    type(f_string), intent(in) :: num
    character(len=*), intent(in) :: a
    character(len=len_trim(adjustl(a))+len_trim(num%msg)) :: c
    c=trim(adjustl(a))//trim(num%msg)
  end function string_and_msg

  pure function msg_and_msg(a,num) result(c)
    implicit none
    type(f_string), intent(in) :: num
    type(f_string), intent(in) :: a
    character(len=len_trim(a%msg)+len_trim(num%msg)) :: c
    c=trim(adjustl(a%msg))//trim(num%msg)
!!$    type(f_string) :: c
!!$    call f_strcpy(c%msg,a//trim(yaml_toa(num)))
  end function msg_and_msg

  pure function msg_and_string(a,num) result(c)
    implicit none
    type(f_string), intent(in) :: a
    character(len=*), intent(in) :: num
    character(len=len_trim(a%msg)+len_trim(num)) :: c
    c=trim(adjustl(a%msg))//trim(num)
!!$    type(f_string) :: c
!!$    call f_strcpy(c%msg,a//trim(yaml_toa(num)))
  end function msg_and_string

  pure subroutine msg_to_string(string,msg)
    implicit none
    character(len=*), intent(out) :: string
    type(f_string), intent(in) :: msg
    call f_strcpy(string,msg%msg)
  end subroutine msg_to_string

  pure subroutine string_to_msg(msg,string)
    implicit none
    character(len=*), intent(in) :: string
    type(f_string), intent(out) :: msg
    call f_strcpy(msg%msg,string//char(0))
  end subroutine string_to_msg

  !function which attach two strings each other
  pure function attach_ci(s,num) result(c)
    implicit none
    integer(f_integer), intent(in) :: num
    character(len=*), intent(in) :: s
    type(f_string) :: c
    call f_strcpy(c%msg,trim(s)//trim(adjustl(yaml_toa(num))))
  end function attach_ci

  pure function attach_cli(s,num) result(c)
    implicit none
    integer(f_long), intent(in) :: num
    character(len=*), intent(in) :: s
    type(f_string) :: c
    call f_strcpy(c%msg,trim(s)//trim(adjustl(yaml_toa(num))))
  end function attach_cli

  pure function attach_lic(num,s) result(c)
    implicit none
    integer(f_long), intent(in) :: num
    character(len=*), intent(in) :: s
    type(f_string) :: c
    call f_strcpy(c%msg,trim(adjustl(yaml_toa(num)))//trim(adjustl(s)))
  end function attach_lic

  pure function attach_cd(s,num) result(c)
    implicit none
    real(f_double), intent(in) :: num
    character(len=*), intent(in) :: s
    type(f_string) :: c
    call f_strcpy(c%msg,trim(s)//trim(adjustl(yaml_toa(num))))
  end function attach_cd

  pure function yaml_itoa_fmt(num,fmt) result(c)
    implicit none
    integer(f_integer), intent(in) :: num
    character(len=*), intent(in) :: fmt
!!$    type(f_string) :: c
!!$    call f_strcpy(c%msg,trim(adjustl(yaml_toa(num,fmt))))
    character(len=max_value_length) :: c
    c=yaml_toa(num,fmt)
  end function yaml_itoa_fmt

  pure function yaml_litoa_fmt(num,fmt) result(c)
    implicit none
    integer(f_long), intent(in) :: num
    character(len=*), intent(in) :: fmt
!!$    type(f_string) :: c
!!$    call f_strcpy(c%msg,trim(adjustl(yaml_toa(num,fmt))))
    character(len=max_value_length) :: c
    c=yaml_toa(num,fmt)
  end function yaml_litoa_fmt

  pure function yaml_dtoa_fmt(num,fmt) result(c)
    implicit none
    real(f_double), intent(in) :: num
    character(len=*), intent(in) :: fmt
!!$    type(f_string) :: c
!!$    call f_strcpy(c%msg,trim(adjustl(yaml_toa(num,fmt))))
    character(len=max_value_length) :: c
    c=yaml_toa(num,fmt)
  end function yaml_dtoa_fmt

  pure function yaml_ctoa_fmt(num,fmt) result(c)
    implicit none
    character(len=*), intent(in) :: num
    character(len=*), intent(in) :: fmt
!!$    type(f_string) :: c
!!$    call f_strcpy(c%msg,trim(adjustl(yaml_toa(num,fmt))))
    character(len=max_value_length) :: c
    c=yaml_toa(num,fmt)
  end function yaml_ctoa_fmt


  !> Shifts characters in in the string 'str' n positions (positive values
  !! denote a right shift and negative values denote a left shift). Characters
  !! that are shifted off the end are lost. Positions opened up by the shift
  !! are replaced by spaces.
  !! This routine has been downloaded from the website http://gbenthien.net/strings/index.html
  pure subroutine shiftstr(str,n)
    implicit none
    integer, intent(in) :: n
    character(len=*), intent(inout) :: str
    !local variables
    integer :: lenstr,nabs

    lenstr=len(str)
    nabs=iabs(n)
    if(nabs>=lenstr) then
       str=repeat(' ',lenstr)
       return
    end if
    if(n<0) str=str(nabs+1:)//repeat(' ',nabs)  ! shift left
    if(n>0) str=repeat(' ',nabs)//str(:lenstr-nabs)  ! shift right

  end subroutine shiftstr
  
  !> convert a fortran string into a stack-allocated array
  function f_char_ptr(str)
    implicit none
    character(len=*), intent(in) :: str
    character, dimension(len_trim(str)+1) :: f_char_ptr
    !local variables
    integer :: i
    do i=1,len_trim(str)
       f_char_ptr(i)=str(i:i)
    end do
    !i=len_trim(str)+1 !not needed by fortran norm
    f_char_ptr(i)=char(0)
  end function f_char_ptr

  subroutine convert_f_char_ptr(src,dest)
    implicit none
    character, dimension(*), intent(in) :: src
    character(len=*), intent(out) :: dest
    !local variables
    integer :: i
    dest(1:len(dest))=' '
    seek_and_copy: do i=1,len(dest)
       if (src(i)==char(0)) exit seek_and_copy
       dest(i:i)=src(i)
    end do seek_and_copy
  end subroutine convert_f_char_ptr

end module yaml_strings
