!> @file
!! Test the error handling part of flib
!! @example errs.f90
!! This is an example how to use the error handling.
!! @author
!!    Copyright (C) 2013-2014 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!> Routine testing the error handling part of flib
subroutine test_error_handling()
  use f_precisions, only: f_loc
  use yaml_output
!!! [Use]
  use dictionaries
  use yaml_strings, only: f_char_ptr
!!! [Use]
  implicit none
  !local variables
  integer :: ival,ierr,ERR_TOTO,ERR_TITI,ERR_GRAVE
  character(len=128) :: msg
  external :: abort_toto,abort_titi,abort1,abort2
  type(dictionary), pointer :: dict

  call yaml_comment('Error Handling Module Test',hfill='~')

  call f_err_severe_override(abort2)

!!! [Error Define]
  call f_err_define(err_name='ERR_TOTO',&
       err_msg='This is the error message for the error of kind 1 and it is written extensively'//&
       ' on purpose to see whether yaml module prints it',&
       err_action='For this error, contact the routine developer at mail at univ dot gov',&
       err_id=ERR_TOTO,callback=abort_toto)

  call f_err_define(err_name='ERR_TITI',err_msg='test2',err_id=ERR_TITI,&
       callback=abort_titi,callback_data=f_loc(ival))

  call f_err_define(err_name='ERR_GRAVE',err_msg='test2',err_id=ERR_GRAVE,&
       callback=f_err_severe)
!!! [Error Define]
  call yaml_map("Raising the TOTO error, errcode",ERR_TOTO)

  if (f_err_raise(.true.,'Extra message added',err_id=ERR_TOTO)) continue ! return

  call yaml_map('Print the error ID',f_get_last_error())
    call yaml_map("Raising the TOTO error, by name, without condition",'ERR_TOTO')
  if (f_err_raise(err_msg='Extra message added again',err_name='ERR_TOTO')) continue ! return

  call yaml_map("Callback done, errcode",ERR_TOTO)

!  call f_err_severe_restore()
  if (f_err_raise(.true.,'Generic error raised, some message here')) continue ! return

  call f_err_clean()

  call f_dump_possible_errors(f_char_ptr('This is the list of the errors'))

  call f_err_set_callback(abort2)

  call yaml_map("Callback done",f_err_raise(.true.,'Now TITI error has been raised',err_id=ERR_TITI))
  call yaml_map("Error check value",f_err_check())
  call yaml_map("Error check code",f_err_check(err_id=ERR_TOTO))
  call yaml_map("Error check code2",f_err_check(err_id=ERR_TITI))
  call yaml_map("Error check code, name",f_err_check(err_name='ERR_TOTO'))
  call yaml_map("Error check code, name",f_err_check(err_name='ERR_TITI'))

  call f_err_clean()

 !Test the nested try
  call yaml_comment("Test open try")
  call yaml_map('Error check value before try',f_err_check())
  call f_err_open_try()
     call f_err_throw('one',err_name='ERR_TOTO')
     call yaml_map("Number of errors(1)",f_get_no_of_errors())
     call f_err_open_try()
        call yaml_map("Number of errors(2)",f_get_no_of_errors())
        call f_err_throw('two',err_name='ERR_TOTO')
        if (f_err_check()) then
           ierr=f_get_last_error(msg)
           call yaml_map("ID",ierr)
           call yaml_map("MSG",msg)
           dict=>f_get_error_dict()
           call yaml_dict_dump(dict)
        end if
        call f_err_open_try()
           call f_err_throw('three',err_name='ERR_TOTO')
           call yaml_map("Number of errors(3)",f_get_no_of_errors())
        call f_err_close_try()
     call f_err_close_try()
  call f_err_close_try()
  call yaml_map('Error check value after try',f_err_check())

  call f_err_unset_callback()
  call f_err_severe_restore()

end subroutine test_error_handling


subroutine abort1()
  use yaml_output
  use exception_callbacks
  implicit none
  call f_dump_last_error()
  call yaml_comment('Ouille',hfill='!')
end subroutine abort1


subroutine abort2()
  use yaml_output
  use exception_callbacks
  implicit none
  call f_dump_last_error()
  call yaml_comment('Aie',hfill='!')
end subroutine abort2


subroutine abort_toto()
  use yaml_output
  use exception_callbacks
  implicit none
  call f_dump_last_error()
  call yaml_comment('TOTO',hfill='!')
end subroutine abort_toto


subroutine abort_titi()
  use yaml_output
  use exception_callbacks
  implicit none
  call f_dump_last_error()
  call yaml_comment('TITI',hfill='!')
end subroutine abort_titi
