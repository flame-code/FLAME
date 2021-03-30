!> @file
!!  Module defining handling of errors
!!  In the spirit of this module each error message is handled separately.
!!  The user of this module is able to define the stopping routine, the new error and the error codes.
!! @author Luigi Genovese
!!    Copyright (C) 2012-2015 BigDFT group <br>
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS


!!$module error_handling
!!$  use exception_callbacks
!!$  use dictionaries, dict_msg_len=>max_field_length
!!$  implicit none
!!$
!!$  private
!!$
!!$  !some parameters
!!$  character(len=*), parameter :: ERRID='Id'
!!$  character(len=*), parameter :: ERRMSG='Message'
!!$  character(len=*), parameter :: ERRACT='Action'
!!$  character(len=*), parameter :: ERRCLBK='Callback Procedure Address'
!!$  character(len=*), parameter :: errclbkadd='Callback Procedure Data Address'
!!$
!!$  character(len=*), parameter :: ERRUNSPEC='UNSPECIFIED'
!!$  character(len=*), parameter :: ERRUNDEF='UNKNOWN'
!!$
!!$  integer :: ERR_GENERIC,ERR_SUCCESS,ERR_NOT_DEFINED
!!$
!!$  type(dictionary), pointer :: dict_errors=>null() !< the global dictionaries of possible errors, nullified if not initialized
!!$  type(dictionary), pointer :: dict_present_error=>null() !< local pointer of present error, nullified if success
!!$
!!$  !> Stack of dict_present_error for nested try (open and close)
!!$  type, private :: error_stack
!!$    type(dictionary), pointer :: current => null()   !< dict_present_error point to here.
!!$    type(error_stack), pointer :: previous => null() !< previous error
!!$  end type error_stack
!!$
!!$  type(error_stack), pointer :: error_pipelines=>null() !< Stack of errors for try clause
!!$
!!$  public :: f_err_initialize,f_err_finalize
!!$  public :: f_err_define,f_err_check,f_err_raise,f_err_clean,f_get_error_dict
!!$
!!$  !public variables of the callback module
!!$  public :: f_err_set_callback,f_err_unset_callback
!!$  public :: f_err_severe,f_err_severe_override,f_err_severe_restore
!!$  public :: f_loc
!!$
!!$contains

  !> Initialize the module, in case it has not be done so
  !! the double initialization should not be possible when library is
  !! correctly set up
  subroutine f_err_initialize()
    implicit none
    !local variables
    call f_err_unset_callback()
    call f_err_unset_last_error_callback()
    call f_err_unset_all_errors_callback()
    if (associated(error_pipelines)) then
       call error_pipelines_clean()
    end if
    allocate(error_pipelines)
    call dict_init(error_pipelines%current)
    dict_present_error=>error_pipelines%current
    call dict_init(dict_errors)
  end subroutine f_err_initialize


  !> Call at the end of the program to finalize the error module (deallocation + report)
  subroutine f_err_finalize()
    implicit none
    call f_err_unset_callback()
    call f_err_unset_last_error_callback()
    call f_err_unset_all_errors_callback()
    call f_err_severe_restore()
    call dict_free(dict_errors)
    call error_pipelines_clean()
  end subroutine f_err_finalize


  !> Clean the stack of dict_present_error for nested try
  subroutine error_pipelines_clean()
    implicit none
    type(error_stack), pointer :: stack
    nullify(dict_present_error)
    do while(associated(error_pipelines))
      call dict_free(error_pipelines%current)
      stack=>error_pipelines%previous
      deallocate(error_pipelines)
      error_pipelines=>stack
    end do
  end subroutine error_pipelines_clean


  !> Define a new error specification and returns the corresponding error code
  !! optionally, an error-specific callback function can be defined
  !! @warning
  !!   gfortran complains on runtime (with option -Wextra) that this routine is recursive.
  recursive subroutine f_err_define(err_name,err_msg,err_id,err_action,callback,callback_data)
    implicit none
    character(len=*), intent(in) :: err_name               !< name of the error
    character(len=*), intent(in) :: err_msg                !< error message
    integer, intent(out) :: err_id                         !< code of the error
    integer(kind=8), intent(in), optional :: callback_data !< ??? not really sure
    character(len=*), intent(in), optional :: err_action   !< ??? not really sure
    external :: callback
    optional :: callback
    !local variable
    type(dictionary), pointer :: dict_error

    !assure initialization of the library in case of misuse
    if (.not. associated(dict_errors)) then !call f_err_initialize()
       write(0,*)'error_handling library not initialized, f_lib_initialized should be called'
       call f_err_severe()
    end if

    err_id=ERR_GENERIC
!    if (associated(dict_errors)) then
       err_id=dict_len(dict_errors)
!    end if

    call dict_init(dict_error)
    call set(dict_error//err_name//ERRID,err_id)
    call set(dict_error//err_name//ERRMSG,err_msg)
    if (present(err_action)) then
       if (len_trim(err_action) /=0)&
            call set(dict_error//err_name//ERRACT,trim(err_action))
    end if
    if (present(callback)) then
       call set(dict_error//err_name//ERRCLBK,f_loc(callback))
       if (present(callback_data)) &
            call set(dict_error//err_name//errclbkadd,callback_data)
    end if

    !if (associated(dict_errors))
    call add(dict_errors,dict_error)

  end subroutine f_err_define


  !> This function returns true if a generic error has been raised
  !! in case of specified errors, it returns true if an error of this kind has been raised
  function f_err_check(err_id,err_name)
    !Returns `.true.` is the error identified by the input arguments is present, `.false.` otherwse. If the arguments are absent 
    !returns `.true.` if any arror has been produced
    implicit none
    integer, intent(in), optional :: err_id            !The code of the error to be checked for
    character(len=*), intent(in), optional :: err_name ! Name of the error to search
    logical :: f_err_check                             ! 
    include 'get_err-inc.f90'

    !check if a particular error has been found
    if (get_error==-1) then
       f_err_check = dict_len(dict_present_error) /= 0
    else
       !otherwise check is some error is present
       f_err_check =  get_error/=0
    end if

  end function f_err_check


  function f_err_raise_str(condition,err_msg,err_id,err_name,callback,callback_data) result(f_err_raise)
    use yaml_strings, only: yaml_toa,f_string
    !use yaml_output, only: yaml_dict_dump,yaml_map
    implicit none
    logical, intent(in) :: condition             !< the condition which raise the error
    integer, intent(in), optional :: err_id                !< the code of the error to be raised.
    !! it should already have been defined by f_err_define
    character(len=*), intent(in), optional :: err_name     !< error name
    type(f_string), intent(in) :: err_msg      !< error message
    integer(kind=8), intent(in), optional :: callback_data !< ??? not really sure
    external :: callback                                   !< action to be performed ???
    optional :: callback
    logical :: f_err_raise
    !local variables
    !    integer :: new_errcode
!    integer(kind=8) :: clbk_data_add

    f_err_raise=condition
    !once the error has been identified add it to the present errors and call callback function if needed
    if (f_err_raise) then
       call f_err_throw(err_msg=err_msg%msg,err_name=err_name,err_id=err_id,&
            callback=callback,callback_data=callback_data)

!!$       clbk_data_add=callback_data_add
!!$       if (present(callback_data)) clbk_data_add=callback_data
!!$
!!$       if (present(callback)) then
!!$          if (present(err_id)) then
!!$             call f_err_throw(err_msg%msg,err_id=err_id,callback=callback,callback_data=clbk_data_add)
!!$          else if (present(err_name)) then
!!$             call f_err_throw(err_msg%msg,err_name=err_name,callback=callback,callback_data=clbk_data_add)
!!$          else
!!$             call f_err_throw(err_msg%msg,callback=callback,callback_data=clbk_data_add)
!!$          end if
!!$       else
!!$          if (present(err_id)) then
!!$             call f_err_throw(err_msg%msg,err_id=err_id,callback_data=clbk_data_add)
!!$          else if (present(err_name)) then
!!$             call f_err_throw(err_msg%msg,err_name=err_name,callback_data=clbk_data_add)
!!$          else
!!$             call f_err_throw(err_msg%msg,callback_data=clbk_data_add)
!!$          end if
!!$       end if
    end if
  end function f_err_raise_str


  !> This routine should be generalized to allow the possiblity of addin customized message at the
  !! raise of the error. Also customized callback should be allowed
  !! @warning This function is detected as recursive by gfortran
  recursive function f_err_raise(condition,err_msg,err_id,err_name,callback,callback_data)
    use yaml_strings, only: yaml_toa
    !use yaml_output, only: yaml_dict_dump,yaml_map
    implicit none
    logical, intent(in), optional :: condition             !< the condition which raise the error
    integer, intent(in), optional :: err_id                !< the code of the error to be raised.
                                                           !! it should already have been defined by f_err_define
    character(len=*), intent(in), optional :: err_name     !< error name
    character(len=*), intent(in), optional :: err_msg      !< error message
    integer(kind=8), intent(in), optional :: callback_data !< ??? not really sure
    external :: callback                                   !< action to be performed ???
    optional :: callback
    logical :: f_err_raise
    !local variables
!    integer :: new_errcode
!    integer(kind=8) :: clbk_data_add
!    integer(kind=8) :: clbk_add
!    character(len=max_field_length), dimension(1) :: keys
!    type(dictionary), pointer :: dict_tmp
!    character(len=max_field_length) :: message

    if (present(condition)) then
       f_err_raise=condition
    else
       f_err_raise=.true.
    end if
    !once the error has been identified add it to the present errors and call callback function if needed
    if (f_err_raise) then
       call f_err_throw(err_msg=err_msg,err_name=err_name,err_id=err_id,&
            callback=callback,callback_data=callback_data)
!!$       !throw the error with the annoying stuff of optional variables
!!$       if (present(err_msg)) then
!!$          message(1:len(message))=err_msg
!!$       else
!!$          message(1:len(message))='UNKNOWN'
!!$       end if
!!$
!!$       clbk_data_add=callback_data_add
!!$       if (present(callback_data)) clbk_data_add=callback_data
!!$
!!$       if (present(callback)) then
!!$          if (present(err_id)) then
!!$             call f_err_throw(message,err_id=err_id,callback=callback,callback_data=clbk_data_add)
!!$          else if (present(err_name)) then
!!$             call f_err_throw(message,err_name=err_name,callback=callback,callback_data=clbk_data_add)
!!$          else
!!$             call f_err_throw(message,callback=callback,callback_data=clbk_data_add)
!!$          end if
!!$       else
!!$          if (present(err_id)) then
!!$             call f_err_throw(message,err_id=err_id,callback_data=clbk_data_add)
!!$          else if (present(err_name)) then
!!$             call f_err_throw(message,err_name=err_name,callback_data=clbk_data_add)
!!$          else
!!$             call f_err_throw(message,callback_data=clbk_data_add)
!!$          end if
!!$       end if
    end if
  end function f_err_raise


  !> Raise the error indicated
  !! @warning:  This routine might formally call itself, i.e. it uses methods that might raise exceptions
  !! it is developer's responsibility to avoid deadlocks
  recursive subroutine f_err_throw_c(err_msg,err_id,err_name,callback,callback_data)
    use yaml_strings, only: yaml_toa,f_char_ptr
    implicit none
    integer, intent(in), optional :: err_id                    !< The code of the error to be raised.
                                                               !! it should already have been defined by f_err_define
    character(len=*), intent(in), optional :: err_name         !< error name
    character(len=*), intent(in), optional :: err_msg          !< error message
    !type(dictionary), pointer, optional :: err_dict           !< Add a dictionary instead of a message
    integer(kind=8), intent(in), optional :: callback_data     !< ??? not really sure
    external :: callback
    optional :: callback
    !local variables
    integer :: new_errcode
    integer(kind=8) :: clbk_add,clbk_data_add
    character(len=max_field_length), dimension(1) :: keys
    type(dictionary), pointer :: dict_tmp
    character(len=max_field_length) :: message

    !to prevent infinite loop due to not association of the error handling
    if (.not. associated(dict_present_error)) then
       write(0,*) 'error_handling library not initialized'
       call f_err_severe()
    end if
    !search the error ID (the generic error is always the first)
    new_errcode=ERR_GENERIC
    if (present(err_name)) then
       new_errcode= max(dict_errors .index. err_name,ERR_GENERIC)
       !add a potentially verbose error list in case the specified error name has not been found
       if ((dict_errors .index. err_name) < ERR_GENERIC) then
          call f_dump_possible_errors(f_char_ptr('Error raised, name entered= '//trim(err_name)//&
               '. Errorcode found='//trim(yaml_toa(new_errcode))//&
               '; index function returned'//trim(yaml_toa(dict_errors .index. err_name))))
       end if
    else if (present(err_id)) then
       new_errcode=ERR_GENERIC
       if (err_id < dict_len(dict_errors)) then
          new_errcode=err_id
       else
          !add a potentially verbose error list in case the specified error name has not been found
          call f_dump_possible_errors(f_char_ptr('Error raised, id='//trim(yaml_toa(err_id))))
       end if
    end if

    !if (present(err_dict)) then
    !   !to prevent infinite loop due to not association of the error handling
    !   if (.not. associated(err_dict)) then
    !      write(0,*) 'dictionary err_dict or initialized in f_err_throw'
    !      call f_err_severe()
    !   else
    !   call add(dict_present_error,&
    !        dict_new((/ ERRID .is. yaml_toa(new_errcode),ERR_ADD_INFO .is. err_dict/)))
    !   end if
    !else
       if (present(err_msg)) then
          message(1:len(message))=err_msg
       else
          message(1:len(message))='UNKNOWN'
       end if
       call add(dict_present_error,&
            dict_new((/ ERRID .is. yaml_toa(new_errcode),ERR_ADD_INFO .is. message/)))
    !end if

    !if we are in a try-catch environment, no callback has
    !to be called after the error is produced
    if (try_environment) then
       !identify callback function
       clbk_add=f_loc(f_err_ignore)
       clbk_data_add=int(0,kind=8)
    else
       !identify callback function
       clbk_add=callback_add
       clbk_data_add=callback_data_add
       if (present(callback_data)) clbk_data_add=callback_data
       if (present(callback)) then
          clbk_add=f_loc(callback)
       else
          !find the callback in the error definition
          !these data can be inserted in a function
          dict_tmp=>f_get_error_dict(new_errcode)
          !that is how dict_keys function should be called
          if (dict_size(dict_tmp) <= size(keys)) keys=dict_keys(dict_tmp)
          dict_tmp=>dict_tmp//trim(keys(1))

          if (has_key(dict_tmp,ERRCLBK)) clbk_add=dict_tmp//ERRCLBK
          if (has_key(dict_tmp,errclbkadd)) clbk_data_add=dict_tmp//errclbkadd
       end if
    end if
    call err_abort(clbk_add,clbk_data_add)
  end subroutine f_err_throw_c


  subroutine f_err_throw_str(message,err_id,err_name,callback,callback_data)
    use yaml_strings, only: f_string
    implicit none
    integer, intent(in), optional :: err_id                    !< The code of the error to be raised.
    !! it should already have been defined by f_err_define
    character(len=*), intent(in), optional :: err_name         !< error name
    type(f_string), intent(in) :: message         !< error message
    !type(dictionary), pointer, optional :: err_dict           !< Add a dictionary instead of a message
    integer(kind=8), intent(in), optional :: callback_data     !< ??? not really sure
    external :: callback
    optional :: callback
    !local variables
!    integer(kind=8) :: clbk_data_add
!!$    integer(kind=8) :: clbk_add


    call f_err_throw_c(err_msg=message%msg,err_id=err_id,err_name=err_name,&
         callback=callback,callback_data=callback_data)

!!$
!!$    clbk_data_add=callback_data_add
!!$    if (present(callback_data)) clbk_data_add=callback_data
!!$
!!$    if (present(callback)) then
!!$       if (present(err_id)) then
!!$          call f_err_throw_c(message%msg,err_id=err_id,callback=callback,callback_data=clbk_data_add)
!!$       else if (present(err_name)) then
!!$          call f_err_throw_c(message%msg,err_name=err_name,callback=callback,callback_data=clbk_data_add)
!!$       else
!!$          call f_err_throw_c(message%msg,callback=callback,callback_data=clbk_data_add)
!!$       end if
!!$    else
!!$       if (present(err_id)) then
!!$          call f_err_throw_c(message%msg,err_id=err_id,callback_data=clbk_data_add)
!!$       else if (present(err_name)) then
!!$          call f_err_throw_c(message%msg,err_name=err_name,callback_data=clbk_data_add)
!!$       else
!!$          call f_err_throw_c(message%msg,callback_data=clbk_data_add)
!!$       end if
!!$    end if
  end subroutine f_err_throw_str


  !> Get the error ierror as a dictionary
  !! use as dict=>f_get_error_dict()
  function f_get_error_dict(ierror)
    implicit none
    integer, intent(in), optional :: ierror
    type(dictionary), pointer :: f_get_error_dict
    !local variables
    integer :: ierr
    if (present(ierror)) then
       ierr=ierror
    else
       !Last error
       ierr=f_get_last_error() !dict_len(dict_present_error)-1
    end if
    f_get_error_dict=>dict_errors//ierr

  end function f_get_error_dict


  !> Get the number of errors (public)
  function f_get_no_of_errors()
    implicit none
    integer :: f_get_no_of_errors

    f_get_no_of_errors=max(dict_len(dict_present_error),0)
  end function f_get_no_of_errors


  function f_get_past_error(ierr_num,add_msg)
    implicit none
    integer, intent(in) :: ierr_num
    character(len=*), optional :: add_msg
    integer :: f_get_past_error

    f_get_past_error=get_error_id(ierr_num)
    if (present(add_msg)) call get_error_msg(ierr_num,add_msg)

  end function f_get_past_error


  !> Get the Id of the error ierr (private)
  function get_error_id(ierr)
    implicit none
    integer, intent(in) :: ierr
    integer :: get_error_id
    if (ierr >= 0) then
       get_error_id=dict_present_error//ierr//ERRID
    else
       get_error_id=0
    end if
  end function get_error_id


  !> Get the message of the error given by its number (private).
  subroutine get_error_msg(ierr,add_msg)
    implicit none
    integer, intent(in) :: ierr
    character(len=*), intent(out) :: add_msg

    if (ierr >=0) then
       add_msg=dict_present_error//ierr//ERR_ADD_INFO
    else
       add_msg=repeat(' ',len(add_msg))
    end if
  end subroutine get_error_msg


  !> Identify id of last error occurred (public)
  function f_get_last_error(add_msg)
    implicit none
    character(len=*), intent(out), optional :: add_msg
    integer :: f_get_last_error
    !local variables
    integer :: ierr
    ierr=dict_len(dict_present_error)-1
    if (ierr >= 0) then
       f_get_last_error=dict_present_error//ierr//ERRID
       if (present(add_msg)) add_msg=dict_present_error//ierr//ERR_ADD_INFO

    else
       f_get_last_error=0
       if (present(add_msg)) add_msg=repeat(' ',len(add_msg))
    end if
  end function f_get_last_error


  !> Get all the errors that are present on the error pipe
  subroutine f_get_errors(dict)
    implicit none
    type(dictionary), pointer, intent(out) :: dict !<output dictionary containing the result
    nullify(dict)
    call dict_copy(src=dict_present_error,dest=dict)

  end subroutine f_get_errors


  !> Clean the dictionary of present errors
   subroutine f_err_clean()
    implicit none
    nullify(dict_present_error)
    call dict_free(error_pipelines%current)
    call dict_init(error_pipelines%current)
    dict_present_error=>error_pipelines%current
  end subroutine f_err_clean


  !> Clean last error, if any and get message.
  function f_err_pop(err_id, err_name, add_msg)
    implicit none
    integer, intent(in), optional :: err_id            !< The code of the error to be checked for
    character(len=*), intent(in), optional :: err_name !< Name of the error to search
    character(len=*), intent(out), optional :: add_msg !< additional message
    integer :: f_err_pop
    !local variables
    include 'get_err-inc.f90'

    if (get_error == -1) then
       ierr=dict_len(dict_present_error)-1
       get_error = 3
    end if
    if (ierr >= 0 .and. get_error > 0) then
       f_err_pop=dict_present_error//ierr//ERRID
       if (present(add_msg)) add_msg=dict_present_error//ierr//ERR_ADD_INFO
       call dict_remove(dict_present_error, ierr)
       if (.not.associated(dict_present_error)) then
         !call dict_init(dict_present_error)
         call dict_init(error_pipelines%current)
         dict_present_error => error_pipelines%current
       end if
    else
       f_err_pop=0
       if (present(add_msg)) add_msg=repeat(' ',len(add_msg))
    end if
  end function f_err_pop


  !> Activate the exception handling for all errors
  !! also the errors which have f_err_severe as callbacks
  !! multiple calls to f_err_open_try have the same effect of one call
  subroutine f_err_open_try()
    implicit none
    type(error_stack), pointer :: stack
    !call f_err_set_callback(f_err_ignore)
    try_environment=.true.
    allocate(stack)
    stack%previous=>error_pipelines
    error_pipelines=>stack
    call dict_init(error_pipelines%current)
    dict_present_error=>error_pipelines%current

  end subroutine f_err_open_try


  !> Close the try environment. At the end of the try environment
  !! the errors are cleaned. To recover an error in a try environment
  !! the correct behaviour is to perform f_err_check before calling
  !! f_err_close_try
  subroutine f_err_close_try(exceptions)
    implicit none
    type(dictionary), pointer, intent(out), optional :: exceptions !<retrieve the exceptions
    !local variables
    type(error_stack), pointer :: stack

    !retrieve the exceptions if needed
    if (present(exceptions)) then
       nullify(exceptions)
       if (f_err_check()) call f_get_errors(exceptions)
    end if

    !call f_err_unset_callback()
    if (associated(error_pipelines%previous)) then
      nullify(dict_present_error)
      call dict_free(error_pipelines%current)
      stack=>error_pipelines%previous
      deallocate(error_pipelines)
      error_pipelines=>stack
      dict_present_error=>error_pipelines%current
      try_environment=associated(error_pipelines%previous)
    else
      call f_err_clean() !no errors anymore for this stack
      try_environment=.false.
    end if
  end subroutine f_err_close_try


  !> Give all define errors in a dictionary
  function f_get_error_definitions()
    implicit none
    type(dictionary), pointer :: f_get_error_definitions
    f_get_error_definitions => dict_errors
  end function f_get_error_definitions

!!$end module error_handling
