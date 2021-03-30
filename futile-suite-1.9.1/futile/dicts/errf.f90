subroutine bind_err_define(id, name, lname, message, lmessage)
  use dictionaries, only: f_err_define
  implicit none
  integer, intent(out) :: id
  integer, intent(in) :: lname
  character(len = lname), intent(in) :: name
  integer, intent(in) :: lmessage
  character(len = lmessage), intent(in) :: message

  call f_err_define(name, message, id)
end subroutine bind_err_define

subroutine bind_err_define_with_action(id, name, lname, message, lmessage, action, laction)
  use dictionaries, only: f_err_define
  implicit none
  integer, intent(out) :: id
  integer, intent(in) :: lname
  character(len = lname), intent(in) :: name
  integer, intent(in) :: lmessage
  character(len = lmessage), intent(in) :: message
  integer, intent(in) :: laction
  character(len = laction), intent(in) :: action

  call f_err_define(name, message, id, action)
end subroutine bind_err_define_with_action

subroutine bind_err_open_try()
  use dictionaries, only: f_err_open_try
  implicit none

  call f_err_open_try()
end subroutine bind_err_open_try

subroutine bind_err_close_try()
  use dictionaries, only: f_err_close_try
  implicit none

  call f_err_close_try()
end subroutine bind_err_close_try

subroutine bind_err_throw_by_id(id)
  use dictionaries, only: f_err_throw
  implicit none
  integer, intent(in) :: id

  call f_err_throw(err_id = id)
end subroutine bind_err_throw_by_id

subroutine bind_err_throw_by_id_with_mess(message, ln, id)
  use dictionaries, only: f_err_throw
  implicit none
  integer, intent(in) :: ln
  character(len = ln), intent(in) :: message
  integer, intent(in) :: id

  call f_err_throw(message, err_id = id)
end subroutine bind_err_throw_by_id_with_mess

subroutine bind_err_throw_by_name(name, lname)
  use dictionaries, only: f_err_throw
  implicit none
  integer, intent(in) :: lname
  character(len = lname), intent(in) :: name

  call f_err_throw(err_name = name)
end subroutine bind_err_throw_by_name

subroutine bind_err_throw_by_name_with_mess(message, ln, name, lname)
  use dictionaries, only: f_err_throw
  implicit none
  integer, intent(in) :: ln
  character(len = ln), intent(in) :: message
  integer, intent(in) :: lname
  character(len = lname), intent(in) :: name

  call f_err_throw(message, err_name = name)
end subroutine bind_err_throw_by_name_with_mess

subroutine bind_err_check(ret)
  use dictionaries, only: f_err_check
  implicit none
  integer, intent(out) :: ret

  ret = 0
  if (f_err_check()) ret = 1
end subroutine bind_err_check

subroutine bind_err_check_by_id(ret, id)
  use dictionaries, only: f_err_check
  implicit none
  integer, intent(out) :: ret
  integer, intent(in) :: id

  ret = 0
  if (f_err_check(err_id = id)) ret = 1
end subroutine bind_err_check_by_id

subroutine bind_err_check_by_name(ret, name, ln)
  use dictionaries, only: f_err_check
  implicit none
  integer, intent(out) :: ret
  integer, intent(in) :: ln
  character(len = ln), intent(in) :: name

  ret = 0
  if (f_err_check(err_name = name)) ret = 1
end subroutine bind_err_check_by_name
