module module_f_bind
  use f_precisions, only: f_address, f_loc

  type f_bind
     integer :: n_args
     integer(f_address) :: callback
     integer(f_address), dimension(2) :: args
  end type f_bind

contains
  
  subroutine f_bind_undefined(bind)
    type(f_bind), intent(inout) :: bind
    bind%n_args = 0
    bind%callback = 0
    bind%args = 0
  end subroutine f_bind_undefined

  subroutine f_bind_prepare(bind)
    type(f_bind), intent(inout) :: bind
    bind%args = 0
  end subroutine f_bind_prepare

  subroutine f_bind_execute(bind)
    type(f_bind), intent(in) :: bind
    external :: obj1

    if (bind%callback == 0) return
    if (bind%n_args /= count(bind%args /= 0)) stop

    select case(bind%n_args)
    case (0)
       call call_external_c_fromadd(bind%callback)
    case (1)
       call call_external_c_fromadd_data(bind%callback, bind%args(1))
    case (2)
       call call_external_c_fromadd_data_data(bind%callback, bind%args(1), bind%args(2))
    end select
  end subroutine f_bind_execute

end module module_f_bind

subroutine f_bind_define(bind, callback, n_args)
  use module_f_bind
  type(f_bind), intent(inout) :: bind
  integer, intent(in) :: n_args
  external :: callback

  bind%n_args = n_args
  bind%callback = f_loc(callback)
  bind%args = 0
end subroutine f_bind_define

subroutine f_bind_add_arg(bind, obj)
  use module_f_bind
  type(f_bind), intent(inout) :: bind
  external :: obj

  integer :: i

  do i = 1, size(bind%args)
     if (bind%args(i) == 0) then
        bind%args(i) = f_loc(obj)
        return
     end if
  end do
end subroutine f_bind_add_arg
