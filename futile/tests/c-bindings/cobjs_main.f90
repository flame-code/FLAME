module my_object
  type object
     integer :: a_value
  end type object
end module my_object

program test_cobjs

  use module_f_objects
  use my_object

  implicit none

  type(object) :: obj
  type(signal_ctx) :: sig

  call f_lib_initialize()

  call f_object_new_("my_object")  
  call f_object_add_signal("my_object", "init", 1)
  call f_object_add_signal("my_object", "set", 1)

  call connect_in_c()

  if (f_object_signal_prepare("my_object", "init", sig)) then
     call f_object_signal_add_arg(sig, obj)
     call f_object_signal_emit(sig)
  else
     obj%a_value = -1
  end if
  call obj_print(obj)

  call f_lib_finalize()

end program test_cobjs

subroutine obj_print(obj)
  use my_object
  use yaml_output
  type(object), intent(in) :: obj

  call yaml_map("value", obj%a_value)
end subroutine obj_print

subroutine obj_set(obj, val)
  use my_object
  use module_f_objects
  type(object), intent(inout) :: obj
  integer, intent(in) :: val

  type(signal_ctx) :: sig

  obj%a_value = val
  if (f_object_signal_prepare("my_object", "set", sig)) then
     call f_object_signal_add_arg(sig, obj)
     call f_object_signal_emit(sig)
  end if
end subroutine obj_set
