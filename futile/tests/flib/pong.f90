!!! [init]
subroutine init(caller, val)
!!! [init]
  use module_f_objects
  use yaml_output

  character(len = *), intent(in) :: caller
  integer, intent(in) :: val

  integer :: sid
  type(kernel_ctx) :: kernel

  interface
     subroutine pong(p)
       use dictionaries
       character(len = max_field_length), intent(out) :: p
     end subroutine pong
  end interface

  call yaml_mapping_open("plugin")
  call yaml_map("initialised", .true.)
  call yaml_map("caller", caller)
  call yaml_map("value", val)
  call yaml_mapping_close()

  call f_object_kernel_new(kernel, pong, 1)
  call f_object_signal_connect("class", "ping", kernel, sid)
end subroutine init

subroutine pong(p)
  use dictionaries
  character(len = max_field_length), intent(out) :: p

  write(p, "(A)") "pong"
end subroutine pong
