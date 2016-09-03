!!! [init]
subroutine init()
!!! [init]
  use module_f_objects
  use yaml_output

  integer :: sid
  type(kernel_ctx) :: kernel

  interface
     subroutine pong(p)
       use dictionaries
       character(len = max_field_length), intent(out) :: p
     end subroutine pong
  end interface

  call yaml_map("plugin", "initialised")

  call f_object_kernel_new(kernel, pong, 1)
  call f_object_signal_connect("class", "ping", kernel, sid)
end subroutine init

subroutine pong(p)
  use dictionaries
  character(len = max_field_length), intent(out) :: p

  write(p, "(A)") "pong"
end subroutine pong
