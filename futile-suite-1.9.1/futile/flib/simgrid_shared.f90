! Fortran bindings for simgrid shared allocations
module smpi_shared
  use iso_c_binding
  public :: smpi_shared_malloc, smpi_shared_free
  interface

  function smpi_shared_malloc(size, file, line) bind(C, name='smpi_shared_malloc')
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char
    implicit none
    type(c_ptr) :: smpi_shared_malloc
    integer(kind=8), value :: size
    character(kind=c_char), intent(in) :: file(*)
    integer(kind=8), value :: line
  end function smpi_shared_malloc

  subroutine smpi_shared_free(p) bind(C, name='smpi_shared_free')
    use f_precisions, only: f_address
    implicit none
    integer(f_address), intent(in), value :: p
  end subroutine smpi_shared_free

  end interface
end module smpi_shared
