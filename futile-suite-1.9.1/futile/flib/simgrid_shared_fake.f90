!> fake module to substitute the simgrid shared allocation module..
!option 1 : crash
!option 2 : default to malloc ?
module smpi_shared
  use iso_c_binding
  public :: smpi_shared_malloc, smpi_shared_free
  interface
    pure function smpi_shared_malloc(size, file, line) bind(C, name='malloc')!'malloc')
    use, intrinsic :: iso_c_binding, only: c_ptr, c_char
      implicit none
    type(c_ptr) :: smpi_shared_malloc
    integer(kind=8), intent(in), value :: size
    character(kind=c_char), intent(in) :: file(*)
    integer(kind=8), intent(in), value :: line
    integer ierror
!    stop "trying to SMPI run shared allocations without --enable-simgrid-shared"
    end function smpi_shared_malloc

    pure subroutine smpi_shared_free(p) bind(C, name='mybindfree')!'malloc')
    use f_precisions, only: f_address
      implicit none
    integer(f_address), intent(in) :: p
!    stop "trying to SMPI run shared allocations without --enable-simgrid-shared"
    end subroutine smpi_shared_free
    
  end interface
end module smpi_shared

