!*****************************************************************************************
module mod_defs
  !Same way used in module f_precisions of futile
  !use, intrinsic :: iso_fortran_env
  implicit none
  public
  !for reals and complex, to be verified if supported
  integer, parameter:: fsp=selected_real_kind(6, 37) !simple precision specification, to be used for real variables
  integer, parameter:: fdp=selected_real_kind(15, 307) !double precision for real 
  !integer, parameter:: fqp=selected_real_kind(33, 4931) !parameter indicating the quadrupole precision, if supported by the fortran processor
  integer, parameter:: fqp=selected_real_kind(15, 307) !TO_BE_CORRECTED
  !integer, parameter:: fsp=REAL32
  !integer, parameter:: fdp=REAL64
  !integer, parameter:: fqp=REAL128
end module mod_defs
!*****************************************************************************************
