!*****************************************************************************************
module mod_defs
  !Same way used in module f_precisions of futile
  implicit none
  public
  !for reals and complex, to be verified if supported
  integer, parameter:: fsp=selected_real_kind(6, 37) !simple precision specification, to be used for real variables
  integer, parameter:: fdp=selected_real_kind(15, 307) !double precision for real 
  integer, parameter:: fqp=selected_real_kind(33, 4931) !parameter indicating the quadrupole precision, if supported by the fortran processor
end module mod_defs
!*****************************************************************************************
