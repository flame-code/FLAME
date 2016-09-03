!> @file
!!  Define a module to extend numerical functions
!! @author
!!    Copyright (C) 2015-2015 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
module numerics
  use f_precisions, only: db=>f_double
  implicit none
  private

  !constants
  real(db), parameter, public :: pi=3.141592653589793238462643383279502884197169399375105820974944592_db
  real(db), parameter, public :: twopi=6.283185307179586476925286766559005768394338798750211641949889185_db
  real(db), parameter, public :: fourpi=12.56637061435917295385057353311801153678867759750042328389977837_db
  real(db), parameter, public :: oneopi=0.3183098861837906715377675267450287240689192914809128974953346881_db
  real(db), parameter, public :: oneotwopi=0.1591549430918953357688837633725143620344596457404564487476673441_db
  real(db), parameter, public :: oneofourpi=0.07957747154594766788444188168625718101722982287022822437383367203_db
  real(db), parameter, public :: oneoeightpi=0.03978873577297383394222094084312859050861491143511411218691683601_db
  !rationals
  real(db), parameter, public :: onehalf=0.5_db
  real(db), parameter, public :: onethird=0.33333333333333333333333333333333333333333333333333333333333333333333_db
  real(db), parameter, public :: onequarter=0.25_db
  

  !>> Physical constants.
  !> 1 AU in angstroem
  real(db), parameter, public :: Bohr_Ang = 0.52917721092_db
  !> 1 Hartree, in cm^-1 (from abinit 5.7.x)
  real(db), parameter, public :: Ha_cmm1=219474.6313705_db 
  !> 1 Hartree, in eV
  real(db), parameter, public :: Ha_eV=27.21138505_db                           !< 1 Hartree in eV
  real(db), parameter, public :: eV_Ha=3.674932379e-2_db                        !< 1 ev, in Hartree
  real(db), parameter, public :: Ha_K=315774.664550534774_db                    !< 1 Hartree, in Kelvin
  real(db), parameter, public :: Ha_THz=6579.683920722_db                       !< 1 Hartree, in THz
  real(db), parameter, public :: Ha_J=4.35974394d-18                            !< 1 Hartree, in J
  real(db), parameter, public :: e_Cb=1.602176487d-19                           !< minus the electron charge, in Coulomb
  real(db), parameter, public :: kb_HaK=8.617343d-5/Ha_eV                       !< Boltzmann constant in Ha/K
  real(db), parameter, public :: amu_emass=1.660538782e-27_db/9.10938215e-31_db !< 1 atomic mass unit, in electronic mass
  real(db), parameter, public :: AU_GPa=29421.010901602753_db                   !< 1 Ha/Bohr^3 in GPa
  real(db), parameter, public :: Radian_Degree = 57.29577951308232087679_db     !< 1 radian in degrees
  real(db), parameter, public :: eVAng_HaBohr = Bohr_Ang*eV_Ha                  !< convert forces from eV/Angstroem to hartree/bohr
  real(db), parameter, public :: Debye_AU = 0.393430307_db                      !< 1 Debye in Hartree*Bohr
  !>  1 AU of force in dyn
  real(db), parameter, public :: dyn_AU=8.238722514e-3_db
  real(db), parameter, public :: kcalMol_Ha = 0.001593601437458137_db        !< from kcal_th/mol to hartree
  !!(thermochemical calorie used in amber: 1cal_th=4.184J)
  !!also see: http://archive.ambermd.org/201009/0039.html
  !!convert forces from kcal_th/mol/angstrom to hartree/bohr
  real(db), parameter, public :: kcalMolAng_HaBohr =0.0008432975639921999_db 


  
  interface safe_exp
     module procedure safe_dexp
  end interface safe_exp

  public :: safe_exp

  contains

    !> fpe-free way of calling exp.
    !! Crop the results to zero in the case of underflow
    pure function safe_dexp(x,extra_crop_order,underflow) result(ex)
      implicit none
      !> argument of the exponential function
      double precision, intent(in) :: x
      !> determine the value under which the result is assumed to be zero
      double precision, intent(in), optional :: underflow
      !> further restrict the valid range of the function
      !! by the order of magnitude indicated.
      !! Useful when the function has to be multiplied by extra terms
      !! The default is log of epsilon**2
      integer, intent(in), optional :: extra_crop_order
      double precision :: ex
      !local variables
      !> if the exponent is bigger than this value, the result is tiny(1.0)
      double precision, parameter :: mn_expo=-708.396418532264d0 ! = log(tiny(1.d0))
      !> if the exponent is lower than this value, the result is huge(1.0)
      double precision, parameter :: mx_expo=709.78271289338397d0 ! = log(huge(1.d0))
      !> the value of the cropping
      double precision, parameter :: crop_expo=72.0873067782343d0 ! = -2*log(epsilon(1.d0))
      double precision :: crop,mn,mx

      if (x==0.d0) then
         ex=1.d0
         return
      end if
      crop=crop_expo
      if (present(extra_crop_order)) crop=real(extra_crop_order,kind=8)
      mn=mn_expo+crop
      mx=mx_expo-crop
      if (present(underflow)) mn=log(abs(underflow))
      if (x > mn .and. x< mx) then
         ex=exp(x)
      else if (x <= mn) then
         ex=0.d0
      else
         ex=exp(mx)
      end if

    end function safe_dexp

    !> give a function which takes into account overflows and underflows even in the gaussian arguments
    pure function safe_gaussian(x0,x,alpha) result(gau)
      implicit none
      double precision, intent(in) :: x0 !< gaussian center
      double precision, intent(in) :: x !< argument
      !double precision, intent(in), optional :: sigma !<standard deviation
      double precision, intent(in) :: alpha !< exponent
      double precision :: gau
      !local variables
      !> if the sqrt is bigger than this value, the result is tiny(1.0)
      double precision, parameter :: mn_sqrt= sqrt(tiny(1.d0))
      !> if the sqrt is lower than this value, the result is huge(1.0)
      double precision, parameter :: mx_sqrt= sqrt(huge(1.d0))

      double precision :: gau_arg,xd

      !evaluate in safe way gau_arg
      xd=abs(x-x0) !assume that this is legal
      if (xd > mn_sqrt .and. xd< mx_sqrt) then
         xd=xd*xd
         gau_arg=-alpha*xd
         !if everything goes fine
         gau=safe_exp(gau_arg)
      else if (x <= mn_sqrt) then
         gau=1.d0
      else
         gau=0.d0
      end if
    end function safe_gaussian

end module numerics
