!{\src2tex{textfont=tt}}
!!****m* ABINIT/defs_basis
!! NAME
!! defs_basis
!!
!! FUNCTION
!! This module contains definitions for a number of named constants and
!! physical constants.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (HM, XG,XW)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!!
!! NOTES
!! Of the named constants,
!! by far the most important are those that define the 'kind' types of
!! virtually all the variables used in a (well-written) FORTRAN 90 code
!! the content of this file is derived from 'Numerical Recipes in Fortran 90'
!! W.H. Press et al., volume 2 of 'Fortran Numerical Recipes', Cambridge
!! University Press, Second Edition (1996), p. 937 and 1361
!!
!! SOURCE

module defs_basis

 implicit none

!Real constants related to the golden number
 real(8), parameter :: gold=1.618033988749894848204586834365638117720309179d0
 real(8), parameter :: goldenratio=2.d0-gold

!Real constants derived from pi
 real(8), parameter :: pi=3.141592653589793238462643383279502884197d0
!The following are not used
!real(8), parameter :: rad_to_deg=180.d0/pi
!real(8), parameter :: deg_to_rad=one/rad_to_deg
!real(8), parameter :: half_pi=pi*half
!real(8), parameter :: third_pi=pi*third
!real(8), parameter :: quarter_pi=pi*quarter
!real(8), parameter :: two_thirds_pi=two_thirds*pi



!Real physical constants
!Revised fundamental constants from http://physics.nist.gov/cuu/Constants/index.html
!(from 2006 least squares adjustment)
 real(8), parameter :: Bohr_Ang=0.52917720859d0    ! 1 Bohr, in Angstrom
 real(8), parameter :: Ha_cmm1=219474.6313705d0  ! 1 Hartree, in cm^-1
 real(8), parameter :: Ha_eV=27.21138386d0 ! 1 Hartree, in eV
 real(8), parameter :: Ha_K=315774.65d0 ! 1Hartree, in Kelvin
 real(8), parameter :: Ha_THz=6579.683920722d0 ! 1 Hartree, in THz
 real(8), parameter :: Ha_J=4.35974394d-18    !1 Hartree, in J
 real(8), parameter :: e_Cb=1.602176487d-19 ! minus the electron charge, in Coulomb
 real(8), parameter :: kb_HaK=8.617343d-5/Ha_eV ! Boltzmann constant in Ha/K
 real(8), parameter :: amu_emass=1.660538782d-27/9.10938215d-31 ! 1 atomic mass unit, in electronic mass
!This value is 1Ha/bohr^3 in 1d9 J/m^3
!real(8), parameter :: HaBohr3_GPa=29421.033d0 ! 1 Ha/Bohr^3, in GPa
 real(8), parameter :: HaBohr3_GPa=Ha_eV/Bohr_Ang**3*e_Cb*1.0d+21 ! 1 Ha/Bohr^3, in GPa
 real(8), parameter :: Avogadro=6.02214179d23 ! per mole
 real(8), parameter :: AmuBohr2_Cm2=e_Cb*1.0d20/(Bohr_Ang*Bohr_Ang)
 real(8), parameter :: InvFineStruct=137.035999679d0  ! Inverse of fine structure constant
 real(8), parameter :: Sp_Lt=2.99792458d8/2.1876912633d6 ! speed of light in atomic units
 real(8), parameter :: Time_Sec=2.418884326505D-17 !  Atomic unit of time, in seconds
 real(8), parameter :: BField_Tesla=0.0 ! Atomic unit of induction field, in Tesla.
 real(8), parameter :: Ha_kcalmol=Ha_J*Avogadro*0.238902957619d0*1.d-3 !A hartree in kcal/mol

end module defs_basis
!!***

