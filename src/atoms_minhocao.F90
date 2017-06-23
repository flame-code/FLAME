!{\src2tex{textfont=tt}}
!!****f* ABINIT/atmdata
!! NAME
!! atmdata
!!
!! FUNCTION
!! Return atomic data : symbol, covalent radius, atomic mass
!! Atomic masses are those recommended by the commission on Atomic Weights and
!! Isotopic Abundances, Inorganic Chemistry Division, IUPAC, in
!! Pure Appl. Chem. 60, 841 (1988). For Tc, Pm, Po to Ac, Pa and beyond U,
!! none of the isotopes has a half-life greater than 3.0d10 years, and
!! the values provided here do not come from that source.
!!
!! COPYRIGHT
!! Copyright (C) 2000-2010 ABINIT group (XG, MJV)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!
!! INPUTS
!! znucl=atomic number (a real(dp) number ! the nearest
!!  integer is selected in the routine ...)
!!
!! OUTPUT
!! amu=atomic mass (Masses beyond element 103 are fixed at 260)
!! rcov=covalent radius   (Elements beyond 86 have an estimated covalent radius)
!! character(len=2) symbol=atomic symbol
!!
!! PARENTS
!!      abi_etsf_geo_put,bonds_lgth_angles,fresid,invars1,invars2,lwf
!!      make_bonds,mlwfovlp_setup,out1dm,prt_cml,prt_cml2,recip_ylm
!!      wvl_init_type_proj,wvl_tail_corrections
!!
!! CHILDREN
!!
!! SOURCE

subroutine atmdata(amu,rcov,symbol,znucl)
use defs_basis
use blj_params
 implicit none

!Arguments ------------------------------------
!scalars
 real(8),intent(in) :: znucl
 real(8),intent(out) :: amu,rcov
 character(len=2),intent(out) :: symbol

!Local variables-------------------------------
!scalars
 integer :: nucl

! *************************************************************************

 nucl=nint(znucl)
 select case (nucl)
   case(0)
     amu= 1.d0       ; rcov=1.d0/Bohr_Ang    ; symbol='no'
   case(1)
     amu=1.00794d0   ; rcov=0.32d0/Bohr_Ang ; symbol=' H'
   case(2)
     amu=4.002602d0  ; rcov=0.93d0/Bohr_Ang ; symbol='He'
   case(3)
     amu=6.941d0     ; rcov=1.23d0/Bohr_Ang ; symbol='Li'
   case(4)
     amu=9.012182d0  ; rcov=0.90d0/Bohr_Ang ; symbol='Be'
   case(5)
     amu=10.811d0    ; rcov=0.80d0/Bohr_Ang ; symbol=' B'
   case(6)
     amu=12.011d0    ; rcov=0.77d0/Bohr_Ang ; symbol=' C'
   case(7)
     amu=14.00674d0  ; rcov=0.74d0/Bohr_Ang ; symbol=' N'
   case(8)
     amu=15.9994d0   ; rcov=0.73d0/Bohr_Ang ; symbol=' O'
   case(9)
     amu=18.9984032d0; rcov=0.72d0/Bohr_Ang ; symbol=' F'

   case(10)
     amu=20.1797d0   ; rcov=0.71d0/Bohr_Ang ; symbol='Ne'
   case(11)
     amu=22.989768d0 ; rcov=1.54d0/Bohr_Ang ; symbol='Na'
   case(12)
     amu=24.3050d0   ; rcov=1.36d0/Bohr_Ang ; symbol='Mg'
   case(13)
     amu=26.981539d0 ; rcov=1.18d0/Bohr_Ang ; symbol='Al'
   case(14)
     amu=28.0855d0   ; rcov=1.11d0/Bohr_Ang ; symbol='Si'
   case(15)
     amu=30.973762d0 ; rcov=1.06d0/Bohr_Ang ; symbol=' P'
   case(16)
     amu=32.066d0    ; rcov=1.02d0/Bohr_Ang ; symbol=' S'
   case(17)
     amu=35.4527d0   ; rcov=0.99d0/Bohr_Ang ; symbol='Cl'
   case(18)
     amu=39.948d0    ; rcov=0.98d0/Bohr_Ang ; symbol='Ar'
   case(19)
     amu=39.0983d0   ; rcov=2.03d0/Bohr_Ang ; symbol=' K'

   case(20)
     amu=40.078d0    ; rcov=1.74d0/Bohr_Ang ; symbol='Ca'
   case(21)
     amu=44.955910d0 ; rcov=1.44d0/Bohr_Ang ; symbol='Sc'
   case(22)
     amu=47.88d0     ; rcov=1.32d0/Bohr_Ang ; symbol='Ti'
   case(23)
     amu=50.9415d0   ; rcov=1.22d0/Bohr_Ang ; symbol=' V'
   case(24)
     amu=51.9961d0   ; rcov=1.18d0/Bohr_Ang ; symbol='Cr'
   case(25)
     amu=54.93805d0  ; rcov=1.17d0/Bohr_Ang ; symbol='Mn'
   case(26)
     amu=55.847d0    ; rcov=1.17d0/Bohr_Ang ; symbol='Fe'
   case(27)
     amu=58.93320d0  ; rcov=1.16d0/Bohr_Ang ; symbol='Co'
   case(28)
     amu=58.69d0     ; rcov=1.15d0/Bohr_Ang ; symbol='Ni'
   case(29)
     amu=63.546d0    ; rcov=1.17d0/Bohr_Ang ; symbol='Cu'

   case(30)
     amu=65.39d0     ; rcov=1.25d0/Bohr_Ang ; symbol='Zn'
   case(31)
     amu=69.723d0    ; rcov=1.26d0/Bohr_Ang ; symbol='Ga'
   case(32)
     amu=72.61d0     ; rcov=1.22d0/Bohr_Ang ; symbol='Ge'
   case(33)
     amu=74.92159d0  ; rcov=1.20d0/Bohr_Ang ; symbol='As'
   case(34)
     amu=78.96d0     ; rcov=1.16d0/Bohr_Ang ; symbol='Se'
   case(35)
     amu=79.904d0    ; rcov=1.14d0/Bohr_Ang ; symbol='Br'
   case(36)
     amu=83.80d0     ; rcov=1.12d0/Bohr_Ang ; symbol='Kr'
   case(37)
     amu=85.4678d0   ; rcov=2.16d0/Bohr_Ang ; symbol='Rb'
   case(38)
     amu=87.62d0     ; rcov=1.91d0/Bohr_Ang ; symbol='Sr'
   case(39)
     amu=88.90585d0  ; rcov=1.62d0/Bohr_Ang ; symbol=' Y'

   case(40)
     amu=91.224d0    ; rcov=1.45d0/Bohr_Ang ; symbol='Zr'
   case(41)
     amu=92.90638d0  ; rcov=1.34d0/Bohr_Ang ; symbol='Nb'
   case(42)
     amu=95.94d0     ; rcov=1.30d0/Bohr_Ang ; symbol='Mo'
   case(43)
     amu=98.9062d0   ; rcov=1.27d0/Bohr_Ang ; symbol='Tc'
   case(44)
     amu=101.07d0    ; rcov=1.25d0/Bohr_Ang ; symbol='Ru'
   case(45)
     amu=102.9055d0  ; rcov=1.25d0/Bohr_Ang ; symbol='Rh'
   case(46)
     amu=106.42d0    ; rcov=1.28d0/Bohr_Ang ; symbol='Pd'
   case(47)
     amu=107.8682d0  ; rcov=1.34d0/Bohr_Ang ; symbol='Ag'
   case(48)
     amu=112.411d0   ; rcov=1.48d0/Bohr_Ang ; symbol='Cd'
   case(49)
     amu=114.82d0    ; rcov=1.44d0/Bohr_Ang ; symbol='In'

   case(50)
     amu=118.710d0   ; rcov=1.41d0/Bohr_Ang ; symbol='Sn'
   case(51)
     amu=121.753d0   ; rcov=1.40d0/Bohr_Ang ; symbol='Sb'
   case(52)
     amu=127.60d0    ; rcov=1.36d0/Bohr_Ang ; symbol='Te'
   case(53)
     amu=126.90447d0 ; rcov=1.33d0/Bohr_Ang ; symbol=' I'
   case(54)
     amu=131.29d0    ; rcov=1.31d0/Bohr_Ang ; symbol='Xe'
   case(55)
     amu=132.90543d0 ; rcov=2.35d0/Bohr_Ang ; symbol='Cs'
   case(56)
     amu=137.327d0   ; rcov=1.98d0/Bohr_Ang ; symbol='Ba'
   case(57)
     amu=138.9055d0  ; rcov=1.69d0/Bohr_Ang ; symbol='La'
   case(58)
     amu=140.115d0   ; rcov=1.65d0/Bohr_Ang ; symbol='Ce'
   case(59)
     amu=140.90765d0 ; rcov=1.65d0/Bohr_Ang ; symbol='Pr'

   case(60)
     amu=144.24d0    ; rcov=1.64d0/Bohr_Ang ; symbol='Nd'
   case(61)
     amu=147.91d0    ; rcov=1.64d0/Bohr_Ang ; symbol='Pm'
   case(62)
     amu=150.36d0    ; rcov=1.62d0/Bohr_Ang ; symbol='Sm'
   case(63)
     amu=151.965d0   ; rcov=1.85d0/Bohr_Ang ; symbol='Eu'
   case(64)
     amu=157.25d0    ; rcov=1.61d0/Bohr_Ang ; symbol='Gd'
   case(65)
     amu=158.92534d0 ; rcov=1.59d0/Bohr_Ang ; symbol='Tb'
   case(66)
     amu=162.50d0    ; rcov=1.59d0/Bohr_Ang ; symbol='Dy'
   case(67)
     amu=164.93032d0 ; rcov=1.57d0/Bohr_Ang ; symbol='Ho'
   case(68)
     amu=167.26d0    ; rcov=1.57d0/Bohr_Ang ; symbol='Er'
   case(69)
     amu=168.93421d0 ; rcov=1.56d0/Bohr_Ang ; symbol='Tm'

   case(70)
     amu=173.04d0    ; rcov=1.70d0/Bohr_Ang ; symbol='Yb'
   case(71)
     amu=174.967d0   ; rcov=1.56d0/Bohr_Ang ; symbol='Lu'
   case(72)
     amu=178.49d0    ; rcov=1.44d0/Bohr_Ang ; symbol='Hf'
   case(73)
     amu=180.9479d0  ; rcov=1.34d0/Bohr_Ang ; symbol='Ta'
   case(74)
     amu=183.85d0    ; rcov=1.30d0/Bohr_Ang ; symbol=' W'
   case(75)
     amu=186.207d0   ; rcov=1.28d0/Bohr_Ang ; symbol='Re'
   case(76)
     amu=190.2d0     ; rcov=1.26d0/Bohr_Ang ; symbol='Os'
   case(77)
     amu=192.22d0    ; rcov=1.27d0/Bohr_Ang ; symbol='Ir'
   case(78)
     amu=195.08d0    ; rcov=1.30d0/Bohr_Ang ; symbol='Pt'
   case(79)
     amu=196.96654d0 ; rcov=1.34d0/Bohr_Ang ; symbol='Au'

   case(80)
     amu=200.59d0    ; rcov=1.49d0/Bohr_Ang ; symbol='Hg'
   case(81)
     amu=204.3833d0  ; rcov=1.48d0/Bohr_Ang ; symbol='Tl'
   case(82)
     amu=207.2d0     ; rcov=1.47d0/Bohr_Ang ; symbol='Pb'
   case(83)
     amu=208.98037d0 ; rcov=1.46d0/Bohr_Ang ; symbol='Bi'
   case(84)
     amu=209.0d0     ; rcov=1.46d0/Bohr_Ang ; symbol='Po'
   case(85)
     amu=210.0d0     ; rcov=1.45d0/Bohr_Ang ; symbol='At'
   case(86)
     amu=222.0d0     ; rcov=1.45d0/Bohr_Ang ; symbol='Rn'
   case(87)
     amu=223.0d0     ; rcov=2.50d0/Bohr_Ang ; symbol='Fr'
   case(88)
     amu=226.0254d0  ; rcov=2.10d0/Bohr_Ang ; symbol='Ra'
   case(89)
     amu=230.0d0     ; rcov=1.85d0/Bohr_Ang ; symbol='Ac'

   case(90)
     amu=232.0381d0  ; rcov=1.65d0/Bohr_Ang ; symbol='Th'
   case(91)
     amu=231.0359d0  ; rcov=1.50d0/Bohr_Ang ; symbol='Pa'
   case(92)
     amu=238.0289d0  ; rcov=1.42d0/Bohr_Ang ; symbol=' U'
   case(93)
     amu=237.0482d0  ; rcov=1.42d0/Bohr_Ang ; symbol='Np'
   case(94)
     amu=242.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Pu'
   case(95)
     amu=243.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Am'
   case(96)
     amu=247.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Cm'
   case(97)
     amu=247.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Bk'
   case(98)
     amu=249.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Cf'
   case(99)
     amu=254.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Es'

   case(100)
     amu=253.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Fm'
   case(101)
     amu=256.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Md'
   case(102)
     amu=254.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='No'
   case(103)
     amu=257.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Lr'
   case(104)
     amu=260.0d0     ; rcov=1.42d0/Bohr_Ang ; symbol='Xx'
!Lennard Jones Particles
   case(201)
     amu=1.d0        ; rcov=2.d0**(1.d0/6.d0)*sigmalj(1,1)/Bohr_Ang   ; symbol='LA'
   case(202)
     amu=1.d0        ; rcov=2.d0**(1.d0/6.d0)*sigmalj(2,2)/Bohr_Ang   ; symbol='LB'
   case default
     amu=1.d0        ; rcov=1.d0/Bohr_Ang   ; symbol='NO'
 end select

end subroutine atmdata
!!***

subroutine mlj_atmdata(amu,sigma,eps,rcov,symbol,znucl)
use defs_basis
 implicit none

!Arguments ------------------------------------
!scalars
 real(8),intent(in) :: znucl
 real(8),intent(out) :: amu,rcov,sigma,eps
 character(len=2),intent(out) :: symbol

!Local variables-------------------------------
!scalars
 integer :: nucl
! *************************************************************************

 nucl=nint(znucl)
 select case (nucl)
   case(200)
     amu=1.d0        ; sigma=0.50d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Aa'
   case(201)                                                                                             
     amu=1.d0        ; sigma=0.75d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ab'
   case(202)                                                                                             
     amu=1.d0        ; sigma=1.00d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ac'
   case(203)                                                                                             
     amu=1.d0        ; sigma=1.25d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ad'
   case(204)                                                                                             
     amu=1.d0        ; sigma=1.50d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ae'
   case(205)                                                                                             
     amu=1.d0        ; sigma=1.75d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Af'
   case(206)                                                                                             
     amu=1.d0        ; sigma=2.00d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ag'
   case(207)                                                                                             
     amu=1.d0        ; sigma=2.25d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ah'
   case(208)                                                                                             
     amu=1.d0        ; sigma=2.50d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ai'
   case(209)                                                                                             
     amu=1.d0        ; sigma=2.75d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Aj'

   case(210)
     amu=1.d0        ; sigma=0.50d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ba'
   case(211)
     amu=1.d0        ; sigma=0.75d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bb'
   case(212)
     amu=1.d0        ; sigma=1.00d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bc'
   case(213)
     amu=1.d0        ; sigma=1.25d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bd'
   case(214)
     amu=1.d0        ; sigma=1.50d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Be'
   case(215)
     amu=1.d0        ; sigma=1.75d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bf'
   case(216)
     amu=1.d0        ; sigma=2.00d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bg'
   case(217)
     amu=1.d0        ; sigma=2.25d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bh'
   case(218)
     amu=1.d0        ; sigma=2.50d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bi'
   case(219)
     amu=1.d0        ; sigma=2.75d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Bj'

   case(220)
     amu=1.d0        ; sigma=0.50d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ca'
   case(221)
     amu=1.d0        ; sigma=0.75d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Cb'
   case(222)
     amu=1.d0        ; sigma=1.00d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Cc'
   case(223)
     amu=1.d0        ; sigma=1.25d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Cd'
   case(224)
     amu=1.d0        ; sigma=1.50d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ce'
   case(225)
     amu=1.d0        ; sigma=1.75d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Cf'
   case(226)
     amu=1.d0        ; sigma=2.00d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Cg'
   case(227)
     amu=1.d0        ; sigma=2.25d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ch'
   case(228)
     amu=1.d0        ; sigma=2.50d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Ci'
   case(229)
     amu=1.d0        ; sigma=2.75d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='Cj'

   case default
     amu=1.d0        ; sigma=1.00d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; symbol='NO'
 end select
end subroutine mlj_atmdata

subroutine symbol2znucl(amu,rcov,symbol,znucl)

 use defs_basis
 use blj_params

 implicit none

!Arguments ------------------------------------
!scalars
 real(8),intent(out) :: znucl
 real(8),intent(out) :: amu,rcov
 character(len=2),intent(in) :: symbol

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 select case (symbol)
  case('no')
   amu=     1.d0         ; rcov=1.d0/Bohr_Ang    ; znucl=0
  case(' H', 'H ')
   amu=     1.00794d0   ; rcov=0.32d0/Bohr_Ang ; znucl=1
  case('He')
   amu=     4.002602d0  ; rcov=0.93d0/Bohr_Ang ; znucl=2
  case('Li')
   amu=     6.941d0     ; rcov=1.23d0/Bohr_Ang ; znucl=3
  case('Be')
   amu=     9.012182d0  ; rcov=0.90d0/Bohr_Ang ; znucl=4
  case(' B', 'B ')
   amu=     10.811d0    ; rcov=0.80d0/Bohr_Ang ; znucl=5
  case(' C', 'C ')
   amu=     12.011d0    ; rcov=0.77d0/Bohr_Ang ; znucl=6
  case(' N', 'N ')
   amu=     14.00674d0  ; rcov=0.74d0/Bohr_Ang ; znucl=7
  case(' O', 'O ')
   amu=     15.9994d0   ; rcov=0.73d0/Bohr_Ang ; znucl=8
  case(' F', 'F ')
   amu=     18.9984032d0; rcov=0.72d0/Bohr_Ang ; znucl=9
            
  case('Ne')
   amu=     20.1797d0   ; rcov=0.71d0/Bohr_Ang ; znucl=10
  case('Na')                                             
   amu=    22.989768d0 ; rcov=1.54d0/Bohr_Ang ; znucl= 11
  case('Mg')                                             
   amu=     24.3050d0   ; rcov=1.36d0/Bohr_Ang ; znucl=12
  case('Al')                                             
   amu=     26.981539d0 ; rcov=1.18d0/Bohr_Ang ; znucl=13
  case('Si')                                             
   amu=     28.0855d0   ; rcov=1.11d0/Bohr_Ang ; znucl=14
  case(' P', 'P ')                                             
   amu=     30.973762d0 ; rcov=1.06d0/Bohr_Ang ; znucl=15
  case(' S', 'S ')                                             
   amu=     32.066d0    ; rcov=1.02d0/Bohr_Ang ; znucl=16
  case('Cl')                                             
   amu=     35.4527d0   ; rcov=0.99d0/Bohr_Ang ; znucl=17
  case('Ar')                                             
   amu=     39.948d0    ; rcov=0.98d0/Bohr_Ang ; znucl=18
  case(' K', 'K ')                                             
   amu=     39.0983d0   ; rcov=2.03d0/Bohr_Ang ; znucl=19
            
  case('Ca')
   amu=     40.078d0    ; rcov=1.74d0/Bohr_Ang ; znucl=20
  case('Sc')                                             
   amu=     44.955910d0 ; rcov=1.44d0/Bohr_Ang ; znucl=21
  case('Ti')                                             
   amu=     47.88d0     ; rcov=1.32d0/Bohr_Ang ; znucl=22
  case(' V', 'V ')                                             
   amu=     50.9415d0   ; rcov=1.22d0/Bohr_Ang ; znucl=23
  case('Cr')                                             
   amu=     51.9961d0   ; rcov=1.18d0/Bohr_Ang ; znucl=24
  case('Mn')                                             
   amu=     54.93805d0  ; rcov=1.17d0/Bohr_Ang ; znucl=25
  case('Fe')                                             
   amu=     55.847d0    ; rcov=1.17d0/Bohr_Ang ; znucl=26
  case('Co')                                             
   amu=     58.93320d0  ; rcov=1.16d0/Bohr_Ang ; znucl=27
  case('Ni')                                             
   amu=     58.69d0     ; rcov=1.15d0/Bohr_Ang ; znucl=28
  case('Cu')                                             
   amu=     63.546d0    ; rcov=1.17d0/Bohr_Ang ; znucl=29
            
  case('Zn')
   amu=     65.39d0     ; rcov=1.25d0/Bohr_Ang ; znucl=30
  case('Ga')                                             
   amu=     69.723d0    ; rcov=1.26d0/Bohr_Ang ; znucl=31
  case('Ge')                                             
   amu=     72.61d0     ; rcov=1.22d0/Bohr_Ang ; znucl=32
  case('As')                                             
   amu=     74.92159d0  ; rcov=1.20d0/Bohr_Ang ; znucl=33
  case('Se')                                             
   amu=     78.96d0     ; rcov=1.16d0/Bohr_Ang ; znucl=34
  case('Br')                                             
   amu=     79.904d0    ; rcov=1.14d0/Bohr_Ang ; znucl=35
  case('Kr')                                             
   amu=     83.80d0     ; rcov=1.12d0/Bohr_Ang ; znucl=36
  case('Rb')                                             
   amu=     85.4678d0   ; rcov=2.16d0/Bohr_Ang ; znucl=37
  case('Sr')                                             
   amu=     87.62d0     ; rcov=1.91d0/Bohr_Ang ; znucl=38
  case(' Y', 'Y ')                                             
   amu=     88.90585d0  ; rcov=1.62d0/Bohr_Ang ; znucl=39
            
  case('Zr')
   amu=     91.224d0    ; rcov=1.45d0/Bohr_Ang ; znucl=40
  case('Nb')                                             
   amu=     92.90638d0  ; rcov=1.34d0/Bohr_Ang ; znucl=41
  case('Mo')                                             
   amu=     95.94d0     ; rcov=1.30d0/Bohr_Ang ; znucl=42
  case('Tc')                                             
   amu=     98.9062d0   ; rcov=1.27d0/Bohr_Ang ; znucl=43
  case('Ru')                                             
   amu=     101.07d0    ; rcov=1.25d0/Bohr_Ang ; znucl=44
  case('Rh')                                             
   amu=     102.9055d0  ; rcov=1.25d0/Bohr_Ang ; znucl=45
  case('Pd')                                             
   amu=     106.42d0    ; rcov=1.28d0/Bohr_Ang ; znucl=46
  case('Ag')                                             
   amu=     107.8682d0  ; rcov=1.34d0/Bohr_Ang ; znucl=47
  case('Cd')                                             
   amu=     112.411d0   ; rcov=1.48d0/Bohr_Ang ; znucl=48
  case('In')                                             
   amu=     114.82d0    ; rcov=1.44d0/Bohr_Ang ; znucl=49
            
  case('Sn')
   amu=     118.710d0   ; rcov=1.41d0/Bohr_Ang ; znucl=50
  case('Sb')                                             
   amu=     121.753d0   ; rcov=1.40d0/Bohr_Ang ; znucl=51
  case('Te')                                             
   amu=     127.60d0    ; rcov=1.36d0/Bohr_Ang ; znucl=52
  case(' I', 'I ')                                             
   amu=     126.90447d0 ; rcov=1.33d0/Bohr_Ang ; znucl=53
  case('Xe')                                             
   amu=     131.29d0    ; rcov=1.31d0/Bohr_Ang ; znucl=54
  case('Cs')                                             
   amu=     132.90543d0 ; rcov=2.35d0/Bohr_Ang ; znucl=55
  case('Ba')                                             
   amu=     137.327d0   ; rcov=1.98d0/Bohr_Ang ; znucl=56
  case('La')                                             
   amu=     138.9055d0  ; rcov=1.69d0/Bohr_Ang ; znucl=57
  case('Ce')                                             
   amu=     140.115d0   ; rcov=1.65d0/Bohr_Ang ; znucl=58
  case('Pr')                                             
   amu=     140.90765d0 ; rcov=1.65d0/Bohr_Ang ; znucl=59
            
  case('Nd')
   amu=     144.24d0    ; rcov=1.64d0/Bohr_Ang ; znucl=60
  case('Pm')                                             
   amu=     147.91d0    ; rcov=1.64d0/Bohr_Ang ; znucl=61
  case('Sm')                                             
   amu=     150.36d0    ; rcov=1.62d0/Bohr_Ang ; znucl=62
  case('Eu')                                             
   amu=     151.965d0   ; rcov=1.85d0/Bohr_Ang ; znucl=63
  case('Gd')                                             
   amu=     157.25d0    ; rcov=1.61d0/Bohr_Ang ; znucl=64
  case('Tb')                                             
   amu=     158.92534d0 ; rcov=1.59d0/Bohr_Ang ; znucl=65
  case('Dy')                                             
   amu=     162.50d0    ; rcov=1.59d0/Bohr_Ang ; znucl=66
  case('Ho')                                             
   amu=     164.93032d0 ; rcov=1.57d0/Bohr_Ang ; znucl=67
  case('Er')                                             
   amu=     167.26d0    ; rcov=1.57d0/Bohr_Ang ; znucl=68
  case('Tm')                                             
   amu=     168.93421d0 ; rcov=1.56d0/Bohr_Ang ; znucl=69
           
  case('Yb')
   amu=     173.04d0    ; rcov=1.70d0/Bohr_Ang ; znucl=70
  case('Lu')                                             
   amu=     174.967d0   ; rcov=1.56d0/Bohr_Ang ; znucl=71
  case('Hf')                                             
   amu=     178.49d0    ; rcov=1.44d0/Bohr_Ang ; znucl=72
  case('Ta')                                             
   amu=     180.9479d0  ; rcov=1.34d0/Bohr_Ang ; znucl=73
  case(' W', 'W ')                                             
   amu=     183.85d0    ; rcov=1.30d0/Bohr_Ang ; znucl=74
  case('Re')                                             
   amu=     186.207d0   ; rcov=1.28d0/Bohr_Ang ; znucl=75
  case('Os')                                             
   amu=     190.2d0     ; rcov=1.26d0/Bohr_Ang ; znucl=76
  case('Ir')                                             
   amu=     192.22d0    ; rcov=1.27d0/Bohr_Ang ; znucl=77
  case('Pt')                                             
   amu=     195.08d0    ; rcov=1.30d0/Bohr_Ang ; znucl=78
  case('Au')                                             
   amu=     196.96654d0 ; rcov=1.34d0/Bohr_Ang ; znucl=79
           
  case('Hg')
   amu=     200.59d0    ; rcov=1.49d0/Bohr_Ang ; znucl=80
  case('Tl')                                             
   amu=     204.3833d0  ; rcov=1.48d0/Bohr_Ang ; znucl=81
  case('Pb')                                             
   amu=     207.2d0     ; rcov=1.47d0/Bohr_Ang ; znucl=82
  case('Bi')                                             
   amu=     208.98037d0 ; rcov=1.46d0/Bohr_Ang ; znucl=83
  case('Po')                                             
   amu=     209.0d0     ; rcov=1.46d0/Bohr_Ang ; znucl=84
  case('At')                                             
   amu=     210.0d0     ; rcov=1.45d0/Bohr_Ang ; znucl=85
  case('Rn')                                             
   amu=     222.0d0     ; rcov=1.45d0/Bohr_Ang ; znucl=86
  case('Fr')                                             
   amu=     223.0d0     ; rcov=2.50d0/Bohr_Ang ; znucl=87
  case('Ra')                                             
   amu=     226.0254d0  ; rcov=2.10d0/Bohr_Ang ; znucl=88
  case('Ac')                                             
   amu=     230.0d0     ; rcov=1.85d0/Bohr_Ang ; znucl=89
           
  case('Th')
   amu=     232.0381d0  ; rcov=1.65d0/Bohr_Ang ; znucl=90
  case('Pa')                                             
   amu=     231.0359d0  ; rcov=1.50d0/Bohr_Ang ; znucl=91
  case(' U', 'U ')                                             
   amu=     238.0289d0  ; rcov=1.42d0/Bohr_Ang ; znucl=92
  case('Np')                                             
   amu=     237.0482d0  ; rcov=1.42d0/Bohr_Ang ; znucl=93
  case('Pu')                                             
   amu=     242.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=94
  case('Am')                                             
   amu=     243.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=95
  case('Cm')                                             
   amu=     247.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=96
  case('Bk')                                             
   amu=     247.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=97
  case('Cf')                                             
   amu=     249.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=98
  case('Es')                                             
   amu=     254.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=89
            
  case('Fm')
   amu=     253.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=100
  case('Md')
   amu=     256.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=101
  case('No')
   amu=     254.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=102
  case('Lr')
   amu=     257.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=103
  case('Xx')
   amu=     260.0d0     ; rcov=1.42d0/Bohr_Ang ; znucl=104
!Lennard Jones Particles
  case('LA')
     amu=1.d0           ; rcov=2.d0**(1.d0/6.d0)*sigmalj(1,1)/Bohr_Ang ; znucl=201
  case('LB')
     amu=1.d0           ; rcov=2.d0**(1.d0/6.d0)*sigmalj(2,2)/Bohr_Ang ; znucl=202
  case default
   stop "Unknown atom!!!"
 end select

end subroutine symbol2znucl





subroutine mlj_symbol2znucl(amu,sigma,eps,rcov,symbol,znucl)
 use defs_basis
 use blj_params

 implicit none

!Arguments ------------------------------------
!scalars
 real(8),intent(out) :: znucl
 real(8),intent(out) :: amu,rcov,sigma,eps
 character(len=2),intent(in) :: symbol

!Local variables-------------------------------
!scalars
 character(len=500) :: msg

! *************************************************************************

 select case (symbol)
   case('Aa')
               amu=1.d0        ; sigma=0.50d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=200
   case('Ab')                                                                                          
               amu=1.d0        ; sigma=0.75d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=201
   case('Ac')                                                                                         
               amu=1.d0        ; sigma=1.00d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=202
   case('Ad')                                                                                        
               amu=1.d0        ; sigma=1.25d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=203
   case('Ae')                                                                                       
               amu=1.d0        ; sigma=1.50d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=204
   case('Af')                                                                                      
               amu=1.d0        ; sigma=1.75d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=205
   case('Ag')                                                                                     
               amu=1.d0        ; sigma=2.00d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=206
   case('Ah')                                                                                    
               amu=1.d0        ; sigma=2.25d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=207
   case('Ai')                                                                                   
               amu=1.d0        ; sigma=2.50d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=208
   case('Aj')                                                                                  
               amu=1.d0        ; sigma=2.75d0   ; eps=0.50d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=209
                                                                                             
   case('Ba')
               amu=1.d0        ; sigma=0.50d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=210
   case('Bb')
               amu=1.d0        ; sigma=0.75d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=211
   case('Bc')
               amu=1.d0        ; sigma=1.00d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=212
   case('Bd')
               amu=1.d0        ; sigma=1.25d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=213
   case('Be')
               amu=1.d0        ; sigma=1.50d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=214
   case('Bf')
               amu=1.d0        ; sigma=1.75d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=215
   case('Bg')
               amu=1.d0        ; sigma=2.00d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=216
   case('Bh')
               amu=1.d0        ; sigma=2.25d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=217
   case('Bi')
               amu=1.d0        ; sigma=2.50d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=218
   case('Bj')
               amu=1.d0        ; sigma=2.75d0   ; eps=0.75d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=219
                                                                                             
   case('Ca')
               amu=1.d0        ; sigma=0.50d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=220
   case('Cb')
               amu=1.d0        ; sigma=0.75d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=221
   case('Cc')
               amu=1.d0        ; sigma=1.00d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=222
   case('Cd')
               amu=1.d0        ; sigma=1.25d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=223
   case('Ce')
               amu=1.d0        ; sigma=1.50d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=224
   case('Cf')
               amu=1.d0        ; sigma=1.75d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=225
   case('Cg')
               amu=1.d0        ; sigma=2.00d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=226
   case('Ch')
               amu=1.d0        ; sigma=2.25d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=227
   case('Ci')
               amu=1.d0        ; sigma=2.50d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=228
   case('Cj')
               amu=1.d0        ; sigma=2.75d0   ; eps=1.00d0  ; rcov=2.d0**(1.d0/6.d0)*sigma/Bohr_Ang  ; znucl=229
   case default
   stop "Unknown atom!!!"
 end select
end subroutine mlj_symbol2znucl

