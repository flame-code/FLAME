subroutine sg_ops(nsym_tot,rsym_all,cryssys_all,brav_all,nsymp_all,nsym_all,ind_rsym_all, &
    LATSGP,CRYSSYS,BRAV,NSYMP,NSYM,RSYM,NSYMMAX)
implicit none
integer, intent(in):: nsym_tot
real(8), intent(in):: rsym_all(4,4,nsym_tot)
integer, intent(in):: cryssys_all(230), brav_all(230)
integer, intent(in):: nsymp_all(230), nsym_all(230), ind_rsym_all(230)
integer:: LATSGP,NSYMP,NSYM,NSYMMAX,CRYSSYS,BRAV
real(8):: RSYM(4,4,NSYMMAX)
integer:: ii
 
!CRYSSYS=1: Triclinic     (  1-  2)
!CRYSSYS=2: Monoclinic    (  3- 15)
!CRYSSYS=3: Orthorhombic  ( 16- 74)
!CRYSSYS=4: Tetragonal    ( 75-142)
!CRYSSYS=5: Trigonal      (143-176)
!CRYSSYS=6: Hexagonal     (144-194)
!CRYSSYS=7: Cubic         (195-230)
 
!BRAV=1: Base Centered A     (A)
!BRAV=2: Base Centered B     (B)
!BRAV=3: Base Centered C     (C)
!BRAV=4: Face centered       (F)
!BRAV=5: Body centered       (I)
!BRAV=6: Hexagonal           (H)
!BRAV=7: Primitive centering (P)
 
if(NSYMMAX.ne.192) stop 'NSYMMAX needs to be 192'
RSYM=0.d0
if(LATSGP<1 .or. LATSGP>230) stop 'Wrong space group index'
NSYMP=nsymp_all(LATSGP)
NSYM=nsym_all(LATSGP)
CRYSSYS=cryssys_all(LATSGP)
BRAV=brav_all(LATSGP)
do ii=ind_rsym_all(LATSGP),ind_rsym_all(LATSGP)+NSYM-1
RSYM(1:4,1:4,ii-ind_rsym_all(LATSGP)+1)=rsym_all(1:4,1:4,ii)
enddo
end subroutine sg_ops
subroutine sg_lims(LATSGP,XYZLIM)
implicit none
integer:: LATSGP
real(8):: XYZLIM(2,3)
 
!CRYSSYS=1: Triclinic     (  1-  2)
!CRYSSYS=2: Monoclinic    (  3- 15)
!CRYSSYS=3: Orthorhombic  ( 16- 74)
!CRYSSYS=4: Tetragonal    ( 75-142)
!CRYSSYS=5: Trigonal      (143-176)
!CRYSSYS=6: Hexagonal     (144-194)
!CRYSSYS=7: Cubic         (195-230)
 
!BRAV=1: Base Centered A     (A)
!BRAV=2: Base Centered B     (B)
!BRAV=3: Base Centered C     (C)
!BRAV=4: Face centered       (F)
!BRAV=5: Body centered       (I)
!BRAV=6: Hexagonal           (H)
!BRAV=7: Primitive centering (P)
 
select case(LATSGP)
  case(  1)
   call sg_lims_init_001(XYZLIM)
  case(  2)
   call sg_lims_init_002(XYZLIM)
  case(  3)
   call sg_lims_init_003(XYZLIM)
  case(  4)
   call sg_lims_init_004(XYZLIM)
  case(  5)
   call sg_lims_init_005(XYZLIM)
  case(  6)
   call sg_lims_init_006(XYZLIM)
  case(  7)
   call sg_lims_init_007(XYZLIM)
  case(  8)
   call sg_lims_init_008(XYZLIM)
  case(  9)
   call sg_lims_init_009(XYZLIM)
  case( 10)
   call sg_lims_init_010(XYZLIM)
  case( 11)
   call sg_lims_init_011(XYZLIM)
  case( 12)
   call sg_lims_init_012(XYZLIM)
  case( 13)
   call sg_lims_init_013(XYZLIM)
  case( 14)
   call sg_lims_init_014(XYZLIM)
  case( 15)
   call sg_lims_init_015(XYZLIM)
  case( 16)
   call sg_lims_init_016(XYZLIM)
  case( 17)
   call sg_lims_init_017(XYZLIM)
  case( 18)
   call sg_lims_init_018(XYZLIM)
  case( 19)
   call sg_lims_init_019(XYZLIM)
  case( 20)
   call sg_lims_init_020(XYZLIM)
  case( 21)
   call sg_lims_init_021(XYZLIM)
  case( 22)
   call sg_lims_init_022(XYZLIM)
  case( 23)
   call sg_lims_init_023(XYZLIM)
  case( 24)
   call sg_lims_init_024(XYZLIM)
  case( 25)
   call sg_lims_init_025(XYZLIM)
  case( 26)
   call sg_lims_init_026(XYZLIM)
  case( 27)
   call sg_lims_init_027(XYZLIM)
  case( 28)
   call sg_lims_init_028(XYZLIM)
  case( 29)
   call sg_lims_init_029(XYZLIM)
  case( 30)
   call sg_lims_init_030(XYZLIM)
  case( 31)
   call sg_lims_init_031(XYZLIM)
  case( 32)
   call sg_lims_init_032(XYZLIM)
  case( 33)
   call sg_lims_init_033(XYZLIM)
  case( 34)
   call sg_lims_init_034(XYZLIM)
  case( 35)
   call sg_lims_init_035(XYZLIM)
  case( 36)
   call sg_lims_init_036(XYZLIM)
  case( 37)
   call sg_lims_init_037(XYZLIM)
  case( 38)
   call sg_lims_init_038(XYZLIM)
  case( 39)
   call sg_lims_init_039(XYZLIM)
  case( 40)
   call sg_lims_init_040(XYZLIM)
  case( 41)
   call sg_lims_init_041(XYZLIM)
  case( 42)
   call sg_lims_init_042(XYZLIM)
  case( 43)
   call sg_lims_init_043(XYZLIM)
  case( 44)
   call sg_lims_init_044(XYZLIM)
  case( 45)
   call sg_lims_init_045(XYZLIM)
  case( 46)
   call sg_lims_init_046(XYZLIM)
  case( 47)
   call sg_lims_init_047(XYZLIM)
  case( 48)
   call sg_lims_init_048(XYZLIM)
  case( 49)
   call sg_lims_init_049(XYZLIM)
  case( 50)
   call sg_lims_init_050(XYZLIM)
  case( 51)
   call sg_lims_init_051(XYZLIM)
  case( 52)
   call sg_lims_init_052(XYZLIM)
  case( 53)
   call sg_lims_init_053(XYZLIM)
  case( 54)
   call sg_lims_init_054(XYZLIM)
  case( 55)
   call sg_lims_init_055(XYZLIM)
  case( 56)
   call sg_lims_init_056(XYZLIM)
  case( 57)
   call sg_lims_init_057(XYZLIM)
  case( 58)
   call sg_lims_init_058(XYZLIM)
  case( 59)
   call sg_lims_init_059(XYZLIM)
  case( 60)
   call sg_lims_init_060(XYZLIM)
  case( 61)
   call sg_lims_init_061(XYZLIM)
  case( 62)
   call sg_lims_init_062(XYZLIM)
  case( 63)
   call sg_lims_init_063(XYZLIM)
  case( 64)
   call sg_lims_init_064(XYZLIM)
  case( 65)
   call sg_lims_init_065(XYZLIM)
  case( 66)
   call sg_lims_init_066(XYZLIM)
  case( 67)
   call sg_lims_init_067(XYZLIM)
  case( 68)
   call sg_lims_init_068(XYZLIM)
  case( 69)
   call sg_lims_init_069(XYZLIM)
  case( 70)
   call sg_lims_init_070(XYZLIM)
  case( 71)
   call sg_lims_init_071(XYZLIM)
  case( 72)
   call sg_lims_init_072(XYZLIM)
  case( 73)
   call sg_lims_init_073(XYZLIM)
  case( 74)
   call sg_lims_init_074(XYZLIM)
  case( 75)
   call sg_lims_init_075(XYZLIM)
  case( 76)
   call sg_lims_init_076(XYZLIM)
  case( 77)
   call sg_lims_init_077(XYZLIM)
  case( 78)
   call sg_lims_init_078(XYZLIM)
  case( 79)
   call sg_lims_init_079(XYZLIM)
  case( 80)
   call sg_lims_init_080(XYZLIM)
  case( 81)
   call sg_lims_init_081(XYZLIM)
  case( 82)
   call sg_lims_init_082(XYZLIM)
  case( 83)
   call sg_lims_init_083(XYZLIM)
  case( 84)
   call sg_lims_init_084(XYZLIM)
  case( 85)
   call sg_lims_init_085(XYZLIM)
  case( 86)
   call sg_lims_init_086(XYZLIM)
  case( 87)
   call sg_lims_init_087(XYZLIM)
  case( 88)
   call sg_lims_init_088(XYZLIM)
  case( 89)
   call sg_lims_init_089(XYZLIM)
  case( 90)
   call sg_lims_init_090(XYZLIM)
  case( 91)
   call sg_lims_init_091(XYZLIM)
  case( 92)
   call sg_lims_init_092(XYZLIM)
  case( 93)
   call sg_lims_init_093(XYZLIM)
  case( 94)
   call sg_lims_init_094(XYZLIM)
  case( 95)
   call sg_lims_init_095(XYZLIM)
  case( 96)
   call sg_lims_init_096(XYZLIM)
  case( 97)
   call sg_lims_init_097(XYZLIM)
  case( 98)
   call sg_lims_init_098(XYZLIM)
  case( 99)
   call sg_lims_init_099(XYZLIM)
  case(100)
   call sg_lims_init_100(XYZLIM)
  case(101)
   call sg_lims_init_101(XYZLIM)
  case(102)
   call sg_lims_init_102(XYZLIM)
  case(103)
   call sg_lims_init_103(XYZLIM)
  case(104)
   call sg_lims_init_104(XYZLIM)
  case(105)
   call sg_lims_init_105(XYZLIM)
  case(106)
   call sg_lims_init_106(XYZLIM)
  case(107)
   call sg_lims_init_107(XYZLIM)
  case(108)
   call sg_lims_init_108(XYZLIM)
  case(109)
   call sg_lims_init_109(XYZLIM)
  case(110)
   call sg_lims_init_110(XYZLIM)
  case(111)
   call sg_lims_init_111(XYZLIM)
  case(112)
   call sg_lims_init_112(XYZLIM)
  case(113)
   call sg_lims_init_113(XYZLIM)
  case(114)
   call sg_lims_init_114(XYZLIM)
  case(115)
   call sg_lims_init_115(XYZLIM)
  case(116)
   call sg_lims_init_116(XYZLIM)
  case(117)
   call sg_lims_init_117(XYZLIM)
  case(118)
   call sg_lims_init_118(XYZLIM)
  case(119)
   call sg_lims_init_119(XYZLIM)
  case(120)
   call sg_lims_init_120(XYZLIM)
  case(121)
   call sg_lims_init_121(XYZLIM)
  case(122)
   call sg_lims_init_122(XYZLIM)
  case(123)
   call sg_lims_init_123(XYZLIM)
  case(124)
   call sg_lims_init_124(XYZLIM)
  case(125)
   call sg_lims_init_125(XYZLIM)
  case(126)
   call sg_lims_init_126(XYZLIM)
  case(127)
   call sg_lims_init_127(XYZLIM)
  case(128)
   call sg_lims_init_128(XYZLIM)
  case(129)
   call sg_lims_init_129(XYZLIM)
  case(130)
   call sg_lims_init_130(XYZLIM)
  case(131)
   call sg_lims_init_131(XYZLIM)
  case(132)
   call sg_lims_init_132(XYZLIM)
  case(133)
   call sg_lims_init_133(XYZLIM)
  case(134)
   call sg_lims_init_134(XYZLIM)
  case(135)
   call sg_lims_init_135(XYZLIM)
  case(136)
   call sg_lims_init_136(XYZLIM)
  case(137)
   call sg_lims_init_137(XYZLIM)
  case(138)
   call sg_lims_init_138(XYZLIM)
  case(139)
   call sg_lims_init_139(XYZLIM)
  case(140)
   call sg_lims_init_140(XYZLIM)
  case(141)
   call sg_lims_init_141(XYZLIM)
  case(142)
   call sg_lims_init_142(XYZLIM)
  case(143)
   call sg_lims_init_143(XYZLIM)
  case(144)
   call sg_lims_init_144(XYZLIM)
  case(145)
   call sg_lims_init_145(XYZLIM)
  case(146)
   call sg_lims_init_146(XYZLIM)
  case(147)
   call sg_lims_init_147(XYZLIM)
  case(148)
   call sg_lims_init_148(XYZLIM)
  case(149)
   call sg_lims_init_149(XYZLIM)
  case(150)
   call sg_lims_init_150(XYZLIM)
  case(151)
   call sg_lims_init_151(XYZLIM)
  case(152)
   call sg_lims_init_152(XYZLIM)
  case(153)
   call sg_lims_init_153(XYZLIM)
  case(154)
   call sg_lims_init_154(XYZLIM)
  case(155)
   call sg_lims_init_155(XYZLIM)
  case(156)
   call sg_lims_init_156(XYZLIM)
  case(157)
   call sg_lims_init_157(XYZLIM)
  case(158)
   call sg_lims_init_158(XYZLIM)
  case(159)
   call sg_lims_init_159(XYZLIM)
  case(160)
   call sg_lims_init_160(XYZLIM)
  case(161)
   call sg_lims_init_161(XYZLIM)
  case(162)
   call sg_lims_init_162(XYZLIM)
  case(163)
   call sg_lims_init_163(XYZLIM)
  case(164)
   call sg_lims_init_164(XYZLIM)
  case(165)
   call sg_lims_init_165(XYZLIM)
  case(166)
   call sg_lims_init_166(XYZLIM)
  case(167)
   call sg_lims_init_167(XYZLIM)
  case(168)
   call sg_lims_init_168(XYZLIM)
  case(169)
   call sg_lims_init_169(XYZLIM)
  case(170)
   call sg_lims_init_170(XYZLIM)
  case(171)
   call sg_lims_init_171(XYZLIM)
  case(172)
   call sg_lims_init_172(XYZLIM)
  case(173)
   call sg_lims_init_173(XYZLIM)
  case(174)
   call sg_lims_init_174(XYZLIM)
  case(175)
   call sg_lims_init_175(XYZLIM)
  case(176)
   call sg_lims_init_176(XYZLIM)
  case(177)
   call sg_lims_init_177(XYZLIM)
  case(178)
   call sg_lims_init_178(XYZLIM)
  case(179)
   call sg_lims_init_179(XYZLIM)
  case(180)
   call sg_lims_init_180(XYZLIM)
  case(181)
   call sg_lims_init_181(XYZLIM)
  case(182)
   call sg_lims_init_182(XYZLIM)
  case(183)
   call sg_lims_init_183(XYZLIM)
  case(184)
   call sg_lims_init_184(XYZLIM)
  case(185)
   call sg_lims_init_185(XYZLIM)
  case(186)
   call sg_lims_init_186(XYZLIM)
  case(187)
   call sg_lims_init_187(XYZLIM)
  case(188)
   call sg_lims_init_188(XYZLIM)
  case(189)
   call sg_lims_init_189(XYZLIM)
  case(190)
   call sg_lims_init_190(XYZLIM)
  case(191)
   call sg_lims_init_191(XYZLIM)
  case(192)
   call sg_lims_init_192(XYZLIM)
  case(193)
   call sg_lims_init_193(XYZLIM)
  case(194)
   call sg_lims_init_194(XYZLIM)
  case(195)
   call sg_lims_init_195(XYZLIM)
  case(196)
   call sg_lims_init_196(XYZLIM)
  case(197)
   call sg_lims_init_197(XYZLIM)
  case(198)
   call sg_lims_init_198(XYZLIM)
  case(199)
   call sg_lims_init_199(XYZLIM)
  case(200)
   call sg_lims_init_200(XYZLIM)
  case(201)
   call sg_lims_init_201(XYZLIM)
  case(202)
   call sg_lims_init_202(XYZLIM)
  case(203)
   call sg_lims_init_203(XYZLIM)
  case(204)
   call sg_lims_init_204(XYZLIM)
  case(205)
   call sg_lims_init_205(XYZLIM)
  case(206)
   call sg_lims_init_206(XYZLIM)
  case(207)
   call sg_lims_init_207(XYZLIM)
  case(208)
   call sg_lims_init_208(XYZLIM)
  case(209)
   call sg_lims_init_209(XYZLIM)
  case(210)
   call sg_lims_init_210(XYZLIM)
  case(211)
   call sg_lims_init_211(XYZLIM)
  case(212)
   call sg_lims_init_212(XYZLIM)
  case(213)
   call sg_lims_init_213(XYZLIM)
  case(214)
   call sg_lims_init_214(XYZLIM)
  case(215)
   call sg_lims_init_215(XYZLIM)
  case(216)
   call sg_lims_init_216(XYZLIM)
  case(217)
   call sg_lims_init_217(XYZLIM)
  case(218)
   call sg_lims_init_218(XYZLIM)
  case(219)
   call sg_lims_init_219(XYZLIM)
  case(220)
   call sg_lims_init_220(XYZLIM)
  case(221)
   call sg_lims_init_221(XYZLIM)
  case(222)
   call sg_lims_init_222(XYZLIM)
  case(223)
   call sg_lims_init_223(XYZLIM)
  case(224)
   call sg_lims_init_224(XYZLIM)
  case(225)
   call sg_lims_init_225(XYZLIM)
  case(226)
   call sg_lims_init_226(XYZLIM)
  case(227)
   call sg_lims_init_227(XYZLIM)
  case(228)
   call sg_lims_init_228(XYZLIM)
  case(229)
   call sg_lims_init_229(XYZLIM)
  case(230)
   call sg_lims_init_230(XYZLIM)
case default
stop 'Wrong space group index'
end select
contains
subroutine sg_lims_init_001(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_002(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_003(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_004(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.499990d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_005(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.499990d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_006(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_007(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_008(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_009(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_010(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_011(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_012(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_013(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_014(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_015(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_016(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_017(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_018(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_019(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_020(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_021(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_022(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_023(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  1.000010d0/)
end subroutine
subroutine sg_lims_init_024(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_025(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_026(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_027(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_028(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.999990d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_029(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.500010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_030(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_031(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_032(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_033(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_034(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_035(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_036(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_037(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_038(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_039(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_040(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.499990d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_041(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_042(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_043(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.499990d0,  0.125010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_044(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_045(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_046(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.999990d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_047(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_048(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_049(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_050(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_051(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_052(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.499990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_053(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_054(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_055(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_056(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_057(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.499990d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_058(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_059(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_060(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_061(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.499990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_062(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.499990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_063(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.499990d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_064(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_065(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_066(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_067(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_068(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_069(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_070(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.125010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_071(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_072(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_073(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_074(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_075(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_076(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.249990d0/)
end subroutine
subroutine sg_lims_init_077(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_078(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.249990d0/)
end subroutine
subroutine sg_lims_init_079(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_080(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.249990d0/)
end subroutine
subroutine sg_lims_init_081(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_082(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_083(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_084(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_085(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_086(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_087(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_088(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_089(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_090(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_091(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_092(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_093(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_094(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_095(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_096(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_097(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_098(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_099(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_100(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_101(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_102(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_103(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_104(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_105(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_106(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_107(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_108(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_109(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_110(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_111(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_112(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_113(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_114(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_115(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_116(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_117(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_118(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_119(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_120(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_121(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_122(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_123(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_124(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_125(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_126(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_127(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_128(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_129(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_130(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_131(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_132(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_133(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_134(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_135(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_136(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_137(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_138(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_139(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_140(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_141(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_142(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.250010d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_143(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_144(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_145(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_146(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_147(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_148(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_149(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_150(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_151(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_152(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_153(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_154(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_155(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_156(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_157(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_158(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_159(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_160(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.333343d0,  0.333323d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_161(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.333343d0,  0.333323d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_162(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_163(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_164(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  1.000010d0/)
end subroutine
subroutine sg_lims_init_165(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_166(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_167(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.166677d0,  0.333343d0/)
end subroutine
subroutine sg_lims_init_168(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_169(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166657d0/)
end subroutine
subroutine sg_lims_init_170(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166657d0/)
end subroutine
subroutine sg_lims_init_171(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_172(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_173(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_174(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_175(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_176(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_177(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_178(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_179(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_180(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_181(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.166677d0/)
end subroutine
subroutine sg_lims_init_182(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_183(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_184(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_185(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_186(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_187(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_188(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_189(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_190(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_191(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_192(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_193(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.333343d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_194(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.666677d0,  0.666677d0,  0.333323d0/)
end subroutine
subroutine sg_lims_init_195(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_196(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_197(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_198(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_199(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_200(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_201(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_202(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_203(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.125010d0,  0.125010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_204(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_205(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.499990d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_206(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_207(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_208(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_209(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_210(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_211(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_212(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_213(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.999990d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_214(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.999990d0,  0.125010d0/)
end subroutine
subroutine sg_lims_init_215(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_216(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_217(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_218(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_219(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_220(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_221(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_222(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_223(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.500010d0/)
end subroutine
subroutine sg_lims_init_224(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_225(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.250010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_226(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.250010d0,  0.250010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_227(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.125010d0,  0.125010d0,  0.999990d0/)
end subroutine
subroutine sg_lims_init_228(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.125010d0,  0.125010d0,  0.499990d0/)
end subroutine
subroutine sg_lims_init_229(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.500010d0,  0.500010d0,  0.250010d0/)
end subroutine
subroutine sg_lims_init_230(XYZLIM)
implicit none
real(8):: XYZLIM(2,3)
   XYZLIM(1,:) =(/ 0.000000d0,  0.000000d0,  0.000000d0/)
   XYZLIM(2,:) =(/ 0.125010d0,  0.125010d0,  0.999990d0/)
end subroutine
end subroutine sg_lims
