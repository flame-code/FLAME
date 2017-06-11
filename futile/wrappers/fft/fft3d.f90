
!>  @file 
!!  @brief 3-dimensional complex-complex FFT routines
!!
!!   When compared to the best vendor implementations on RISC architectures 
!!   it gives close to optimal performance (perhaps loosing 20 percent in speed)
!!   and it is significanly faster than many not so good vendor implementations 
!!   as well as other portable FFT's. 
!!   On all vector machines tested so far (Cray, NEC, Fujitsu) is 
!!   was significantly faster than the vendor routines
!!   The theoretical background is described in :\n
!!   1) S. Goedecker: Rotating a three-dimensional array in optimal
!!   positions for vector processing: Case study for a three-dimensional Fast
!!   Fourier Transform, Comp. Phys. Commun. 76, 294 (1993)
!!   Citing of this reference is greatly appreciated if the routines are used 
!!   for scientific work.
!!
!! @author
!!   Copyright (C) Stefan Goedecker, Lausanne, Switzerland, August 1, 1991
!!   Copyright (C) Stefan Goedecker, Cornell University, Ithaca, USA, 1994
!!   Copyright (C) Stefan Goedecker, MPI Stuttgart, Germany, 1999
!!   Copyright (C) Stefan Goedecker, CEA Grenoble, 2002
!!   Copyright (C) 2002-2010 BigDFT group 
!!   The part for radix 7 added by Alexey Neelov, Basel University,  2008
!!   This file is distributed under the terms of the
!!   GNU General Public License, see ~/COPYING file
!!   or http://www.gnu.org/copyleft/gpl.txt .
!!   For the list of contributors, see ~/AUTHORS 
!!
!! @warning
!!   Presumably good compiler flags:
!!   IBM, serial power 2: xlf -qarch=pwr2 -O2 -qmaxmem=-1
!!   with OpenMP: IBM: xlf_r -qfree -O4 -qarch=pwr3 -qtune=pwr3 -qsmp=omp -qmaxmem=-1 ; 
!!                     a.out
!!   DEC: f90 -O3 -arch ev67 -pipeline
!!   with OpenMP: DEC: f90 -O3 -arch ev67 -pipeline -omp -lelan ; 
!!                     prun -N1 -c4 a.out
!!
!! PERFORMANCE AND THE NCACHE
!!  The most important feature for performance is the right choice of 
!!  the parameter ncache. On a vector machine ncache has to be put to 0.
!!  On a RISC machine with cache, it is very important to find the optimal 
!!  value of NCACHE. NCACHE determines the size of the work array zw, that
!!  has to fit into cache. It has therefore to be chosen to equal roughly 
!!  half the size of the physical cache in units of real(kind=8) numbers.
!!  If the machine has 2 cache levels it can not be predicted which 
!!  cache level will be the most relevant one for choosing ncache. 
!!  The optimal value of ncache can easily be determined by numerical 
!!  experimentation. A too large value of ncache leads to a dramatic 
!!  and sudden decrease of performance, a too small value to a to a 
!!  slow and less dramatic decrease of performance. If NCACHE is set 
!!  to a value so small, that not even a single one dimensional transform 
!!  can be done in the workarray zw, the program stops with an error 
!!  message.
!!

!> Module which contains parameters for FFT3D
module module_fft_sg

implicit none

!< Maximum number of points for FFT (should be same number in fft3d routine)
integer, parameter :: nfft_max=2097152
!< Number of factors in the decomposition
integer, parameter :: n_factors = 7
integer :: i_d,j_d

integer, parameter :: MODE_IO=-1000
integer, parameter :: MODE_IW=-1001
integer, parameter :: MODE_WO=-1002
integer, parameter :: MODE_WW=-1003

!! @warning
!!   some reasonable values of ncache: 
!!   IBM/RS6000/590: 16*1024 ; IBM/RS6000/390: 3*1024 ; 
!!   IBM/PwPC: 1*1024 ; SGI/MIPS/R8000: 16*1024 ; DEC/Alpha/EV5 and EV6 6*1024
!!   But if you care about performance find the optimal value of ncache yourself!
!!       On all vector machines: ncache=0

integer :: ncache =  8*1024  !< To have all available dimensions, ncache should be a multiple of 4*nfft_max (T.D.)
!integer, parameter :: ncache = (4*nfft_max)
!Vectorial computer
!integer, parameter :: ncache = 0


integer, parameter :: ndata = 304  !< Multiple of 3,5,4,6,7,8 (and 2) with certain restrictions
integer, dimension(ndata), parameter :: i_data = (/   &
3,      4,      5,      6,      7,      8,      9,     12,  &
14,     15,     16,     18,     20,     21,     24,     25,  &
27,     28,     30,     32,     35,     36,     40,     42,  &
45,     48,     50,     54,     56,     60,     63,     64,  &
70,     72,     75,     80,     81,     84,     90,     96,  &
100,    105,    108,    112,    120,    125,    126,    128,  &
135,    140,    144,    150,    160,    162,    168,    175,  &
180,    189,    192,    200,    210,    216,    224,    225,  &
240,    243,    252,    256,    270,    280,    288,    300,  &
315,    320,    324,    336,    360,    375,    378,    384,  &
400,    405,    420,    432,    448,    450,    480,    486,  &
500,    504,    512,    525,    540,    560,    567,    576,  &
600,    625,    630,    640,    648,    672,    675,    700,  &
720,    729,    750,    756,    768,    800,    810,    840,  &
864,    875,    896,    900,    945,    960,    972,   1000,  &
1008,   1024,   1050,   1080,   1120,   1125,   1134,   1152,  &
1200,   1215,   1260,   1280,   1296,   1344,   1350,   1400,  &
1440,   1458,   1500,   1512,   1536,   1575,   1600,   1620,  &
1680,   1701,   1728,   1750,   1792,   1800,   1875,   1890,   1920,  &
1944,   2000,   2016,   2025,   2048,   2100,   2160,   2240,  &
2250,   2268,   2304,   2400,   2430,   2500,   2520,   2560,  &
2592,   2625,   2688,   2700,   2800,   2835,   2880,   3000,  &
3024,   3072,   3125,   3150,   3200,   3240,   3360,   3375,  &
3402,   3456,   3500,   3584,   3600,   3750,   3780,   3840,  &
3888,   4000,   4032,   4050,   4096,   4200,   4320,   4375,  &
4480,   4500,   4536,   4608,   4725,   4800,   5000,   5040,  &
5120,   5184,   5250,   5376,   5400,   5600,   5625,   5670,  &
5760,   6000,   6048,   6144,   6250,   6300,   6400,   6480,   6720,  &
6750,   6912,   7000,   7168,   7200,   7500,   7560,   7680,  &
7875,   8000,   8064,   8192,   8400,   8640,   8750,   8960,   9000,  &
9072,   9216,   9375,   9450,   9600,  10000,  10080,  10240,  &
10368,  10500,  10752,  10800,  11200,  11250,  11520,  12000,  &
12096,  12288,  12500,  12600,  12800,  13125,  13440,  13824,  &
14000,  14336,  14400,  15000,  15120,  15360,  15625,  15750,  &
16000,  16128,  16384,  16800,  17280,  17500,  17920,  18000,  &
18432,  18750,  19200,  20000,  20160,  20480,  20736,  21000,  21504,  &
21875,  22400,  22500,  23040,  24000,  32768,  65536,   131072,  262144, &
524288,  1048576, 2097152  /)
! The factors 6 and 7 are only allowed in the first place!
integer, dimension(1+n_factors,ndata) :: ij_data
data ((ij_data(i_d,j_d),i_d=1,1+n_factors),j_d=1,ndata) /  &
3,   3,   1,   1,   1,   1,   1,   1,        4,   4,   1,   1,   1,   1,   1,   1,  &
5,   5,   1,   1,   1,   1,   1,   1,        6,   6,   1,   1,   1,   1,   1,   1,  &
7,   7,   1,   1,   1,   1,   1,   1,        8,   8,   1,   1,   1,   1,   1,   1,  &
9,   3,   3,   1,   1,   1,   1,   1,       12,   4,   3,   1,   1,   1,   1,   1,  &
14,   7,   2,   1,   1,   1,   1,   1,       15,   5,   3,   1,   1,   1,   1,   1,  &
16,   4,   4,   1,   1,   1,   1,   1,       18,   6,   3,   1,   1,   1,   1,   1,  &
20,   5,   4,   1,   1,   1,   1,   1,       21,   7,   3,   1,   1,   1,   1,   1,  &
24,   8,   3,   1,   1,   1,   1,   1,       25,   5,   5,   1,   1,   1,   1,   1,  &
27,   3,   3,   3,   1,   1,   1,   1,       28,   7,   4,   1,   1,   1,   1,   1,  &
30,   6,   5,   1,   1,   1,   1,   1,       32,   8,   4,   1,   1,   1,   1,   1,  &
35,   7,   5,   1,   1,   1,   1,   1,       36,   4,   3,   3,   1,   1,   1,   1,  &
40,   8,   5,   1,   1,   1,   1,   1,       42,   7,   3,   2,   1,   1,   1,   1,  &
45,   5,   3,   3,   1,   1,   1,   1,       48,   4,   4,   3,   1,   1,   1,   1,  &
50,   5,   5,   2,   1,   1,   1,   1,       54,   6,   3,   3,   1,   1,   1,   1,  &
56,   7,   8,   1,   1,   1,   1,   1,       60,   5,   4,   3,   1,   1,   1,   1,  &
63,   7,   3,   3,   1,   1,   1,   1,       64,   8,   8,   1,   1,   1,   1,   1,  &
70,   7,   5,   2,   1,   1,   1,   1,       72,   8,   3,   3,   1,   1,   1,   1,  &
75,   5,   5,   3,   1,   1,   1,   1,       80,   5,   4,   4,   1,   1,   1,   1,  &
81,   3,   3,   3,   3,   1,   1,   1,       84,   7,   4,   3,   1,   1,   1,   1,  &
90,   6,   5,   3,   1,   1,   1,   1,       96,   8,   4,   3,   1,   1,   1,   1,  &
100,   5,   5,   4,   1,   1,   1,   1,      105,   7,   5,   3,   1,   1,   1,   1,  &
108,   4,   3,   3,   3,   1,   1,   1,      112,   7,   4,   4,   1,   1,   1,   1,  &
120,   8,   5,   3,   1,   1,   1,   1,      125,   5,   5,   5,   1,   1,   1,   1,  &
126,   7,   3,   3,   2,   1,   1,   1,      128,   8,   4,   4,   1,   1,   1,   1,  &
135,   5,   3,   3,   3,   1,   1,   1,      140,   7,   5,   4,   1,   1,   1,   1,  &
144,   6,   8,   3,   1,   1,   1,   1,      150,   6,   5,   5,   1,   1,   1,   1,  &
160,   8,   5,   4,   1,   1,   1,   1,      162,   6,   3,   3,   3,   1,   1,   1,  &
168,   7,   8,   3,   1,   1,   1,   1,      175,   7,   5,   5,   1,   1,   1,   1,  &
180,   5,   4,   3,   3,   1,   1,   1,      189,   7,   3,   3,   3,   1,   1,   1,  &
192,   6,   8,   4,   1,   1,   1,   1,      200,   8,   5,   5,   1,   1,   1,   1,  &
210,   7,   3,   5,   2,   1,   1,   1,      216,   8,   3,   3,   3,   1,   1,   1,  &
224,   7,   8,   4,   1,   1,   1,   1,      225,   5,   5,   3,   3,   1,   1,   1,  &
240,   6,   8,   5,   1,   1,   1,   1,      243,   3,   3,   3,   3,   3,   1,   1,  &
252,   7,   4,   3,   3,   1,   1,   1,      256,   8,   8,   4,   1,   1,   1,   1,  &
270,   6,   5,   3,   3,   1,   1,   1,      280,   7,   8,   5,   1,   1,   1,   1,  &
288,   8,   4,   3,   3,   1,   1,   1,      300,   5,   5,   4,   3,   1,   1,   1,  &
315,   7,   5,   3,   3,   1,   1,   1,      320,   5,   4,   4,   4,   1,   1,   1,  &
324,   4,   3,   3,   3,   3,   1,   1,      336,   7,   4,   4,   3,   1,   1,   1,  &
360,   8,   5,   3,   3,   1,   1,   1,      375,   5,   5,   5,   3,   1,   1,   1,  &
378,   7,   3,   3,   3,   2,   1,   1,      384,   8,   4,   4,   3,   1,   1,   1,  &
400,   5,   5,   4,   4,   1,   1,   1,      405,   5,   3,   3,   3,   3,   1,   1,  &
420,   7,   5,   4,   3,   1,   1,   1,      432,   4,   4,   3,   3,   3,   1,   1,  &
448,   7,   8,   8,   1,   1,   1,   1,      450,   6,   5,   5,   3,   1,   1,   1,  &
480,   8,   5,   4,   3,   1,   1,   1,      486,   6,   3,   3,   3,   3,   1,   1,  &
500,   5,   5,   5,   4,   1,   1,   1,      504,   7,   8,   3,   3,   1,   1,   1,  &
512,   8,   8,   8,   1,   1,   1,   1,      525,   7,   5,   5,   3,   1,   1,   1,  &
540,   5,   4,   3,   3,   3,   1,   1,      560,   7,   5,   4,   4,   1,   1,   1,  &
567,   7,   3,   3,   3,   3,   1,   1,      576,   4,   4,   4,   3,   3,   1,   1,  &
600,   8,   5,   5,   3,   1,   1,   1,      625,   5,   5,   5,   5,   1,   1,   1,  &
630,   7,   5,   3,   3,   2,   1,   1,      640,   8,   5,   4,   4,   1,   1,   1,  &
648,   8,   3,   3,   3,   3,   1,   1,      672,   7,   8,   4,   3,   1,   1,   1,  &
675,   5,   5,   3,   3,   3,   1,   1,      700,   7,   5,   5,   4,   1,   1,   1,  &
720,   5,   4,   4,   3,   3,   1,   1,      729,   3,   3,   3,   3,   3,   3,   1,  &
750,   6,   5,   5,   5,   1,   1,   1,      756,   7,   4,   3,   3,   3,   1,   1,  &
768,   4,   4,   4,   4,   3,   1,   1,      800,   8,   5,   5,   4,   1,   1,   1,  &
810,   6,   5,   3,   3,   3,   1,   1,      840,   7,   8,   5,   3,   1,   1,   1,  &
864,   8,   4,   3,   3,   3,   1,   1,      875,   7,   5,   5,   5,   1,   1,   1,  &
896,   7,   8,   4,   4,   1,   1,   1,      900,   5,   5,   4,   3,   3,   1,   1,  &
945,   7,   5,   3,   3,   3,   1,   1,      960,   5,   4,   4,   4,   3,   1,   1,  &
972,   4,   3,   3,   3,   3,   3,   1,     1000,   8,   5,   5,   5,   1,   1,   1,  &
1008,   7,   8,   3,   3,   2,   1,   1,     1024,   4,   4,   4,   4,   4,   1,   1,  &
1050,   7,   5,   5,   3,   2,   1,   1,     1080,   6,   5,   4,   3,   3,   1,   1,  &
1120,   7,   8,   5,   4,   1,   1,   1,     1125,   5,   5,   5,   3,   3,   1,   1,  &
1134,   7,   3,   3,   3,   3,   2,   1,     1152,   6,   4,   4,   4,   3,   1,   1,  &
1200,   6,   8,   5,   5,   1,   1,   1,     1215,   5,   3,   3,   3,   3,   3,   1,  &
1260,   7,   5,   4,   3,   3,   1,   1,     1280,   8,   8,   5,   4,   1,   1,   1,  &
1296,   6,   8,   3,   3,   3,   1,   1,     1344,   7,   8,   4,   3,   2,   1,   1,  &
1350,   6,   5,   5,   3,   3,   1,   1,     1400,   7,   8,   5,   5,   1,   1,   1,  &
1440,   6,   5,   4,   4,   3,   1,   1,     1458,   6,   3,   3,   3,   3,   3,   1,  &
1500,   5,   5,   5,   4,   3,   1,   1,     1512,   7,   8,   3,   3,   3,   1,   1,  &
1536,   6,   8,   8,   4,   1,   1,   1,     1575,   7,   5,   5,   3,   3,   1,   1,  &
1600,   8,   8,   5,   5,   1,   1,   1,     1620,   5,   4,   3,   3,   3,   3,   1,  &
1680,   7,   8,   5,   3,   2,   1,   1,     1701,   7,   3,   3,   3,   3,   3,   1,  &
1728,   6,   8,   4,   3,   3,   1,   1,     1750,   7,   5,   5,   5,   2,   1,   1,  &
1792,   7,   8,   8,   4,   1,   1,   1,  &
1800,   6,   5,   5,   4,   3,   1,   1,     1875,   5,   5,   5,   5,   3,   1,   1,  &
1890,   7,   5,   3,   3,   3,   2,   1,     1920,   6,   5,   4,   4,   4,   1,   1,  &
1944,   6,   4,   3,   3,   3,   3,   1,     2000,   5,   5,   5,   4,   4,   1,   1,  &
2016,   7,   8,   4,   3,   3,   1,   1,     2025,   5,   5,   3,   3,   3,   3,   1,  &
2048,   8,   4,   4,   4,   4,   1,   1,     2100,   7,   5,   5,   4,   3,   1,   1,  &
2160,   6,   8,   5,   3,   3,   1,   1,     2240,   7,   5,   4,   4,   4,   1,   1,  &
2250,   6,   5,   5,   5,   3,   1,   1,     2268,   7,   4,   3,   3,   3,   3,   1,  &
2304,   6,   8,   4,   4,   3,   1,   1,     2400,   6,   5,   5,   4,   4,   1,   1,  &
2430,   6,   5,   3,   3,   3,   3,   1,     2500,   5,   5,   5,   5,   4,   1,   1,  &
2520,   7,   8,   5,   3,   3,   1,   1,     2560,   8,   5,   4,   4,   4,   1,   1,  &
2592,   6,   4,   4,   3,   3,   3,   1,     2625,   7,   5,   5,   5,   3,   1,   1,  &
2688,   7,   8,   4,   4,   3,   1,   1,     2700,   5,   5,   4,   3,   3,   3,   1,  &
2800,   7,   5,   5,   4,   4,   1,   1,     2835,   7,   5,   3,   3,   3,   3,   1,  &
2880,   6,   8,   5,   4,   3,   1,   1,     3000,   6,   5,   5,   5,   4,   1,   1,  &
3024,   7,   4,   4,   3,   3,   3,   1,     3072,   6,   8,   4,   4,   4,   1,   1,  &
3125,   5,   5,   5,   5,   5,   1,   1,     3150,   7,   5,   5,   3,   3,   2,   1,  &
3200,   8,   5,   5,   4,   4,   1,   1,     3240,   6,   5,   4,   3,   3,   3,   1,  &
3360,   7,   8,   5,   4,   3,   1,   1,     3375,   5,   5,   5,   3,   3,   3,   1,  &
3402,   7,   3,   3,   3,   3,   3,   2,     3456,   6,   4,   4,   4,   3,   3,   1,  &
3500,   7,   5,   5,   5,   4,   1,   1,     3584,   7,   8,   8,   8,   1,   1,   1,  &
3600,   6,   8,   5,   5,   3,   1,   1,     3750,   6,   5,   5,   5,   5,   1,   1,  &
3780,   7,   5,   4,   3,   3,   3,   1,     3840,   6,   8,   5,   4,   4,   1,   1,  &
3888,   6,   8,   3,   3,   3,   3,   1,     4000,   8,   5,   5,   5,   4,   1,   1,  &
4032,   7,   4,   4,   4,   3,   3,   1,     4050,   6,   5,   5,   3,   3,   3,   1,  &
4096,   8,   8,   4,   4,   4,   1,   1,     4200,   7,   8,   5,   5,   3,   1,   1,  &
4320,   6,   5,   4,   4,   3,   3,   1,     4375,   7,   5,   5,   5,   5,   1,   1,  &
4480,   7,   8,   5,   4,   4,   1,   1,     4500,   5,   5,   5,   4,   3,   3,   1,  &
4536,   7,   8,   3,   3,   3,   3,   1,     4608,   6,   8,   8,   4,   3,   1,   1,  &
4725,   7,   5,   5,   3,   3,   3,   1,     4800,   6,   8,   5,   5,   4,   1,   1,  &
5000,   8,   5,   5,   5,   5,   1,   1,     5040,   7,   5,   4,   4,   3,   3,   1,  &
5120,   8,   8,   5,   4,   4,   1,   1,     5184,   6,   8,   4,   3,   3,   3,   1,  &
5250,   7,   5,   5,   5,   3,   2,   1,     5376,   7,   4,   4,   4,   4,   3,   1,  &
5400,   6,   5,   5,   4,   3,   3,   1,     5600,   7,   8,   5,   5,   4,   1,   1,  &
5625,   5,   5,   5,   5,   3,   3,   1,     5670,   7,   5,   3,   3,   3,   3,   2,  &
5760,   6,   8,   8,   5,   3,   1,   1,     6000,   6,   8,   5,   5,   5,   1,   1,  &
6048,   7,   8,   4,   3,   3,   3,   1,     6144,   6,   8,   8,   4,   4,   1,   1,  &
6250,   5,   5,   5,   5,   5,   2,   1,   &
6300,   7,   5,   5,   4,   3,   3,   1,     6400,   8,   8,   5,   5,   4,   1,   1,  &
6480,   6,   8,   5,   3,   3,   3,   1,     6720,   7,   5,   4,   4,   4,   3,   1,  &
6750,   6,   5,   5,   5,   3,   3,   1,     6912,   6,   8,   4,   4,   3,   3,   1,  &
7000,   7,   8,   5,   5,   5,   1,   1,     7168,   7,   4,   4,   4,   4,   4,   1,  &
7200,   6,   5,   5,   4,   4,   3,   1,     7500,   5,   5,   5,   5,   4,   3,   1,  &
7560,   7,   8,   5,   3,   3,   3,   1,     7680,   6,   8,   8,   5,   4,   1,   1,  &
7875,   7,   5,   5,   5,   3,   3,   1,     8000,   8,   8,   5,   5,   5,   1,   1,  &
8064,   7,   8,   4,   4,   3,   3,   1,     8192,   8,   8,   8,   4,   4,   1,   1,  &
8400,   7,   8,   5,   5,   3,   2,   1,     8640,   8,   8,   5,   3,   3,   3,   1,  &
8750,   7,   5,   5,   5,   5,   2,   1,   &
8960,   7,   8,   8,   5,   4,   1,   1,     9000,   8,   5,   5,   5,   3,   3,   1,  &
9072,   7,   8,   3,   3,   3,   3,   2,     9216,   6,   8,   8,   8,   3,   1,   1,  &
9375,   5,   5,   5,   5,   5,   3,   1,     9450,   7,   5,   5,   3,   3,   3,   2,  &
9600,   8,   5,   5,   4,   4,   3,   1,    10000,   5,   5,   5,   5,   4,   4,   1,  &
10080,   7,   8,   5,   4,   3,   3,   1,    10240,   8,   8,   8,   5,   4,   1,   1,  &
10368,   6,   8,   8,   3,   3,   3,   1,    10500,   7,   5,   5,   5,   4,   3,   1,  &
10752,   7,   8,   8,   8,   3,   1,   1,    10800,   6,   8,   5,   5,   3,   3,   1,  &
11200,   7,   8,   8,   5,   5,   1,   1,    11250,   6,   5,   5,   5,   5,   3,   1,  &
11520,   8,   8,   5,   4,   3,   3,   1,    12000,   8,   5,   5,   5,   4,   3,   1,  &
12096,   7,   8,   8,   3,   3,   3,   1,    12288,   8,   8,   8,   8,   3,   1,   1,  &
12500,   5,   5,   5,   5,   5,   4,   1,    12600,   7,   8,   5,   5,   3,   3,   1,  &
12800,   8,   8,   8,   5,   5,   1,   1,    13125,   7,   5,   5,   5,   5,   3,   1,  &
13440,   7,   8,   5,   4,   4,   3,   1,    13824,   8,   8,   8,   3,   3,   3,   1,  &
14000,   7,   5,   5,   5,   4,   4,   1,    14336,   7,   8,   4,   4,   4,   4,   1,  &
14400,   8,   8,   5,   5,   3,   3,   1,    15000,   8,   5,   5,   5,   5,   3,   1,  &
15120,   7,   8,   5,   3,   3,   3,   2,    15360,   6,   8,   8,   8,   5,   1,   1,  &
15625,   5,   5,   5,   5,   5,   5,   1,    15750,   7,   5,   5,   5,   3,   3,   2,  &
16000,   8,   5,   5,   5,   4,   4,   1,    16128,   7,   8,   8,   4,   3,   3,   1,  &
16384,   8,   8,   8,   8,   4,   1,   1,    16800,   7,   8,   5,   5,   4,   3,   1,  &
17280,   6,   8,   8,   5,   3,   3,   1,    17500,   7,   5,   5,   5,   5,   4,   1,  &
17920,   7,   8,   5,   4,   4,   4,   1,    18000,   6,   8,   5,   5,   5,   3,   1,  &
18432,   8,   8,   8,   4,   3,   3,   1,    18750,   6,   5,   5,   5,   5,   5,   1,  &
19200,   8,   8,   5,   5,   4,   3,   1,    20000,   8,   5,   5,   5,   5,   4,   1,  &
20160,   7,   8,   8,   5,   3,   3,   1,    20480,   8,   8,   8,   8,   5,   1,   1,  &
20736,   8,   8,   3,   3,   3,   3,   4, &
21000,   7,   8,   5,   5,   5,   3,   1,    21504,   7,   8,   8,   3,   4,   4,   1,  &
21875,   7,   5,   5,   5,   5,   5,   1,    22400,   7,   8,   5,   5,   4,   4,   1,  &
22500,   6,   5,   5,   5,   5,   3,   2, &
23040,   8,   8,   8,   5,   3,   3,   1,    24000,   8,   8,   5,   5,   5,   3,   1,  &
32768,   8,   8,   8,   8,   8,   1,   1,    65536,   8,   8,   8,   8,   8,   2,   1,  &
131072,  8,   8,   8,   8,   8,   4,   1,    262144,  8,   8,   8,   8,   8,   8,   1,  &
524288,  8,   8,   8,   8,   8,   8,   2,    1048576, 8,   8,   8,   8,   8,   8,   4,  &
2097152, 8,   8,   8,   8,   8,   8,   8    /  

contains
  !>impulse coordinate, from 0,...,n/2+1,-n/2+1,...-1
  pure function p_index(i,n) result(p)
    implicit none
    integer, intent(in) :: i,n
    integer :: p
    p=i-(i/(n/2+2))*n-1
  end function p_index

  !>real space coordinate, from 0,...,n-1
  pure function i_index(p,n) result(j)
    implicit none
    integer, intent(in) :: p,n
    integer :: j
    j=p-((p+n)/n-1)*n
  end function i_index

end module module_fft_sg


!>  Set the cache size for the FFT.
!!
subroutine set_cache_size(nc)
   use module_fft_sg
   implicit none
   !Arguments
   integer, intent(in) :: nc
   !Local variables
   ncache = nc
END SUBROUTINE set_cache_size


!>  Give a number n_next > n compatible for the FFT
!!
subroutine fourier_dim(n,n_next)
   use module_fft_sg
   implicit none
   !Arguments
   integer, intent(in) :: n
   integer, intent(out) :: n_next
   !Local variables
   integer :: i
   loop_data: do i=1,ndata
      if (n <= i_data(i)) then
         n_next = i_data(i)
         return
      end if
   end do loop_data
   write(unit=*,fmt=*) "fourier_dim: ",n," is bigger than ",i_data(ndata)
   stop
END SUBROUTINE fourier_dim


!>  Give the dimensions of fft arrays
!!
subroutine dimensions_fft(n1,n2,n3,nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b)
   implicit none
   integer, intent(in)::n1,n2,n3
   integer,intent(out)::nd1,nd2,nd3,n1f,n3f,n1b,n3b,nd1f,nd3f,nd1b,nd3b
   ! Array sizes for the real-to-complex FFT: note that n1(there)=n1(here)+1
   ! and the same for n2,n3.
   nd1=n1+2
   nd2=n2+2
   nd3=n3+2
   ! n1b>=n1f;   n3f>=n3b
   n1f=(n1+2)/2 
   n3f=(n3+1)/2+1

   n1b=(n1+1)/2+1
   n3b=(n3+2)/2

   nd1f=n1f+1
   nd3f=n3f+1

   nd1b=n1b+1
   nd3b=n3b+1
END SUBROUTINE dimensions_fft


!1-dim complex-complex FFT routine
!the array in input is the first part, the output is the inzee
!the input is destroyed
subroutine fft_1d_ctoc(isign,nfft,n,zinout,inzee)
   use module_fft_sg
   implicit none
   integer, intent(in) :: n,nfft,isign
   integer, intent(out) :: inzee
   real(kind=8), dimension(2,nfft*n,2), intent(inout) :: zinout
   !local variables
   integer :: ic,i,ntrig
   !automatic arrays for the FFT
   integer, dimension(n_factors) :: after,now,before
   real(kind=8), dimension(:,:), allocatable :: trig

   ntrig=n
   allocate(trig(2,ntrig))
   !arrays for the FFT (to be halved)
   call ctrig_sg(n,ntrig,trig,after,before,now,isign,ic)
   !perform the FFT 
   inzee=1
   !write(15,*)halfft_cache(:,:,inzee)
   do i=1,ic
      call fftstp_sg(nfft,nfft,n,nfft,n,&
      zinout(1,1,inzee),zinout(1,1,3-inzee),&
      ntrig,trig,after(i),now(i),before(i),isign)
      inzee=3-inzee
   enddo

   deallocate(trig)
END SUBROUTINE fft_1d_ctoc

!> routine constituting the Building-block of other routines
subroutine FFT_1d(n,nfft,zinout,isign,inzee,transpose,real_input)
  use module_fft_sg
  implicit none
  logical, intent(in) :: transpose,real_input
  integer, intent(in) :: n,nfft,isign
  integer, intent(out) :: inzee
  real(kind=8), dimension(2,nfft*n,2), intent(inout) :: zinout
  !local variables
  integer :: ic,ntrig,npr,iam
  !automatic arrays for the FFT
  integer, dimension(n_factors) :: after,now,before
  real(kind=8), dimension(:,:), allocatable :: trig
  real(kind=8), allocatable, dimension(:,:,:) :: zw  

  ntrig=n
  allocate(trig(2,ntrig))
  !arrays for the FFT (to be halved)
  call ctrig_sg(n,ntrig,trig,after,before,now,isign,ic)
  !perform the FFT 

!!!!$omp critical
   allocate(zw(2,ncache/4,2))
!!!!$omp end critical

  inzee=1
  npr=1
!!!!$       npr=omp_get_num_threads()
  iam=0
  if (real_input) then
     call fft_1d_base(nfft,n,nfft,n,nfft,nfft*n,n,&
          ncache,ntrig,trig,after,now,before,ic,&
          isign,inzee,inzee,transpose,iam,npr,zinout,zw)
  else
     call fft_1d_base(nfft,n,nfft,n,nfft,nfft*n,n,&
          ncache,ntrig,trig,after,now,before,ic,&
          isign,inzee,inzee,transpose,iam,npr,zinout,zw)
  end if

  deallocate(zw)
!!!!!!!!!!$omp end parallel  
  deallocate(trig)
end subroutine FFT_1d

subroutine FFT_3d(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)

   use module_fft_sg
   use f_precisions
   implicit none !real(kind=8) (a-h,o-z)
   !Arguments
   integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,i_sign
   integer, intent(inout) :: inzee
   real(kind=8), intent(inout) :: z(2,nd1*nd2*nd3,2)
   !Local variables
   integer, dimension(n_factors) :: after,now,before
   real(kind=8), dimension(:), save, allocatable :: zw  
   real(kind=8), dimension(:,:), save, allocatable :: trig
   !local variables
   integer :: iam,ic,mm,npr,ntrig,nfft,outzee
   !$omp threadprivate(zw,trig)
   !$ integer :: omp_get_num_threads,omp_get_thread_num

   if (max(n1,n2,n3).gt.nfft_max) then
      write(*,*) 'One of the dimensions:',n1,n2,n3,' is bigger than ',nfft_max
      stop
   end if

   ! check whether input values are reasonable
   if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
   if (i_sign.ne.1 .and. i_sign.ne.-1) stop 'wrong i_sign'
   if (n1.gt.nd1) stop 'n1>nd1'
   if (n2.gt.nd2) stop 'n2>nd2'
   if (n3.gt.nd3) stop 'n3>nd3'

   ntrig=max(n1,n2,n3)

   npr=1
   iam=0

   !$omp parallel default(none) &
   !$omp shared(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee,outzee,ncache,ntrig,npr,nfft,mm) &
   !$omp private(ic,before,after,now,iam)
   !$ npr=omp_get_num_threads()
   !$ iam=omp_get_thread_num()
   !$omp critical (allocate_critical_fft3d)
   allocate(zw(ncache))
   allocate(trig(2,ntrig))
   !$omp end critical (allocate_critical_fft3d)
   call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)
   nfft=nd1*n2
   mm=nd1*nd2
   !$omp barrier !make sure that here everybody is ready
   call fft_1d_base(mm,nd3,mm,nd3,nfft,nd1*nd2*nd3,n3,&
        ncache,ntrig,trig,after,now,before,ic,&
        i_sign,inzee,outzee,.true.,iam,npr,z,zw)

   !$omp barrier
   if (n2 /= n3) then
      call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
   end if
   nfft=nd3*n1
   mm=nd3*nd1
   call fft_1d_base(mm,nd2,mm,nd2,nfft,nd1*nd2*nd3,n2,&
        ncache,ntrig,trig,after,now,before,ic,&
        i_sign,outzee,inzee,.true.,iam,npr,z,zw)
   !$omp barrier
   if (n1 /= n2) then
      call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
   end if
   nfft=nd2*n3
   mm=nd2*nd3
   call fft_1d_base(mm,nd1,mm,nd1,nfft,nd1*nd2*nd3,n1,&
        ncache,ntrig,trig,after,now,before,ic,&
        i_sign,inzee,outzee,.true.,iam,npr,z,zw)

   !$omp critical (deallocate_critical_fft3d)
   deallocate(zw)
   deallocate(trig)
   !$omp end critical (deallocate_critical_fft3d)
   !$omp end parallel  
   inzee=outzee
 end subroutine FFT_3d

!!$ subroutine FFT_3d_rtoc(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)
!!$
!!$   use module_fft_sg
!!$   implicit none !real(kind=8) (a-h,o-z)
!!$   !Arguments
!!$   integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,i_sign
!!$   integer, intent(inout) :: inzee
!!$   real(kind=8), intent(inout) :: z(2,nd1*nd2*nd3,2)
!!$   !Local variables
!!$   integer, dimension(n_factors) :: after,now,before
!!$   real(kind=8), allocatable, dimension(:,:,:) :: zw  
!!$   real(kind=8), dimension(:,:), allocatable :: trig
!!$   !local variables
!!$   integer :: iam,ic,mm,npr,ntrig,nfft
!!$
!!$   if (max(n1,n2,n3).gt.nfft_max) then
!!$      write(*,*) 'One of the dimensions:',n1,n2,n3,' is bigger than ',nfft_max
!!$      stop
!!$   end if
!!$
!!$   ! check whether input values are reasonable
!!$   if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
!!$   if (i_sign.ne.1 .and. i_sign.ne.-1) stop 'wrong i_sign'
!!$   if (n1.gt.nd1) stop 'n1>nd1'
!!$   if (n2.gt.nd2) stop 'n2>nd2'
!!$   if (n3.gt.nd3) stop 'n3>nd3'
!!$
!!$   ntrig=max(n1,n2,n3)
!!$   allocate(trig(2,ntrig))
!!$
!!$   ! Intel IFC does not understand default(private)
!!$!!!!!$omp parallel  default(private) &
!!$!!!!$omp parallel & 
!!$!!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
!!$!!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee,ncache) 
!!$
!!$   npr=1
!!$!!!!$       npr=omp_get_num_threads()
!!$   iam=0
!!$!!!!$       iam=omp_get_thread_num()
!!$
!!$!!!!$omp critical
!!$   allocate(zw(2,ncache/4,2))
!!$!!!!$omp end critical
!!$
!!$   call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)
!!$   nfft=nd1f*n2
!!$   mm=nd1f*nd2
!!$   !   call x0_to_z1_simple(x0,z1,inzee)
!!$   call x0_to_z1(x0,z1,inzee)
!!$   call fft_1d_base(mm,nd3,mm,nd3,nfft,nd1*nd2*nd3,n3,&
!!$        ncache,ntrig,trig,after,now,before,ic,&
!!$        i_sign,inzee,.true.,iam,npr,z1,zw)
!!$   call z1_to_z3(z1,z3,inzee)
!!$
!!$!!!!!!!!!$omp barrier
!!$   if (n2.ne.n3) then
!!$      call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
!!$   end if
!!$   nfft=nd3f*n1
!!$   mm=nd3f*nd1
!!$   call fft_1d_base(mm,nd2,mm,nd2,nfft,nd1*nd2*nd3,n2,&
!!$        ncache,ntrig,trig,after,now,before,ic,&
!!$        i_sign,inzee,.true.,iam,npr,z3,zw)
!!$!!!!!!!!!$omp barrier
!!$   if (n1.ne.n2) then
!!$      call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
!!$   end if
!!$   nfft=nd2*n3f
!!$   mm=nd2*nd3f
!!$   call fft_1d_base(mm,nd1,mm,nd1,nfft,nd1*nd2*nd3,n1,&
!!$        ncache,ntrig,trig,after,now,before,ic,&
!!$        i_sign,inzee,.true.,iam,npr,z3,zw)
!!$   deallocate(zw)
!!$!!!!!!!!!!$omp end parallel  
!!$   deallocate(trig)
!!$ END SUBROUTINE FFT_3D_RTOC

!!$subroutine fft_1d_base_rtoc(ndat_in,ld_in,ndat_out,ld_out,nfft,nin,n,&
!!$     nout,
!!$     ncache,ntrig,trig,after,now,before,ic,&
!!$     i_sign,inzee,transpose,iam,nthread,x,z,zw)
!!$  use f_precisions
!!$  use module_fft_sg, only: n_factors
!!$  implicit none
!!$  logical, intent(in) :: transpose
!!$  integer, intent(in) :: ndat_in,ld_in,ndat_out,ld_out,nfft,nin,nout,n
!!$  integer, intent(in) :: ncache,ntrig,ic
!!$  integer, intent(in) :: i_sign,iam,nthread
!!$  integer, intent(inout) :: inzee
!!$  integer, dimension(n_factors) :: after,now,before
!!$  real(f_double), dimension(2,ntrig), intent(in) :: trig
!!$  real(f_double), dimension(nin,2), intent(in) :: x
!!$  real(f_double), dimension(2,nout,2), intent(out) :: z
!!$  real(f_double), dimension(2,ncache/4,2), intent(inout) :: zw
!!$
!!$  !   call x0_to_z1_simple(x0,z1,inzee)
!!$  call x0_to_z1(x0,z1,inzee)
!!$  call fft_1d_base(mm,nd3,mm,nd3,nfft,nd1*nd2*nd3,n3,&
!!$       ncache,ntrig,trig,after,now,before,ic,&
!!$       i_sign,inzee,.true.,iam,npr,z1,zw)
!!$  call z1_to_z3(z1,z3,inzee)
!!$
!!$nout=max(nd1*nd3f,nd1f*nd3)*nd2
!!$
!!$!omp parallelization to be added
!!$lotomp=n2/nthread+1
!!$i2s=iam*lotomp
!!$i2e=min((iam+1)*lotomp,n2)
!!$n2omp=i2e-i2s
!!$ins=i2s*nd3*nd1f+1
!!$outs=i2s*nd1*nd3f+1
!!$call f_zero(2*nd1*nd3f*n2omp,z(1,outs,3-inzee))
!!$  call zhalf_to_z(n3f,n3,nd3,nd3f,&
!!$       n1f,nd1f,nd1,n2,n1,&
!!$       z(1,ins,inzee),z(1,outs,3-inzee))
!!$  inzee=3-inzee
!!$
!!$  contains
!!$
!!$    subroutine x0_to_z1(x0,z1,inzee)
!!$      ! Transform the real array x0 into a complex z1
!!$      ! real      part of z1: elements of x0 with odd  i1
!!$      ! imaginary part of z1: elements of x0 with even i1
!!$      implicit none
!!$      !Arguments
!!$      integer,intent(in)::inzee
!!$      real(kind=8),intent(in):: x0(n1,n2,n3)
!!$      real(kind=8),intent(out)::z1(2,nd1f,nd2,nd3,2)
!!$      !Local variables
!!$      integer :: i2,i3
!!$      z1=0.d0
!!$      do i3=1,n3
!!$         do i2=1,n2
!!$            ! 2*n1f=n1 for even n1
!!$            ! 2*n1f=n1+1 for odd n1. Then, we copy one more element than
!!$            ! necessary, but that's no problem.
!!$            call my_copy(z1(1,1,i2,i3,inzee),x0(1,i2,i3))
!!$         end do
!!$      end do
!!$    END SUBROUTINE x0_to_z1
!!$
!!$    subroutine my_copy(x,y)
!!$      ! copies complex array y into complex array x
!!$      implicit none
!!$      !Arguments
!!$      real(kind=8) :: x(2,n1f),y(2,n1f)
!!$      x=y
!!$    END SUBROUTINE my_copy
!!$
!!$    subroutine x0_to_z1_simple(x0,z1,inzee)
!!$      ! Transform the real array x0 into a complex z1
!!$      ! real      part of z1: elements of x0 with odd  i1
!!$      ! imaginary part of z1: elements of x0 with even i1
!!$      implicit none
!!$      !Arguments
!!$      integer, intent(in) :: inzee
!!$      real(kind=8), intent(in)  :: x0(n1,n2,n3)
!!$      real(kind=8), intent(out) :: z1(2,nd1f,nd2,nd3,2)
!!$      !Local variables
!!$      integer :: i1,i2,i3
!!$      if (n1f*2.eq.n1) then
!!$         do i3=1,n3
!!$            do i2=1,n2
!!$               do i1=1,n1f
!!$                  z1(1,i1,i2,i3,inzee)=x0(2*i1-1,i2,i3)
!!$                  z1(2,i1,i2,i3,inzee)=x0(2*i1  ,i2,i3)
!!$               end do
!!$            end do
!!$         end do
!!$      else ! n1=2*n1f-1
!!$         do i3=1,n3
!!$            do i2=1,n2
!!$               do i1=1,n1f-1
!!$                  z1(1,i1,i2,i3,inzee)=x0(2*i1-1,i2,i3)
!!$                  z1(2,i1,i2,i3,inzee)=x0(2*i1  ,i2,i3)
!!$               end do
!!$               z1(1,n1f,i2,i3,inzee)=x0(n1,i2,i3)
!!$            end do
!!$         end do
!!$      end if
!!$    END SUBROUTINE x0_to_z1_simple

!!$end subroutine fft_1d_base_rtoc

!!$n_left_in=n3f
!!$n_left_out=n3
!!$ldl_in=nd3
!!$ldl_out=nd3f
!!$n_half=n1f
!!$ld_in=nd1f
!!$ld_out=nd1
!!$n_right=n2
!!$n=n1
subroutine zhalf_to_z(n_left_in,n_left_out,ldl_in,ldl_out,&
     n_half,ld_in,ld_out,n_right,n,&
     z1,z3)
  ! transforms the array z1 that stores elements of z corresponding to even 
  ! and odd values of i1, as symmetric and antisymmetric combinations w.r.t.
  ! flip of i3,
  ! into the array z3 that stores only elements of z with i3=<ldl_out
  implicit none
  !Arguments
  integer, intent(in) :: n_left_in,n_left_out,ldl_in,ldl_out,n_half
  integer, intent(in) :: ld_in,ld_out,n_right,n
  real(kind=8), intent(in) :: z1(2,ldl_in,ld_in,*)
  !> this array has to be put to zero before
  real(kind=8), intent(inout) :: z3( 2,ldl_out,ld_out ,*)
  !Local variables
  integer :: i1,i2,i3

  if (n_half*2.eq.n) then
     ! i3=1
     do i2=1,n_right
        do i1=1,n_half
           z3(1,1,2*i1-1,i2)= 2.d0*z1(1,1,i1,i2)
           z3(2,1,2*i1-1,i2)= 0.d0
           z3(1,1,2*i1  ,i2)= 2.d0*z1(2,1,i1,i2)
           z3(2,1,2*i1  ,i2)= 0.d0
        end do
     end do

     do i2=1,n_right
        do i1=1,n_half
           do i3=2,n_left_in
              z3(1,i3,2*i1-1,i2)= z1(1,i3,i1,i2)+z1(1,n_left_out+2-i3,i1,i2)
              z3(1,i3,2*i1  ,i2)= z1(2,i3,i1,i2)+z1(2,n_left_out+2-i3,i1,i2)
              !!If statement added to cope with ifort 2011 optimization bug
              if (i3 /= n_left_out+2-i3) then
                 z3(2,i3,2*i1-1,i2)= z1(2,i3,i1,i2)-z1(2,n_left_out+2-i3,i1,i2)
                 z3(2,i3,2*i1  ,i2)=-z1(1,i3,i1,i2)+z1(1,n_left_out+2-i3,i1,i2)
              else
                 z3(2,i3,2*i1-1,i2)=0.0d0! z1(2,i3,i1,i2)-z1(2,n_left_out+2-i3,i1,i2)
                 z3(2,i3,2*i1  ,i2)=0.0d0!-z1(1,i3,i1,i2)+z1(1,n_left_out+2-i3,i1,i2)
              end if
           end do
        end do
     end do
  else ! n=2*n_half-1
     ! i3=1
     do i2=1,n_right
        do i1=1,n_half-1
           z3(1,1,2*i1-1,i2)= 2.d0*z1(1,1,i1,i2)
           z3(2,1,2*i1-1,i2)= 0.d0
           z3(1,1,2*i1  ,i2)= 2.d0*z1(2,1,i1,i2)
           z3(2,1,2*i1  ,i2)= 0.d0
        end do
     end do

     do i2=1,n_right
        do i1=1,n_half-1
           do i3=2,n_left_in
              z3(1,i3,2*i1-1,i2)= z1(1,i3,i1,i2)+z1(1,n_left_out+2-i3,i1,i2)
              z3(2,i3,2*i1-1,i2)= z1(2,i3,i1,i2)-z1(2,n_left_out+2-i3,i1,i2)
              z3(1,i3,2*i1  ,i2)= z1(2,i3,i1,i2)+z1(2,n_left_out+2-i3,i1,i2)
              z3(2,i3,2*i1  ,i2)=-z1(1,i3,i1,i2)+z1(1,n_left_out+2-i3,i1,i2)
           end do
        end do
     end do

     ! i1=n_half is treated separately: 2*n_half-1=n, but terms with 2*n_half are
     ! omitted

     do i2=1,n_right
        z3(1,1,n,i2)= 2.d0*z1(1,1,n_half,i2)
        z3(2,1,n,i2)= 0.d0
     end do

     do i2=1,n_right
        do i3=2,n_left_in
           z3(1,i3,n,i2)= z1(1,i3,n_half,i2)+z1(1,n_left_out+2-i3,n_half,i2)
           z3(2,i3,n,i2)= z1(2,i3,n_half,i2)-z1(2,n_left_out+2-i3,n_half,i2)
        end do
     end do
  end if
END SUBROUTINE zhalf_to_z

!!$subroutine fft_parallel_block(n,n2,n2p,n1d,n2d,n3d,&
!!$     nfft,lot,ldo,i2s,i2ps,i3s,&
!!$     ntrig,trig,after,now,before,ic,isign,&
!!$     zin,zcache,zout,inzee,tocache)
!!$  use f_precisions, only: dp => f_double
!!$  use module_fft_sg, only: n_factors
!!$  implicit none
!!$  logical, intent(in) :: tocache
!!$  integer, intent(in) :: n,n2,n2p,n1d,n2d,n3d,nfft,lot,ldo,i3s
!!$  integer, intent(in) :: ntrig,ic,isign
!!$  integer, intent(inout) :: inzee
!!$  integer, dimension(n_factors), intent(in) :: after,now,before
!!$  real(dp), dimension(2,ntrig), intent(in) :: trig
!!$  real(dp), dimension(2,n1d,n2d,n3d,*), intent(in) :: zin
!!$  integer, intent(inout) :: i2s,i2ps
!!$  real(dp), dimension(2*lot*n,2), intent(inout) :: zcache
!!$  real(dp), dimension(*), intent(inout) :: zout
!!$
!!$  if (inzee==FFT_SCRATCH) then
!!$     !input: i1,j2,j3,jp2,(jp3)
!!$     call transpose_and_pad_input([n,n2,n2p],[n1d,n2d,n3d],n1d,lot,nfft,i3s,&
!!$          i2s,i2ps,zin,zcache)
!!$     !output: j2,jp2,i1,j3,(jp3)
!!$     inzee=1
!!$  end if
!!$
!!$  !performing FFT
!!$  !input: i2,i1,j3,(jp3)
!!$  call fft_inorder(tocache,n,nfft,lot,ldo,inzee,&
!!$       ntrig,trig,after,now,before,ic,isign,&
!!$       zcache,zout)
!!$  !output: i2,I1,j3,(jp3)  
!!$
!!$  !in case it should be retrieved extract the cache results in the zout array
!!$  
!!$  !input: i2,I1,j3,(jp3)  
!!$  call transpose_output(ndims,lds,n1dim,lot,nfft,j3,&
!!$     i2s,i2ps,zcache(1,inzee),zout)
!!$  !output: I1,j2,j3,jp2,(jp3)
!!$
!!$end subroutine fft_parallel_block

!!$subroutine fft_inorder(cacheonly,n,ndat,ldz,ldo,inzee,&
!!$     ntrig,trig,after,now,before,ic,isign,&
!!$     zcache,zout)
!!$  use f_precisions, only: dp => f_double
!!$  use module_fft_sg, only: n_factors
!!$  implicit none
!!$  logical, intent(in) :: cacheonly !<temporary, if .true. zout is ignored
!!$  integer, intent(in) :: n,ndat,ldz,ldo
!!$  integer, intent(in) :: ntrig,ic,isign
!!$  integer, intent(inout) :: inzee
!!$  integer, dimension(n_factors), intent(in) :: after,now,before
!!$  real(dp), dimension(2,ntrig), intent(in) :: trig
!!$  real(dp), dimension(2*ldz*n,2), intent(in) :: zcache
!!$  real(dp), dimension(2,*), intent(out) :: zout
!!$  !local variable
!!$  integer :: i,iccache
!!$  !output: J2,Jp2,I1,j3,(jp3)
!!$  !performing FFT
!!$  !input: I2,I1,j3,(jp3)
!!$  if (cacheonly) then
!!$     iccache=ic
!!$  else
!!$     iccache=ic-1
!!$  end if
!!$
!!$  do i=1,iccache ! ic-1
!!$     call fftstp_sg(ldz,ndat,n,ldz,n,zcache(1,inzee),zcache(1,3-inzee),&
!!$          ntrig,trig,after(i),now(i),before(i),isign)
!!$     inzee=3-inzee
!!$  enddo
!!$  !storing the last step into zt array
!!$  if (iccache==ic) return
!!$  i=ic
!!$  call fftstp_sg(ldz,ndat,n,ldo,n,zcache(1,inzee),zout,&
!!$       ntrig,trig,after(i),now(i),before(i),isign)
!!$  !output: I2,i1,j3,(jp3)
!!$end subroutine fft_inorder


subroutine fft_1d_internal(mode,n,ndat_in,ld_in,ndat_out,ld_out,ldw,&
     nlines,nchuncks,&
     ntrig,trig,after,now,before,ic,i_sign,&
     inzee,transpose,iam,nthread,zi,zw,zo)
  use module_fft_sg
  use f_precisions
  implicit none
  logical, intent(in) :: transpose
  integer, intent(in) :: mode
  integer, intent(in) :: n,ndat_in,ld_in,ndat_out,ld_out,ldw,nlines,nchuncks
  integer, intent(in) :: ntrig,ic
  integer, intent(in) :: i_sign,iam,nthread
  integer, intent(inout) :: inzee
  integer, dimension(n_factors) :: after,now,before
  real(f_double), dimension(2,ntrig), intent(in) :: trig
  real(f_double), dimension(2,*), intent(in) :: zi
  real(f_double), dimension(2,*), intent(out) :: zo
  real(f_double), dimension(2,ldw,2), intent(inout) :: zw
  !local variables
  integer :: i,lotomp,ma,mb,nfft_th,jline,jline_transposed,jompa,jompb,inzet

  lotomp=(nlines)/nthread+1
  jompa=iam*lotomp+1
  jompb=min((iam+1)*lotomp,nlines)
  inzet=inzee
  i=1
  if (ic==1) then
     nfft_th=jompb-jompa+1
     jline=jompa
     call fft_kernel(mode,transpose)
  else
     !from here onwards ic>1
     do jline=jompa,jompb,nchuncks
        ma=jline
        mb=min(jline+(nchuncks-1),jompb)
        nfft_th=mb-ma+1
        i=1
        select case(mode)
        case(MODE_IO,MODE_IW)
           call fft_kernel(MODE_IW,.false.)
        case(MODE_WW,MODE_WO)
           call fft_kernel(MODE_WW,.false.)
        end select
        do i=2,ic-1
           call fft_kernel(MODE_WW,.false.)
        end do
        i=ic
        select case(mode)
        case(MODE_IO,MODE_WO)
           call fft_kernel(MODE_WO,transpose)
        case(MODE_WW,MODE_IW)
           call fft_kernel(MODE_WW,transpose)
        end select
     end do
  end if
  if (iam==0) inzee=inzet !as it is a shared variable

contains

  subroutine fft_kernel(internal_mode,transposed)
    implicit none
    logical, intent(in) :: transposed
    integer, intent(in) :: internal_mode

    select case(internal_mode)
    case(MODE_IO)
       if (transposed) then
          jline_transposed=(jline-1)*ld_out+1
          call fftrot_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
               zi(1,jline),zo(1,jline_transposed),&
               ntrig,trig,after(i),now(i),before(i),i_sign)
       else
          call fftstp_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
               zi(1,jline),zo(1,jline),&
               ntrig,trig,after(i),now(i),before(i),i_sign)

       end if
    case(MODE_IW)
       if (transposed) then
          call fftrot_sg(ndat_in,nfft_th,ld_in,nchuncks,n,&
               zi(1,jline),zw(1,1,inzet),&
               ntrig,trig,after(i),now(i),before(i),i_sign)
       else
          call fftstp_sg(ndat_in,nfft_th,ld_in,nchuncks,n,&
               zi(1,jline),zw(1,1,inzet), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
       end if
    case(MODE_WO)
       if (transposed) then
          jline_transposed=(jline-1)*ld_out+1
          call fftrot_sg(nchuncks,nfft_th,n,ndat_out,ld_out,&
               zw(1,1,inzet),zo(1,jline_transposed), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
       else
          call fftstp_sg(nchuncks,nfft_th,n,ndat_out,ld_out,&
               zw(1,1,inzet),zo(1,jline), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
       end if
    case(MODE_WW)
       if (transposed) then
          call fftrot_sg(nchuncks,nfft_th,n,nchuncks,n,&
               zw(1,1,inzet),zw(1,1,3-inzet), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
       else
          call fftstp_sg(nchuncks,nfft_th,n,nchuncks,n,&
               zw(1,1,inzet),zw(1,1,3-inzet), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
       end if
       inzet=3-inzet
    end select
  end subroutine fft_kernel
end subroutine fft_1d_internal


subroutine fft_1d_base(ndat_in,ld_in,ndat_out,ld_out,nfft,ninout,n,&
     ncache,ntrig,trig,after,now,before,ic,&
     i_sign,inzee,outzee,transpose,iam,nthread,z,zw)
  use f_precisions
  use module_fft_sg, only: n_factors
  implicit none
  logical, intent(in) :: transpose
  integer, intent(in) :: ndat_in,ld_in,ndat_out,ld_out,nfft,ninout,n
  integer, intent(in) :: ncache,ntrig,ic
  integer, intent(in) :: i_sign,iam,nthread
  integer, intent(in) :: inzee
  integer, intent(out) :: outzee
  integer, dimension(n_factors) :: after,now,before
  real(f_double), dimension(2,ntrig), intent(in) :: trig
  real(f_double), dimension(2,ninout,2), intent(inout) :: z
  real(f_double), dimension(2,ncache/4,2), intent(inout) :: zw
  !local variables
  logical :: real_input = .false.
  integer :: i,lotomp,ma,mb,nfft_th,j,jj,jompa,jompb,inzeep,inzet,lot,nn

  lot=max(1,ncache/(4*n))
  nn=lot

  inzet=inzee
  if (ic.eq.1 .or. ncache==0) then
     i=ic
     lotomp=nfft/nthread+1
     ma=iam*lotomp+1
     mb=min((iam+1)*lotomp,nfft)
     nfft_th=mb-ma+1
     j=ma
     if (transpose) then
        jj=j*ld_out-ld_out+1
        !input z(2,ndat_in,ld_in,inzet)
        !output z(2,ld_out,ndat_out,3-inzet)          
        call fftrot_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
             z(1,j,inzet),z(1,jj,3-inzet), &
             ntrig,trig,after(i),now(i),before(i),i_sign)
        if (real_input) then !only works when ld_out==n
           inzet=3-inzet
           call unpack_rfft(ndat_out,n,&
                z(1,jj,inzet),z(1,jj,3-inzet))
        end if
     else
        !input z(2,ndat_in,ld_in,inzet)
        !output z(2,ndat_out,ld_out,3-inzet)
        call fftstp_sg(ndat_in,nfft_th,ld_in,ndat_out,ld_out,&
             z(1,j,inzet),z(1,j,3-inzet), &
             ntrig,trig,after(i),now(i),before(i),i_sign)
        if (real_input) then !only works when ld_out==n
           inzet=3-inzet
           !here maybe the ndat_out has to be rethought
           call unpack_rfft_t((ndat_out-1)/2+1,ndat_out,n,&
                z(1,j,inzet),z(1,j,3-inzet))
        end if
     end if
  else
     lotomp=(nfft)/nthread+1
     jompa=iam*lotomp+1
     jompb=min((iam+1)*lotomp,nfft)
     do j=jompa,jompb,lot
        ma=j
        mb=min(j+(lot-1),jompb)
        nfft_th=mb-ma+1
        i=1
        inzeep=2
        call fftstp_sg(ndat_in,nfft_th,ld_in,nn,n,&
             z(1,j,inzet),zw(1,1,3-inzeep), &
             ntrig,trig,after(i),now(i),before(i),i_sign)
        inzeep=1
        do i=2,ic-1
           call fftstp_sg(nn,nfft_th,n,nn,n,&
                zw(1,1,inzeep),zw(1,1,3-inzeep), &
                ntrig,trig,after(i),now(i),before(i),i_sign)
           inzeep=3-inzeep
        end do
        i=ic
        if (transpose) then
           jj=j*ld_out-ld_out+1
           call fftrot_sg(nn,nfft_th,n,ndat_out,ld_out,&
                zw(1,1,inzeep),z(1,jj,3-inzet), &
                ntrig,trig,after(i),now(i),before(i),i_sign)
           if (real_input) then !only works when ld_out==n
              inzet=3-inzet
              call unpack_rfft(ndat_out,n,&
                   z(1,jj,inzet),z(1,jj,3-inzet))
           end if
        else
           call fftstp_sg(nn,nfft_th,n,ndat_out,ld_out,&
                zw(1,1,inzeep),z(1,j,3-inzet), &
                ntrig,trig,after(i),now(i),before(i),i_sign)
           if (real_input) then !only works when ld_out==n
              inzet=3-inzet
              call unpack_rfft_t((ndat_out-1)/2+1,ndat_out,n,&
                   z(1,j,inzet),z(1,j,3-inzet))
           end if
        end if
     end do
  end if
  inzet=3-inzet
  if (iam==0) outzee=inzet !as it is a shared variable

end subroutine fft_1d_base

!> routine to create the actual result of a complex fft from
!! a set of fft which have been packed from a real input
!! @warning: this only works if n is even
subroutine unpack_rfft(ndat,n,zin,zout)
  use f_precisions, only: f_double
  implicit none
  integer, intent(in) :: ndat,n
  real(f_double), intent(in) :: zin(2,n,*) !according ndat is odd or even(ndat-1)/2+1)
  real(f_double), intent(out) :: zout(2,n/2,2*ndat)
  !local variables
  integer :: idat,i,jp
  do idat=1,ndat
     zout(1,1,2*idat-1)=2.0_f_double*zin(1,1,idat)
     zout(2,1,2*idat-1)=0.0_f_double
     do i=2,n/2
        jp=n+2-i !j(-pj,n)+1
        zout(1,i,2*idat-1)=zin(1,i,idat)+zin(1,jp,idat)
        zout(2,i,2*idat-1)=zin(2,i,idat)-zin(2,jp,idat)
     end do
     zout(1,1,2*idat  )=2.0_f_double*zin(2,1,idat)
     zout(2,1,2*idat  )=0.0_f_double
     do i=2,n/2
        jp=n+2-i !j(-pj,n)+1
        zout(1,i,2*idat  )=zin(2,i,idat)+zin(2,jp,idat)
        zout(2,i,2*idat  )=-zin(1,i,idat)+zin(1,jp,idat)
     end do
  end do
  return
  if (2*(ndat/2)==ndat) return
  !last case, idat=(ndat-1)/2+1
  idat=(ndat-1)/2+1 
  zout(1,1,2*idat-1)=2.0_f_double*zin(1,1,idat)
  zout(2,1,2*idat-1)=0.0_f_double
  do i=2,n/2
     jp=n+2-i !j(-pj,n)+1
     zout(1,i,2*idat-1)=zin(1,i,idat)+zin(1,jp,idat)
     zout(2,i,2*idat-1)=zin(2,i,idat)-zin(2,jp,idat)
  end do

end subroutine unpack_rfft

!> routine to create the actual result of a complex fft from
!! a set of fft which have been packed from a real input
!! @warning: this only works if n is even
subroutine unpack_rfft_t(ndath,ndat,n,zin,zout)
  use f_precisions, only: f_double
  implicit none
  integer, intent(in) :: ndat,n,ndath
  real(f_double), intent(in) :: zin(2,ndath,n) !according ndat is odd or even(ndat-1)/2+1)
  real(f_double), intent(out) :: zout(2,ndat,n/2)
  !local variables
  integer :: idat,i,jp

  do idat=1,ndat/2
     zout(1,2*idat-1,1)=2.0_f_double*zin(1,idat,1)
     zout(2,2*idat-1,1)=0.0_f_double
     zout(1,2*idat  ,1)=2.0_f_double*zin(2,idat,1)
     zout(2,2*idat  ,1)=0.0_f_double
     do i=2,n/2
        jp=n+2-i !j(-pj,n)+1
        zout(1,2*idat-1,i)=zin(1,idat,i)+zin(1,idat,jp)
        zout(2,2*idat-1,i)=zin(2,idat,i)-zin(2,idat,jp)
        zout(1,2*idat  ,i)=zin(2,idat,i)+zin(2,idat,jp)
        zout(2,2*idat  ,i)=-zin(1,idat,i)+zin(1,idat,jp)
     end do
  end do
  if (2*(ndat/2)==ndat) return
  !last case, idat=(ndat-1)/2+1
  idat=(ndat-1)/2+1 
  zout(1,2*idat-1,1)=2.0_f_double*zin(1,idat,1)
  zout(2,2*idat-1,1)=0.0_f_double
  do i=2,n/2
     jp=n+2-i !j(-pj,n)+1
     zout(1,2*idat-1,i)=zin(1,idat,i)+zin(1,idat,jp)
     zout(2,2*idat-1,i)=zin(2,idat,i)-zin(2,idat,jp)
  end do

end subroutine unpack_rfft_t


!> @brief Calculates the discrete Fourier transform
!!
!! F(I1,I2,I3)= S_(j1,j2,j3) EXP(i_sign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
!!
!!  with optimal performance on vector computer, workstations and 
!!  multiprocessor shared memory computers using OpenMP compiler directives
!!   INPUT:
!!       n1,n2,n3:physical dimension of the transform. It must be a 
!!                product of the prime factors 2,3,5, but greater than 3. 
!!               If two ni's are equal it is recommended to place them 
!!               behind each other.
!!       nd1,nd2,nd3:memory dimension of Z. ndi must always be greater or 
!!                   equal than ni. On a vector machine, it is recomended 
!!                  to chose ndi=ni if ni is odd and ndi=ni+1 if ni is 
!!                  even to obtain optimal execution speed. On RISC 
!!                  machines ndi=ni is usually fine for odd ni, for even 
!!                  ni one should try ndi=ni+1, ni+2, ni+4 to find the 
!!                  optimal performance. 
!!      inzee=1: first part of Z is data (input) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           Z(1,i1,i2,i3,inzee)=real(R(i1,i2,i3))
!!           Z(2,i1,i2,i3,inzee)=imag(R(i1,i2,i3))
!!   OUTPUT:
!!      inzee=1: first part of Z is data (output) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
!!           imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
!!      inzee on output is in general different from inzee on input
!!   The input data are always overwritten independently of the 
!!  value of inzee.
!! 
!!
!!
subroutine FFT(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee)

   use module_fft_sg
   implicit real(kind=8) (a-h,o-z)
   !!!$      interface
   !!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
   !!!!$        end function omp_get_num_threads
   !!!$      end interface
   !!!!$      interface
   !!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
   !!!!$        end function omp_get_thread_num
   !!!!$      end interface

   !Arguments
   integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3,i_sign
   integer, intent(inout) :: inzee
   real(kind=8), intent(inout) :: z(2,nd1*nd2*nd3,2)
   !Local variables
   integer, dimension(n_factors) :: after,now,before
   real(kind=8), allocatable, dimension(:,:,:) :: zw  
   real(kind=8), dimension(:,:), allocatable :: trig

   if (max(n1,n2,n3).gt.nfft_max) then
      write(*,*) 'One of the dimensions:',n1,n2,n3,' is bigger than ',nfft_max
      stop
   end if

   ! check whether input values are reasonable
   if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
   if (i_sign.ne.1 .and. i_sign.ne.-1) stop 'wrong i_sign'
   if (n1.gt.nd1) stop 'n1>nd1'
   if (n2.gt.nd2) stop 'n2>nd2'
   if (n3.gt.nd3) stop 'n3>nd3'

   ntrig=max(n1,n2,n3)
   allocate(trig(2,ntrig))

   ! vector computer with memory banks:
   if (ncache.eq.0) then
      call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)
      nfft=nd1*n2
      mm=nd1*nd2
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd3,mm,nd3,z(1,1,inzee),z(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee

      if (n2.ne.n3) then
         call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
      end if
      nfft=nd3*n1
      mm=nd3*nd1
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd2,mm,nd2,z(1,1,inzee),z(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee

      if (n1.ne.n2) then
         call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
      end if
      nfft=nd2*n3
      mm=nd2*nd3
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd1,mm,nd1,z(1,1,inzee),z(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee

      ! RISC machine with cache:
   else
      ! Intel IFC does not understand default(private)
      !!!!!$omp parallel  default(private) &
      !!!!$omp parallel & 
      !!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
      !!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee,ncache) 
      npr=1
      !!!!$       npr=omp_get_num_threads()
      iam=0
      !!!!$       iam=omp_get_thread_num()
      !      write(6,*) 'npr,iam',npr,iam
      ! Critical section only necessary on Intel
      !!!!$omp critical
      allocate(zw(2,ncache/4,2))
      !!!!$omp end critical

      inzet=inzee
      ! TRANSFORM ALONG Z AXIS

      mm=nd1*nd2
      m=nd3
      lot=max(1,ncache/(4*n3))
      nn=lot
      n=n3
      if (2*n*lot*2.gt.ncache) then
         stop 'ncache1'
      end if

      call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)

      if (ic.eq.1) then
         i=ic
         lotomp=(nd1*n2)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd1*n2)
         nfft=mb-ma+1
         j=ma
         jj=j*nd3-nd3+1
         call fftrot_sg(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)

      else

         lotomp=(nd1*n2)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd1*n2)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd3-nd3+1

            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1

            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if

      inzet=3-inzet

      !!!!!!!!!$omp barrier

      ! TRANSFORM ALONG Y AXIS
      mm=nd3*nd1
      m=nd2
      lot=max(1,ncache/(4*n2))
      nn=lot
      n=n2
      if (2*n*lot*2.gt.ncache) then
         stop 'ncache2'
      end if

      if (n2.ne.n3) then
         call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
      end if

      if (ic.eq.1) then
         i=ic
         lotomp=(nd3*n1)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd3*n1)
         nfft=mb-ma+1
         j=ma
         jj=j*nd2-nd2+1
         call fftrot_sg(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
      else
         lotomp=(nd3*n1)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd3*n1)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd2-nd2+1
            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1
            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if
      inzet=3-inzet

      !!!!!!!!$omp barrier

      ! TRANSFORM ALONG X AXIS
      mm=nd2*nd3
      m=nd1
      lot=max(1,ncache/(4*n1))
      nn=lot
      n=n1
      if (2*n*lot*2.gt.ncache) then
         write(*,*) ncache,lot,n,2*n*lot*2
         stop 'ncache3'
      end if

      if (n1.ne.n2) then
         call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
      end if

      if (ic.eq.1) then
         i=ic
         lotomp=(nd2*n3)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd2*n3)
         nfft=mb-ma+1
         j=ma
         jj=j*nd1-nd1+1
         call fftrot_sg(mm,nfft,m,mm,m,z(1,j,inzet),z(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
      else

         lotomp=(nd2*n3)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd2*n3)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd1-nd1+1
            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1
            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if
      inzet=3-inzet

      deallocate(zw)
      if (iam.eq.0) inzee=inzet
      !!!!!!!!!!$omp end parallel  

   end if

   deallocate(trig)
END SUBROUTINE FFT


!> @brief Calculates the discrete Fourier transform 
!!
!! F(I1,I2,I3) = S_(j1,j2,j3) EXP(i_sign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
!!  with optimal performance on vector computer, workstations and 
!!  multiprocessor shared memory computers using OpenMP compiler directives
!!   INPUT:
!!       n1,n2,n3:physical dimension of the transform. It must be a 
!!                product of the prime factors 2,3,5, but greater than 3. 
!!               If two ni's are equal it is recommended to place them 
!!               behind each other.
!!       nd1,nd2,nd3:memory dimension of Z. ndi must always be greater or 
!!                   equal than ni. On a vector machine, it is recomended 
!!                  to chose ndi=ni if ni is odd and ndi=ni+1 if ni is 
!!                  even to obtain optimal execution speed. On RISC 
!!                  machines ndi=ni is usually fine for odd ni, for even 
!!                  ni one should try ndi=ni+1, ni+2, ni+4 to find the 
!!                  optimal performance. 
!!      inzee=1: first part of Z is data (input) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           Z(1,i1,i2,i3,inzee)=real(R(i1,i2,i3))
!!           Z(2,i1,i2,i3,inzee)=imag(R(i1,i2,i3))
!!   OUTPUT:
!!      inzee=1: first part of Z is data (output) array, 
!!               second part work array
!!      inzee=2: first part of Z is work array, second part data array
!!           real(F(i1,i2,i3))=Z(1,i1,i2,i3,inzee)
!!           imag(F(i1,i2,i3))=Z(2,i1,i2,i3,inzee)
!!      inzee on output is in general different from inzee on input
!!   The input data are always overwritten independently of the 
!!  value of inzee.
!!
!!
subroutine FFT_for(n1,n2,n3,n1f,n3f,nd1,nd2,nd3,nd1f,nd3f,x0,z1,z3,inzee)

   use module_fft_sg
   implicit real(kind=8) (a-h,o-z)

   !!!!!!!$      interface
   !!!!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
   !!!!!!$        end function omp_get_num_threads
   !!!!!!$      end interface
   !!!!!!!$      interface
   !!!!!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
   !!!!!!!!$        end function omp_get_thread_num
   !!!!!!!!$      end interface

   !Arguments
   integer, intent(in) :: n1,n2,n3,n1f,n3f
   integer, intent(in) :: nd1,nd2,nd3,nd1f,nd3f
   integer, intent(inout) :: inzee
   real(kind=8), intent(in) :: x0(n1,n2,n3)
   real(kind=8), intent(out) :: z3(2,nd1*nd2*nd3f,2)
   real(kind=8), intent(inout) :: z1(2,nd1f*nd2*nd3,2) ! work array
   !Local variables
   real(kind=8), allocatable, dimension(:,:,:) :: zw  
   integer, dimension(n_factors) :: after,now,before
   real(kind=8), dimension(:,:), allocatable :: trig
   integer::mm,nffta,i_sign=1

   if (max(n1,n2,n3).gt.nfft_max) then
      write(*,*) 'Dimension bigger than ', nfft_max
      stop
   end if

   inzee=1

   !*******************Alexey*********************************************************************
   !         ncache=0
   !**********************************************************************************************

   ! check whether input values are reasonable
   if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
   if (i_sign.ne.1 .and. i_sign.ne.-1) stop 'wrong i_sign'
   if (n1.gt.nd1) stop 'n1>nd1'
   if (n2.gt.nd2) stop 'n2>nd2'
   if (n3.gt.nd3) stop 'n3>nd3'

   ntrig=max(n1,n2,n3)
   allocate(trig(2,ntrig))

   !   call x0_to_z1_simple(x0,z1,inzee)
   call x0_to_z1(x0,z1,inzee)

   if (ncache.eq.0) then
      ! vector computer with memory banks:
      call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)
      mm=nd1f*nd2 
      nffta=nd1f*n2
      do i=1,ic-1
         call fftstp_sg(mm,nffta,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nffta,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee
      call z1_to_z3(z1,z3,inzee)
      if (n2.ne.n3) then
         call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
      end if
      nfft=nd3f*n1
      mm=nd3f*nd1
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd2,mm,nd2,z3(1,1,inzee),z3(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd2,mm,nd2,z3(1,1,inzee),z3(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee
      if (n1.ne.n2) then
         call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
      end if
      nfft=nd2*n3f
      mm=nd2*nd3f
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do 
      i=ic
      call fftrot_sg(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee

      ! RISC machine with cache:
   else
      ! Intel IFC does not understand default(private)
      !!!!$omp parallel  default(private) &
      !!!!$omp parallel & 
      !!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
      !!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee,ncache) 
      npr=1
      !!!!!!!$       npr=omp_get_num_threads()
      iam=0
      !!!!!!!$       iam=omp_get_thread_num()
      !!!!!        write(6,*) 'npr,iam',npr,iam
      ! Critical section only necessary on Intel
      !!!!!$omp critical
      allocate(zw(2,ncache/4,2))
      !!!!!$omp end critical

      inzet=inzee
      ! TRANSFORM ALONG Z AXIS

      mm=nd1f*nd2
      m=nd3
      lot=max(1,ncache/(4*n3))
      nn=lot
      n=n3
      if (2*n*lot*2.gt.ncache) then 
         write(*,*) ncache,lot,n,2*n*lot*2
         stop 'ncache1'
      end if
      call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)
      if (ic.eq.1) then
         i=ic
         lotomp=(nd1f*n2)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd1f*n2)
         nfft=mb-ma+1
         j=ma
         jj=j*nd3-nd3+1
         call fftrot_sg(mm,nfft,m,mm,m,z1(1,j,inzet),z1(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
      else
         lotomp=(nd1f*n2)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd1f*n2)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd3-nd3+1
            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z1(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1
            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z1(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if

      inzet=3-inzet

      call z1_to_z3(z1,z3,inzet)
      !!!!!!!!!$omp barrier

      ! TRANSFORM ALONG Y AXIS
      mm=nd3f*nd1
      m=nd2
      lot=max(1,ncache/(4*n2))
      nn=lot
      n=n2
      if (2*n*lot*2.gt.ncache) then
         write(*,*) ncache,lot,n,2*n*lot*2
         stop 'ncache2'
      end if
      if (n2.ne.n3) then
         call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
      end if
      if (ic.eq.1) then
         i=ic
         lotomp=(nd3f*n1)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd3f*n1)
         nfft=mb-ma+1
         j=ma
         jj=j*nd2-nd2+1
         call fftrot_sg(mm,nfft,m,mm,m,z3(1,j,inzet),z3(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
      else
         lotomp=(nd3f*n1)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd3f*n1)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd2-nd2+1
            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z3(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1
            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z3(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if
      inzet=3-inzet

      !!!!!!!!!$omp barrier

      ! TRANSFORM ALONG X AXIS
      mm=nd2*nd3f
      m=nd1
      lot=max(1,ncache/(4*n1))
      nn=lot
      n=n1
      if (2*n*lot*2.gt.ncache) then
         write(*,*) ncache,lot,n,2*n*lot*2
         stop 'ncache3'
      end if
      if (n1.ne.n2) then
         call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
      end if
      if (ic.eq.1) then
         i=ic
         lotomp=(nd2*n3f)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd2*n3f)
         nfft=mb-ma+1
         j=ma
         jj=j*nd1-nd1+1
         call fftrot_sg(mm,nfft,m,mm,m,z3(1,j,inzet),z3(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
      else
         lotomp=(nd2*n3f)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd2*n3f)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd1-nd1+1
            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z3(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1
            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z3(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if
      inzet=3-inzet
      deallocate(zw)
      if (iam.eq.0) inzee=inzet
      !!!!!!!!!$omp end parallel  

   end if

   deallocate(trig)

   contains

   subroutine x0_to_z1(x0,z1,inzee)
      ! Transform the real array x0 into a complex z1
      ! real      part of z1: elements of x0 with odd  i1
      ! imaginary part of z1: elements of x0 with even i1
      implicit none
      !Arguments
      integer,intent(in)::inzee
      real(kind=8),intent(in):: x0(n1,n2,n3)
      real(kind=8),intent(out)::z1(2,nd1f,nd2,nd3,2)
      !Local variables
      integer :: i2,i3
z1=0.d0
      do i3=1,n3
         do i2=1,n2
            ! 2*n1f=n1 for even n1
            ! 2*n1f=n1+1 for odd n1. Then, we copy one more element than
            ! necessary, but that's no problem.
            call my_copy(z1(1,1,i2,i3,inzee),x0(1,i2,i3))
         end do
      end do
   END SUBROUTINE x0_to_z1

   subroutine my_copy(x,y)
      ! copies complex array y into complex array x
      implicit none
      !Arguments
      real(kind=8) :: x(2,n1f),y(2,n1f)
      x=y
   END SUBROUTINE my_copy

   subroutine x0_to_z1_simple(x0,z1,inzee)
      ! Transform the real array x0 into a complex z1
      ! real      part of z1: elements of x0 with odd  i1
      ! imaginary part of z1: elements of x0 with even i1
      implicit none
      !Arguments
      integer, intent(in) :: inzee
      real(kind=8), intent(in)  :: x0(n1,n2,n3)
      real(kind=8), intent(out) :: z1(2,nd1f,nd2,nd3,2)
      !Local variables
      integer :: i1,i2,i3
      if (n1f*2.eq.n1) then
         do i3=1,n3
            do i2=1,n2
               do i1=1,n1f
                  z1(1,i1,i2,i3,inzee)=x0(2*i1-1,i2,i3)
                  z1(2,i1,i2,i3,inzee)=x0(2*i1  ,i2,i3)
               end do
            end do
         end do
      else ! n1=2*n1f-1
         do i3=1,n3
            do i2=1,n2
               do i1=1,n1f-1
                  z1(1,i1,i2,i3,inzee)=x0(2*i1-1,i2,i3)
                  z1(2,i1,i2,i3,inzee)=x0(2*i1  ,i2,i3)
               end do
               z1(1,n1f,i2,i3,inzee)=x0(n1,i2,i3)
            end do
         end do
      end if
   END SUBROUTINE x0_to_z1_simple

   subroutine z1_to_z3(z1,z3,inzee)
      ! transforms the array z1 that stores elements of z corresponding to even 
      ! and odd values of i1, as symmetric and antisymmetric combinations w.r.t.
      ! flip of i3,
      ! into the array z3 that stores only elements of z with i3=<nd3f
      implicit none
      !Arguments
      integer, intent(in) :: inzee
      real(kind=8), intent(in) :: z1(2,nd3,nd1f,nd2,2)
      real(kind=8), intent(out) :: z3( 2,nd3f,nd1 ,nd2,2)
      !Local variables
      integer :: i1,i2,i3
z3=0.d0
      if (n1f*2.eq.n1) then
         ! i3=1
         do i2=1,n2
            do i1=1,n1f
               z3(1,1,2*i1-1,i2,inzee)= 2.d0*z1(1,1,i1,i2,inzee)
               z3(2,1,2*i1-1,i2,inzee)= 0.d0
               z3(1,1,2*i1  ,i2,inzee)= 2.d0*z1(2,1,i1,i2,inzee)
               z3(2,1,2*i1  ,i2,inzee)= 0.d0
            end do
         end do

         do i2=1,n2
            do i1=1,n1f
               do i3=2,n3f
                  z3(1,i3,2*i1-1,i2,inzee)= z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                  z3(1,i3,2*i1  ,i2,inzee)= z1(2,i3,i1,i2,inzee)+z1(2,n3+2-i3,i1,i2,inzee)
                  !!If statement added to cope with ifort 2011 optimization bug
                  if (i3 /= n3+2-i3) then
                     z3(2,i3,2*i1-1,i2,inzee)= z1(2,i3,i1,i2,inzee)-z1(2,n3+2-i3,i1,i2,inzee)
                     z3(2,i3,2*i1  ,i2,inzee)=-z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                  else
                     z3(2,i3,2*i1-1,i2,inzee)=0.0d0! z1(2,i3,i1,i2,inzee)-z1(2,n3+2-i3,i1,i2,inzee)
                     z3(2,i3,2*i1  ,i2,inzee)=0.0d0!-z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                  end if
               end do
            end do
         end do
      else ! n1=2*n1f-1
         ! i3=1
         do i2=1,n2
            do i1=1,n1f-1
               z3(1,1,2*i1-1,i2,inzee)= 2.d0*z1(1,1,i1,i2,inzee)
               z3(2,1,2*i1-1,i2,inzee)= 0.d0
               z3(1,1,2*i1  ,i2,inzee)= 2.d0*z1(2,1,i1,i2,inzee)
               z3(2,1,2*i1  ,i2,inzee)= 0.d0
            end do
         end do

         do i2=1,n2
            do i1=1,n1f-1
               do i3=2,n3f
                  z3(1,i3,2*i1-1,i2,inzee)= z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
                  z3(2,i3,2*i1-1,i2,inzee)= z1(2,i3,i1,i2,inzee)-z1(2,n3+2-i3,i1,i2,inzee)
                  z3(1,i3,2*i1  ,i2,inzee)= z1(2,i3,i1,i2,inzee)+z1(2,n3+2-i3,i1,i2,inzee)
                  z3(2,i3,2*i1  ,i2,inzee)=-z1(1,i3,i1,i2,inzee)+z1(1,n3+2-i3,i1,i2,inzee)
               end do
            end do
         end do

         ! i1=n1f is treated separately: 2*n1f-1=n1, but terms with 2*n1f are
         ! omitted

         do i2=1,n2
            z3(1,1,n1,i2,inzee)= 2.d0*z1(1,1,n1f,i2,inzee)
            z3(2,1,n1,i2,inzee)= 0.d0
         end do

         do i2=1,n2
            do i3=2,n3f
               z3(1,i3,n1,i2,inzee)= z1(1,i3,n1f,i2,inzee)+z1(1,n3+2-i3,n1f,i2,inzee)
               z3(2,i3,n1,i2,inzee)= z1(2,i3,n1f,i2,inzee)-z1(2,n3+2-i3,n1f,i2,inzee)
            end do
         end do
      end if
   END SUBROUTINE z1_to_z3

END SUBROUTINE fft_for


!> Calculates the discrete fourier transform 
!! F(I1,I2,I3) = S_(j1,j2,j3) EXP(i_sign*i*2*pi*(j1*i1/n1+j2*i2/n2+j3*i3/n3)) R(j1,j2,j3)
!!
subroutine FFT_back(n1,n2,n3,n1b,n3f,n3b,nd1,nd2,nd3,nd1b,nd3f,nd3b,y,z1,z3,inzee)

   use module_fft_sg
   implicit real(kind=8) (a-h,o-z)

   !!!!!!!$      interface
   !!!!!!!$        integer ( kind=4 ) function omp_get_num_threads ( )
   !!!!!!$        end function omp_get_num_threads
   !!!!!!$      end interface
   !!!!!!!$      interface
   !!!!!!!$        integer ( kind=4 ) function omp_get_thread_num ( )
   !!!!!!!$        end function omp_get_thread_num
   !!!!!!!!$      end interface

   !Arguments
   integer, intent(in) :: n1,n2,n3,nd1,nd2,nd3
   integer, intent(in) :: n1b,n3f,n3b,nd1b,nd3f,nd3b
   real(kind=8),intent(out) :: y(n1,n2,n3)
   real(kind=8),intent(inout) :: z3(2,nd1*nd2*nd3b,2)
   real(kind=8)                :: z1(2,nd1b*nd2*nd3,2) ! work array
   !Local variables
   real(kind=8), allocatable, dimension(:,:,:) :: zw  
   real(kind=8), dimension(:,:), allocatable :: trig
   integer, dimension(n_factors) :: after,now,before

   i_sign=-1
   if (max(n1,n2,n3).gt.1024) stop '1024'

   ! check whether input values are reasonable
   if (inzee.le.0 .or. inzee.ge.3) stop 'wrong inzee'
   if (i_sign.ne.1 .and. i_sign.ne.-1) stop 'wrong i_sign'
   if (n1.gt.nd1) stop 'n1>nd1'
   if (n2.gt.nd2) stop 'n2>nd2'
   if (n3.gt.nd3) stop 'n3>nd3'

   ntrig=max(n1,n2,n3)
   allocate(trig(2,ntrig))

   !call z3_to_z1(z3,z1,inzee)

   if (ncache.eq.0) then
      ! vector computer with memory banks:
      call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)
      nfft=nd1b*n2
      mm=nd1b*nd2
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd3,mm,nd3,z1(1,1,inzee),z1(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)

      inzee=3-inzee

      if (n2.ne.n3) then
         call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
      end if
      nfft=nd3*n1b
      mm=nd3*nd1b
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd2,mm,nd2,z1(1,1,inzee),z1(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd2,mm,nd2,z1(1,1,inzee),z1(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee

      ! here we transform back from z1 to z3
      call z1_to_z3(z1,z3,inzee)

      if (n1.ne.n2) then
         call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
      end if
      nfft=nd2*n3b
      mm=nd2*nd3b
      do i=1,ic-1
         call fftstp_sg(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
         inzee=3-inzee
      end do
      i=ic
      call fftrot_sg(mm,nfft,nd1,mm,nd1,z3(1,1,inzee),z3(1,1,3-inzee), &
      ntrig,trig,after(i),now(i),before(i),i_sign)
      inzee=3-inzee

      call z3_to_y(z3,y,inzee)

      ! RISC machine with cache:
   else
      ! Intel IFC does not understand default(private)
      !!!!!!!!!$omp parallel  default(private) &
      !!!!!$omp parallel & 
      !!!!!$omp private(zw,trig,before,after,now,i,j,iam,npr,jj,ma,mb,mm,ic,n,m,jompa,jompb,lot,lotomp,inzeep,inzet,nn,nfft) &
      !!!!!$omp shared(n1,n2,n3,nd1,nd2,nd3,z,i_sign,inzee,ncache) 
      npr=1
      !!!!!$     npr=omp_get_num_threads()
      iam=0
      !!!$       iam=omp_get_thread_num()
      !      write(6,*) 'npr,iam',npr,iam
      ! Critical section only necessary on Intel
      !!!!!!$omp critical
      allocate(zw(2,ncache/4,2))
      !!!!!!!$omp end critical

      inzet=inzee
      ! TRANSFORM ALONG Z AXIS

      mm=nd1b*nd2
      m=nd3
      lot=max(1,ncache/(4*n3))
      nn=lot
      n=n3
      if (2*n*lot*2.gt.ncache) then
         write(*,*) ncache,lot,n,2*n*lot*2
         stop 'ncache1'
      end if

      call ctrig_sg(n3,ntrig,trig,after,before,now,i_sign,ic)

      if (ic.eq.1) then
         i=ic
         lotomp=(nd1b*n2)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd1b*n2)
         nfft=mb-ma+1
         j=ma
         jj=j*nd3-nd3+1
         call fftrot_sg(mm,nfft,m,mm,m,z1(1,j,inzet),z1(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
      else

         lotomp=(nd1b*n2)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd1b*n2)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd3-nd3+1

            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z1(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1

            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z1(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if

      inzet=3-inzet

      !!!!!!!!$omp barrier

      ! TRANSFORM ALONG Y AXIS
      mm=nd3*nd1b
      m=nd2
      lot=max(1,ncache/(4*n2))
      nn=lot
      n=n2
      if (2*n*lot*2.gt.ncache) then
         write(*,*) ncache,lot,n,2*n*lot*2
         stop 'ncache2'
      end if
      if (n2.ne.n3) then
         call ctrig_sg(n2,ntrig,trig,after,before,now,i_sign,ic)
      end if

      if (ic.eq.1) then
         i=ic
         lotomp=(nd3*n1b)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd3*n1b)
         nfft=mb-ma+1
         j=ma
         jj=j*nd2-nd2+1
         call fftrot_sg(mm,nfft,m,mm,m,z1(1,j,inzet),z1(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)
      else

         lotomp=(nd3*n1b)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd3*n1b)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd2-nd2+1

            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z1(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1

            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do

            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z1(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if
      inzet=3-inzet

      call z1_to_z3(z1,z3,inzet)
      !!!!!!$omp barrier

      ! TRANSFORM ALONG X AXIS
      mm=nd2*nd3b
      m=nd1
      lot=max(1,ncache/(4*n1))
      nn=lot
      n=n1
      if (2*n*lot*2.gt.ncache) then
         write(*,*) ncache,lot,n,2*n*lot*2
         stop 'ncache3'
      end if
      if (n1.ne.n2) then
         call ctrig_sg(n1,ntrig,trig,after,before,now,i_sign,ic)
      end if

      if (ic.eq.1) then
         i=ic
         lotomp=(nd2*n3b)/npr+1
         ma=iam*lotomp+1
         mb=min((iam+1)*lotomp,nd2*n3b)
         nfft=mb-ma+1
         j=ma
         jj=j*nd1-nd1+1
         call fftrot_sg(mm,nfft,m,mm,m,z3(1,j,inzet),z3(1,jj,3-inzet), &
         ntrig,trig,after(i),now(i),before(i),i_sign)

      else

         lotomp=(nd2*n3b)/npr+1
         jompa=iam*lotomp+1
         jompb=min((iam+1)*lotomp,nd2*n3b)
         do j=jompa,jompb,lot
            ma=j
            mb=min(j+(lot-1),jompb)
            nfft=mb-ma+1
            jj=j*nd1-nd1+1

            i=1
            inzeep=2
            call fftstp_sg(mm,nfft,m,nn,n,z3(1,j,inzet),zw(1,1,3-inzeep), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
            inzeep=1

            do i=2,ic-1
               call fftstp_sg(nn,nfft,n,nn,n,zw(1,1,inzeep),zw(1,1,3-inzeep), &
               ntrig,trig,after(i),now(i),before(i),i_sign)
               inzeep=3-inzeep
            end do
            i=ic
            call fftrot_sg(nn,nfft,n,mm,m,zw(1,1,inzeep),z3(1,jj,3-inzet), &
            ntrig,trig,after(i),now(i),before(i),i_sign)
         end do
      end if
      inzet=3-inzet

      call z3_to_y(z3,y,inzet)
      deallocate(zw)
      if (iam.eq.0) inzee=inzet
      !!!!!!!!!!!$omp end parallel  


   end if

   deallocate(trig)

   contains

   subroutine z3_to_z1(z3,z1,inzee)
      ! transforms the data from the format z3:
      ! output of fft_for: stores only elements of z with i3=<nd3f
      ! 
      ! to the format z1:
      ! input of fft_back: stores only elements of z with i1=<nd1b
      implicit none
      integer,intent(in)::inzee
      real(kind=8),intent(in):: z3(2,nd1 ,nd2,nd3f,2)
      real(kind=8),intent(out)::z1(2,nd1b,nd2,nd3,2)
      integer i1,i2,i3
z1=0.d0
      ! i3=1: then z1 is contained in z3 
      do i2=1,n2
         do i1=1,n1b
            z1(1,i1,i2,1,inzee)=z3(1,i1,i2,1,inzee)
            z1(2,i1,i2,1,inzee)=z3(2,i1,i2,1,inzee)
         end do
      end do    

      do i3=2,n3f
         ! i2=1
         ! i1=1
         z1(1,1,1,i3,inzee)=z3(1,1,1,i3,inzee)
         z1(2,1,1,i3,inzee)=z3(2,1,1,i3,inzee)

         z1(1,1,1,n3+2-i3,inzee)=z3(1,1,1,i3,inzee)
         z1(2,1,1,n3+2-i3,inzee)=-z3(2,1,1,i3,inzee)

         ! i2=1
         do i1=2,n1b
            z1(1,i1,1,i3,inzee)=z3(1,i1,1,i3,inzee)
            z1(2,i1,1,i3,inzee)=z3(2,i1,1,i3,inzee)

            z1(1,i1,1,n3+2-i3,inzee)= z3(1,n1+2-i1,1,i3,inzee)
            z1(2,i1,1,n3+2-i3,inzee)=-z3(2,n1+2-i1,1,i3,inzee)
         end do

         do i2=2,n2
            ! i1=1
            z1(1,1,i2,i3,inzee)=z3(1,1,i2,i3,inzee)
            z1(2,1,i2,i3,inzee)=z3(2,1,i2,i3,inzee)

            z1(1,1,i2,n3+2-i3,inzee)= z3(1,1,n2+2-i2,i3,inzee)
            z1(2,1,i2,n3+2-i3,inzee)=-z3(2,1,n2+2-i2,i3,inzee)

            do i1=2,n1b
               z1(1,i1,i2,i3,inzee)=z3(1,i1,i2,i3,inzee)
               z1(2,i1,i2,i3,inzee)=z3(2,i1,i2,i3,inzee)

               z1(1,i1,i2,n3+2-i3,inzee)= z3(1,n1+2-i1,n2+2-i2,i3,inzee)
               z1(2,i1,i2,n3+2-i3,inzee)=-z3(2,n1+2-i1,n2+2-i2,i3,inzee)
            end do
         end do
      end do
   END SUBROUTINE z3_to_z1

   subroutine z1_to_z3(z1,z3,inzee)
      ! transforms the data from the format z1:
      ! stores the part of z with i1=<nd1b 
      ! to the format z3:
      ! stores the elements of z with even and odd values of i3
      ! as its even and odd parts w.r.t. flip of n1
      implicit none
      integer,intent(in)::inzee
      real(kind=8),intent(in):: z1(2,nd2,nd3,nd1b,2)
      real(kind=8),intent(out)::z3(2,nd2,nd3b,nd1,2)
      integer i1,i2,i3
z3=0.d0
      if (2*n3b.eq.n3) then 
         ! i1=1
         do i3=1,n3b
            do i2=1,n2
               z3(1,i2,i3,1,inzee)= z1(1,i2,2*i3-1,1,inzee)-z1(2,i2,2*i3,1,inzee)
               z3(2,i2,i3,1,inzee)= z1(2,i2,2*i3-1,1,inzee)+z1(1,i2,2*i3,1,inzee)
            end do
         end do

         do i1=2,n1b
            do i3=1,n3b
               do i2=1,n2
                  z3(1,i2,i3,i1     ,inzee)= z1(1,i2,2*i3-1,i1,inzee)-z1(2,i2,2*i3,i1,inzee)
                  z3(2,i2,i3,i1     ,inzee)= z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
                  z3(1,i2,i3,n1+2-i1,inzee)= z1(1,i2,2*i3-1,i1,inzee)+z1(2,i2,2*i3,i1,inzee)
                  z3(2,i2,i3,n1+2-i1,inzee)=-z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
               end do
            end do
         end do
      else  ! 2*n3b=n3+1
         ! i1=1
         do i3=1,n3b-1
            do i2=1,n2
               z3(1,i2,i3,1,inzee)= z1(1,i2,2*i3-1,1,inzee)-z1(2,i2,2*i3,1,inzee)
               z3(2,i2,i3,1,inzee)= z1(2,i2,2*i3-1,1,inzee)+z1(1,i2,2*i3,1,inzee)
            end do
         end do

         do i1=2,n1b
            do i3=1,n3b-1
               do i2=1,n2
                  z3(1,i2,i3,i1     ,inzee)= z1(1,i2,2*i3-1,i1,inzee)-z1(2,i2,2*i3,i1,inzee)
                  z3(2,i2,i3,i1     ,inzee)= z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
                  z3(1,i2,i3,n1+2-i1,inzee)= z1(1,i2,2*i3-1,i1,inzee)+z1(2,i2,2*i3,i1,inzee)
                  z3(2,i2,i3,n1+2-i1,inzee)=-z1(2,i2,2*i3-1,i1,inzee)+z1(1,i2,2*i3,i1,inzee)
               end do
            end do
         end do

         ! i3=n3b is treated separately: 2*n3b-1=n3, but the terms with 2*n3b are
         ! omitted

         ! i1=1
         do i2=1,n2
            z3(1,i2,n3b,1,inzee)= z1(1,i2,n3,1,inzee)
            z3(2,i2,n3b,1,inzee)= z1(2,i2,n3,1,inzee)
         end do

         do i1=2,n1b
            do i2=1,n2
               z3(1,i2,n3b,i1     ,inzee)= z1(1,i2,n3,i1,inzee)
               z3(2,i2,n3b,i1     ,inzee)= z1(2,i2,n3,i1,inzee)
               z3(1,i2,n3b,n1+2-i1,inzee)= z1(1,i2,n3,i1,inzee)
               z3(2,i2,n3b,n1+2-i1,inzee)=-z1(2,i2,n3,i1,inzee)
            end do
         end do

      end if
   END SUBROUTINE z1_to_z3

   subroutine z3_to_y(z3,y,inzee)
      ! transforms the output of FFT: z3, for which:
      ! real part of      z3 contains elements of y with odd  i3
      ! imaginary part of z3 contains elements of y with even i3

      ! into the final output: real array y.
      implicit none
      integer,intent(in)::inzee
      real(kind=8),intent(in)::z3(2,nd1,nd2,nd3b,2)
      real(kind=8),intent(out)::y(n1,n2,n3)
      integer i1,i2,i3
      real(kind=8) fac

      fac=.5d0/(n1*n2*n3)

      if (2*n3b.eq.n3) then 
         do i3=1,n3b
            do i2=1,n2
               do i1=1,n1
                  y(i1,i2,2*i3-1)=z3(1,i1,i2,i3,inzee)*fac
                  y(i1,i2,2*i3  )=z3(2,i1,i2,i3,inzee)*fac
               end do
            end do
         end do
      else ! 2*n3b=n3+1
         do i3=1,n3b-1
            do i2=1,n2
               do i1=1,n1
                  y(i1,i2,2*i3-1)=z3(1,i1,i2,i3,inzee)*fac
                  y(i1,i2,2*i3  )=z3(2,i1,i2,i3,inzee)*fac
               end do
            end do
         end do

         ! i3=n3b is treated separately: 2*n3b-1=n3, but the terms with 2*n3b are
         ! omitted

         do i2=1,n2
            do i1=1,n1
               y(i1,i2,n3)=z3(1,i1,i2,n3b,inzee)*fac
            end do
         end do
      end if
   END SUBROUTINE z3_to_y

END SUBROUTINE fft_back


!> ctrig_sg
!! @warning
!!  Different factorizations affect the performance
!!  Factoring 64 as 4*4*4 might for example be faster on some machines than 8*8.
!!
subroutine ctrig_sg(n,ntrig,trig,after,before,now,i_sign,ic)

   use module_fft_sg
   implicit real(kind=8) (a-h,o-z)

   ! Arguments
   integer :: n,i_sign,ic
   real(kind=8), dimension(2,ntrig) :: trig
   integer, dimension(n_factors) :: after,now,before
   ! Local variables

   do i=1,ndata
      if (n.eq.ij_data(1,i)) then
         ic=0
         do j=1,n_factors
            itt=ij_data(1+j,i)
            if (itt.gt.1) then
               ic=ic+1
               now(j)=ij_data(1+j,i)
            else
               goto 1000
            end if
         end do
         goto 1000
      end if
   end do
   print *,'VALUE OF',n,'NOT ALLOWED FOR FFT, ALLOWED VALUES ARE:'
   write(6,"(15(i8))") (ij_data(1,i),i=1,ndata)
   stop
   1000 continue

   after(1)=1
   before(ic)=1
   do i=2,ic
      after(i)=after(i-1)*now(i-1)
      before(ic-i+1)=before(ic-i+2)*now(ic-i+2)
   end do

   !   write(6,"(6(i3))") (after(i),i=1,ic)
   !   write(6,"(6(i3))") (now(i),i=1,ic)
   !   write(6,"(6(i3))") (before(i),i=1,ic)

   twopi=6.283185307179586d0
   angle=real(i_sign,kind=8)*twopi/real(n,kind=8)
   if (mod(n,2).eq.0) then
      nh=n/2
      trig(1,1)=1.d0
      trig(2,1)=0.d0
      trig(1,nh+1)=-1.d0
      trig(2,nh+1)=0.d0
      do i=1,nh-1
         trigc=cos(real(i,kind=8)*angle)
         trigs=sin(real(i,kind=8)*angle)
         trig(1,i+1)=trigc
         trig(2,i+1)=trigs
         trig(1,n-i+1)=trigc
         trig(2,n-i+1)=-trigs
      end do
   else
      nh=(n-1)/2
      trig(1,1)=1.d0
      trig(2,1)=0.d0
      do i=1,nh
         trigc=cos(real(i,kind=8)*angle)
         trigs=sin(real(i,kind=8)*angle)
         trig(1,i+1)=trigc
         trig(2,i+1)=trigs
         trig(1,n-i+1)=trigc
         trig(2,n-i+1)=-trigs
      end do
   end if

END SUBROUTINE ctrig_sg


!> fftstp_sg
!!
subroutine fftstp_sg(mm,nfft,m,nn,n,zin,zout,ntrig,trig,after,now,before,i_sign)

   !use module_fft_sg, except_this_one_A => fftstp_sg

   implicit real(kind=8) (a-h,o-z)

   ! Arguments
   integer :: mm,nfft,m,nn,n,after,before,i_sign
   real(kind=8), intent(in) :: trig(2,ntrig)
   real(kind=8) :: zin(2,mm,m),zout(2,nn,n)
   ! Local variables
   integer :: atn,atb
   atn=after*now
   atb=after*before

   ! sqrt(.5d0)
   rt2i=0.7071067811865475d0

   ! Radix 2
   if (now.eq.2) then
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nout1=nout1+atn
         nout2=nout1+after
         do j=1,nfft
            r1=zin(1,j,nin1)  !c1=(r1,s1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2) !c2=(r2,s2)
            s2=zin(2,j,nin2)
            zout(1,j,nout1)= r2 + r1  !c1+c2
            zout(2,j,nout1)= s2 + s1
            zout(1,j,nout2)= r1 - r2 !c1-c2
            zout(2,j,nout2)= s1 - s2
         end do
      end do

      !     Big loop
      Big_Loop: do ia=2,after
         ias=ia-1
         if (2*ias.eq.after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     zout(1,j,nout1)= r1 - r2
                     zout(2,j,nout1)= s2 + s1
                     zout(1,j,nout2)= r2 + r1
                     zout(2,j,nout2)= s1 - s2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     zout(1,j,nout1)= r2 + r1
                     zout(2,j,nout1)= s1 - s2
                     zout(1,j,nout2)= r1 - r2
                     zout(2,j,nout2)= s2 + s1
                  end do
               end do
            end if
         else if (4*ias.eq.after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r - s)*rt2i
                     s2=(r + s)*rt2i
                     zout(1,j,nout1)= r2 + r1
                     zout(2,j,nout1)= s2 + s1
                     zout(1,j,nout2)= r1 - r2
                     zout(2,j,nout2)= s1 - s2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     zout(1,j,nout1)= r2 + r1
                     zout(2,j,nout1)= s2 + s1
                     zout(1,j,nout2)= r1 - r2
                     zout(2,j,nout2)= s1 - s2
                  end do
               end do
            end if
         else if (4*ias.eq.3*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(r - s)*rt2i
                     zout(1,j,nout1)= r1 - r2
                     zout(2,j,nout1)= s2 + s1
                     zout(1,j,nout2)= r2 + r1
                     zout(2,j,nout2)= s1 - s2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(s - r)*rt2i
                     s2=(r + s)*rt2i
                     zout(1,j,nout1)= r2 + r1
                     zout(2,j,nout1)= s1 - s2
                     zout(1,j,nout2)= r1 - r2
                     zout(2,j,nout2)= s2 + s1
                  end do
               end do
            end if
         else
            itrig=ias*before+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nout1=nout1+atn
               nout2=nout1+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  zout(1,j,nout1)= r2 + r1
                  zout(2,j,nout1)= s2 + s1
                  zout(1,j,nout2)= r1 - r2
                  zout(2,j,nout2)= s1 - s2
               end do
            end do
         end if

      end do Big_Loop
      !     End of Big_loop
      ! End of if (now.eq.2) (radix 2)

      ! Radix 4
   else if (now.eq.4) then
      if (i_sign.eq.1) then 
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r=r1 + r3
               s=r2 + r4
               zout(1,j,nout1) = r + s
               zout(1,j,nout3) = r - s
               r=r1 - r3
               s=s2 - s4
               zout(1,j,nout2) = r - s 
               zout(1,j,nout4) = r + s
               r=s1 + s3
               s=s2 + s4
               zout(2,j,nout1) = r + s 
               zout(2,j,nout3) = r - s
               r=s1 - s3
               s=r2 - r4
               zout(2,j,nout2) = r + s 
               zout(2,j,nout4) = r - s
            end do
         end do
         do ia=2,after
            ias=ia-1
            if (2*ias.eq.after) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r-s)*rt2i
                     s2=(r+s)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(r + s)*rt2i
                     s4=(r - s)*rt2i
                     r=r1 - r3
                     s=r2 - r4
                     zout(1,j,nout1) = r + s
                     zout(1,j,nout3) = r - s
                     r=r1 + r3
                     s=s2 - s4
                     zout(1,j,nout2) = r - s 
                     zout(1,j,nout4) = r + s
                     r=s1 + s3
                     s=s2 + s4
                     zout(2,j,nout1) = r + s 
                     zout(2,j,nout3) = r - s
                     r=s1 - s3
                     s=r2 + r4
                     zout(2,j,nout2) = r + s 
                     zout(2,j,nout4) = r - s
                  end do
               end do
            else
               itt=ias*before
               itrig=itt+1
               cr2=trig(1,itrig)
               ci2=trig(2,itrig)
               itrig=itrig+itt
               cr3=trig(1,itrig)
               ci3=trig(2,itrig)
               itrig=itrig+itt
               cr4=trig(1,itrig)
               ci4=trig(2,itrig)
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=r*cr2 - s*ci2
                     s2=r*ci2 + s*cr2
                     r=zin(1,j,nin3)
                     s=zin(2,j,nin3)
                     r3=r*cr3 - s*ci3
                     s3=r*ci3 + s*cr3
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=r*cr4 - s*ci4
                     s4=r*ci4 + s*cr4
                     r=r1 + r3
                     s=r2 + r4
                     zout(1,j,nout1) = r + s
                     zout(1,j,nout3) = r - s
                     r=r1 - r3
                     s=s2 - s4
                     zout(1,j,nout2) = r - s 
                     zout(1,j,nout4) = r + s
                     r=s1 + s3
                     s=s2 + s4
                     zout(2,j,nout1) = r + s 
                     zout(2,j,nout3) = r - s
                     r=s1 - s3
                     s=r2 - r4
                     zout(2,j,nout2) = r + s 
                     zout(2,j,nout4) = r - s
                  end do
               end do
            end if
         end do
      else
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r=r1 + r3
               s=r2 + r4
               zout(1,j,nout1) = r + s
               zout(1,j,nout3) = r - s
               r=r1 - r3
               s=s2 - s4
               zout(1,j,nout2) = r + s
               zout(1,j,nout4) = r - s
               r=s1 + s3
               s=s2 + s4
               zout(2,j,nout1) = r + s
               zout(2,j,nout3) = r - s
               r=s1 - s3
               s=r2 - r4
               zout(2,j,nout2) = r - s
               zout(2,j,nout4) = r + s
            end do
         end do
         do ia=2,after
            ias=ia-1
            if (2*ias.eq.after) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(s - r)*rt2i
                     s4=(r + s)*rt2i
                     r=r1 + r3
                     s=r2 + r4
                     zout(1,j,nout1) = r + s
                     zout(1,j,nout3) = r - s
                     r=r1 - r3
                     s=s2 + s4
                     zout(1,j,nout2) = r + s
                     zout(1,j,nout4) = r - s
                     r=s1 - s3
                     s=s2 - s4
                     zout(2,j,nout1) = r + s
                     zout(2,j,nout3) = r - s
                     r=s1 + s3
                     s=r2 - r4
                     zout(2,j,nout2) = r - s
                     zout(2,j,nout4) = r + s
                  end do
               end do
            else
               itt=ias*before
               itrig=itt+1
               cr2=trig(1,itrig)
               ci2=trig(2,itrig)
               itrig=itrig+itt
               cr3=trig(1,itrig)
               ci3=trig(2,itrig)
               itrig=itrig+itt
               cr4=trig(1,itrig)
               ci4=trig(2,itrig)
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=r*cr2 - s*ci2
                     s2=r*ci2 + s*cr2
                     r=zin(1,j,nin3)
                     s=zin(2,j,nin3)
                     r3=r*cr3 - s*ci3
                     s3=r*ci3 + s*cr3
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=r*cr4 - s*ci4
                     s4=r*ci4 + s*cr4
                     r=r1 + r3
                     s=r2 + r4
                     zout(1,j,nout1) = r + s
                     zout(1,j,nout3) = r - s
                     r=r1 - r3
                     s=s2 - s4
                     zout(1,j,nout2) = r + s
                     zout(1,j,nout4) = r - s
                     r=s1 + s3
                     s=s2 + s4
                     zout(2,j,nout1) = r + s
                     zout(2,j,nout3) = r - s
                     r=s1 - s3
                     s=r2 - r4
                     zout(2,j,nout2) = r - s
                     zout(2,j,nout4) = r + s
                  end do
               end do
            end if
         end do
      end if
      ! End of else if (now.eq.4) (radix 4)

      ! Radix 8
   else if (now.eq.8) then
      if (i_sign.eq.-1) then 
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nin5=nin4+atb
            nin6=nin5+atb
            nin7=nin6+atb
            nin8=nin7+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            nout5=nout4+after
            nout6=nout5+after
            nout7=nout6+after
            nout8=nout7+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r5=zin(1,j,nin5)
               s5=zin(2,j,nin5)
               r6=zin(1,j,nin6)
               s6=zin(2,j,nin6)
               r7=zin(1,j,nin7)
               s7=zin(2,j,nin7)
               r8=zin(1,j,nin8)
               s8=zin(2,j,nin8)
               r=r1 + r5
               s=r3 + r7
               ap=r + s
               am=r - s
               r=r2 + r6
               s=r4 + r8
               bp=r + s
               bm=r - s
               r=s1 + s5
               s=s3 + s7
               cp=r + s
               cm=r - s
               r=s2 + s6
               s=s4 + s8
               dp1=r + s
               dm=r - s
               zout(1,j,nout1) = ap + bp
               zout(2,j,nout1) = cp + dp1
               zout(1,j,nout5) = ap - bp
               zout(2,j,nout5) = cp - dp1
               zout(1,j,nout3) = am + dm
               zout(2,j,nout3) = cm - bm
               zout(1,j,nout7) = am - dm
               zout(2,j,nout7) = cm + bm
               r=r1 - r5
               s=s3 - s7
               ap=r + s
               am=r - s
               r=s1 - s5
               s=r3 - r7
               bp=r + s
               bm=r - s
               r=s4 - s8
               s=r2 - r6
               cp=r + s
               cm=r - s
               r=s2 - s6
               s=r4 - r8
               dp1=r + s
               dm=r - s
               r = ( cp + dm)*rt2i
               s = ( dm - cp)*rt2i
               cp= ( cm + dp1)*rt2i
               dp1 = ( cm - dp1)*rt2i
               zout(1,j,nout2) = ap + r
               zout(2,j,nout2) = bm + s
               zout(1,j,nout6) = ap - r
               zout(2,j,nout6) = bm - s
               zout(1,j,nout4) = am + cp
               zout(2,j,nout4) = bp + dp1
               zout(1,j,nout8) = am - cp
               zout(2,j,nout8) = bp - dp1
            end do
         end do
         do ia=2,after
            ias=ia-1
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            itrig=itrig+itt
            cr5=trig(1,itrig)
            ci5=trig(2,itrig)
            itrig=itrig+itt
            cr6=trig(1,itrig)
            ci6=trig(2,itrig)
            itrig=itrig+itt
            cr7=trig(1,itrig)
            ci7=trig(2,itrig)
            itrig=itrig+itt
            cr8=trig(1,itrig)
            ci8=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nin4=nin3+atb
               nin5=nin4+atb
               nin6=nin5+atb
               nin7=nin6+atb
               nin8=nin7+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               nout4=nout3+after
               nout5=nout4+after
               nout6=nout5+after
               nout7=nout6+after
               nout8=nout7+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=zin(1,j,nin4)
                  s=zin(2,j,nin4)
                  r4=r*cr4 - s*ci4
                  s4=r*ci4 + s*cr4
                  r=zin(1,j,nin5)
                  s=zin(2,j,nin5)
                  r5=r*cr5 - s*ci5
                  s5=r*ci5 + s*cr5
                  r=zin(1,j,nin6)
                  s=zin(2,j,nin6)
                  r6=r*cr6 - s*ci6
                  s6=r*ci6 + s*cr6
                  r=zin(1,j,nin7)
                  s=zin(2,j,nin7)
                  r7=r*cr7 - s*ci7
                  s7=r*ci7 + s*cr7
                  r=zin(1,j,nin8)
                  s=zin(2,j,nin8)
                  r8=r*cr8 - s*ci8
                  s8=r*ci8 + s*cr8
                  r=r1 + r5
                  s=r3 + r7
                  ap=r + s
                  am=r - s
                  r=r2 + r6
                  s=r4 + r8
                  bp=r + s
                  bm=r - s
                  r=s1 + s5
                  s=s3 + s7
                  cp=r + s
                  cm=r - s
                  r=s2 + s6
                  s=s4 + s8
                  dp1=r + s
                  dm=r - s
                  zout(1,j,nout1) = ap + bp
                  zout(2,j,nout1) = cp + dp1
                  zout(1,j,nout5) = ap - bp
                  zout(2,j,nout5) = cp - dp1
                  zout(1,j,nout3) = am + dm
                  zout(2,j,nout3) = cm - bm
                  zout(1,j,nout7) = am - dm
                  zout(2,j,nout7) = cm + bm
                  r=r1 - r5
                  s=s3 - s7
                  ap=r + s
                  am=r - s
                  r=s1 - s5
                  s=r3 - r7
                  bp=r + s
                  bm=r - s
                  r=s4 - s8
                  s=r2 - r6
                  cp=r + s
                  cm=r - s
                  r=s2 - s6
                  s=r4 - r8
                  dp1=r + s
                  dm=r - s
                  r = ( cp + dm)*rt2i
                  s = ( dm - cp)*rt2i
                  cp= ( cm + dp1)*rt2i
                  dp1 = ( cm - dp1)*rt2i
                  zout(1,j,nout2) = ap + r
                  zout(2,j,nout2) = bm + s
                  zout(1,j,nout6) = ap - r
                  zout(2,j,nout6) = bm - s
                  zout(1,j,nout4) = am + cp
                  zout(2,j,nout4) = bp + dp1
                  zout(1,j,nout8) = am - cp
                  zout(2,j,nout8) = bp - dp1
               end do
            end do
         end do
      else ! else for i_sign.eq.1
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nin5=nin4+atb
            nin6=nin5+atb
            nin7=nin6+atb
            nin8=nin7+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            nout5=nout4+after
            nout6=nout5+after
            nout7=nout6+after
            nout8=nout7+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r5=zin(1,j,nin5)
               s5=zin(2,j,nin5)
               r6=zin(1,j,nin6)
               s6=zin(2,j,nin6)
               r7=zin(1,j,nin7)
               s7=zin(2,j,nin7)
               r8=zin(1,j,nin8)
               s8=zin(2,j,nin8)
               r=r1 + r5
               s=r3 + r7
               ap=r + s
               am=r - s
               r=r2 + r6
               s=r4 + r8
               bp=r + s
               bm=r - s
               r=s1 + s5
               s=s3 + s7
               cp=r + s
               cm=r - s
               r=s2 + s6
               s=s4 + s8
               dp1=r + s
               dm=r - s
               zout(1,j,nout1) = ap + bp
               zout(2,j,nout1) = cp + dp1
               zout(1,j,nout5) = ap - bp
               zout(2,j,nout5) = cp - dp1
               zout(1,j,nout3) = am - dm
               zout(2,j,nout3) = cm + bm
               zout(1,j,nout7) = am + dm
               zout(2,j,nout7) = cm - bm
               r= r1 - r5
               s=-s3 + s7
               ap=r + s
               am=r - s
               r=s1 - s5
               s=r7 - r3
               bp=r + s
               bm=r - s
               r=-s4 + s8
               s= r2 - r6
               cp=r + s
               cm=r - s
               r=-s2 + s6
               s= r4 - r8
               dp1=r + s
               dm=r - s
               r = ( cp + dm)*rt2i
               s = ( cp - dm)*rt2i
               cp= ( cm + dp1)*rt2i
               dp1= ( dp1 - cm)*rt2i
               zout(1,j,nout2) = ap + r
               zout(2,j,nout2) = bm + s
               zout(1,j,nout6) = ap - r
               zout(2,j,nout6) = bm - s
               zout(1,j,nout4) = am + cp
               zout(2,j,nout4) = bp + dp1
               zout(1,j,nout8) = am - cp
               zout(2,j,nout8) = bp - dp1
            end do
         end do
         do ia=2,after
            ias=ia-1
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            itrig=itrig+itt
            cr5=trig(1,itrig)
            ci5=trig(2,itrig)
            itrig=itrig+itt
            cr6=trig(1,itrig)
            ci6=trig(2,itrig)
            itrig=itrig+itt
            cr7=trig(1,itrig)
            ci7=trig(2,itrig)
            itrig=itrig+itt
            cr8=trig(1,itrig)
            ci8=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nin4=nin3+atb
               nin5=nin4+atb
               nin6=nin5+atb
               nin7=nin6+atb
               nin8=nin7+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               nout4=nout3+after
               nout5=nout4+after
               nout6=nout5+after
               nout7=nout6+after
               nout8=nout7+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=zin(1,j,nin4)
                  s=zin(2,j,nin4)
                  r4=r*cr4 - s*ci4
                  s4=r*ci4 + s*cr4
                  r=zin(1,j,nin5)
                  s=zin(2,j,nin5)
                  r5=r*cr5 - s*ci5
                  s5=r*ci5 + s*cr5
                  r=zin(1,j,nin6)
                  s=zin(2,j,nin6)
                  r6=r*cr6 - s*ci6
                  s6=r*ci6 + s*cr6
                  r=zin(1,j,nin7)
                  s=zin(2,j,nin7)
                  r7=r*cr7 - s*ci7
                  s7=r*ci7 + s*cr7
                  r=zin(1,j,nin8)
                  s=zin(2,j,nin8)
                  r8=r*cr8 - s*ci8
                  s8=r*ci8 + s*cr8
                  r=r1 + r5
                  s=r3 + r7
                  ap=r + s
                  am=r - s
                  r=r2 + r6
                  s=r4 + r8
                  bp=r + s
                  bm=r - s
                  r=s1 + s5
                  s=s3 + s7
                  cp=r + s
                  cm=r - s
                  r=s2 + s6
                  s=s4 + s8
                  dp1=r + s
                  dm=r - s
                  zout(1,j,nout1) = ap + bp
                  zout(2,j,nout1) = cp + dp1
                  zout(1,j,nout5) = ap - bp
                  zout(2,j,nout5) = cp - dp1
                  zout(1,j,nout3) = am - dm
                  zout(2,j,nout3) = cm + bm
                  zout(1,j,nout7) = am + dm
                  zout(2,j,nout7) = cm - bm
                  r= r1 - r5
                  s=-s3 + s7
                  ap=r + s
                  am=r - s
                  r=s1 - s5
                  s=r7 - r3
                  bp=r + s
                  bm=r - s
                  r=-s4 + s8
                  s= r2 - r6
                  cp=r + s
                  cm=r - s
                  r=-s2 + s6
                  s= r4 - r8
                  dp1=r + s
                  dm=r - s
                  r = ( cp + dm)*rt2i
                  s = ( cp - dm)*rt2i
                  cp= ( cm + dp1)*rt2i
                  dp1= ( dp1 - cm)*rt2i
                  zout(1,j,nout2) = ap + r
                  zout(2,j,nout2) = bm + s
                  zout(1,j,nout6) = ap - r
                  zout(2,j,nout6) = bm - s
                  zout(1,j,nout4) = am + cp
                  zout(2,j,nout4) = bp + dp1
                  zout(1,j,nout8) = am - cp
                  zout(2,j,nout8) = bp - dp1
               end do
            end do
         end do
      end if  !end if of i_sign
      ! end of radix 8

      ! Radix 3
   else if (now.eq.3) then
      ! .5d0*sqrt(3.d0)
      bb=real(i_sign,kind=8)*0.8660254037844387d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         do j=1,nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r=r2 + r3
            s=s2 + s3
            zout(1,j,nout1) = r + r1
            zout(2,j,nout1) = s + s1
            r1=r1 - .5d0*r
            s1=s1 - .5d0*s
            r2=bb*(r2-r3)
            s2=bb*(s2-s3)
            zout(1,j,nout2) = r1 - s2 
            zout(2,j,nout2) = s1 + r2
            zout(1,j,nout3) = r1 + s2 
            zout(2,j,nout3) = s1 - r2
         end do
      end do
      loop_3000: do ia=2,after
         ias=ia-1
         if (4*ias.eq.3*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     r3=zin(1,j,nin3)
                     s3=zin(2,j,nin3)
                     r=r3 + r2
                     s=s2 - s3
                     zout(1,j,nout1) = r1 - r
                     zout(2,j,nout1) = s + s1
                     r1=r1 + .5d0*r
                     s1=s1 - .5d0*s        
                     r2=bb*(r2-r3)        
                     s2=bb*(s2+s3)
                     zout(1,j,nout2) = r1 - s2 
                     zout(2,j,nout2) = s1 - r2
                     zout(1,j,nout3) = r1 + s2 
                     zout(2,j,nout3) = s1 + r2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     r3=zin(1,j,nin3)
                     s3=zin(2,j,nin3)
                     r=r2 - r3
                     s=s2 + s3
                     zout(1,j,nout1) = r + r1
                     zout(2,j,nout1) = s1 - s
                     r1=r1 - .5d0*r
                     s1=s1 + .5d0*s        
                     r2=bb*(r2+r3)        
                     s2=bb*(s2-s3)
                     zout(1,j,nout2) = r1 + s2 
                     zout(2,j,nout2) = s1 + r2
                     zout(1,j,nout3) = r1 - s2 
                     zout(2,j,nout3) = s1 - r2
                  end do
               end do
            end if
         else if (8*ias.eq.3*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r - s)*rt2i
                     s2=(r + s)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3) 
                     r=r2 - r3
                     s=s2 + s3
                     zout(1,j,nout1) = r + r1
                     zout(2,j,nout1) = s + s1
                     r1=r1 - .5d0*r
                     s1=s1 - .5d0*s        
                     r2=bb*(r2+r3)        
                     s2=bb*(s2-s3)
                     zout(1,j,nout2) = r1 - s2 
                     zout(2,j,nout2) = s1 + r2
                     zout(1,j,nout3) = r1 + s2 
                     zout(2,j,nout3) = s1 - r2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=r2 + r3
                     s=s2 - s3
                     zout(1,j,nout1) = r + r1
                     zout(2,j,nout1) = s + s1
                     r1=r1 - .5d0*r
                     s1=s1 - .5d0*s        
                     r2=bb*(r2-r3)        
                     s2=bb*(s2+s3)
                     zout(1,j,nout2) = r1 - s2 
                     zout(2,j,nout2) = s1 + r2
                     zout(1,j,nout3) = r1 + s2 
                     zout(2,j,nout3) = s1 - r2
                  end do
               end do
            end if
         else
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=r2 + r3
                  s=s2 + s3
                  zout(1,j,nout1) = r + r1
                  zout(2,j,nout1) = s + s1
                  r1=r1 - .5d0*r
                  s1=s1 - .5d0*s
                  r2=bb*(r2-r3)
                  s2=bb*(s2-s3)
                  zout(1,j,nout2) = r1 - s2 
                  zout(2,j,nout2) = s1 + r2
                  zout(1,j,nout3) = r1 + s2 
                  zout(2,j,nout3) = s1 - r2
               end do
            end do
         end if
      end do loop_3000
      ! End of if (now.eq.3) (radix 3)

      ! Radix 5
   else if (now.eq.5) then
      ! cos(2.d0*pi/5.d0)
      cos2=0.3090169943749474d0
      ! cos(4.d0*pi/5.d0)
      cos4=-0.8090169943749474d0
      ! sin(2.d0*pi/5.d0)
      sin2=real(i_sign,kind=8)*0.9510565162951536d0
      ! sin(4.d0*pi/5.d0)
      sin4=real(i_sign,kind=8)*0.5877852522924731d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nin4=nin3+atb
         nin5=nin4+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         nout4=nout3+after
         nout5=nout4+after
         do j=1,nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r4=zin(1,j,nin4)
            s4=zin(2,j,nin4)
            r5=zin(1,j,nin5)
            s5=zin(2,j,nin5)
            r25 = r2 + r5
            r34 = r3 + r4
            s25 = s2 - s5
            s34 = s3 - s4
            zout(1,j,nout1) = r1 + r25 + r34
            r = r1 + cos2*r25 + cos4*r34
            s = sin2*s25 + sin4*s34
            zout(1,j,nout2) = r - s
            zout(1,j,nout5) = r + s
            r = r1 + cos4*r25 + cos2*r34
            s = sin4*s25 - sin2*s34
            zout(1,j,nout3) = r - s
            zout(1,j,nout4) = r + s
            r25 = r2 - r5
            r34 = r3 - r4
            s25 = s2 + s5
            s34 = s3 + s4
            zout(2,j,nout1) = s1 + s25 + s34
            r = s1 + cos2*s25 + cos4*s34
            s = sin2*r25 + sin4*r34
            zout(2,j,nout2) = r + s
            zout(2,j,nout5) = r - s
            r = s1 + cos4*s25 + cos2*s34
            s = sin4*r25 - sin2*r34
            zout(2,j,nout3) = r + s
            zout(2,j,nout4) = r - s
         end do
      end do
      loop_5000: do ia=2,after
         ias=ia-1
         if (8*ias.eq.5*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nin5=nin4+atb        
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  nout5=nout4+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r - s)*rt2i
                     s2=(r + s)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3) 
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(r + s)*rt2i
                     s4=(r - s)*rt2i
                     r5=zin(1,j,nin5)
                     s5=zin(2,j,nin5)
                     r25 = r2 - r5
                     r34 = r3 + r4
                     s25 = s2 + s5
                     s34 = s3 - s4
                     zout(1,j,nout1) = r1 + r25 - r34
                     r = r1 + cos2*r25 - cos4*r34 
                     s = sin2*s25 + sin4*s34
                     zout(1,j,nout2) = r - s
                     zout(1,j,nout5) = r + s
                     r = r1 + cos4*r25 - cos2*r34 
                     s = sin4*s25 - sin2*s34
                     zout(1,j,nout3) = r - s
                     zout(1,j,nout4) = r + s
                     r25 = r2 + r5
                     r34 = r4 - r3
                     s25 = s2 - s5
                     s34 = s3 + s4
                     zout(2,j,nout1) = s1 + s25 + s34
                     r = s1 + cos2*s25 + cos4*s34
                     s = sin2*r25 + sin4*r34
                     zout(2,j,nout2) = r + s
                     zout(2,j,nout5) = r - s
                     r = s1 + cos4*s25 + cos2*s34
                     s = sin4*r25 - sin2*r34
                     zout(2,j,nout3) = r + s
                     zout(2,j,nout4) = r - s
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nin5=nin4+atb        
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  nout5=nout4+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(s - r)*rt2i
                     s4=(r + s)*rt2i
                     r5=zin(1,j,nin5)
                     s5=zin(2,j,nin5)
                     r25 = r2 - r5
                     r34 = r3 + r4
                     s25 = s2 + s5
                     s34 = s4 - s3
                     zout(1,j,nout1) = r1 + r25 + r34
                     r = r1 + cos2*r25 + cos4*r34
                     s = sin2*s25 + sin4*s34
                     zout(1,j,nout2) = r - s
                     zout(1,j,nout5) = r + s
                     r = r1 + cos4*r25 + cos2*r34
                     s = sin4*s25 - sin2*s34
                     zout(1,j,nout3) = r - s
                     zout(1,j,nout4) = r + s
                     r25 = r2 + r5
                     r34 = r3 - r4
                     s25 = s2 - s5
                     s34 = s3 + s4
                     zout(2,j,nout1) = s1 + s25 - s34
                     r = s1 + cos2*s25 - cos4*s34
                     s = sin2*r25 + sin4*r34
                     zout(2,j,nout2) = r + s
                     zout(2,j,nout5) = r - s
                     r = s1 + cos4*s25 - cos2*s34
                     s = sin4*r25 - sin2*r34
                     zout(2,j,nout3) = r + s
                     zout(2,j,nout4) = r - s
                  end do
               end do
            end if
         else !if of (ias...
            ias=ia-1
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            itrig=itrig+itt
            cr5=trig(1,itrig)
            ci5=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nin4=nin3+atb
               nin5=nin4+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               nout4=nout3+after
               nout5=nout4+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=zin(1,j,nin4)
                  s=zin(2,j,nin4)
                  r4=r*cr4 - s*ci4
                  s4=r*ci4 + s*cr4
                  r=zin(1,j,nin5)
                  s=zin(2,j,nin5)
                  r5=r*cr5 - s*ci5
                  s5=r*ci5 + s*cr5
                  r25 = r2 + r5
                  r34 = r3 + r4
                  s25 = s2 - s5
                  s34 = s3 - s4
                  zout(1,j,nout1) = r1 + r25 + r34
                  r = r1 + cos2*r25 + cos4*r34
                  s = sin2*s25 + sin4*s34
                  zout(1,j,nout2) = r - s
                  zout(1,j,nout5) = r + s
                  r = r1 + cos4*r25 + cos2*r34
                  s = sin4*s25 - sin2*s34
                  zout(1,j,nout3) = r - s
                  zout(1,j,nout4) = r + s
                  r25 = r2 - r5
                  r34 = r3 - r4
                  s25 = s2 + s5
                  s34 = s3 + s4
                  zout(2,j,nout1) = s1 + s25 + s34
                  r = s1 + cos2*s25 + cos4*s34
                  s = sin2*r25 + sin4*r34
                  zout(2,j,nout2) = r + s
                  zout(2,j,nout5) = r - s
                  r = s1 + cos4*s25 + cos2*s34
                  s = sin4*r25 - sin2*r34
                  zout(2,j,nout3) = r + s
                  zout(2,j,nout4) = r - s
               end do
            end do
         end if
      end do loop_5000
      ! End of if now.eq.5 (radix 5)

      ! Radix 6
   else if (now.eq.6) then
      ! .5d0*sqrt(3.d0)
      bb=real(i_sign,kind=8)*0.8660254037844387d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nin4=nin3+atb
         nin5=nin4+atb
         nin6=nin5+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         nout4=nout3+after
         nout5=nout4+after
         nout6=nout5+after
         do j=1,nfft
            r2=zin(1,j,nin3)
            s2=zin(2,j,nin3)
            r3=zin(1,j,nin5)
            s3=zin(2,j,nin5)
            r=r2 + r3
            s=s2 + s3
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            ur1 = r + r1
            ui1 = s + s1
            r1=r1 - .5d0*r
            s1=s1 - .5d0*s
            r=r2-r3
            s=s2-s3
            ur2 = r1 - s*bb
            ui2 = s1 + r*bb
            ur3 = r1 + s*bb
            ui3 = s1 - r*bb

            r2=zin(1,j,nin6)
            s2=zin(2,j,nin6)
            r3=zin(1,j,nin2)
            s3=zin(2,j,nin2)
            r=r2 + r3
            s=s2 + s3
            r1=zin(1,j,nin4)
            s1=zin(2,j,nin4)
            vr1 = r + r1
            vi1 = s + s1
            r1=r1 - .5d0*r
            s1=s1 - .5d0*s
            r=r2-r3
            s=s2-s3
            vr2 = r1 - s*bb
            vi2 = s1 + r*bb
            vr3 = r1 + s*bb
            vi3 = s1 - r*bb

            zout(1,j,nout1)=ur1+vr1
            zout(2,j,nout1)=ui1+vi1
            zout(1,j,nout5)=ur2+vr2
            zout(2,j,nout5)=ui2+vi2
            zout(1,j,nout3)=ur3+vr3
            zout(2,j,nout3)=ui3+vi3
            zout(1,j,nout4)=ur1-vr1
            zout(2,j,nout4)=ui1-vi1
            zout(1,j,nout2)=ur2-vr2
            zout(2,j,nout2)=ui2-vi2
            zout(1,j,nout6)=ur3-vr3
            zout(2,j,nout6)=ui3-vi3
         end do
      end do
      ! End of radix 6

      !Radix 7
   else if (now.eq.7) then
      ia=1
      nin1=ia-after
      nout1=ia-atn

      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nin4=nin3+atb
         nin5=nin4+atb
         nin6=nin5+atb
         nin7=nin6+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         nout4=nout3+after
         nout5=nout4+after
         nout6=nout5+after
         nout7=nout6+after
         do j=1,nfft

            cos2= 0.62348980185873353053d0
            cos4=-0.22252093395631440429d0
            cos6=-0.90096886790241912624d0

            sin2= 0.78183148246802980871d0*real(i_sign,kind=8)
            sin4= 0.97492791218182360702d0*real(i_sign,kind=8)
            sin6= 0.43388373911755812048d0*real(i_sign,kind=8)

            r1=zin(1,j,nin1);   s1=zin(2,j,nin1)
            r2=zin(1,j,nin2);   s2=zin(2,j,nin2)
            r3=zin(1,j,nin3);   s3=zin(2,j,nin3)
            r4=zin(1,j,nin4);   s4=zin(2,j,nin4)
            r5=zin(1,j,nin5);   s5=zin(2,j,nin5)
            r6=zin(1,j,nin6);   s6=zin(2,j,nin6)
            r7=zin(1,j,nin7);   s7=zin(2,j,nin7)

            rp2=r2+r7; rp3=r3+r6; rp4=r4+r5
            rm2=r2-r7; rm3=r3-r6; rm4=r4-r5
            sp2=s2+s7; sp3=s3+s6; sp4=s4+s5
            sm2=s2-s7; sm3=s3-s6; sm4=s4-s5 !3*4=12 flops

            ur2=r1+cos2*rp2+cos4*rp3+cos6*rp4
            vr2=   sin2*sm2+sin4*sm3+sin6*sm4

            ui2=s1+cos2*sp2+cos4*sp3+cos6*sp4
            vi2=   sin2*rm2+sin4*rm3+sin6*rm4

            ur3=r1+cos4*rp2+cos6*rp3+cos2*rp4
            vr3=   sin4*sm2-sin6*sm3-sin2*sm4

            ui3=s1+cos4*sp2+cos6*sp3+cos2*sp4
            vi3=   sin4*rm2-sin6*rm3-sin2*rm4

            ur4=r1+cos6*rp2+cos2*rp3+cos4*rp4
            vr4=   sin6*sm2-sin2*sm3+sin4*sm4

            ui4=s1+cos6*sp2+cos2*sp3+cos4*sp4
            vi4=   sin6*rm2-sin2*rm3+sin4*rm4   ! 11*6=66 flops, of them
            !6*6=36 multiplications

            zout(1,j,nout1)=r1+rp2+rp3+rp4; zout(2,j,nout1)=s1+sp2+sp3+sp4
            zout(1,j,nout2)=ur2-vr2;        zout(2,j,nout2)=ui2+vi2
            zout(1,j,nout3)=ur3-vr3;        zout(2,j,nout3)=ui3+vi3
            zout(1,j,nout4)=ur4-vr4;        zout(2,j,nout4)=ui4+vi4
            zout(1,j,nout5)=ur4+vr4;        zout(2,j,nout5)=ui4-vi4
            zout(1,j,nout6)=ur3+vr3;        zout(2,j,nout6)=ui3-vi3
            zout(1,j,nout7)=ur2+vr2;        zout(2,j,nout7)=ui2-vi2 
            !(3+6)*2 =18 flops
         end do

         ! In total: 12+66+18=96flops
      end do
      ! End of radix 7

   else 
      write(*,'(a,i6)') 'Error fftstp_sg: radix not defined ',now
      stop
   end if !end of now

END SUBROUTINE fftstp_sg


!> fftrot_sg
!!
subroutine fftrot_sg(mm,nfft,m,nn,n,zin,zout,ntrig,trig,after,now,before,i_sign)

   !n(c) use module_fft_sg
   implicit real(kind=8) (a-h,o-z)

   integer :: after,before,atn,atb
   dimension trig(2,ntrig),zin(2,mm,m),zout(2,n,nn)

   atn=after*now
   atb=after*before

   !  sqrt(.5d0)
   rt2i=0.7071067811865475d0

   ! Radix 2
   if (now.eq.2) then
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nout1=nout1+atn
         nout2=nout1+after
         do j=1,nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            zout(1,nout1,j)= r2 + r1
            zout(2,nout1,j)= s2 + s1
            zout(1,nout2,j)= r1 - r2
            zout(2,nout2,j)= s1 - s2
         end do
      end do
      loop_2000: do ia=2,after
         ias=ia-1
         if (2*ias.eq.after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     zout(1,nout1,j)= r1 - r2
                     zout(2,nout1,j)= s2 + s1
                     zout(1,nout2,j)= r2 + r1
                     zout(2,nout2,j)= s1 - s2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     zout(1,nout1,j)= r2 + r1
                     zout(2,nout1,j)= s1 - s2
                     zout(1,nout2,j)= r1 - r2
                     zout(2,nout2,j)= s2 + s1
                  end do
               end do
            end if
         else if (4*ias.eq.after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r - s)*rt2i
                     s2=(r + s)*rt2i
                     zout(1,nout1,j)= r2 + r1
                     zout(2,nout1,j)= s2 + s1
                     zout(1,nout2,j)= r1 - r2
                     zout(2,nout2,j)= s1 - s2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     zout(1,nout1,j)= r2 + r1
                     zout(2,nout1,j)= s2 + s1
                     zout(1,nout2,j)= r1 - r2
                     zout(2,nout2,j)= s1 - s2
                  end do
               end do
            end if
         else if (4*ias.eq.3*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(r - s)*rt2i
                     zout(1,nout1,j)= r1 - r2
                     zout(2,nout1,j)= s2 + s1
                     zout(1,nout2,j)= r2 + r1
                     zout(2,nout2,j)= s1 - s2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(s - r)*rt2i
                     s2=(r + s)*rt2i
                     zout(1,nout1,j)= r2 + r1
                     zout(2,nout1,j)= s1 - s2
                     zout(1,nout2,j)= r1 - r2
                     zout(2,nout2,j)= s2 + s1
                  end do
               end do
            end if
         else
            itrig=ias*before+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nout1=nout1+atn
               nout2=nout1+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  zout(1,nout1,j)= r2 + r1
                  zout(2,nout1,j)= s2 + s1
                  zout(1,nout2,j)= r1 - r2
                  zout(2,nout2,j)= s1 - s2
               end do
            end do
         end if
      end do loop_2000
      !End of radix 2

      ! Radix 4
   else if (now.eq.4) then
      if (i_sign.eq.1) then 
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r=r1 + r3
               s=r2 + r4
               zout(1,nout1,j) = r + s
               zout(1,nout3,j) = r - s
               r=r1 - r3
               s=s2 - s4
               zout(1,nout2,j) = r - s 
               zout(1,nout4,j) = r + s
               r=s1 + s3
               s=s2 + s4
               zout(2,nout1,j) = r + s 
               zout(2,nout3,j) = r - s
               r=s1 - s3
               s=r2 - r4
               zout(2,nout2,j) = r + s 
               zout(2,nout4,j) = r - s
            end do
         end do
         loop_4000: do ia=2,after
            ias=ia-1
            if (2*ias.eq.after) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r-s)*rt2i
                     s2=(r+s)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(r + s)*rt2i
                     s4=(r - s)*rt2i
                     r=r1 - r3
                     s=r2 - r4
                     zout(1,nout1,j) = r + s
                     zout(1,nout3,j) = r - s
                     r=r1 + r3
                     s=s2 - s4
                     zout(1,nout2,j) = r - s 
                     zout(1,nout4,j) = r + s
                     r=s1 + s3
                     s=s2 + s4
                     zout(2,nout1,j) = r + s 
                     zout(2,nout3,j) = r - s
                     r=s1 - s3
                     s=r2 + r4
                     zout(2,nout2,j) = r + s 
                     zout(2,nout4,j) = r - s
                  end do
               end do
            else
               itt=ias*before
               itrig=itt+1
               cr2=trig(1,itrig)
               ci2=trig(2,itrig)
               itrig=itrig+itt
               cr3=trig(1,itrig)
               ci3=trig(2,itrig)
               itrig=itrig+itt
               cr4=trig(1,itrig)
               ci4=trig(2,itrig)
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=r*cr2 - s*ci2
                     s2=r*ci2 + s*cr2
                     r=zin(1,j,nin3)
                     s=zin(2,j,nin3)
                     r3=r*cr3 - s*ci3
                     s3=r*ci3 + s*cr3
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=r*cr4 - s*ci4
                     s4=r*ci4 + s*cr4
                     r=r1 + r3
                     s=r2 + r4
                     zout(1,nout1,j) = r + s
                     zout(1,nout3,j) = r - s
                     r=r1 - r3
                     s=s2 - s4
                     zout(1,nout2,j) = r - s 
                     zout(1,nout4,j) = r + s
                     r=s1 + s3
                     s=s2 + s4
                     zout(2,nout1,j) = r + s 
                     zout(2,nout3,j) = r - s
                     r=s1 - s3
                     s=r2 - r4
                     zout(2,nout2,j) = r + s 
                     zout(2,nout4,j) = r - s
                  end do
               end do
            end if
         end do loop_4000

      else
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r=r1 + r3
               s=r2 + r4
               zout(1,nout1,j) = r + s
               zout(1,nout3,j) = r - s
               r=r1 - r3
               s=s2 - s4
               zout(1,nout2,j) = r + s
               zout(1,nout4,j) = r - s
               r=s1 + s3
               s=s2 + s4
               zout(2,nout1,j) = r + s
               zout(2,nout3,j) = r - s
               r=s1 - s3
               s=r2 - r4
               zout(2,nout2,j) = r - s
               zout(2,nout4,j) = r + s
            end do
         end do
         loop_4100: do ia=2,after
            ias=ia-1
            if (2*ias.eq.after) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(s - r)*rt2i
                     s4=(r + s)*rt2i
                     r=r1 + r3
                     s=r2 + r4
                     zout(1,nout1,j) = r + s
                     zout(1,nout3,j) = r - s
                     r=r1 - r3
                     s=s2 + s4
                     zout(1,nout2,j) = r + s
                     zout(1,nout4,j) = r - s
                     r=s1 - s3
                     s=s2 - s4
                     zout(2,nout1,j) = r + s
                     zout(2,nout3,j) = r - s
                     r=s1 + s3
                     s=r2 - r4
                     zout(2,nout2,j) = r - s
                     zout(2,nout4,j) = r + s
                  end do
               end do
            else
               itt=ias*before
               itrig=itt+1
               cr2=trig(1,itrig)
               ci2=trig(2,itrig)
               itrig=itrig+itt
               cr3=trig(1,itrig)
               ci3=trig(2,itrig)
               itrig=itrig+itt
               cr4=trig(1,itrig)
               ci4=trig(2,itrig)
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=r*cr2 - s*ci2
                     s2=r*ci2 + s*cr2
                     r=zin(1,j,nin3)
                     s=zin(2,j,nin3)
                     r3=r*cr3 - s*ci3
                     s3=r*ci3 + s*cr3
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=r*cr4 - s*ci4
                     s4=r*ci4 + s*cr4
                     r=r1 + r3
                     s=r2 + r4
                     zout(1,nout1,j) = r + s
                     zout(1,nout3,j) = r - s
                     r=r1 - r3
                     s=s2 - s4
                     zout(1,nout2,j) = r + s
                     zout(1,nout4,j) = r - s
                     r=s1 + s3
                     s=s2 + s4
                     zout(2,nout1,j) = r + s
                     zout(2,nout3,j) = r - s
                     r=s1 - s3
                     s=r2 - r4
                     zout(2,nout2,j) = r - s
                     zout(2,nout4,j) = r + s
                  end do
               end do
            end if
         end do loop_4100
      end if
      !End of radix 4

      ! Radix 8
   else if (now.eq.8) then
      if (i_sign.eq.-1) then 
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nin5=nin4+atb
            nin6=nin5+atb
            nin7=nin6+atb
            nin8=nin7+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            nout5=nout4+after
            nout6=nout5+after
            nout7=nout6+after
            nout8=nout7+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r5=zin(1,j,nin5)
               s5=zin(2,j,nin5)
               r6=zin(1,j,nin6)
               s6=zin(2,j,nin6)
               r7=zin(1,j,nin7)
               s7=zin(2,j,nin7)
               r8=zin(1,j,nin8)
               s8=zin(2,j,nin8)
               r=r1 + r5
               s=r3 + r7
               ap=r + s
               am=r - s
               r=r2 + r6
               s=r4 + r8
               bp=r + s
               bm=r - s
               r=s1 + s5
               s=s3 + s7
               cp=r + s
               cm=r - s
               r=s2 + s6
               s=s4 + s8
               dp1=r + s
               dm=r - s
               zout(1,nout1,j) = ap + bp
               zout(2,nout1,j) = cp + dp1
               zout(1,nout5,j) = ap - bp
               zout(2,nout5,j) = cp - dp1
               zout(1,nout3,j) = am + dm
               zout(2,nout3,j) = cm - bm
               zout(1,nout7,j) = am - dm
               zout(2,nout7,j) = cm + bm
               r=r1 - r5
               s=s3 - s7
               ap=r + s
               am=r - s
               r=s1 - s5
               s=r3 - r7
               bp=r + s
               bm=r - s
               r=s4 - s8
               s=r2 - r6
               cp=r + s
               cm=r - s
               r=s2 - s6
               s=r4 - r8
               dp1=r + s
               dm=r - s
               r = ( cp + dm)*rt2i
               s = ( dm - cp)*rt2i
               cp= ( cm + dp1)*rt2i
               dp1= ( cm - dp1)*rt2i
               zout(1,nout2,j) = ap + r
               zout(2,nout2,j) = bm + s
               zout(1,nout6,j) = ap - r
               zout(2,nout6,j) = bm - s
               zout(1,nout4,j) = am + cp
               zout(2,nout4,j) = bp + dp1
               zout(1,nout8,j) = am - cp
               zout(2,nout8,j) = bp - dp1
            end do
         end do

         do ia=2,after
            ias=ia-1
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            itrig=itrig+itt
            cr5=trig(1,itrig)
            ci5=trig(2,itrig)
            itrig=itrig+itt
            cr6=trig(1,itrig)
            ci6=trig(2,itrig)
            itrig=itrig+itt
            cr7=trig(1,itrig)
            ci7=trig(2,itrig)
            itrig=itrig+itt
            cr8=trig(1,itrig)
            ci8=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nin4=nin3+atb
               nin5=nin4+atb
               nin6=nin5+atb
               nin7=nin6+atb
               nin8=nin7+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               nout4=nout3+after
               nout5=nout4+after
               nout6=nout5+after
               nout7=nout6+after
               nout8=nout7+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=zin(1,j,nin4)
                  s=zin(2,j,nin4)
                  r4=r*cr4 - s*ci4
                  s4=r*ci4 + s*cr4
                  r=zin(1,j,nin5)
                  s=zin(2,j,nin5)
                  r5=r*cr5 - s*ci5
                  s5=r*ci5 + s*cr5
                  r=zin(1,j,nin6)
                  s=zin(2,j,nin6)
                  r6=r*cr6 - s*ci6
                  s6=r*ci6 + s*cr6
                  r=zin(1,j,nin7)
                  s=zin(2,j,nin7)
                  r7=r*cr7 - s*ci7
                  s7=r*ci7 + s*cr7
                  r=zin(1,j,nin8)
                  s=zin(2,j,nin8)
                  r8=r*cr8 - s*ci8
                  s8=r*ci8 + s*cr8
                  r=r1 + r5
                  s=r3 + r7
                  ap=r + s
                  am=r - s
                  r=r2 + r6
                  s=r4 + r8
                  bp=r + s
                  bm=r - s
                  r=s1 + s5
                  s=s3 + s7
                  cp=r + s
                  cm=r - s
                  r=s2 + s6
                  s=s4 + s8
                  dp1=r + s
                  dm=r - s
                  zout(1,nout1,j) = ap + bp
                  zout(2,nout1,j) = cp + dp1
                  zout(1,nout5,j) = ap - bp
                  zout(2,nout5,j) = cp - dp1
                  zout(1,nout3,j) = am + dm
                  zout(2,nout3,j) = cm - bm
                  zout(1,nout7,j) = am - dm
                  zout(2,nout7,j) = cm + bm
                  r=r1 - r5
                  s=s3 - s7
                  ap=r + s
                  am=r - s
                  r=s1 - s5
                  s=r3 - r7
                  bp=r + s
                  bm=r - s
                  r=s4 - s8
                  s=r2 - r6
                  cp=r + s
                  cm=r - s
                  r=s2 - s6
                  s=r4 - r8
                  dp1=r + s
                  dm=r - s
                  r = ( cp + dm)*rt2i
                  s = ( dm - cp)*rt2i
                  cp= ( cm + dp1)*rt2i
                  dp1= ( cm - dp1)*rt2i
                  zout(1,nout2,j) = ap + r
                  zout(2,nout2,j) = bm + s
                  zout(1,nout6,j) = ap - r
                  zout(2,nout6,j) = bm - s
                  zout(1,nout4,j) = am + cp
                  zout(2,nout4,j) = bp + dp1
                  zout(1,nout8,j) = am - cp
                  zout(2,nout8,j) = bp - dp1
               end do
            end do
         end do

      else !i_sgn == 1
         ia=1
         nin1=ia-after
         nout1=ia-atn
         do ib=1,before
            nin1=nin1+after
            nin2=nin1+atb
            nin3=nin2+atb
            nin4=nin3+atb
            nin5=nin4+atb
            nin6=nin5+atb
            nin7=nin6+atb
            nin8=nin7+atb
            nout1=nout1+atn
            nout2=nout1+after
            nout3=nout2+after
            nout4=nout3+after
            nout5=nout4+after
            nout6=nout5+after
            nout7=nout6+after
            nout8=nout7+after
            do j=1,nfft
               r1=zin(1,j,nin1)
               s1=zin(2,j,nin1)
               r2=zin(1,j,nin2)
               s2=zin(2,j,nin2)
               r3=zin(1,j,nin3)
               s3=zin(2,j,nin3)
               r4=zin(1,j,nin4)
               s4=zin(2,j,nin4)
               r5=zin(1,j,nin5)
               s5=zin(2,j,nin5)
               r6=zin(1,j,nin6)
               s6=zin(2,j,nin6)
               r7=zin(1,j,nin7)
               s7=zin(2,j,nin7)
               r8=zin(1,j,nin8)
               s8=zin(2,j,nin8)
               r=r1 + r5
               s=r3 + r7
               ap=r + s
               am=r - s
               r=r2 + r6
               s=r4 + r8
               bp=r + s
               bm=r - s
               r=s1 + s5
               s=s3 + s7
               cp=r + s
               cm=r - s
               r=s2 + s6
               s=s4 + s8
               dp1=r + s
               dm=r - s
               zout(1,nout1,j) = ap + bp
               zout(2,nout1,j) = cp + dp1
               zout(1,nout5,j) = ap - bp
               zout(2,nout5,j) = cp - dp1
               zout(1,nout3,j) = am - dm
               zout(2,nout3,j) = cm + bm
               zout(1,nout7,j) = am + dm
               zout(2,nout7,j) = cm - bm
               r= r1 - r5
               s=-s3 + s7
               ap=r + s
               am=r - s
               r=s1 - s5
               s=r7 - r3
               bp=r + s
               bm=r - s
               r=-s4 + s8
               s= r2 - r6
               cp=r + s
               cm=r - s
               r=-s2 + s6
               s= r4 - r8
               dp1=r + s
               dm=r - s
               r = ( cp + dm)*rt2i
               s = ( cp - dm)*rt2i
               cp= ( cm + dp1)*rt2i
               dp1= ( dp1 - cm)*rt2i
               zout(1,nout2,j) = ap + r
               zout(2,nout2,j) = bm + s
               zout(1,nout6,j) = ap - r
               zout(2,nout6,j) = bm - s
               zout(1,nout4,j) = am + cp
               zout(2,nout4,j) = bp + dp1
               zout(1,nout8,j) = am - cp
               zout(2,nout8,j) = bp - dp1
            end do
         end do

         do ia=2,after
            ias=ia-1
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            itrig=itrig+itt
            cr5=trig(1,itrig)
            ci5=trig(2,itrig)
            itrig=itrig+itt
            cr6=trig(1,itrig)
            ci6=trig(2,itrig)
            itrig=itrig+itt
            cr7=trig(1,itrig)
            ci7=trig(2,itrig)
            itrig=itrig+itt
            cr8=trig(1,itrig)
            ci8=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nin4=nin3+atb
               nin5=nin4+atb
               nin6=nin5+atb
               nin7=nin6+atb
               nin8=nin7+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               nout4=nout3+after
               nout5=nout4+after
               nout6=nout5+after
               nout7=nout6+after
               nout8=nout7+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=zin(1,j,nin4)
                  s=zin(2,j,nin4)
                  r4=r*cr4 - s*ci4
                  s4=r*ci4 + s*cr4
                  r=zin(1,j,nin5)
                  s=zin(2,j,nin5)
                  r5=r*cr5 - s*ci5
                  s5=r*ci5 + s*cr5
                  r=zin(1,j,nin6)
                  s=zin(2,j,nin6)
                  r6=r*cr6 - s*ci6
                  s6=r*ci6 + s*cr6
                  r=zin(1,j,nin7)
                  s=zin(2,j,nin7)
                  r7=r*cr7 - s*ci7
                  s7=r*ci7 + s*cr7
                  r=zin(1,j,nin8)
                  s=zin(2,j,nin8)
                  r8=r*cr8 - s*ci8
                  s8=r*ci8 + s*cr8
                  r=r1 + r5
                  s=r3 + r7
                  ap=r + s
                  am=r - s
                  r=r2 + r6
                  s=r4 + r8
                  bp=r + s
                  bm=r - s
                  r=s1 + s5
                  s=s3 + s7
                  cp=r + s
                  cm=r - s
                  r=s2 + s6
                  s=s4 + s8
                  dp1=r + s
                  dm=r - s
                  zout(1,nout1,j) = ap + bp
                  zout(2,nout1,j) = cp + dp1
                  zout(1,nout5,j) = ap - bp
                  zout(2,nout5,j) = cp - dp1
                  zout(1,nout3,j) = am - dm
                  zout(2,nout3,j) = cm + bm
                  zout(1,nout7,j) = am + dm
                  zout(2,nout7,j) = cm - bm
                  r= r1 - r5
                  s=-s3 + s7
                  ap=r + s
                  am=r - s
                  r=s1 - s5
                  s=r7 - r3
                  bp=r + s
                  bm=r - s
                  r=-s4 + s8
                  s= r2 - r6
                  cp=r + s
                  cm=r - s
                  r=-s2 + s6
                  s= r4 - r8
                  dp1=r + s
                  dm=r - s
                  r = ( cp + dm)*rt2i
                  s = ( cp - dm)*rt2i
                  cp= ( cm + dp1)*rt2i
                  dp1= ( dp1 - cm)*rt2i
                  zout(1,nout2,j) = ap + r
                  zout(2,nout2,j) = bm + s
                  zout(1,nout6,j) = ap - r
                  zout(2,nout6,j) = bm - s
                  zout(1,nout4,j) = am + cp
                  zout(2,nout4,j) = bp + dp1
                  zout(1,nout8,j) = am - cp
                  zout(2,nout8,j) = bp - dp1
               end do
            end do
         end do !ia=2,after

      end if !i_sign
      !End of radix 8

      ! Radix 3
   else if (now.eq.3) then 
      ! .5d0*sqrt(3.d0)
      bb=i_sign*0.8660254037844387d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         do j=1,nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r=r2 + r3
            s=s2 + s3
            zout(1,nout1,j) = r + r1
            zout(2,nout1,j) = s + s1
            r1=r1 - .5d0*r
            s1=s1 - .5d0*s
            r2=bb*(r2-r3)
            s2=bb*(s2-s3)
            zout(1,nout2,j) = r1 - s2 
            zout(2,nout2,j) = s1 + r2
            zout(1,nout3,j) = r1 + s2 
            zout(2,nout3,j) = s1 - r2
         end do
      end do

      do ia=2,after
         ias=ia-1
         if (4*ias.eq.3*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     r3=zin(1,j,nin3)
                     s3=zin(2,j,nin3)
                     r=r2 + r3
                     s=s2 - s3
                     zout(1,nout1,j) = r1 - r
                     zout(2,nout1,j) = s + s1
                     r1=r1 + .5d0*r
                     s1=s1 - .5d0*s        
                     r2=bb*(r2-r3)        
                     s2=bb*(s2+s3)
                     zout(1,nout2,j) = r1 - s2 
                     zout(2,nout2,j) = s1 - r2
                     zout(1,nout3,j) = r1 + s2 
                     zout(2,nout3,j) = s1 + r2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r2=zin(2,j,nin2)
                     s2=zin(1,j,nin2)
                     r3=zin(1,j,nin3)
                     s3=zin(2,j,nin3)
                     r=r2 - r3
                     s=s2 + s3
                     zout(1,nout1,j) = r + r1
                     zout(2,nout1,j) = s1 - s
                     r1=r1 - .5d0*r
                     s1=s1 + .5d0*s        
                     r2=bb*(r2+r3)        
                     s2=bb*(s2-s3)
                     zout(1,nout2,j) = r1 + s2 
                     zout(2,nout2,j) = s1 + r2
                     zout(1,nout3,j) = r1 - s2 
                     zout(2,nout3,j) = s1 - r2
                  end do
               end do
            end if

         else if (8*ias.eq.3*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r - s)*rt2i
                     s2=(r + s)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3) 
                     r=r2 - r3
                     s=s2 + s3
                     zout(1,nout1,j) = r + r1
                     zout(2,nout1,j) = s + s1
                     r1=r1 - .5d0*r
                     s1=s1 - .5d0*s        
                     r2=bb*(r2+r3)        
                     s2=bb*(s2-s3)
                     zout(1,nout2,j) = r1 - s2 
                     zout(2,nout2,j) = s1 + r2
                     zout(1,nout3,j) = r1 + s2 
                     zout(2,nout3,j) = s1 - r2
                  end do
               end do
            else
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=r2 + r3
                     s=s2 - s3
                     zout(1,nout1,j) = r + r1
                     zout(2,nout1,j) = s + s1
                     r1=r1 - .5d0*r
                     s1=s1 - .5d0*s        
                     r2=bb*(r2-r3)        
                     s2=bb*(s2+s3)
                     zout(1,nout2,j) = r1 - s2 
                     zout(2,nout2,j) = s1 + r2
                     zout(1,nout3,j) = r1 + s2 
                     zout(2,nout3,j) = s1 - r2
                  end do
               end do
            end if
         else !ia.eq...
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=r2 + r3
                  s=s2 + s3
                  zout(1,nout1,j) = r + r1
                  zout(2,nout1,j) = s + s1
                  r1=r1 - .5d0*r
                  s1=s1 - .5d0*s
                  r2=bb*(r2-r3)
                  s2=bb*(s2-s3)
                  zout(1,nout2,j) = r1 - s2 
                  zout(2,nout2,j) = s1 + r2
                  zout(1,nout3,j) = r1 + s2 
                  zout(2,nout3,j) = s1 - r2
               end do
            end do
         end if
      end do !ia=2,after
      ! End of radix 3

      ! Radix 5
   else if (now.eq.5) then
      ! cos(2.d0*pi/5.d0)
      cos2=0.3090169943749474d0
      ! cos(4.d0*pi/5.d0)
      cos4=-0.8090169943749474d0
      ! sin(2.d0*pi/5.d0)
      sin2=i_sign*0.9510565162951536d0
      ! sin(4.d0*pi/5.d0)
      sin4=i_sign*0.5877852522924731d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nin4=nin3+atb
         nin5=nin4+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         nout4=nout3+after
         nout5=nout4+after
         do j=1,nfft
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            r2=zin(1,j,nin2)
            s2=zin(2,j,nin2)
            r3=zin(1,j,nin3)
            s3=zin(2,j,nin3)
            r4=zin(1,j,nin4)
            s4=zin(2,j,nin4)
            r5=zin(1,j,nin5)
            s5=zin(2,j,nin5)
            r25 = r2 + r5
            r34 = r3 + r4
            s25 = s2 - s5
            s34 = s3 - s4
            zout(1,nout1,j) = r1 + r25 + r34
            r = r1 + cos2*r25 + cos4*r34
            s = sin2*s25 + sin4*s34
            zout(1,nout2,j) = r - s
            zout(1,nout5,j) = r + s
            r = r1 + cos4*r25 + cos2*r34
            s = sin4*s25 - sin2*s34
            zout(1,nout3,j) = r - s
            zout(1,nout4,j) = r + s
            r25 = r2 - r5
            r34 = r3 - r4
            s25 = s2 + s5
            s34 = s3 + s4
            zout(2,nout1,j) = s1 + s25 + s34
            r = s1 + cos2*s25 + cos4*s34
            s = sin2*r25 + sin4*r34
            zout(2,nout2,j) = r + s
            zout(2,nout5,j) = r - s
            r = s1 + cos4*s25 + cos2*s34
            s = sin4*r25 - sin2*r34
            zout(2,nout3,j) = r + s
            zout(2,nout4,j) = r - s
         end do
      end do
      loop_5000: do ia=2,after
         ias=ia-1
         if (8*ias.eq.5*after) then
            if (i_sign.eq.1) then
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nin5=nin4+atb        
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  nout5=nout4+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r - s)*rt2i
                     s2=(r + s)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3) 
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(r + s)*rt2i
                     s4=(r - s)*rt2i
                     r5=zin(1,j,nin5)
                     s5=zin(2,j,nin5)
                     r25 = r2 - r5
                     r34 = r3 + r4
                     s25 = s2 + s5
                     s34 = s3 - s4
                     zout(1,nout1,j) = r1 + r25 - r34
                     r = r1 + cos2*r25 - cos4*r34
                     s = sin2*s25 + sin4*s34
                     zout(1,nout2,j) = r - s
                     zout(1,nout5,j) = r + s
                     r = r1 + cos4*r25 - cos2*r34
                     s = sin4*s25 - sin2*s34
                     zout(1,nout3,j) = r - s
                     zout(1,nout4,j) = r + s
                     r25 = r2 + r5
                     r34 = r4 - r3
                     s25 = s2 - s5
                     s34 = s3 + s4
                     zout(2,nout1,j) = s1 + s25 + s34
                     r = s1 + cos2*s25 + cos4*s34
                     s = sin2*r25 + sin4*r34
                     zout(2,nout2,j) = r + s
                     zout(2,nout5,j) = r - s
                     r = s1 + cos4*s25 + cos2*s34
                     s = sin4*r25 - sin2*r34
                     zout(2,nout3,j) = r + s
                     zout(2,nout4,j) = r - s
                  end do
               end do
            else !i_sign.eq.-1
               nin1=ia-after
               nout1=ia-atn
               do ib=1,before
                  nin1=nin1+after
                  nin2=nin1+atb
                  nin3=nin2+atb
                  nin4=nin3+atb
                  nin5=nin4+atb        
                  nout1=nout1+atn
                  nout2=nout1+after
                  nout3=nout2+after
                  nout4=nout3+after
                  nout5=nout4+after
                  do j=1,nfft
                     r1=zin(1,j,nin1)
                     s1=zin(2,j,nin1)
                     r=zin(1,j,nin2)
                     s=zin(2,j,nin2)
                     r2=(r + s)*rt2i
                     s2=(s - r)*rt2i
                     r3=zin(2,j,nin3)
                     s3=zin(1,j,nin3)
                     r=zin(1,j,nin4)
                     s=zin(2,j,nin4)
                     r4=(s - r)*rt2i
                     s4=(r + s)*rt2i
                     r5=zin(1,j,nin5)
                     s5=zin(2,j,nin5)
                     r25 = r2 - r5
                     r34 = r3 + r4
                     s25 = s2 + s5
                     s34 = s4 - s3
                     zout(1,nout1,j) = r1 + r25 + r34
                     r = r1 + cos2*r25 + cos4*r34
                     s = sin2*s25 + sin4*s34
                     zout(1,nout2,j) = r - s
                     zout(1,nout5,j) = r + s
                     r = r1 + cos4*r25 + cos2*r34
                     s = sin4*s25 - sin2*s34
                     zout(1,nout3,j) = r - s
                     zout(1,nout4,j) = r + s
                     r25 = r2 + r5
                     r34 = r3 - r4
                     s25 = s2 - s5
                     s34 = s3 + s4
                     zout(2,nout1,j) = s1 + s25 - s34
                     r = s1 + cos2*s25 - cos4*s34
                     s = sin2*r25 + sin4*r34
                     zout(2,nout2,j) = r + s
                     zout(2,nout5,j) = r - s
                     r = s1 + cos4*s25 - cos2*s34
                     s = sin4*r25 - sin2*r34
                     zout(2,nout3,j) = r + s
                     zout(2,nout4,j) = r - s
                  end do
               end do
            end if
         else
            ias=ia-1
            itt=ias*before
            itrig=itt+1
            cr2=trig(1,itrig)
            ci2=trig(2,itrig)
            itrig=itrig+itt
            cr3=trig(1,itrig)
            ci3=trig(2,itrig)
            itrig=itrig+itt
            cr4=trig(1,itrig)
            ci4=trig(2,itrig)
            itrig=itrig+itt
            cr5=trig(1,itrig)
            ci5=trig(2,itrig)
            nin1=ia-after
            nout1=ia-atn
            do ib=1,before
               nin1=nin1+after
               nin2=nin1+atb
               nin3=nin2+atb
               nin4=nin3+atb
               nin5=nin4+atb
               nout1=nout1+atn
               nout2=nout1+after
               nout3=nout2+after
               nout4=nout3+after
               nout5=nout4+after
               do j=1,nfft
                  r1=zin(1,j,nin1)
                  s1=zin(2,j,nin1)
                  r=zin(1,j,nin2)
                  s=zin(2,j,nin2)
                  r2=r*cr2 - s*ci2
                  s2=r*ci2 + s*cr2
                  r=zin(1,j,nin3)
                  s=zin(2,j,nin3)
                  r3=r*cr3 - s*ci3
                  s3=r*ci3 + s*cr3
                  r=zin(1,j,nin4)
                  s=zin(2,j,nin4)
                  r4=r*cr4 - s*ci4
                  s4=r*ci4 + s*cr4
                  r=zin(1,j,nin5)
                  s=zin(2,j,nin5)
                  r5=r*cr5 - s*ci5
                  s5=r*ci5 + s*cr5
                  r25 = r2 + r5
                  r34 = r3 + r4
                  s25 = s2 - s5
                  s34 = s3 - s4
                  zout(1,nout1,j) = r1 + r25 + r34
                  r = r1 + cos2*r25 + cos4*r34
                  s = sin2*s25 + sin4*s34
                  zout(1,nout2,j) = r - s
                  zout(1,nout5,j) = r + s
                  r = r1 + cos4*r25 + cos2*r34
                  s = sin4*s25 - sin2*s34
                  zout(1,nout3,j) = r - s
                  zout(1,nout4,j) = r + s
                  r25 = r2 - r5
                  r34 = r3 - r4
                  s25 = s2 + s5
                  s34 = s3 + s4
                  zout(2,nout1,j) = s1 + s25 + s34
                  r = s1 + cos2*s25 + cos4*s34
                  s = sin2*r25 + sin4*r34
                  zout(2,nout2,j) = r + s
                  zout(2,nout5,j) = r - s
                  r = s1 + cos4*s25 + cos2*s34
                  s = sin4*r25 - sin2*r34
                  zout(2,nout3,j) = r + s
                  zout(2,nout4,j) = r - s
               end do
            end do
         end if
      end do loop_5000
      !End of radix 5

      ! Radix 6
   else if (now.eq.6) then
      ! .5d0*sqrt(3.d0)
      bb=i_sign*0.8660254037844387d0
      ia=1
      nin1=ia-after
      nout1=ia-atn
      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nin4=nin3+atb
         nin5=nin4+atb
         nin6=nin5+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         nout4=nout3+after
         nout5=nout4+after
         nout6=nout5+after
         do j=1,nfft
            r2=zin(1,j,nin3)
            s2=zin(2,j,nin3)
            r3=zin(1,j,nin5)
            s3=zin(2,j,nin5)
            r=r2 + r3
            s=s2 + s3
            r1=zin(1,j,nin1)
            s1=zin(2,j,nin1)
            ur1 = r + r1
            ui1 = s + s1
            r1=r1 - .5d0*r
            s1=s1 - .5d0*s
            r=r2-r3
            s=s2-s3
            ur2 = r1 - s*bb
            ui2 = s1 + r*bb
            ur3 = r1 + s*bb
            ui3 = s1 - r*bb

            r2=zin(1,j,nin6)
            s2=zin(2,j,nin6)
            r3=zin(1,j,nin2)
            s3=zin(2,j,nin2)
            r=r2 + r3
            s=s2 + s3
            r1=zin(1,j,nin4)
            s1=zin(2,j,nin4)
            vr1 = r + r1
            vi1 = s + s1
            r1=r1 - .5d0*r
            s1=s1 - .5d0*s
            r=r2-r3
            s=s2-s3
            vr2 = r1 - s*bb
            vi2 = s1 + r*bb
            vr3 = r1 + s*bb
            vi3 = s1 - r*bb

            zout(1,nout1,j)=ur1+vr1
            zout(2,nout1,j)=ui1+vi1
            zout(1,nout5,j)=ur2+vr2
            zout(2,nout5,j)=ui2+vi2
            zout(1,nout3,j)=ur3+vr3
            zout(2,nout3,j)=ui3+vi3
            zout(1,nout4,j)=ur1-vr1
            zout(2,nout4,j)=ui1-vi1
            zout(1,nout2,j)=ur2-vr2
            zout(2,nout2,j)=ui2-vi2
            zout(1,nout6,j)=ur3-vr3
            zout(2,nout6,j)=ui3-vi3
         end do
      end do
      !End of radix 6

      ! Radix 7
   else if (now.eq.7) then

      ia=1
      nin1=ia-after
      nout1=ia-atn

      do ib=1,before
         nin1=nin1+after
         nin2=nin1+atb
         nin3=nin2+atb
         nin4=nin3+atb
         nin5=nin4+atb
         nin6=nin5+atb
         nin7=nin6+atb
         nout1=nout1+atn
         nout2=nout1+after
         nout3=nout2+after
         nout4=nout3+after
         nout5=nout4+after
         nout6=nout5+after
         nout7=nout6+after
         do j=1,nfft

            cos2= 0.62348980185873353053d0
            cos4=-0.22252093395631440429d0
            cos6=-0.90096886790241912624d0

            sin2= 0.78183148246802980871d0*real(i_sign,kind=8)
            sin4= 0.97492791218182360702d0*real(i_sign,kind=8)
            sin6= 0.43388373911755812048d0*real(i_sign,kind=8)

            r1=zin(1,j,nin1);   s1=zin(2,j,nin1)
            r2=zin(1,j,nin2);   s2=zin(2,j,nin2)
            r3=zin(1,j,nin3);   s3=zin(2,j,nin3)
            r4=zin(1,j,nin4);   s4=zin(2,j,nin4)
            r5=zin(1,j,nin5);   s5=zin(2,j,nin5)
            r6=zin(1,j,nin6);   s6=zin(2,j,nin6)
            r7=zin(1,j,nin7);   s7=zin(2,j,nin7)

            rp2=r2+r7; rp3=r3+r6; rp4=r4+r5
            rm2=r2-r7; rm3=r3-r6; rm4=r4-r5
            sp2=s2+s7; sp3=s3+s6; sp4=s4+s5
            sm2=s2-s7; sm3=s3-s6; sm4=s4-s5

            ur2=r1+cos2*rp2+cos4*rp3+cos6*rp4
            vr2=   sin2*sm2+sin4*sm3+sin6*sm4

            ui2=s1+cos2*sp2+cos4*sp3+cos6*sp4
            vi2=   sin2*rm2+sin4*rm3+sin6*rm4

            ur3=r1+cos4*rp2+cos6*rp3+cos2*rp4
            vr3=   sin4*sm2-sin6*sm3-sin2*sm4

            ui3=s1+cos4*sp2+cos6*sp3+cos2*sp4
            vi3=   sin4*rm2-sin6*rm3-sin2*rm4

            ur4=r1+cos6*rp2+cos2*rp3+cos4*rp4
            vr4=   sin6*sm2-sin2*sm3+sin4*sm4

            ui4=s1+cos6*sp2+cos2*sp3+cos4*sp4
            vi4=   sin6*rm2-sin2*rm3+sin4*rm4

            zout(1,nout1,j)=r1+rp2+rp3+rp4; zout(2,nout1,j)=s1+sp2+sp3+sp4
            zout(1,nout2,j)=ur2-vr2;        zout(2,nout2,j)=ui2+vi2
            zout(1,nout3,j)=ur3-vr3;        zout(2,nout3,j)=ui3+vi3
            zout(1,nout4,j)=ur4-vr4;        zout(2,nout4,j)=ui4+vi4
            zout(1,nout5,j)=ur4+vr4;        zout(2,nout5,j)=ui4-vi4
            zout(1,nout6,j)=ur3+vr3;        zout(2,nout6,j)=ui3-vi3
            zout(1,nout7,j)=ur2+vr2;        zout(2,nout7,j)=ui2-vi2

         end do
      end do
      !End of radix 7

   else
      stop 'error fftrot_sg'
   end if

END SUBROUTINE fftrot_sg

