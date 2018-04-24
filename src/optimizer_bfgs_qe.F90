!
! Copyright (C) 2002-2006 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------

MODULE save_bfgs
implicit none
REAL(8),allocatable :: spos_p(:)
REAL(8),allocatable :: sgrad_p(:)
INTEGER             :: sscf_iter
INTEGER             :: sbfgs_iter
INTEGER             :: sgdiis_iter
real(8)             :: senergy_p
REAL(8),allocatable :: spos_old(:,:)
REAL(8),allocatable :: sgrad_old(:,:)
REAL(8),allocatable :: sinv_hess(:,:)
INTEGER             :: str_min_hit
REAL(8)             :: snr_step_length
LOGICAL             :: prev_bfgs
END MODULE


MODULE constants
  !----------------------------------------------------------------------------
  !
!  USE kinds, ONLY : DP
  !
  ! ... The constants needed everywhere
  !
  IMPLICIT NONE
  !
  SAVE
  !
  ! ... Mathematical constants
  ! 
  REAL(8), PARAMETER :: pi     = 3.14159265358979323846d0 
  REAL(8), PARAMETER :: tpi    = 2.0d0 * pi
  REAL(8), PARAMETER :: fpi    = 4.0d0 * pi
  REAL(8), PARAMETER :: sqrtpi = 1.77245385090551602729d0 
  REAL(8), PARAMETER :: sqrtpm1= 1.0d0 / sqrtpi
  REAL(8), PARAMETER :: sqrt2  = 1.41421356237309504880d0
  !
  ! ... Physical constants, SI (NIST CODATA 2006), Web Version 5.1
  !     http://physics.nist.gov/constants
  REAL(8), PARAMETER :: H_PLANCK_SI      = 6.62606896d-34   ! J s
  REAL(8), PARAMETER :: K_BOLTZMANN_SI   = 1.3806504d-23    ! J K^-1 
  REAL(8), PARAMETER :: ELECTRON_SI      = 1.602176487d-19  ! C
  REAL(8), PARAMETER :: ELECTRONVOLT_SI  = 1.602176487d-19  ! J  
  REAL(8), PARAMETER :: ELECTRONMASS_SI  = 9.10938215d-31   ! Kg
  REAL(8), PARAMETER :: HARTREE_SI       = 4.35974394d-18   ! J
  REAL(8), PARAMETER :: RYDBERG_SI       = HARTREE_SI/2.0d0   ! J
  REAL(8), PARAMETER :: BOHR_RADIUS_SI   = 0.52917720859d-10 ! m
  REAL(8), PARAMETER :: AMU_SI           = 1.660538782d-27  ! Kg
  REAL(8), PARAMETER :: C_SI             = 2.99792458d+8    ! m sec^-1
  REAL(8), PARAMETER :: MUNOUGHT_SI      = fpi*1.0d-7       ! N A^-2
  REAL(8), PARAMETER :: EPSNOUGHT_SI     = 1.0d0 / (MUNOUGHT_SI * &
                                                       C_SI**2) ! F m^-1
  !
  ! ... Physical constants, atomic units:
  ! ... AU for "Hartree" atomic units (e = m = hbar = 1)
  ! ... RY for "Rydberg" atomic units (e^2=2, m=1/2, hbar=1)
  !
  REAL(8), PARAMETER :: K_BOLTZMANN_AU   = K_BOLTZMANN_SI / HARTREE_SI
  REAL(8), PARAMETER :: K_BOLTZMANN_RY   = K_BOLTZMANN_SI / RYDBERG_SI
  !
  ! ... Unit conversion factors: energy and masses
  !
  REAL(8), PARAMETER :: AUTOEV           = HARTREE_SI / ELECTRONVOLT_SI
  REAL(8), PARAMETER :: RYTOEV           = AUTOEV / 2.0d0
  REAL(8), PARAMETER :: AMU_AU           = AMU_SI / ELECTRONMASS_SI
  REAL(8), PARAMETER :: AMU_RY           = AMU_AU / 2.0d0
  !
  ! ... Unit conversion factors: atomic unit of time, in s and ps
  !
  REAL(8), PARAMETER :: AU_SEC           = H_PLANCK_SI/tpi/HARTREE_SI
  REAL(8), PARAMETER :: AU_PS            = AU_SEC * 1.0d+12
  !
  ! ... Unit conversion factors: pressure (1 Pa = 1 J/m^3, 1GPa = 10 Kbar )
  !
  REAL(8), PARAMETER :: AU_GPA           = HARTREE_SI / BOHR_RADIUS_SI ** 3 &
                                            / 1.0d+9 
  REAL(8), PARAMETER :: RY_KBAR          = 10.0d0 * AU_GPA / 2.0d0
  !
  ! ... Unit conversion factors: 1 debye = 10^-18 esu*cm 
  ! ...                                  = 3.3356409519*10^-30 C*m 
  ! ...                                  = 0.208194346 e*A
  ! ... ( 1 esu = (0.1/c) Am, c=299792458 m/s)
  !
  REAL(8), PARAMETER :: DEBYE_SI         = 3.3356409519d0 * 1.0d-30 ! C*m 
  REAL(8), PARAMETER :: AU_DEBYE         = ELECTRON_SI * BOHR_RADIUS_SI / &
                                            DEBYE_SI
  !
  REAL(8), PARAMETER :: eV_to_kelvin = ELECTRONVOLT_SI / K_BOLTZMANN_SI
  REAL(8), PARAMETER :: ry_to_kelvin = RYDBERG_SI / K_BOLTZMANN_SI
  !
  ! .. Unit conversion factors: Energy to wavelength
  !
  REAL(8), PARAMETER :: EVTONM = 1d+9 * H_PLANCK_SI * C_SI / &
                                  &ELECTRONVOLT_SI
  REAL(8), PARAMETER :: RYTONM = 1d+9 * H_PLANCK_SI * C_SI / RYDBERG_SI
  !
  !  Speed of light in atomic units
  !
  REAL(8), PARAMETER :: C_AU             = C_SI / BOHR_RADIUS_SI * AU_SEC
  !
  ! ... zero up to a given accuracy
  !
  REAL(8), PARAMETER :: eps4  = 1.0d-4
  REAL(8), PARAMETER :: eps6  = 1.0d-6
  REAL(8), PARAMETER :: eps8  = 1.0d-8
  REAL(8), PARAMETER :: eps12 = 1.0d-12
  REAL(8), PARAMETER :: eps14 = 1.0d-14
  REAL(8), PARAMETER :: eps16 = 1.0d-16
  REAL(8), PARAMETER :: eps24 = 1.0d-24
  REAL(8), PARAMETER :: eps32 = 1.0d-32
  !
  REAL(8), PARAMETER :: gsmall = 1.0d-12
  !
  REAL(8), PARAMETER :: e2 = 2.0d0      ! the square of the electron charge
  REAL(8), PARAMETER :: degspin = 2.0d0 ! the number of spins per level
  !
  !!!!!! COMPATIBIILITY
  !
  REAL(8), PARAMETER :: BOHR_RADIUS_CM = BOHR_RADIUS_SI * 100.0d0
  REAL(8), PARAMETER :: BOHR_RADIUS_ANGS = BOHR_RADIUS_CM * 1.0d8
  REAL(8), PARAMETER :: ANGSTROM_AU = 1.0d0/BOHR_RADIUS_ANGS
  REAL(8), PARAMETER :: DIP_DEBYE = AU_DEBYE
  REAL(8), PARAMETER :: AU_TERAHERTZ  = AU_PS
  REAL(8), PARAMETER :: AU_TO_OHMCMM1 = 46000.0d0 ! (ohm cm)^-1
  REAL(8), PARAMETER :: RY_TO_THZ = 1.0d0 / AU_TERAHERTZ / FPI
  REAL(8), PARAMETER :: RY_TO_GHZ = RY_TO_THZ*1000.0d0
  REAL(8), PARAMETER :: RY_TO_CMM1 = 1.d+10 * RY_TO_THZ / C_SI
  !

END MODULE constants
!
! Copyright (C) 2003-2005 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!----------------------------------------------------------------------------
MODULE basic_algebra_routines
  !----------------------------------------------------------------------------
  !
  ! ... Written by Carlo Sbraccia ( 16/12/2003 )
  !
  ! ... This module contains a limited number of functions and operators
  ! ... for vectorial algebra. Wherever possible the appropriate BLAS routine
  ! ... ( always the double precision version ) is used.
  !
  ! ... List of public methods :
  !
  !  x .dot. y          dot product between vectors ( <x|y> )
  !  x .ext. y          external (vector) product between vectors ( <x|y> )
  !  norm( x )          norm of a vector ( SQRT(<x|x>) )
  !  A .times. x        matrix-vector multiplication ( A|x> )
  !  x .times. A        vector-matrix multiplication ( <x|A )
  !  matrix( x, y )     vector-vector multiplication ( |x><y| )
  !  identity( N )      identity matrix of rank N
  !
!  USE kinds,  ONLY : DP
  !
  IMPLICIT NONE
  !
  INTERFACE OPERATOR( .dot. )
     !
     MODULE PROCEDURE dot_product_
     !
  END INTERFACE
  !
  INTERFACE OPERATOR( .ext. )
     !
     MODULE PROCEDURE external_product_
     !
  END INTERFACE
  !  
  INTERFACE OPERATOR( .times. )
     !
     MODULE PROCEDURE matrix_times_vector, vector_times_matrix
     !
  END INTERFACE
  !
  CONTAINS
     !
     !-----------------------------------------------------------------------
     FUNCTION dot_product_( vec1, vec2 ) result(dotprod)
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(8), INTENT(IN) :: vec1(:), vec2(:)
       REAL(8)             :: dotprod
       !
       REAL(8) :: ddot
       EXTERNAL    ddot
       !
       dotprod = ddot( SIZE( vec1 ), vec1, 1, vec2, 1 )
       !
     END FUNCTION dot_product_
     !
     !-----------------------------------------------------------------------
     FUNCTION external_product_( vec1, vec2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(8), INTENT(IN) :: vec1(:), vec2(:)
       REAL(8)             :: external_product_(SIZE( vec1 ))
       !
       !
       external_product_(1) = + vec1(2)*vec2(3) - vec1(3)*vec2(2)
       external_product_(2) = - vec1(1)*vec2(3) - vec1(3)*vec2(1)
       external_product_(3) = + vec1(1)*vec2(2) - vec1(2)*vec2(1)
       !
     END FUNCTION external_product_
     !     
     !----------------------------------------------------------------------- 
     FUNCTION norm( vec ) result(dnrm)
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(8), INTENT(IN) :: vec(:)
       REAL(8)             :: dnrm
       !
       REAL(8) :: dnrm2
       EXTERNAL    dnrm2   
       !
       dnrm = dnrm2( SIZE( vec ), vec, 1 )
       !
     END FUNCTION norm
     !
     !-----------------------------------------------------------------------
     FUNCTION matrix_times_vector( mat, vec )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(8), INTENT(IN) :: vec(:)
       REAL(8), INTENT(IN) :: mat(:,:)
       REAL(8)             :: matrix_times_vector(SIZE( vec ))
! gfortran hack
       REAL(8)             :: aux(SIZE( vec ))
       INTEGER              :: dim1
       !
       dim1 = SIZE( vec )
       !
       CALL DGEMV( 'N', dim1, dim1, 1.0d0, mat, dim1, vec, 1, 0.0d0, &
                   aux, 1 ) 
! gfortran hack 
!                  matrix_times_vector, 1 )
       matrix_times_vector = aux
       !
     END FUNCTION  matrix_times_vector
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION vector_times_matrix( vec, mat )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(8), INTENT(IN) :: vec(:)
       REAL(8), INTENT(IN) :: mat(:,:)
       REAL(8)             :: vector_times_matrix(SIZE( vec ))
! gfortran hack
       REAL(8)             :: aux(SIZE( vec ))
       INTEGER              :: dim1
       !
       dim1 = SIZE( vec )
       !
       CALL DGEMV( 'T', dim1, dim1, 1.0d0, mat, dim1, vec, 1, 0.0d0, &
                   aux, 1) 
! gfortran hack 
!                  vector_times_matrix, 1 )
       vector_times_matrix = aux
       !
     END FUNCTION vector_times_matrix
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION matrix( vec1, vec2 )
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       REAL(8), INTENT(IN) :: vec1(:), vec2(:)
       REAL(8)             :: matrix(SIZE( vec1 ),SIZE( vec2 ))
!#ifdef __GFORTRAN
! gfortran hack - explicit preprocessing is used because this hack
! costs an additional matrix allocation, which may not be a good idea
       REAL(8)             :: aux(SIZE( vec1 ),SIZE( vec2 ))
!#endif
       INTEGER              :: dim1, dim2
       !
       dim1 = SIZE( vec1 )
       dim2 = SIZE( vec2 )
       !
!#ifdef __GFORTRAN
       !
       aux = 0.0d0
       CALL DGER( dim1, dim2, 1.0d0, vec1, 1, vec2, 1, aux, dim1 )
       matrix = aux
!#else
!       !
!       matrix = 0.0d0
!       CALL DGER( dim1, dim2, 1.0d0, vec1, 1, vec2, 1, matrix, dim1 )
!       !
!#endif
       !
     END FUNCTION matrix
     !
     !
     !-----------------------------------------------------------------------
     FUNCTION identity( dim ) result(iden)
       !-----------------------------------------------------------------------
       !
       IMPLICIT NONE
       !
       INTEGER, INTENT(IN) :: dim
       REAL(8)            :: iden(dim,dim)
       INTEGER             :: i
       !
       iden = 0.0d0
       !
       FORALL( i = 1:dim ) iden(i,i) = 1.0d0
       !
     END FUNCTION identity
     !    
END MODULE basic_algebra_routines
! Copyright (C) 2003-2010 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
MODULE bfgs_module
   !----------------------------------------------------------------------------
   !
   ! ... Ionic relaxation through the Newton-Raphson optimization scheme
   ! ... based on the Broyden-Fletcher-Goldfarb-Shanno algorithm for the
   ! ... estimate of the inverse Hessian matrix.
   ! ... The ionic relaxation is performed converting cartesian (and cell) 
   ! ... positions into internal coordinates.
   ! ... The algorithm uses a "trust radius" line search based on Wolfe 
   ! ... conditions. Steps are rejected until the first Wolfe condition
   ! ... (sufficient energy decrease) is satisfied. Updated step length
   ! ... is estimated from quadratic interpolation. 
   ! ... When the step is accepted inverse hessian is updated according to 
   ! ... BFGS scheme and a new search direction is obtained from NR or GDIIS
   ! ... method. The corresponding step length is limited by trust_radius_max 
   ! ... and can't be larger than the previous step multiplied by a certain 
   ! ... factor determined by Wolfe and other convergence conditions.
   !
   ! ... Originally written ( 5/12/2003 ) and maintained ( 2003-2007 ) by 
   ! ... Carlo Sbraccia
   ! ... Modified for variable-cell-shape relaxation ( 2007-2008 ) by 
   ! ...   Javier Antonio Montoya, Lorenzo Paulatto and Stefano de Gironcoli
   ! ... Re-analyzed by Stefano de Gironcoli ( 2010 )
   !
   ! ... references :
   !
   ! ... 1) Roger Fletcher, Practical Methods of Optimization, John Wiley and
   ! ...    Sons, Chichester, 2nd edn, 1987.
   ! ... 2) Salomon R. Billeter, Alexander J. Turner, Walter Thiel,
   ! ...    Phys. Chem. Chem. Phys. 2, 2177 (2000).
   ! ... 3) Salomon R. Billeter, Alessandro Curioni, Wanda Andreoni,
   ! ...    Comput. Mat. Science 27, 437, (2003).
   ! ... 4) Ren Weiqing, PhD Thesis: Numerical Methods for the Study of Energy
   ! ...    Landscapes and Rare Events.
   !
   !
!   USE kinds,     ONLY : DP
!   USE io_files,  ONLY : iunbfgs, prefix
!   USE constants, ONLY : eps16
 !  USE cell_base, ONLY : iforceh  !This is fixlat parameters
   !
   USE basic_algebra_routines
   !
   IMPLICIT NONE
   !
   PRIVATE
   !
   ! ... public methods
   !
   PUBLIC :: bfgs, terminate_bfgs
   !
   ! ... public variables
   !
   PUBLIC :: bfgs_ndim,        &
             trust_radius_ini, trust_radius_min, trust_radius_max, &
             w_1,              w_2
   !
   ! ... global module variables
   !
   SAVE
   !
   CHARACTER (len=8) :: fname="energy" ! name of the function to be minimized
   !
   REAL(8), ALLOCATABLE :: &
      pos(:),            &! positions + cell
      grad(:),           &! gradients + cell_force
      pos_p(:),          &! positions at the previous accepted iteration
      grad_p(:),         &! gradients at the previous accepted iteration
      inv_hess(:,:),     &! inverse hessian matrix (updated using BFGS formula)
      metric(:,:),       &
      h_block(:,:),      &
      hinv_block(:,:),   &
      step(:),           &! the (new) search direction (normalized NR step)
      step_old(:),       &! the previous search direction (normalized NR step)
      pos_old(:,:),      &! list of m old positions - used only by gdiis
      grad_old(:,:),     &! list of m old gradients - used only by gdiis
      pos_best(:)         ! best extrapolated positions - used only by gdiis
   REAL(8) :: &
      nr_step_length,    &! length of (new) Newton-Raphson step
      nr_step_length_old,&! length of previous Newton-Raphson step
      trust_radius,      &! new displacement along the search direction
      trust_radius_old,  &! old displacement along the search direction
      energy_p            ! energy at previous accepted iteration
   INTEGER :: &
      scf_iter,          &! number of scf iterations
      bfgs_iter,         &! number of bfgs iterations
      gdiis_iter,        &! number of gdiis iterations
      tr_min_hit = 0      ! set to 1 if the trust_radius has already been
                          ! set to the minimum value at the previous step
                          ! set to 2 if trust_radius is reset again: exit
   LOGICAL :: &
      conv_bfgs           ! .TRUE. when bfgs convergence has been achieved
   !
   ! ... default values for the following variables are set in
   ! ... Modules/read_namelist.f90 (SUBROUTINE ions_defaults)
   !
   ! ... Note that trust_radius_max, trust_radius_min, trust_radius_ini,
   ! ... w_1, w_2, bfgs_ndim have a default value, but can also be assigned
   ! ... in the input.
   !
   INTEGER :: &
      bfgs_ndim           ! dimension of the subspace for GDIIS
                          ! fixed to 1 for standard BFGS algorithm
   REAL(8)  :: &
      trust_radius_ini,  &! suggested initial displacement
      trust_radius_min,  &! minimum allowed displacement
      trust_radius_max    ! maximum allowed displacement

   REAL(8)  ::          &! parameters for Wolfe conditions
      w_1,               &! 1st Wolfe condition: sufficient energy decrease
      w_2                 ! 2nd Wolfe condition: sufficient gradient decrease
   REAL(8), PUBLIC  :: &
    upscale            ! maximum reduction of convergence threshold
   
!My own
   REAL(8):: iforceh(3,3)=1.d0
   REAL(8), PARAMETER :: eps16 = 1.0d-16
   CHARACTER(10):: prefix="QBFGS"
   INTEGER:: iunbfgs=121
   !
CONTAINS
   !
   !------------------------------------------------------------------------
   SUBROUTINE bfgs( pos_in, h, energy, grad_in, fcell, fixion, scratch, stdout,&
                 energy_thr, grad_thr, cell_thr, energy_error, grad_error,     &
                 cell_error, istep, nstep, step_accepted, stop_bfgs, lmovecell)
      !------------------------------------------------------------------------
      !
      ! ... list of input/output arguments :
      !
      !  pos            : vector containing 3N coordinates of the system ( x )
      !  energy         : energy of the system ( V(x) )
      !  grad           : vector containing 3N components of grad( V(x) )
      !  fixion         : vector used to freeze a deg. of freedom
      !  scratch        : scratch directory
      !  stdout         : unit for standard output
      !  energy_thr     : treshold on energy difference for BFGS convergence
      !  grad_thr       : treshold on grad difference for BFGS convergence
      !                    the largest component of grad( V(x) ) is considered
      !  energy_error   : energy difference | V(x_i) - V(x_i-1) |
      !  grad_error     : the largest component of
      !                         | grad(V(x_i)) - grad(V(x_i-1)) |
      !  cell_error     : the largest component of: omega*(stress-press*I)
      !  nstep          : the maximun nuber of scf-steps
      !  step_accepted  : .TRUE. if a new BFGS step is done
      !  stop_bfgs      : .TRUE. if BFGS convergence has been achieved
      !
      IMPLICIT NONE
      !
      REAL(8),         INTENT(INOUT) :: pos_in(:)
      REAL(8),         INTENT(INOUT) :: h(3,3)
      REAL(8),         INTENT(INOUT) :: energy
      REAL(8),         INTENT(INOUT) :: grad_in(:)
      REAL(8),         INTENT(INOUT) :: fcell(3,3)
      INTEGER,          INTENT(IN)    :: fixion(:)
      CHARACTER(LEN=*), INTENT(IN)    :: scratch
      INTEGER,          INTENT(IN)    :: stdout
      REAL(8),         INTENT(IN)    :: energy_thr, grad_thr, cell_thr
      INTEGER,          INTENT(OUT)   :: istep
      INTEGER,          INTENT(IN)    :: nstep
      REAL(8),         INTENT(OUT)   :: energy_error, grad_error, cell_error
      LOGICAL,          INTENT(OUT)   :: step_accepted, stop_bfgs
      LOGICAL,          INTENT(IN)    :: lmovecell
      !
      INTEGER  :: n, i, j, k, nat
      LOGICAL  :: lwolfe
      REAL(8) :: dE0s, den
      ! ... for scaled coordinates
      REAL(8) :: hinv(3,3),g(3,3),ginv(3,3),garbage, omega
      !
      !
      lwolfe=.false.
      n = SIZE( pos_in ) + 9
      nat = size (pos_in) / 3
!      if (nat*3 /= size (pos_in)) call errore('bfgs',' strange dimension',1)
      if (nat*3 /= size (pos_in)) write(*,*)'bfgs',' strange dimension',1
      !
      ! ... work-space allocation
      !
      ALLOCATE( pos(    n ) )
      ALLOCATE( grad(   n ) )
      !
      ALLOCATE( grad_old( n, bfgs_ndim ) )
      ALLOCATE( pos_old(  n, bfgs_ndim ) )
      !
      ALLOCATE( inv_hess( n, n ) )
      !
      ALLOCATE( pos_p(    n ) )
      ALLOCATE( grad_p(   n ) )
      ALLOCATE( step(     n ) )
      ALLOCATE( step_old( n ) )
      ALLOCATE( pos_best( n ) )
      ! ... scaled coordinates work-space
      ALLOCATE( hinv_block( n-9, n-9 ) )
      ! ... cell related work-space
      ALLOCATE( metric( n , n  ) )
      !
      ! ... the BFGS file read (pos & grad) in scaled coordinates
      !
      call invmat(3, h, hinv, omega)
      ! volume is defined to be positve even for left-handed vector triplet
      omega = abs(omega) 
      !
      hinv_block = 0.d0
      FORALL ( k=0:nat-1, i=1:3, j=1:3 ) hinv_block(i+3*k,j+3*k) = hinv(i,j)
      !
      ! ... generate metric to work with scaled ionic coordinates
      g = MATMUL(TRANSPOSE(h),h)
      call invmat(3,g,ginv,garbage)
      metric = 0.d0
      FORALL ( k=0:nat-1,   i=1:3, j=1:3 ) metric(i+3*k,j+3*k) = g(i,j)
      FORALL ( k=nat:nat+2, i=1:3, j=1:3 ) metric(i+3*k,j+3*k) = 0.04 * omega * ginv(i,j)
      !
      ! ... generate bfgs vectors for the degrees of freedom and their gradients
      pos = 0.0
      pos(1:n-9) = pos_in
      if (lmovecell) FORALL( i=1:3, j=1:3)  pos( n-9 + j+3*(i-1) ) = h(i,j)
      grad = 0.0
      grad(1:n-9) = grad_in
      if (lmovecell) FORALL( i=1:3, j=1:3) grad( n-9 + j+3*(i-1) ) = fcell(i,j)*iforceh(i,j)
      !
      ! if the cell moves the quantity to be minimized is the enthalpy
      IF ( lmovecell ) fname="enthalpy"
      !
CALL read_bfgs_file( pos, grad, fixion, energy, scratch, n, stdout )
!
scf_iter = scf_iter + 1
istep    = scf_iter
!
! ... convergence is checked here
!
energy_error = ABS( energy_p - energy )
grad_error = MAXVAL( ABS( MATMUL( TRANSPOSE(hinv_block), grad(1:n-9)) ) )
conv_bfgs = energy_error < energy_thr
conv_bfgs = conv_bfgs .AND. ( grad_error < grad_thr )
!
IF( lmovecell) THEN
  cell_error = MAXVAL( ABS( MATMUL ( TRANSPOSE ( RESHAPE( grad(n-8:n), (/ 3, 3 /) ) ),&
				     TRANSPOSE(h) ) ) ) / omega
  conv_bfgs = conv_bfgs .AND. ( cell_error < cell_thr ) 
!#undef DEBUG
!#ifdef DEBUG
   write (*,'(3f15.10)') TRANSPOSE ( RESHAPE( grad(n-8:n), (/ 3, 3 /) ) )
   write (*,*)
   write (*,'(3f15.10)') TRANSPOSE(h)
   write (*,*)
   write (*,'(3f15.10)') MATMUL (TRANSPOSE( RESHAPE( grad(n-8:n), (/ 3, 3 /) ) ),&
				     TRANSPOSE(h) ) / omega
   write (*,*)
!   write (*,*) cell_error/cell_thr*0.5d0
!#endif
END IF
!
! ... converged (or useless to go on): quick return
!
conv_bfgs = conv_bfgs .OR. ( tr_min_hit > 1 )
IF ( conv_bfgs ) GOTO 1000
!
! ... some output is written
!
WRITE( UNIT = stdout, &
   & FMT = '(/,5X,"number of scf cycles",T30,"= ",I3)' ) scf_iter
WRITE( UNIT = stdout, &
   & FMT = '(5X,"number of bfgs steps",T30,"= ",I3,/)' ) bfgs_iter
IF ( scf_iter > 1 ) WRITE( UNIT = stdout, &
   & FMT = '(5X,A," old",T30,"= ",F18.10," Ry")' ) fname,energy_p
WRITE( UNIT = stdout, &
   & FMT = '(5X,A," new",T30,"= ",F18.10," Ry",/)' ) fname,energy
!
! ... the bfgs algorithm starts here
!
IF ( .NOT. energy_wolfe_condition( energy ) .AND. (scf_iter > 1) ) THEN
 !
 ! ... the previous step is rejected, line search goes on
 !
 step_accepted = .FALSE.
 !
 WRITE( UNIT = stdout, &
      & FMT = '(5X,"CASE: ",A,"_new > ",A,"_old",/)' ) fname,fname
 !
 ! ... the new trust radius is obtained by quadratic interpolation
 !
 ! ... E(s) = a*s*s + b*s + c      ( we use E(0), dE(0), E(s') )
 !
 ! ... s_min = - 0.5*( dE(0)*s'*s' ) / ( E(s') - E(0) - dE(0)*s' )
 !
 !if (abs(scnorm(step_old(:))-1.d0) > 1.d-10) call errore('bfgs', &
 if (abs(scnorm(step_old(:))-1.d0) > 1.d-10) write(*,*)'bfgs', &
	  ' step_old is NOT normalized ',1
 ! (normalized) search direction is the same as in previous step
 step(:) = step_old(:)
 !
 dE0s = ( grad_p(:) .dot. step(:) ) * trust_radius_old
 !IF (dE0s > 0.d0 ) CALL errore( 'bfgs', &
 IF (dE0s > 0.d0 ) write(*,*)  'bfgs', &
	  'dE0s is positive which should never happen', 1 
 den = energy - energy_p - dE0s
 !
 ! estimate new trust radius by interpolation
 trust_radius = - 0.5d0*dE0s*trust_radius_old / den
 !
 WRITE( UNIT = stdout, &
      & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
      trust_radius
 !
 ! ... values from the last succeseful bfgs step are restored
 !
 pos(:)  = pos_p(:)
 energy  = energy_p
 grad(:) = grad_p(:)
 !
 IF ( trust_radius < trust_radius_min ) THEN
    !
    ! ... the history is reset ( the history can be reset at most two
    ! ... consecutive times )
    !
    WRITE( UNIT = stdout, &
	   FMT = '(/,5X,"trust_radius < trust_radius_min")' )
    WRITE( UNIT = stdout, FMT = '(/,5X,"resetting bfgs history",/)' )
    !
    ! ... if tr_min_hit=1 the history has already been reset at the 
    ! ... previous step : something is going wrong
    !
    IF ( tr_min_hit == 1 ) THEN
!               CALL infomsg( 'bfgs', &
       write(*,*)  'bfgs', &
		    'history already reset at previous step: stopping' 
       tr_min_hit = 2 
    ELSE
       tr_min_hit = 1
    END IF
    !
    CALL reset_bfgs( n )
    !
    step(:) = - ( inv_hess(:,:) .times. grad(:) )
    ! normalize step but remember its length
    nr_step_length = scnorm(step)
    step(:) = step(:) / nr_step_length
    !
    trust_radius = min(trust_radius_ini, nr_step_length)
    !
 ELSE
    !
    tr_min_hit = 0
    !
 END IF
 !
ELSE
 !
 ! ... a new bfgs step is done
 !
 bfgs_iter = bfgs_iter + 1
 !
 IF ( bfgs_iter == 1 ) THEN
    !
    ! ... first iteration
    !
    step_accepted = .FALSE.
    !
 ELSE
    !
    step_accepted = .TRUE.
    !
    nr_step_length_old = nr_step_length
    !
    WRITE( UNIT = stdout, &
	 & FMT = '(5X,"CASE: ",A,"_new < ",A,"_old",/)' ) fname,fname
    !
    CALL check_wolfe_conditions( lwolfe, energy, grad )
    !
    CALL update_inverse_hessian( pos, grad, n, stdout )
    !
 END IF
 ! compute new search direction and store NR step length
 IF ( bfgs_ndim > 1 ) THEN
    !
    ! ... GDIIS extrapolation
    !
    CALL gdiis_step()
    !
 ELSE
    !
    ! ... standard Newton-Raphson step
    !
    step(:) = - ( inv_hess(:,:) .times. grad(:) )
    !
 END IF
 IF ( ( grad(:) .dot. step(:) ) > 0.0d0 ) THEN
    !
    WRITE( UNIT = stdout, &
	   FMT = '(5X,"uphill step: resetting bfgs history",/)' )
    !
    CALL reset_bfgs( n )
    step(:) = - ( inv_hess(:,:) .times. grad(:) )
    !
 END IF
 !
 ! normalize the step and save the step length
 nr_step_length = scnorm(step)
 step(:) = step(:) / nr_step_length
 !
 ! ... the new trust radius is computed
 !
 IF ( bfgs_iter == 1 ) THEN
    !
    trust_radius = min(trust_radius_ini, nr_step_length)
    tr_min_hit = 0
    !
 ELSE
    !
    CALL compute_trust_radius( lwolfe, energy, grad, n, stdout )
    !
 END IF
 !
 WRITE( UNIT = stdout, &
      & FMT = '(5X,"new trust radius",T30,"= ",F18.10," bohr")' ) &
      trust_radius
 !
END IF
!
! ... step along the bfgs direction
!
IF ( nr_step_length < eps16 ) &
 !CALL errore( 'bfgs', 'NR step-length unreasonably short', 1 )
 write(*,*) 'bfgs', 'NR step-length unreasonably short', 1 
!
! ... information required by next iteration is saved here ( this must
! ... be done before positions are updated )
!
CALL write_bfgs_file( pos, energy, grad, scratch, n )
!
! ... positions and cell are updated
!
pos(:) = pos(:) + trust_radius * step(:)
!
1000  stop_bfgs = conv_bfgs .OR. ( scf_iter >= nstep ) 
! ... input ions+cell variables
IF ( lmovecell ) FORALL( i=1:3, j=1:3) h(i,j) = pos( n-9 + j+3*(i-1) )
pos_in = pos(1:n-9)
! ... update forces
grad_in = grad(1:n-9)
!
! ... work-space deallocation
!
DEALLOCATE( pos )
DEALLOCATE( grad )
DEALLOCATE( pos_p )
DEALLOCATE( grad_p )
DEALLOCATE( pos_old )
DEALLOCATE( grad_old )
DEALLOCATE( inv_hess )
DEALLOCATE( step )
DEALLOCATE( step_old )
DEALLOCATE( pos_best )
DEALLOCATE( hinv_block )
DEALLOCATE( metric )
!
RETURN
!
CONTAINS
!
!--------------------------------------------------------------------
SUBROUTINE gdiis_step()
 !--------------------------------------------------------------------
 USE basic_algebra_routines
 IMPLICIT NONE
 !
 REAL(8), ALLOCATABLE :: res(:,:), overlap(:,:), work(:)
 INTEGER,  ALLOCATABLE :: iwork(:)
 INTEGER               :: k, k_m, info
 REAL(8)              :: gamma0
 !
 !
 gdiis_iter = gdiis_iter + 1
 !
 k   = MIN( gdiis_iter, bfgs_ndim )
 k_m = k + 1
 !
 ALLOCATE( res( n, k ) )
 ALLOCATE( overlap( k_m, k_m ) )
 ALLOCATE( work( k_m ), iwork( k_m ) )
 !
 work(:)  = 0.0d0
 iwork(:) = 0
 !
 ! ... the new direction is added to the workspace
 !
 DO i = bfgs_ndim, 2, -1
    !
    pos_old(:,i)  = pos_old(:,i-1)
    grad_old(:,i) = grad_old(:,i-1)
    !
 END DO
 !
 pos_old(:,1)  = pos(:)
 grad_old(:,1) = grad(:)
 !
 ! ... |res_i> = H^-1 \times |g_i>
 !
 CALL DGEMM( 'N', 'N', n, k, n, 1.0d0, &
	     inv_hess, n, grad_old, n, 0.0d0, res, n )
 !
 ! ... overlap_ij = <grad_i|res_j>
 !
 CALL DGEMM( 'T', 'N', k, k, n, 1.0d0, &
	     res, n, res, n, 0.0d0, overlap, k_m )
 !
 overlap( :, k_m) = 1.0d0
 overlap(k_m, : ) = 1.0d0
 overlap(k_m,k_m) = 0.0d0
 !
 ! ... overlap is inverted via Bunch-Kaufman diagonal pivoting method
 !
 CALL DSYTRF( 'U', k_m, overlap, k_m, iwork, work, k_m, info )
 CALL DSYTRI( 'U', k_m, overlap, k_m, iwork, work, info )
 !CALL errore( 'gdiis_step', 'error in Bunch-Kaufman inversion', info )
 write(*,*) 'gdiis_step', 'error in Bunch-Kaufman inversion', info 
 !
 ! ... overlap is symmetrised
 !
 FORALL( i = 1:k_m, j = 1:k_m, j > i ) overlap(j,i) = overlap(i,j)
 !
 pos_best(:) = 0.0d0
 step(:)     = 0.0d0
 !
 DO i = 1, k
    !
    gamma0 = overlap(k_m,i)
    !
    pos_best(:) = pos_best(:) + gamma0*pos_old(:,i)
    !
    step(:) = step(:) - gamma0*res(:,i)
    !
 END DO
 !
 ! ... the step must be consistent with the last positions
 !
 step(:) = step(:) + ( pos_best(:) - pos(:) )
 !
 IF ( ( grad(:) .dot. step(:) ) > 0.0d0 ) THEN
    !
    ! ... if the extrapolated direction is uphill use only the
    ! ... last gradient and reset gdiis history
    !
    step(:) = - ( inv_hess(:,:) .times. grad(:) )
    !
    gdiis_iter = 0
    !
 END IF
 !
 DEALLOCATE( res, overlap, work, iwork )
 !
END SUBROUTINE gdiis_step
!
END SUBROUTINE bfgs
!
!------------------------------------------------------------------------
SUBROUTINE reset_bfgs( n )
!------------------------------------------------------------------------
! ... inv_hess in re-initalized to the initial guess 
! ... defined as the inverse metric 
!
use save_bfgs, only: prev_bfgs
INTEGER, INTENT(IN) :: n
!
REAL(8) :: garbage
!
call invmat(n, metric, inv_hess, garbage)
!
gdiis_iter = 0
!MYYYYYYYYYYYYYy
    prev_bfgs=.false.
!MYYYYYYYYYYYYYy
!
END SUBROUTINE reset_bfgs
!
!------------------------------------------------------------------------
SUBROUTINE read_bfgs_file( pos, grad, fixion, energy, scratch, n, stdout )
!------------------------------------------------------------------------
!
use save_bfgs
IMPLICIT NONE
!
REAL(8),         INTENT(INOUT) :: pos(:)
REAL(8),         INTENT(INOUT) :: grad(:)
INTEGER,          INTENT(IN)    :: fixion(:)
CHARACTER(LEN=*), INTENT(IN)    :: scratch
INTEGER,          INTENT(IN)    :: n
INTEGER,          INTENT(IN)    :: stdout
REAL(8),         INTENT(INOUT) :: energy
!
CHARACTER(LEN=256) :: bfgs_file
REAL(8) :: garbage
!
!
!!!bfgs_file = TRIM( scratch ) // TRIM( prefix ) // '.bfgs'
!!!!
!!!INQUIRE( FILE = TRIM( bfgs_file ) , EXIST = file_exists )
!!!!
IF ( prev_bfgs ) THEN
 !
 ! ... bfgs is restarted from file
 !
!!!! OPEN( UNIT = iunbfgs, FILE = TRIM( bfgs_file ), &
!!!!       STATUS = 'UNKNOWN', ACTION = 'READ' )
!!!! !
!!!! READ( iunbfgs, * ) pos_p
!!!! READ( iunbfgs, * ) grad_p
!!!! READ( iunbfgs, * ) scf_iter
!!!! READ( iunbfgs, * ) bfgs_iter
!!!! READ( iunbfgs, * ) gdiis_iter
!!!! READ( iunbfgs, * ) energy_p
!!!! READ( iunbfgs, * ) pos_old
!!!! READ( iunbfgs, * ) grad_old
!!!! READ( iunbfgs, * ) inv_hess
!!!! READ( iunbfgs, * ) tr_min_hit
!!!! READ( iunbfgs, * ) nr_step_length
!!!! !
!!!! CLOSE( UNIT = iunbfgs )
 pos_p            =spos_p
 grad_p           =sgrad_p
 scf_iter         =sscf_iter
 bfgs_iter        =sbfgs_iter
 gdiis_iter       =sgdiis_iter
 energy_p         =senergy_p
 pos_old          =spos_old
 grad_old         =sgrad_old
 inv_hess         =sinv_hess
 tr_min_hit       =str_min_hit
 nr_step_length   =snr_step_length
 !
 step_old = ( pos(:) - pos_p(:) ) 
 trust_radius_old = scnorm( step_old )
 step_old = step_old / trust_radius_old
 !
ELSE
 !
 ! ... bfgs initialization
 !
 WRITE( UNIT = stdout, FMT = '(/,5X,"BFGS Geometry Optimization")' )
 !
 ! initialize the inv_hess to the inverse of the metric
 call invmat(n, metric, inv_hess, garbage)
 !
 pos_p      = 0.0d0
 grad_p     = 0.0d0
 scf_iter   = 0
 bfgs_iter  = 0
 gdiis_iter = 0
 energy_p   = energy
 step_old   = 0.0d0
 nr_step_length = 0.0d0
 !
 trust_radius_old = trust_radius_ini
 !
 pos_old(:,:)  = 0.0d0
 grad_old(:,:) = 0.0d0
 !
 tr_min_hit = 0
 !
END IF
prev_bfgs=.true.
!
END SUBROUTINE read_bfgs_file
!
!------------------------------------------------------------------------
SUBROUTINE write_bfgs_file( pos, energy, grad, scratch, n)
      !------------------------------------------------------------------------
      !
use save_bfgs
      IMPLICIT NONE
      !
      INTEGER,         INTENT(IN) :: n
      REAL(8),         INTENT(IN) :: pos(:)
      REAL(8),         INTENT(IN) :: energy
      REAL(8),         INTENT(IN) :: grad(:)
      CHARACTER(LEN=*), INTENT(IN) :: scratch
      !
      !
!!!!      OPEN( UNIT = iunbfgs, FILE = TRIM( scratch )//TRIM( prefix )//'.bfgs', &
!!!!            STATUS = 'UNKNOWN', ACTION = 'WRITE' )
!!!!      !
!!!!      WRITE( iunbfgs, * ) pos
!!!!      WRITE( iunbfgs, * ) grad
!!!!      WRITE( iunbfgs, * ) scf_iter
!!!!      WRITE( iunbfgs, * ) bfgs_iter
!!!!      WRITE( iunbfgs, * ) gdiis_iter
!!!!      WRITE( iunbfgs, * ) energy
!!!!      WRITE( iunbfgs, * ) pos_old
!!!!      WRITE( iunbfgs, * ) grad_old
!!!!      WRITE( iunbfgs, * ) inv_hess
!!!!      WRITE( iunbfgs, * ) tr_min_hit
!!!!      WRITE( iunbfgs, * ) nr_step_length
!!!!      !
!!!!      CLOSE( UNIT = iunbfgs )
 if(.not.allocated(sinv_hess  ))   ALLOCATE( sinv_hess( n, n ) )
 if(.not.allocated(spos_p     ))   ALLOCATE( spos_p(    n ) )
 if(.not.allocated(sgrad_p    ))   ALLOCATE( sgrad_p(   n ) )
 if(.not.allocated(sgrad_old  ))   ALLOCATE( sgrad_old( n, bfgs_ndim ) )
 if(.not.allocated(spos_old   ))   ALLOCATE( spos_old(  n, bfgs_ndim ) )
 spos_p            =pos
 sgrad_p           =grad
 sscf_iter         =scf_iter
 sbfgs_iter        =bfgs_iter
 sgdiis_iter       =gdiis_iter
 senergy_p         =energy
 spos_old          =pos_old
 sgrad_old         =grad_old
 sinv_hess         =inv_hess
 str_min_hit       =tr_min_hit
 snr_step_length   =nr_step_length
      !
   END SUBROUTINE write_bfgs_file
   !
   !------------------------------------------------------------------------
   SUBROUTINE update_inverse_hessian( pos, grad, n, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      REAL(8), INTENT(IN)  :: pos(:)
      REAL(8), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      INTEGER,  INTENT(IN)  :: stdout
      INTEGER               :: info
      !
      REAL(8), ALLOCATABLE :: y(:), s(:)
      REAL(8), ALLOCATABLE :: Hy(:), yH(:)
      REAL(8)              :: sdoty, sBs, Theta
      REAL(8), ALLOCATABLE :: B(:,:)
      !
      ALLOCATE( y( n ), s( n ), Hy( n ), yH( n ) )
      !
      s(:) = pos(:)  - pos_p(:)
      y(:) = grad(:) - grad_p(:)
      !
      sdoty = ( s(:) .dot. y(:) )
      !
      IF ( ABS( sdoty ) < eps16 ) THEN
         !
         ! ... the history is reset
         !
         WRITE( stdout, '(/,5X,"WARNING: unexpected ", &
                         &     "behaviour in update_inverse_hessian")' )
         WRITE( stdout, '(  5X,"         resetting bfgs history",/)' )
         !
         CALL reset_bfgs( n )
         !
         RETURN
         !
      ELSE
!       Conventional Curvature Trap here
!       See section 18.2 (p538-539 ) of Nocedal and Wright "Numerical
!       Optimization"for instance
!       LDM Addition, April 2011
!
!       While with the Wolfe conditions the Hessian in most cases
!       remains positive definite, if one is far from the minimum
!       and/or "bonds" are being made/broken the curvature condition
!              Hy = s ; or s = By
!       cannot be satisfied if s.y < 0. In addition, if s.y is small
!       compared to s.B.s too greedy a step is taken.
!
!       The trap below is conventional and "OK", and has been around
!       for ~ 30 years but, unfortunately, is rarely mentioned in
!       introductory texts and hence often neglected.
!
!       First, solve for inv_hess*t = s ; i.e. t = B*s
!       Use yH as workspace here

        ALLOCATE (B(n,n) )
        B = inv_hess
        yH= s
        call DPOSV('U',n,1,B,n, yH, n, info)
!       Info .ne. 0 should be trapped ...
        if(info .ne. 0)write( stdout, '(/,5X,"WARNING: info=",i3," for Hessian")' )info
        DEALLOCATE ( B )
!
!       Calculate s.B.s
        sBs = ( s(:) .dot. yH(:) )
!
!       Now the trap itself
        if ( sdoty < 0.20D0*sBs ) then
!               Conventional damping
                Theta = 0.8D0*sBs/(sBs-sdoty)
                WRITE( stdout, '(/,5X,"WARNING: bfgs curvature condition ", &
                &     "failed, Theta=",F6.3)' )theta
                y = Theta*y + (1.D0 - Theta)*yH
        endif
      END IF
      !
      Hy(:) = ( inv_hess .times. y(:) )
      yH(:) = ( y(:) .times. inv_hess )
      !
      ! ... BFGS update
      !
      inv_hess = inv_hess + 1.0d0 / sdoty * &
                 ( ( 1.0d0 + ( y .dot. Hy ) / sdoty ) * matrix( s, s ) - &
                  ( matrix( s, yH ) +  matrix( Hy, s ) ) )
      !
      DEALLOCATE( y, s, Hy, yH )
      !
      RETURN
      !
   END SUBROUTINE update_inverse_hessian
   !
   !------------------------------------------------------------------------
   SUBROUTINE check_wolfe_conditions( lwolfe, energy, grad )
      !------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: energy
      REAL(8), INTENT(IN)  :: grad(:)
      LOGICAL,  INTENT(OUT) :: lwolfe
      !
      lwolfe =  energy_wolfe_condition ( energy ) .AND. &
                gradient_wolfe_condition ( grad )
      !
   END SUBROUTINE check_wolfe_conditions
   !
   !------------------------------------------------------------------------
   FUNCTION energy_wolfe_condition ( energy ) result(res)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: energy
      LOGICAL:: res
      !
      res = &
          ( energy-energy_p ) < w_1 * ( grad_p.dot.step_old ) * trust_radius_old
!          ( energy-energy_p ) < abs(w_1 * ( grad_p.dot.step_old ) * trust_radius_old)
      write(*,*) "EWOLFE", energy-energy_p,w_1 * ( grad_p.dot.step_old ) * trust_radius_old
      !
   END FUNCTION energy_wolfe_condition
   !
   !------------------------------------------------------------------------
   FUNCTION gradient_wolfe_condition ( grad ) result(res)
      !------------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8), INTENT(IN)  :: grad(:)
      LOGICAL:: res
      !
      res = &
          ABS( grad .dot. step_old ) < - w_2 * ( grad_p .dot. step_old )
      !
   END FUNCTION gradient_wolfe_condition
   !
   !------------------------------------------------------------------------
   SUBROUTINE compute_trust_radius( lwolfe, energy, grad, n, stdout )
      !------------------------------------------------------------------------
      !
      IMPLICIT NONE
      !
      LOGICAL,  INTENT(IN)  :: lwolfe
      REAL(8), INTENT(IN)  :: energy
      REAL(8), INTENT(IN)  :: grad(:)
      INTEGER,  INTENT(IN)  :: n
      INTEGER,  INTENT(IN)  :: stdout
      !
      REAL(8) :: a
      LOGICAL  :: ltest
      !
      ltest = ( energy - energy_p ) < w_1 * ( grad_p .dot. step_old ) * trust_radius_old
      ltest = ltest .AND. ( nr_step_length_old > trust_radius_old )
      !
      IF ( ltest ) THEN
         a = 1.5d0
      ELSE
         a = 1.1d0
      END IF
      IF ( lwolfe ) a = 2.d0 * a
      !
      trust_radius = MIN( trust_radius_max, a*trust_radius_old, nr_step_length )
      !
      IF ( trust_radius < trust_radius_min ) THEN
         !
         ! ... the history is reset
         !
         ! ... if tr_min_hit the history has already been reset at the 
         ! ... previous step : something is going wrong
         !
         IF ( tr_min_hit == 1 ) THEN
!            CALL infomsg( 'bfgs', &
             write(*,*) 'bfgs', &
                          'history already reset at previous step: stopping' 
            tr_min_hit = 2 
         ELSE
            tr_min_hit = 1
         END IF
         !
         WRITE( UNIT = stdout, &
                FMT = '(5X,"small trust_radius: resetting bfgs history",/)' )
         !
         CALL reset_bfgs( n )
         step(:) = - ( inv_hess(:,:) .times. grad(:) )
         !
         nr_step_length = scnorm(step)
         step(:) = step(:) / nr_step_length
         !
         trust_radius = min(trust_radius_min, nr_step_length )
         !
      ELSE
         !
         tr_min_hit = 0
         !
      END IF
      !
   END SUBROUTINE compute_trust_radius
   !
   !----------------------------------------------------------------------- 
   FUNCTION scnorm1( vect ) result(res)
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: vect(:)
      REAL(8):: res
      !
      res = SQRT( DOT_PRODUCT( vect  ,  MATMUL( metric, vect ) ) )
      !
   END FUNCTION scnorm1
   !
   !----------------------------------------------------------------------- 
   FUNCTION scnorm( vect ) result(res)
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL(8), INTENT(IN) :: vect(:)
      REAL(8) :: ss
      INTEGER :: i,k,l,n
      REAL(8):: res
      !
      res = 0.d0
      n = SIZE (vect) / 3
      do i=1,n
         ss = 0.d0
         do k=1,3
            do l=1,3
               ss = ss + &
                    vect(k+(i-1)*3)*metric(k+(i-1)*3,l+(i-1)*3)*vect(l+(i-1)*3)
            end do
         end do
         res = MAX (res, SQRT (ss) )
      end do
      !
   END FUNCTION scnorm
   !
   !------------------------------------------------------------------------
   SUBROUTINE terminate_bfgs( energy, energy_thr, grad_thr, cell_thr, &
                              lmovecell, stdout, scratch )
      !------------------------------------------------------------------------
      !
!      USE io_files, ONLY : prefix, delete_if_present
use save_bfgs, only: prev_bfgs
      !
      IMPLICIT NONE
      REAL(8),         INTENT(IN) :: energy, energy_thr, grad_thr, cell_thr
      LOGICAL,          INTENT(IN) :: lmovecell
      INTEGER,          INTENT(IN) :: stdout
      CHARACTER(LEN=*), INTENT(IN) :: scratch
      !
      IF ( conv_bfgs ) THEN
         !
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"bfgs converged in ",I3," scf cycles and ", &
              &         I3," bfgs steps")' ) scf_iter, bfgs_iter
         IF ( lmovecell ) THEN
            WRITE( UNIT = stdout, &
              & FMT = '(5X,"(criteria: energy < ",ES8.1,", force < ",ES8.1, &
              &       ", cell < ",ES8.1,")")') energy_thr, grad_thr, cell_thr
         ELSE
            WRITE( UNIT = stdout, &
              & FMT = '(5X,"(criteria: energy < ",ES8.1,", force < ",ES8.1, &
              &                        ")")') energy_thr, grad_thr
         END IF
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         WRITE( UNIT = stdout, &
              & FMT = '(/,5X,"Final ",A," = ",F18.10," Ry")' ) fname, energy
         !
!         CALL delete_if_present( TRIM( scratch ) // TRIM( prefix ) // '.bfgs' )
         prev_bfgs=.false.
         !
      ELSE
         !
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"The maximum number of steps has been reached.")' )
         WRITE( UNIT = stdout, &
                FMT = '(/,5X,"End of BFGS Geometry Optimization")' )
         !
      END IF
      !
   END SUBROUTINE terminate_bfgs
   !
END MODULE bfgs_module
!
! Copyright (C) 2002-2011 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!----------------------------------------------------------------------------
!SUBROUTINE move_ions()
subroutine GEOPT_qbfgs(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,folder)
 use global, only: units,max_kpt,ka1,kb1,kc1,confine
 use defs_basis
 use interface_code
 use modsocket, only: sock_extra_string
 use save_bfgs, only: prev_bfgs,sbfgs_iter
  !----------------------------------------------------------------------------
  !
  ! ... This routine moves the ions according to the requested scheme:
  !
  ! ... lbfgs               bfgs minimizations
  ! ... lmd                 molecular dynamics ( verlet of vcsmd )
  ! ... lmd + lconstrain    constrained molecular dynamics,
  !
  ! ... coefficients for potential and wavefunctions extrapolation are
  ! ... also computed here
  !
  USE constants,              ONLY : e2, eps8, ry_kbar
!  USE io_global,              ONLY : stdout
!  USE io_files,               ONLY : tmp_dir, iunupdate, seqopn
!  USE kinds,                  ONLY : DP
!  USE cell_base,              ONLY : alat, at, bg, omega, cell_force, fix_volume, fix_area
!  USE cellmd,                 ONLY : omega_old, at_old, press, lmovecell, calc
!  USE ions_base,              ONLY : nat, ityp, tau, if_pos
!  USE fft_base,               ONLY : dfftp
!  USE fft_base,               ONLY : dffts
!  USE grid_subroutines,       ONLY : realspace_grids_init
!  USE gvect,                  ONLY : gcutm
!  USE gvecs,                  ONLY : gcutms
!  USE grid_subroutines,       ONLY : realspace_grids_init
!  USE symm_base,              ONLY : checkallsym
!  USE ener,                   ONLY : etot
!  USE force_mod,              ONLY : force, sigma
!  USE control_flags,          ONLY : istep, nstep, upscale, lbfgs, ldamped, &
!                                     lconstrain, conv_ions, use_SMC, &
!                                     lmd, llang, history, tr2
!  USE relax,                  ONLY : epse, epsf, epsp, starting_scf_threshold
!  USE lsda_mod,               ONLY : lsda, absmag
!  USE mp_images,              ONLY : intra_image_comm
!  USE io_global,              ONLY : ionode_id, ionode
!  USE mp,                     ONLY : mp_bcast
  USE bfgs_module,            ONLY : bfgs,&
                                     terminate_bfgs,&
                                     upscale,&
                                     bfgs_ndim,&
                                     trust_radius_max,&
                                     trust_radius_min,&
                                     trust_radius_ini,&
                                     w_1,& 
                                     w_2
!  USE basic_algebra_routines, ONLY : norm
!  USE dynamics_module,        ONLY : verlet, langevin_md, proj_verlet
!  USE dynamics_module,        ONLY : smart_MC
!  USE dfunct,                 only : newd
  !
  use mod_parini, only: typ_parini
  IMPLICIT NONE
  type(typ_parini), intent(in):: parini
  type(typ_parini), intent(inout):: parres
  !
  LOGICAL, SAVE         :: lcheck_mag = .FALSE., &
                           restart_with_starting_magnetiz = .FALSE., &
                           lcheck_cell= .TRUE., &
                           final_cell_calculation=.FALSE.
!  REAL(8), ALLOCATABLE :: tauold(:,:,:)
  REAL(8)              :: energy_error, gradient_error, cell_error
  LOGICAL               :: step_accepted, exst
  REAL(8), ALLOCATABLE :: pos(:), grad(:)
  REAL(8)              :: h(3,3), fcell(3,3)=0.d0, epsp1
  INTEGER,  ALLOCATABLE :: fixion(:)
  real(8) :: tr
  !

!My inputs
  real(8):: latvec_in(3,3),xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,counter,xred(3,parini%nat),fcart(3,parini%nat),latvec(3,3)
  integer:: iprec
  character(40):: folder
  INTEGER :: stdout = 6    ! unit connected to standard output
  logical :: conv_ions
  logical :: lmovecell=.true.,lconstrain=.false.  
  real(8) :: alat, at(3,3), bg(3,3), at_old(3,3), omega, omega_old
  real(8) :: etot, press, sigma(3,3)
  real(8) ::     epse = 1.0d-10,  &! threshold on total energy
       epsf=1.d-10,               &! threshold on forces
       epsp=1.d-10,               &! threshold on pressure
       starting_scf_threshold=1.d-10 ! elf-explanatory
  character(40):: tmp_dir="./",filename
  integer:: istep, nstep, itime,iexit
  real(8):: tr2
  character(4):: fn4
  logical:: getwfk
!For multiprec
  logical:: multiprec,cellfix_done
  real(8):: tolmxf_switch,cellfix_switch
!For writing files
  real(8):: pressure, fmax, fmax_at,fmax_lat,enthalpy,vol,en0000
!If the pure QE implementation with their units and forces should be used
  logical:: qe_units
!Latvec correction io
  integer:: latvec_io
latvec_io=0
!multiprec is hardcoded and, if true, starts a geopt with iprec==2, and then switches 
!to iprec==1 when the fmax==tolmxf_switch. The switch only occurs once
 max_kpt=.false.
 multiprec=.true.
 tolmxf_switch=1.d-3
 cellfix_switch=1.d-3
 cellfix_done=.false.

qe_units=.true.


!  IF (use_SMC) CALL smart_MC()  ! for smart monte carlo method
  !
  ! ... only one node does the calculation in the parallel case
  !
!  IF ( ionode ) THEN
     !
     conv_ions = .FALSE.
     prev_bfgs = .false.
     !

!Define QE stuff
     alat=1.d0
     at = latvec_in / alat                            ! and so the cell
     CALL recips( at(1,1),at(1,2),at(1,3), bg(1,1),bg(1,2),bg(1,3) )
     xred=xred_in
     fcart=fcart_in
     latvec=latvec_in
     call getvol(latvec,omega)
     omega_old=omega
     at_old=at
     press=parini%target_pressure_habohr*2.d0!In ry/bohr^3
     pressure=parini%target_pressure_habohr
!Default values of definable parameters
     bfgs_ndim=parini%qbfgs_bfgs_ndim
     trust_radius_max=parini%qbfgs_trust_radius_max
     trust_radius_min=parini%qbfgs_trust_radius_min
     trust_radius_ini=parini%qbfgs_trust_radius_ini
     w_1=parini%qbfgs_w_1
     w_2=parini%qbfgs_w_2
     nstep=parini%paropt_geopt%nit
     sbfgs_iter=0
     upscale=100.D0

!     ALLOCATE( tauold( 3, nat, 3 ) )
     !
     ! ... the file containing old positions is opened 
     ! ... ( needed for extrapolation )
     !
  !   CALL seqopn( iunupdate, 'update', 'FORMATTED', exst ) 
  !   !
  !   IF ( exst ) THEN
  !      !
  !      READ( UNIT = iunupdate, FMT = * ) history
  !      READ( UNIT = iunupdate, FMT = * ) tauold
  !      !
  !   ELSE
  !      !
  !      history = 0
  !      tauold  = 0.D0
  !      !
  !      WRITE( UNIT = iunupdate, FMT = * ) history
  !      WRITE( UNIT = iunupdate, FMT = * ) tauold
  !      !
  !   END IF
  !   !
  !   CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
     !
     ! ... save the previous two steps ( a total of three steps is saved )
     !
!     tauold(:,:,3) = tauold(:,:,2)
!     tauold(:,:,2) = tauold(:,:,1)
!     tauold(:,:,1) = tau(:,:)
     !
     ! ... history is updated (a new ionic step has been done)
     !
!     history = MIN( 3, ( history + 1 ) )
     !
     ! ... old positions are written on file
     !
!     CALL seqopn( iunupdate, 'update', 'FORMATTED', exst ) 
!     !
!     WRITE( UNIT = iunupdate, FMT = * ) history
!     WRITE( UNIT = iunupdate, FMT = * ) tauold
!     !
!     CLOSE( UNIT = iunupdate, STATUS = 'KEEP' )
!     !  
!     DEALLOCATE( tauold )
     !
     ! ... do the minimization / dynamics step
     !
!!!     IF ( lmovecell .AND. lconstrain ) THEN
!!!        !
!!!        IF ( lbfgs) CALL errore('move_ions', &
!!!            & 'variable-cell bfgs and constraints not implemented yet', 1 )
!!!        WRITE(*, '(5x,"-------------------------------------------")')
!!!        WRITE(*, '(5x,"NEW FEATURE: constraints with variable cell")')
!!!        WRITE(*, '(5x,"-------------------------------------------")')
!!!        !
!!!     END IF
     !
     ! ... BFGS algorithm is used to minimize ionic configuration
     !
!     bfgs_minimization : &
!     IF ( lbfgs ) THEN
        !
        ! ... the bfgs procedure is used
        !  
        ALLOCATE( pos( 3*parini%nat ), grad( 3*parini%nat ), fixion( 3*parini%nat ) )
        !

!****************************************************************************************************************        
       itime=0
       getwfk=.false.
       write(fn4,'(i4.4)') itime
       sock_extra_string="BFGS"//trim(fn4)
       latvec_in = at * alat
       xred_in=xred
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
       fcart=fcart_in
       etot=etot_in
       sigma(1,1)=-strten_in(1)
       sigma(2,2)=-strten_in(2)
       sigma(3,3)=-strten_in(3)
       sigma(1,2)=-strten_in(6);sigma(2,1)=sigma(1,2)
       sigma(1,3)=-strten_in(5);sigma(3,1)=sigma(1,3)
       sigma(2,3)=-strten_in(4);sigma(3,2)=sigma(2,3)
!****************************************************************************************************************        
!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
if(parini%verb.gt.0) then
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in QBFGS: ",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       endif
endif
call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
       write(*,'(a,i4,2x,i4,4(1x,es17.8),1x,i4)') " # GEOPT QBFGS   ",&
              &itime,sbfgs_iter,enthalpy, fmax, fmax_at,fmax_lat,iprec
!Initial iprec after running the first force call
   if(multiprec) iprec=2
       if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
          iprec=1
       endif
   if(iexit==1) then
     write(*,'(a,i4,2(1x,es25.15))') " # QBFGS converged before entering iterations", itime,enthalpy,fmax
     max_kpt=.false.
     return 
   endif
!write(*,*) "NTIME_GEOPT",parini%paropt_geopt%nit
do itime=1,parini%paropt_geopt%nit
!****************************************************************************************************************        
       write(fn4,'(i4.4)') itime
       sock_extra_string="BFGS"//trim(fn4)
       latvec_in = at * alat
       xred_in=xred
       call get_energyandforces_single(parini,parres,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
       fcart=fcart_in
       etot=etot_in
       sigma(1,1)=-strten_in(1)
       sigma(2,2)=-strten_in(2)
       sigma(3,3)=-strten_in(3)
       sigma(1,2)=-strten_in(6);sigma(2,1)=sigma(1,2)
       sigma(1,3)=-strten_in(5);sigma(3,1)=sigma(1,3)
       sigma(2,3)=-strten_in(4);sigma(3,2)=sigma(2,3)
!Single point calculation finished. Now returning to MD part. All output variables are now in *_in
!****************************************************************************************************************        
!MHM: Write output to file in every step***********************************
       counter=real(itime,8)
       call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
       call getvol(latvec_in,vol)
if(parini%verb.gt.0) then
       write(fn4,'(i4.4)') itime
       filename=trim(folder)//"posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in QBFGS: ",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       if(parini%verb.ge.3) then
       filename=trim(folder)//"posgeopt."//fn4//".vasp"
       call write_atomic_file_poscar(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &parini%char_type(1:parini%ntypat_global),parini%ntypat_global,parini%typat_global,parini%fixat,parini%fixlat,etot_in,pressure,enthalpy,en0000)
       endif
endif
call convcheck(parini,parini%nat,latvec_in,fcart_in,strten_in,parini%target_pressure_habohr,parini%paropt_geopt%strfact,fmax,fmax_at,fmax_lat,parini%paropt_geopt%fmaxtol,iexit)
       write(*,'(a,i4,2x,i4,4(1x,es17.8),1x,i4)') " # GEOPT QBFGS   ",&
              &itime,sbfgs_iter,enthalpy, fmax, fmax_at,fmax_lat,iprec
!*************************************************************************************************************        

        h = at * alat
        pos    =   RESHAPE( xred, (/ 3 * parini%nat /) ) !RESHAPE( tau,    (/ 3 * nat /) )
!        CALL cryst_to_cart( nat, pos, bg, -1 )
if(qe_units) then
        grad   = - RESHAPE( fcart, (/ 3 * parini%nat /) )*2.d0 !- RESHAPE( force,  (/ 3 * nat /) ) * alat
        CALL cryst_to_cart( parini%nat, grad, at, -1 )
!        call fxyz_cart2int(nat,fcart_in,grad,latvec_in)
!        grad=-grad
        !
!        IF ( lmovecell ) THEN
           at_old = at
           omega_old = omega
!           etot = etot + press * omega
           etot = enthalpy*2.d0
           sigma = sigma*2.d0
           CALL cell_force( fcell, - transpose(bg)/alat, sigma, omega, press )
           epsp1 = epsp / ry_kbar
!        END IF
        !
else
!        grad   = - RESHAPE( fcart, (/ 3 * nat /) ) !- RESHAPE( force,  (/ 3 * nat /) ) * alat
!        CALL cryst_to_cart( nat, grad, at, -1 )
        call fxyz_cart2int(parini%nat,fcart_in,grad,latvec_in)
        grad=-grad                                        
!        IF ( lmovecell ) THEN
           at_old = at
           omega_old = omega
!           etot = etot + press * omega
           etot = enthalpy                                
           sigma = sigma                                  
           CALL cell_force( fcell, - transpose(bg)/alat, sigma, omega, parini%target_pressure_habohr )
           epsp1 = epsp / ry_kbar
!        END IF
endif
!!!!NOT YET DONE
        fixion =  1 ! RESHAPE( if_pos, (/ 3 * nat /) )
        CALL bfgs( pos, h, etot, grad, fcell, fixion, tmp_dir, stdout, epse,&
                   epsf, epsp1,  energy_error, gradient_error, cell_error,  &
                   istep, nstep, step_accepted, conv_ions, lmovecell )
        
        
        !
        IF ( lmovecell ) THEN
           ! changes needed only if cell moves
!           if (fix_volume) call impose_deviatoric_strain(alat*at, h)
!           if (fix_area)   call impose_deviatoric_strain_2d(alat*at, h)
           at = h /alat  
           CALL recips( at(1,1),at(1,2),at(1,3), bg(1,1),bg(1,2),bg(1,3) )
           CALL qe_volume( alat, at(1,1), at(1,2), at(1,3), omega )
        END IF
        !
!        CALL cryst_to_cart( nat, pos, at, 1 )
        xred    =   RESHAPE( pos, (/ 3 , parini%nat /) )
        CALL cryst_to_cart( parini%nat, grad, bg, 1 )
        fcart = - RESHAPE( grad, (/ 3, parini%nat /) )
        !
!!!!!        IF ( conv_ions ) THEN
!!!!!!!!           !
!!!!!!!!           IF ( ( lsda .AND. ( absmag < eps8 ) .AND. lcheck_mag ) ) THEN
!!!!!!!!              !
!!!!!!!!              ! ... lsda relaxation :  a final configuration with zero 
!!!!!!!!              ! ...                    absolute magnetization has been found.
!!!!!!!!              !                        A check on this configuration is needed
!!!!!!!!              restart_with_starting_magnetiz = .true.
!!!!!!!!              ! 
!!!!!           IF (lmovecell.and.lcheck_cell) THEN
!!!!!              !
!!!!!              !  After the cell relaxation we make a final calculation
!!!!!              !  with the correct g vectors corresponding to the relaxed
!!!!!              !  cell.
!!!!!              !
!!!!!              final_cell_calculation=.TRUE.
!!!!!              CALL terminate_bfgs ( etot, epse, epsf, epsp, lmovecell, &
!!!!!                                    stdout, tmp_dir )
!!!!!              !
!!!!!           ELSE
!!!!!              !
!!!!!              CALL terminate_bfgs ( etot, epse, epsf, epsp, lmovecell, &
!!!!!                                    stdout, tmp_dir )
!!!!!              !
!!!!!           END IF
!!!!!           !
!!!!!        ELSE
!!!!!           !
           ! ... if a new bfgs step is done, new threshold is computed
           !
           IF ( step_accepted ) THEN
              !
              tr2  = starting_scf_threshold * &
                     MIN( 1.D0, ( energy_error / ( epse * upscale ) ), &
                                ( gradient_error / ( epsf * upscale ) ) )
              tr2  = MAX( ( starting_scf_threshold / upscale ), tr2 ) 
              !
           END IF
          !
           IF ( tr2 > 1.D-10 ) THEN
              WRITE( *, &
                     '(5X,"new conv_thr",T30,"= ",0PF18.10," Ry",/)' ) tr2
           ELSE
              WRITE( *, &
                     '(5X,"new conv_thr",T30,"= ",1PE18.1 ," Ry",/)' ) tr2
           END IF
!!!!!           !
!!!!!           ! ... the logical flag lcheck_mag is set again to .TRUE. (needed if 
!!!!!           ! ... a new configuration with zero absolute magnetization is 
!!!!!           ! ... identified in the following steps of the relaxation)
!!!!!           !
!!!!!!           lcheck_mag = .TRUE.
!!!!!           IF (lmovecell) lcheck_cell = .TRUE.
!!!!!           !
!!!!!        END IF
!!!!!        !
!!!!!!        CALL output_tau( lmovecell, conv_ions )
        !
if(iexit==1) then
    write(*,'(a,i5)') " # QBFGS converged in iterations", itime
    max_kpt=.false.
    exit
endif
!Set precision if necessary
       if(parini%usewf_geopt) then
           getwfk=.true.
       else
           getwfk=.false.
       endif
       if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
           getwfk=.false.
           iprec=1
       endif
!Reset everything, recompute cell and stuff
         if((multiprec.and.itime.ge.parini%paropt_geopt%nit/2).or.&
          &(fmax.lt.1.0d0*tolmxf_switch)) max_kpt=.true.
         if(fmax.lt.cellfix_switch.and..not.cellfix_done.and.(.not.(any(parini%fixlat).or.any(parini%fixat).or.confine.ge.1))) then
!Only perform the cell correction once, presumably close to the end of the optimization run
             call correct_latvec(h,pos,parini%nat,parini%correctalg,latvec_io)
             write(*,*) "New cell found", latvec_io
             cellfix_done=.true.
             if(latvec_io.ne.0) then
                max_kpt=.false.
                getwfk=.false.
                ka1=0;kb1=0;kc1=0
!Reset all
                conv_ions = .FALSE.
                prev_bfgs = .false.
             endif
         endif 
enddo
max_kpt=.false.
        DEALLOCATE( pos, grad, fixion )
        !
!     END IF bfgs_minimization
     !
     ! ... molecular dynamics schemes are used
     !
!     IF ( lmd ) THEN
!        !
!        IF ( calc == ' ' ) THEN
!           !
!           ! ... dynamics algorithms
!           !
!           IF ( ldamped ) THEN
!              !
!              CALL proj_verlet()
!              !
!           ELSE IF ( llang ) THEN
!              !
!              CALL langevin_md()
!              !
!           ELSE
!              !
!              CALL verlet()
!              !
!           END IF
!           !
!        ELSE IF ( calc /= ' ' ) THEN
!           !
!           ! ... variable cell shape md
!           !
!           CALL vcsmd()
!           !
!        END IF
!        !
!     END IF
!     !
!     ! ... before leaving check that the new positions still transform
!     ! ... according to the symmetry of the system.
!     !
!     CALL checkallsym( nat, tau, ityp, dfftp%nr1, dfftp%nr2, dfftp%nr3 )
     !
!  END IF

!  CALL mp_bcast(restart_with_starting_magnetiz,ionode_id,intra_image_comm)
!  CALL mp_bcast(final_cell_calculation,ionode_id,intra_image_comm)
  !
!!  IF ( final_cell_calculation ) THEN
!!     ! 
!!     ! ... Variable-cell optimization: once convergence is achieved, 
!!     ! ... make a final calculation with G-vectors and plane waves
!!     ! ... calculated for the final cell (may differ from the curent
!!     ! ... result, using G_vectors and PWs for the starting cell)
!!     !
!!     WRITE( UNIT = stdout, FMT = 9110 )
!!     WRITE( UNIT = stdout, FMT = 9120 )
!!     !
!!     CALL clean_pw( .FALSE. )
!!     CALL close_files(.TRUE.)
!!     lmovecell=.FALSE.
!!     lcheck_cell=.FALSE.
!!     final_cell_calculation=.FALSE.
!!     lbfgs=.FALSE.
!!     lmd=.FALSE.
!!     lcheck_mag = .FALSE.
!!     restart_with_starting_magnetiz = .FALSE.
!!     ! ... conv_ions is set to .FALSE. to perform a final scf cycle
!!     conv_ions = .FALSE.
!!     ! ... allow re-calculation of FFT grid
!!     !
!!     dfftp%nr1=0; dfftp%nr2=0; dfftp%nr3=0; dffts%nr1=0; dffts%nr2=0; dffts%nr3=0
!!     CALL realspace_grids_init (dfftp, dffts,at, bg, gcutm, gcutms )
!!     CALL init_run()
!!     !
!!  ELSE IF (restart_with_starting_magnetiz) THEN
!!     !
!!     ! ... lsda optimization :  a final configuration with zero 
!!     ! ... absolute magnetization has been found and we check 
!!     ! ... if it is really the minimum energy structure by 
!!     ! ... performing a new scf iteration without any "electronic" history
!!     !
!!     WRITE( UNIT = stdout, FMT = 9010 )
!!     WRITE( UNIT = stdout, FMT = 9020 )
!!     !
!!     lcheck_mag = .FALSE.
!!     restart_with_starting_magnetiz = .FALSE.
!!     ! ... conv_ions is set to .FALSE. to perform a final scf cycle
!!     conv_ions = .FALSE.
!!     !
!!     ! ... re-initialize the potential and wavefunctions
!!     !
!!     CALL potinit()
!!     CALL newd()
!!     CALL wfcinit()
!!     !
!!  END IF
  !
  ! ... broadcast calculated quantities to all nodes
  !
!  CALL mp_bcast( istep,     ionode_id, intra_image_comm )
!  CALL mp_bcast( tau,       ionode_id, intra_image_comm )
!  CALL mp_bcast( force,     ionode_id, intra_image_comm )
!  CALL mp_bcast( tr2,       ionode_id, intra_image_comm )
!  CALL mp_bcast( conv_ions, ionode_id, intra_image_comm )
!  CALL mp_bcast( history,   ionode_id, intra_image_comm )
  !
!  IF ( lmovecell ) THEN
!     !
!     CALL mp_bcast( at,        ionode_id, intra_image_comm )
!     CALL mp_bcast( at_old,    ionode_id, intra_image_comm )
!     CALL mp_bcast( omega,     ionode_id, intra_image_comm )
!     CALL mp_bcast( omega_old, ionode_id, intra_image_comm )
!     CALL mp_bcast( bg,        ionode_id, intra_image_comm )
!     !
!  END IF
  !
  RETURN

9010 FORMAT( /5X,'lsda relaxation :  a final configuration with zero', &
           & /5X,'                   absolute magnetization has been found' )
9020 FORMAT( /5X,'the program is checking if it is really ', &
           &     'the minimum energy structure',             &
           & /5X,'by performing a new scf iteration ',       & 
           &     'without any "electronic" history' )               
  !
9110 FORMAT( /5X,'A final scf calculation at the relaxed structure.' )
9120 FORMAT(  5X,'The G-vectors are recalculated for the final unit cell'/ &
              5X,'Results may differ from those at the preceding step.' )
  !
END SUBROUTINE 
  !
!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------

subroutine recips (a1, a2, a3, b1, b2, b3)
  !---------------------------------------------------------------------
  !
  !   This routine generates the reciprocal lattice vectors b1,b2,b3
  !   given the real space vectors a1,a2,a3. The b's are units of 2 pi/a.
  !
  !     first the input variables
  !
!  use kinds, ONLY: DP
  implicit none
  real(8) :: a1 (3), a2 (3), a3 (3), b1 (3), b2 (3), b3 (3)
  ! input: first direct lattice vector
  ! input: second direct lattice vector
  ! input: third direct lattice vector
  ! output: first reciprocal lattice vector
  ! output: second reciprocal lattice vector
  ! output: third reciprocal lattice vector
  !
  !   then the local variables
  !
  real(8) :: den, s
  ! the denominator
  ! the sign of the permutations
  integer :: iperm, i, j, k, l, ipol
  ! counter on the permutations
  !\
  !  Auxiliary variables
  !/
  !
  ! Counter on the polarizations
  !
  !    first we compute the denominator
  !
  den = 0
  i = 1
  j = 2
  k = 3
  s = 1.d0
100 do iperm = 1, 3
     den = den + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s
  if (s.lt.0.d0) goto 100
  !
  !    here we compute the reciprocal vectors
  !
  i = 1
  j = 2
  k = 3
  do ipol = 1, 3
     b1 (ipol) = (a2 (j) * a3 (k) - a2 (k) * a3 (j) ) / den
     b2 (ipol) = (a3 (j) * a1 (k) - a3 (k) * a1 (j) ) / den
     b3 (ipol) = (a1 (j) * a2 (k) - a1 (k) * a2 (j) ) / den
     l = i
     i = j
     j = k
     k = l
  enddo
  return
end subroutine recips
!
! Copyright (C) 2001-2003 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!-----------------------------------------------------------------------
subroutine cryst_to_cart (nvec, vec, trmat, iflag)
  !-----------------------------------------------------------------------
  !
  !     This routine transforms the atomic positions or the k-point
  !     components from crystallographic to cartesian coordinates 
  !     ( iflag=1 ) and viceversa ( iflag=-1 ).
  !     Output cartesian coordinates are stored in the input ('vec') array
  !
  !
!  USE kinds, ONLY : DP
  implicit none
  !
  integer, intent(in) :: nvec, iflag
  ! nvec:  number of vectors (atomic positions or k-points)
  !        to be transformed from crystal to cartesian and vice versa
  ! iflag: gives the direction of the transformation
  real(8), intent(in) :: trmat (3, 3)
  ! trmat: transformation matrix
  ! if iflag=1:
  !    trmat = at ,  basis of the real-space lattice,       for atoms   or
  !          = bg ,  basis of the reciprocal-space lattice, for k-points
  ! if iflag=-1: the opposite
  real(8), intent(inout) :: vec (3, nvec)
  ! coordinates of the vector (atomic positions or k-points) to be
  ! transformed - overwritten on output
  !
  !    local variables
  !
  integer :: nv, kpol
  ! counter on vectors
  ! counter on polarizations
  real(8) :: vau (3)
  ! workspace
  !
  !     Compute the cartesian coordinates of each vectors
  !     (atomic positions or k-points components)
  !
  do nv = 1, nvec
     if (iflag.eq.1) then
        do kpol = 1, 3
           vau (kpol) = trmat (kpol, 1) * vec (1, nv) + trmat (kpol, 2) &
                * vec (2, nv) + trmat (kpol, 3) * vec (3, nv)
        enddo
     else
        do kpol = 1, 3
           vau (kpol) = trmat (1, kpol) * vec (1, nv) + trmat (2, kpol) &
                * vec (2, nv) + trmat (3, kpol) * vec (3, nv)
        enddo
     endif
     do kpol = 1, 3
        vec (kpol, nv) = vau (kpol)
     enddo
  enddo
  !
  return
end subroutine cryst_to_cart


  subroutine cell_force( fcell, ainv, stress, omega, press)!, wmassIN )
    USE constants, ONLY : eps8
    REAL(8), intent(out) :: fcell(3,3)
    REAL(8), intent(in) :: stress(3,3), ainv(3,3)
    REAL(8), intent(in) :: omega, press
!    REAL(8), intent(in), optional :: wmassIN
    integer        :: i, j
    REAL(8) :: wmass
!    IF (.not. present(wmassIN)) THEN
      wmass = 1.0
!    ELSE
!      wmass = wmassIN
!    END IF
    do j=1,3
      do i=1,3
        fcell(i,j) = ainv(j,1)*stress(i,1) + ainv(j,2)*stress(i,2) + ainv(j,3)*stress(i,3)
      end do
    end do
    do j=1,3
      do i=1,3
        fcell(i,j) = fcell(i,j) - ainv(j,i) * press
      end do
    end do
    IF( wmass < eps8 ) &
!       CALL errore( ' movecell ',' cell mass is less than 0 ! ', 1 )
       write(*,*) ' movecell ',' cell mass is less than 0 ! ', 1 
    fcell = omega * fcell / wmass
    return
  end subroutine cell_force

!
! Copyright (C) 2004 Quantum ESPRESSO group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
subroutine invmat (n, a, a_inv, da)
  !-----------------------------------------------------------------------
  ! computes the inverse "a_inv" of matrix "a", both dimensioned (n,n)
  ! if the matrix is dimensioned 3x3, it also computes determinant "da"
  ! matrix "a" is unchanged on output - LAPACK
  !
!  USE kinds, ONLY : DP
  implicit none
  integer :: n
  real(8), DIMENSION (n,n) :: a, a_inv
  real(8) :: da
  !
  integer :: info, lda, lwork, ipiv (n)
  ! info=0: inversion was successful
  ! lda   : leading dimension (the same as n)
  ! ipiv  : work space for pivoting (assumed of length lwork=n)
  real(8) :: work (n) 
  ! more work space
  !
  lda = n
  lwork=n
  !
  a_inv(:,:) = a(:,:)
  !
  call dgetrf (n, n, a_inv, lda, ipiv, info)
!  call errore ('invmat', 'error in DGETRF', abs (info) )
  if(info.ne.0) write(*,*) 'invmat', 'error in DGETRF', abs (info) 
  call dgetri (n, a_inv, lda, ipiv, work, lwork, info)
!  call errore ('invmat', 'error in DGETRI', abs (info) )
  if(info.ne.0) write(*,*) 'invmat', 'error in DGETRI', abs (info) 
  !
  if (n == 3) then
     da = a(1,1)*(a(2,2)*a(3,3)-a(2,3)*a(3,2)) + &
          a(1,2)*(a(2,3)*a(3,1)-a(2,1)*a(3,3)) + &
          a(1,3)*(a(2,1)*a(3,2)-a(3,1)*a(2,2))
!     IF (ABS(da) < 1.d-10) CALL errore(' invmat ',' singular matrix ', 1)
     IF (ABS(da) < 1.d-10) write(*,*)' invmat ',' singular matrix ', 1
  else
     da = 0.d0
  end if

  return
end subroutine invmat

!
! Copyright (C) 2001 PWSCF group
! This file is distributed under the terms of the
! GNU General Public License. See the file `License'
! in the root directory of the present distribution,
! or http://www.gnu.org/copyleft/gpl.txt .
!
!
!---------------------------------------------------------------------
subroutine qe_volume (alat, a1, a2, a3, omega)
  !---------------------------------------------------------------------
  !
  !     Compute the volume of the unit cell
  !
!  use kinds, ONLY: DP
  implicit none
  !
  !     First the I/O variables
  !
  real(8) :: alat, a1 (3), a2 (3), a3 (3), omega
  ! input:  lattice parameter (unit length)
  ! input: the first lattice vector
  ! input: the second lattice vector
  ! input: the third lattice vector
  ! input: the volume of the unit cell
  !
  !    Here the local variables required by the routine
  !

  real(8) :: s
  ! the sign of a permutation
  integer :: i, j, k, l, iperm
  !\
  ! \
  ! /   auxiliary indices
  !/
  ! counter on permutations
  !
  !   Compute the volume
  !
  omega = 0.d0
  s = 1.d0
  i = 1
  j = 2
  k = 3
101 do iperm = 1, 3
     omega = omega + s * a1 (i) * a2 (j) * a3 (k)
     l = i
     i = j
     j = k
     k = l
  enddo
  i = 2
  j = 1
  k = 3
  s = - s

  if (s.lt.0.d0) goto 101

  omega = abs (omega) * alat**3
  return
end subroutine qe_volume
