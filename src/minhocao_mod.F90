module global
implicit none
!Set of parms.in
!  real(8):: target_pressure_gpa !Target pressures
!  real(8):: target_pressure_habohr  !Target pressures
!  integer:: nat                   !Number of atoms
  integer:: nmol                  !Number of molecules
!  integer:: ntypat                !Number of atom types
!  integer,allocatable:: fragarr(:)
!  integer,allocatable:: fragsize(:)
!  integer,allocatable:: lhead(:)
!  integer,allocatable:: llist(:) 
!  real(8),allocatable:: znucl(:)
!  real(8),allocatable:: amu(:)
!  real(8),allocatable:: rcov(:) 
!  integer,allocatable:: typat(:)
!  logical,allocatable:: fixat(:)
!  logical:: fixlat(7)             !Contains the information of the cell constraints: a,b,c,alpha,beta,gamma,cellshape
!  integer:: ntime_md              !Maximum number of iterations during MD
!  integer:: ntime_geopt           !Maximum number of iterations during GEOPT
!  real(8):: bmass                 !Cell mass during MD and FIRE
!  integer:: mdmin                 !Number of enthalpy minima crossed unit stop MD
!  integer:: mdmin_min,mdmin_max   !min,max number of enthalpy minima crossed unit stop MD, only if automatically determined
!  real(8):: dtion_md              !Initial timestep for MD 
  !real(8):: dtion_fire,dtion_fire_min,dtion_fire_max      !Initial timestep for FIRE, minimal_timestep, maximal_timestep
  real(8):: dtion_fire_min,dtion_fire_max      !Initial timestep for FIRE, minimal_timestep, maximal_timestep
!  real(8):: tolmxf                !Force tolerance for GEOPT convergance 
!  real(8):: strfact               !Factor to multiply stress 
  character(40):: units           !Either angstroem or bohr
!  integer:: ka,kb,kc              !The number of kpoints in each dimension
!  integer:: siesta_kpt_mode       !If 1, the kpoint mesh is defined by cutoff length, else a monkhorst pack mesh is generated (only siesta)
!  integer:: vasp_kpt_mode         !If 1, the kpoint mesh is defined by mesh length, else a monkhorst pack mesh is generated (only vasp)
!  integer:: abinit_kpt_mode       !If 1, the kpoint mesh is defined by kptrlen length, else a monkhorst pack mesh is generated (only abinit)
!  character(2),allocatable:: char_type(:) 
!  real(8):: dkpt1,dkpt2           !Precisions of the kpt mesh if generated automatically
!  character(20):: code            !What code should be used: abinit or siesta
!  logical::  geopt_ext            !At the moment only used for siesta: if true, the external geometry optimizer is used
!  real(8):: alpha_lat, alpha_at   !Stepsize for softening the atomic and lattice coordinates
!  real(8):: alphax_lat , alphax_at !Stepsize for BFGS of the atomic and lattice coordinates
!  integer:: nsoften               !Number of softening steps
!  logical:: usewf_md,usewf_geopt,usewf_soften !Defines when the wavefunctions should be reused in the next step
!  character(5):: geopt_method
  integer:: ka1,kb1,kc1           !The previously used kpt mesh are stored in these variables, only abinit 
  logical:: max_kpt               !If true, the single point in abinit will first evaluate if the new set is better and choose the better one. Default=false
  logical:: reuse_kpt             !If true, the single point in abinit will reuse previous kpt mesh. Default=false
  logical:: reduced               !If true, all output files will be written in reduced coordinates. Initiallized when reading poscur.ascii
!  logical:: findsym               !If true, findsym will be used to get symmetry informations on the fly
!  logical:: finddos               !If true, the DOS at the Fermi level will be evaluated at the end of every geometry optimization 
!  logical:: auto_soft             !If true, the softening stepsize will be adjusted during run
!  logical:: auto_mdmin            !If true, the mdmin parameter will be adjusted during run
!  logical:: auto_dtion_md         !If true, the timestep during MD will be adjusted during run
!  logical:: auto_kpt              !Currently a dummy variable
!  integer:: nit_per_min           !Target number of md steps per md minimum crossing
!  integer:: md_algo               !Algorithm for VCMD: 1=PR, 2=Cleveland, 3=Wentzcovitch
!  integer:: md_integrator         !Integrator for VCMD: 1=Verlet, 2=Velocity-Verlet, 3=Beeman
!  real(8):: md_presscomp          !Pressure compensation during MD by substracting the kinetic energy pressure from the external pressure
!  logical:: mol_soften            !Switch on molecular softening
!  integer:: correctalg            !Method to perform cell corrections
!  integer:: bc                    !1: periodic, 2:free, 3:surface/slab
!  integer:: verb                  !0: very little output, 1: normal output, 2: folders for geopt and md, 3: output stress and forces
  integer:: confine               !0: No confinement, 1: confinement used, but not currently, 2: confinement in action, 3: confinement always on
!  logical:: use_confine           !if true, confinement is enable, otherwise disabled
!  logical:: energy_conservation   !Only used in fixed cell MD
!  logical:: voids                 !If or if not to use void creating LJ particles in the cell
!  logical:: core_rep              !If or if not to add a purely repulsive force on top of the atoms
end module global

module mod_fire
implicit none
   integer,parameter:: Nmin=5
   real(8),parameter:: finc=1.1d0
   real(8),parameter:: fdec=0.5d0
   real(8),parameter:: alphastart=2.d-1
   real(8),parameter:: falpha=0.99d0
   real(8),parameter:: latmass_rel_fire=5.d-1
   real(8),parameter:: latmass_at=1.d0
   real(8):: dtmax    !=5.d+1, will be provided trough params.in
   real(8):: dtmin    !=1.d0,  will be provided trough params.in
end module mod_fire
!> Module which define the type parameterminimization
module minpar
   implicit none

   type parameterminimization
      !> general parameters for all methods
      character (len=10) :: approach
      integer :: iter
      integer :: iflag
      integer :: history
      logical ::converged
      !>parameters for print information
      integer :: verbosity
      integer :: MSAVE
      integer :: MP
      integer :: LP
      integer :: MAXFEV
      integer :: FINSTEP
      integer :: ncount_cluster_x
      integer :: maxiter_lat
      double precision :: ALPHA 
      double precision :: GTOL
      double precision :: XTOL
      double precision :: FTOL
      double precision :: STPMIN
      double precision :: STPMAX
      double precision :: BETAX 
      double precision :: BETAX_LAT
      double precision :: frac_fluct
      double precision :: forcemax
      logical :: DIAGCO
      logical :: IWRITE
   end type parameterminimization

   type(parameterminimization) :: parmin_bfgs

end module minpar

module fingerprint 
   implicit none
   save
!All
!   real(8):: fp_rcut   
!   integer:: fp_nl
   integer:: fp_method
!   integer:: fp_at_nmax
   integer:: fp_all_nmax
!   character(20):: fp_method_ch
!   real(8):: fp_sigma   
!   real(8):: fp_dbin
!Oganov parameters
   real(8):: fp_11_rcut
   real(8):: fp_11_sigma   
   real(8):: fp_11_dbin
   integer:: fp_11_fp_size
   integer:: fp_11_fp_dim
   integer,allocatable:: fp_11_nkinds_sum(:)
!CALYPSO parameters
   integer:: fp_12_nl
   integer:: fp_12_fp_dim
   real(8),allocatable:: fp_12_r_cut(:)
!Modified CALYPSO parameters
   integer:: fp_13_nl
   integer:: fp_13_fp_dim
   real(8),allocatable:: fp_13_r_cut(:)
!xyz2sm parameters
!   integer:: fp_14_m
!   real(8):: fp_14_w1
!   real(8):: fp_14_w2
!Continuous Oganov parameters
   real(8):: fp_15_rcut
   real(8):: fp_15_sigma  
   integer:: fp_15_fp_size
   integer:: fp_15_fp_dim
   integer,allocatable:: fp_15_nkinds_sum(:)
!Continuous Atomic Oganov parameters
   real(8):: fp_16_rcut
   real(8):: fp_16_sigma  
   integer:: fp_16_fp_size
   integer:: fp_16_fp_dim
   integer,allocatable:: fp_16_nkinds_sum(:)
! gom parameters
!   real(8):: fp_17_width_cutoff
!   real(8):: fp_17_nex_cutoff
!   integer:: fp_17_natx_sphere
!   integer:: fp_17_lseg
!   character(len=2) :: fp_17_orbital
!Molecular GOM
!
!
!end module
!!Need to put it into the FP parameters module
!module parameter_molgom
!   implicit none
!   save
!   character:: fp_18_orbital
!   logical, parameter :: write_files = .false.
!   logical, parameter :: clustering = .false.
!   integer :: fp_18_cluster_number = 20
!   integer, parameter :: nat=20
!   integer, parameter :: ntypat=4
!   integer :: fp_18_principleev = 6
!   integer :: fp_18_lseg!=1
!   integer, parameter :: nconf=177
!   integer :: fp_18_molecules=4
!   integer :: fp_18_expaparameter = 4
!   integer :: fp_18_nex_cutoff = 3
!   integer :: fp_18_molecules_sphere = 50
!   real*8  :: fp_18_width_cutoff = 1.d0
!   real*8  :: fp_18_width_overlap = 1.d0
   real*8  :: fp_18_large_vanradius = 1.7d0/0.52917720859d0
!   real(8),allocatable  :: rvan(:) !nat*molecules)
!   character(len=2),allocatable:: finalchar(:) ! dimension(nat*molecules) 
end module !parameter_molgom

!module confinement
!implicit none
!Confinement parameters
!integer:: nconfine       !number of different confinements
!integer,allocatable:: conf_dim(:)    !1,2 or 3 for each of the 3 dimensions latvec(:,1), latvec(:,2), latvec(:,3)
!integer,allocatable:: conf_av(:)     !0: no confinement (should never occur)
!                                     !1: confinement with respect to a fixed value along latvec(:,i)
!                                     !2: confinement with respect to the average
!integer,allocatable:: conf_exp(:)    !The polynomial order for each confinement
!real(8),allocatable:: conf_prefac(:) !The polynomial predactor for each confinement
!real(8),allocatable:: conf_cut(:)    !The cutoff distance from each confinement equilibrium
!real(8),allocatable:: conf_eq(:)     !Equlibrium position of confinement along the confinement direction, will be filled to average or fixed value
!integer,allocatable:: conf_nat(:)    !How many atoms per confinement
!integer,allocatable:: conf_list(:,:) !List of atoms per confinement
!character(1),allocatable:: conf_cartred(:)!Cartesian or reduced coordinates, only if conf_eq is provided. d,D,r,R for reduced, C,c,K,k for cartesian
!end module

module blj_params
     implicit none
     save
     real(8):: sigmalj(2,2)=1.d0 !sigma(i,j) is the sigma between particle i and particle j
     real(8):: epslj(2,2)=1.d0   !eps(i,j) is the epsilon between particle i and particle j
     real(8):: rcut(2,2)         !rcut(i,j) is the cutoff distance between particle i and j
     real(8):: alpha_lj          !cutoff discatnce
end module

module mlj_params
     implicit none
     save
     real(8),allocatable:: sigmamlj(:,:) !sigma(i,j) is the sigma between particle i and particle j
     real(8),allocatable:: epsmlj(:,:)   !eps(i,j) is the epsilon between particle i and particle j
     real(8),allocatable:: rcutmlj(:,:)    !rcut(i,j) is the cutoff distance between particle i and j
     real(8),allocatable:: alphamlj(:,:) !cutoff discatnce
end module

module tersoff_params
     implicit none
     save
     real(8):: rcut1(2,2)         !cutoff
     real(8):: rcut2(2,2)         !cutoff
     integer,allocatable:: Kinds_tersoff(:)   !Kinds, 1 for carbon and 2 for silicon
     logical:: only_c
end module

module void_lj_params
     implicit none
     save
     real(8),allocatable:: sigmavoidlj(:) !sigma(i) is the sigma between then i-th LJ particle and an atom
     real(8),allocatable:: epsvoidlj(:)   !eps(i) is the epsilon between then i-th LJ particle and an atom
     real(8),allocatable:: rcutvoidlj(:)  !rcut(i) is the cutoff distance for the i-th particle
     real(8):: alpha_lj                 !cutoff distance
     integer:: nat_atoms,nat_lj,ntypat_atoms,ntypat_lj
end module

module tb_lj_params
     implicit none
     save
     real(8),allocatable:: sigmatblj(:) !sigma(i) is the sigma between then i-th LJ particle a TB particle
     real(8),allocatable:: epstblj(:)   !eps(i) is the epsilon between then i-th LJ particle a TB particle
     real(8),allocatable:: rcuttblj(:)  !rcut(i) is the cutoff distance for the i-th particle
     real(8):: alpha_lj                 !cutoff distance
     integer:: n_silicon,n_h,n_lj,n_tb
end module

MODULE String_Utility 
   IMPLICIT NONE 
   PRIVATE 
   PUBLIC :: StrUpCase 
   PUBLIC :: StrLowCase 
   CHARACTER( * ), PRIVATE, PARAMETER :: LOWER_CASE = 'abcdefghijklmnopqrstuvwxyz' 
   CHARACTER( * ), PRIVATE, PARAMETER :: UPPER_CASE = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' 
CONTAINS 
   FUNCTION StrUpCase ( Input_String ) RESULT ( Output_String ) 
     ! -- Argument and result 
     CHARACTER( * ), INTENT( IN )     :: Input_String 
     CHARACTER( LEN( Input_String ) ) :: Output_String 
     ! -- Local variables 
     INTEGER :: i, n 
     ! -- Copy input string 
     Output_String = Input_String 
     ! -- Loop over string elements 
     DO i = 1, LEN( Output_String ) 
       ! -- Find location of letter in lower case constant string 
       n = INDEX( LOWER_CASE, Output_String( i:i ) ) 
       ! -- If current substring is a lower case letter, make it upper case 
       IF ( n /= 0 ) Output_String( i:i ) = UPPER_CASE( n:n ) 
     END DO 
   END FUNCTION StrUpCase 
   FUNCTION StrLowCase ( Input_String ) RESULT ( Output_String ) 
     ! -- Argument and result 
     CHARACTER( * ), INTENT( IN )     :: Input_String 
     CHARACTER( LEN( Input_String ) ) :: Output_String 
     ! -- Local variables 
     INTEGER :: i, n 
     ! -- Copy input string 
     Output_String = Input_String 
     ! -- Loop over string elements 
     DO i = 1, LEN( Output_String ) 
       ! -- Find location of letter in upper case constant string 
       n = INDEX( UPPER_CASE, Output_String( i:i ) ) 
       ! -- If current substring is an upper case letter, make it lower case 
       IF ( n /= 0 ) Output_String( i:i ) = LOWER_CASE( n:n ) 
     END DO 
   END FUNCTION StrLowCase 
END MODULE String_Utility 

!!module modipi
!!implicit none
!!      INTEGER:: ipi_socket, ipi_inet, ipi_port        ! socket ID & address of the server
!!      CHARACTER(LEN=1024) :: ipi_host
!!      INTEGER, PARAMETER :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
!!      CHARACTER(LEN=60)::ipi_sock_extra_string="                                                            "
!!      real(8):: ipi_ecutwf(2)
!!end module modipi

module modsocket
implicit none
!      INTEGER:: sock_socket, sock_inet, sock_port        ! socket ID & address of the socket
      INTEGER:: sock_socket
!      CHARACTER(LEN=1024) :: sock_host
      INTEGER, PARAMETER  :: MSGLEN=12   ! length of the headers of the driver/wrapper communication protocol
      CHARACTER(LEN=60)   :: sock_extra_string="                                                            "
!      real(8)             :: sock_ecutwf(2)
end module modsocket

!module mod_sqnm
!implicit none
!      real(8):: sqnm_beta_lat
!      real(8):: sqnm_beta_at
!      !integer:: sqnm_nhist
!      real(8):: sqnm_maxrise
!      real(8):: sqnm_cutoffRatio
!      real(8):: sqnm_steepthresh             
!      real(8):: sqnm_trustr
!end module mod_sqnm

module steepest_descent
implicit none
      real(8):: sd_beta_lat
      real(8):: sd_beta_at
end module steepest_descent

!module qbfgs
!implicit none
!     !integer::qbfgs_bfgs_ndim!=1
!     real(8)::qbfgs_trust_radius_max!=0.5d0
!     real(8)::qbfgs_trust_radius_min!=1.d-5
!     !real(8)::qbfgs_trust_radius_ini!=0.02D0
!     real(8)::qbfgs_w_1!=0.05D0
!     real(8)::qbfgs_w_2!=0.5D0
!end module


