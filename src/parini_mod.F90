!*****************************************************************************************
module mod_parini
    use dictionaries
    use mod_opt, only: typ_paropt
    implicit none
    type typ_parini
        logical:: exists_yaml_file
        integer:: iunit
        !-----------------------------------------------------------------------
        !parameters of [main]
        character(50):: task='unknown'
        character(256):: cwd='unknown'
        character(50):: rng_type='intrinsic'
        logical:: two_level_geopt=.false.
        integer:: iverbose=0
        integer:: iseed=-2
        character(50):: types_main='unknown'
        integer:: ntypat=-1
        integer:: ltypat(20)=-1
        integer:: iatomnum(20)=-1
        integer:: nrun_lammps=0
        character(5):: stypat(20)='unknown'
        logical:: params_new=.false.
        character(100):: str_typat_global
        !character(50):: stypat_genconf=''
        !-----------------------------------------------------------------------
        !parameters of [minhopp]
        logical:: avail_minhopp=.false.
        integer:: nstep_minhopp=0
        integer:: nsoften_minhopp=7
        integer:: mdmin_minhopp=3
        integer:: minter_minhopp=1
        integer:: nrandoff_minhopp=0
        integer:: npminx_minhopp=1000
        real(8):: etoler_minhopp=1.d-2
        real(8):: eref_minhopp=-1.d50
        real(8):: ekinmax_minhopp=1.d20
        real(8):: alpha1_minhopp=1.d0/1.02d0
        real(8):: alpha2_minhopp=1.02d0
        real(8):: beta1_minhopp=1.05d0
        real(8):: beta2_minhopp=1.05d0
        real(8):: beta3_minhopp=1.d0/1.05d0
        logical:: trajectory_minhopp=.false.
        logical:: print_force_minhopp=.false.
        !-----------------------------------------------------------------------
        !parameters of [geopt]
        logical:: avail_geopt=.false.
        type(typ_paropt):: paropt_geopt
        !-----------------------------------------------------------------------
        !parameters of [geopt_prec]
        logical:: avail_geopt_prec=.false.
        type(typ_paropt):: paropt_geopt_prec
        !-----------------------------------------------------------------------
        !parameters of [saddle_opt]
        logical:: avail_saddle_opt=.false.
        type(typ_paropt):: paropt_saddle_opt
        !-----------------------------------------------------------------------
        !parameters of [bader]
        logical:: avail_bader=.false.
        character(50):: filename_bader='total_density.cube'
        character(50):: approach_bader='unknown'
        character(50):: vacuum_bader='yes'
        !-----------------------------------------------------------------------
        !parameters of [potential]
        logical:: avail_potential=.false.
        character(50):: potential_potential='unknown'
        character(50):: potential_potential_sec='unknown'
        character(50):: potential_ann_boundcheck='none'
        character(256):: component_ff='no'
        logical:: drift_potential=.false.
        logical:: cal_charge= .false.
        logical:: add_repulsive= .true.
        !-----------------------------------------------------------------------
        !parameters of [ann]
        logical:: avail_ann=.false.
        logical:: bondbased_ann=.false.
        character(50):: subtask_ann='unknown'
        character(50):: optimizer_ann='unknown'
        character(50):: approach_ann='atombased'
        character(50):: syslinsolver_ann='direct'
        character(50):: symfunc_type_ann='behler'
        character(50):: symfunc='only_calculate'
        integer:: nstep_opt_ann=100
        integer:: nstep_cep=200
        integer:: nconf_rmse=0
        real(8):: alphax_q=1.d0
        real(8):: alphax_r=5.d-2
        real(8):: ampl_rand=1.d0
        real(8):: rgnrmtol=-1.d0
        real(8):: qgnrmtol=-1.d0
        real(8):: etol_ann !the tolerance difference of energies of two configuration
        real(8):: dtol_ann !distance between two FP
        real(8):: ftol_ann !tolerance for simplex method
        real(8):: weight_hardness
        logical:: normalization_ann=.false.
        logical:: prefit_ann=.false.
        logical:: prefit_centt_ann=.false.
        logical:: read_forces_ann
        logical:: restart_param=.false. 
        integer:: restart_iter=0  
        logical:: print_energy=.false. 
        logical:: fit_hoppint=.false. 
        logical:: save_symfunc_force_ann=.false.
        logical:: save_symfunc_behnam=.false.
        logical:: free_bc_direct=.false.
        !-----------------------------------------------------------------------
        !parameters of [saddle]
        character(50):: method_saddle='unknown'
        logical:: avail_saddle=.false.
        character(256):: str_moving_atoms_rand_saddle
        real(8):: dimsep_saddle=-1.d0
        real(8):: ampl_saddle=-1.d0
        real(8):: dbar
        real(8):: alphax_bs
        real(8):: fnrmtol_coarse
        real(8):: contr_dbar
        real(8):: fnrmtol_contracted
        logical:: bar_contract
        integer:: nstep_contract
        integer:: nstep_bs
        integer:: np_splsad !np-1 is the number of movable anchor points
        integer:: np_neb
        integer:: ns2_splsad !number of extra points along the path, beginning of maximization
        integer:: max_fcalls !maximum number of force calls
        real(8):: vdtol_splsad !tolerance for the derivative of potential at maximum point
        real(8):: dt_saddle
        real(8):: htol_splsad !minimal distance to accept new point during maximization
        real(8):: alphax_saddle !reference step size for geometry optimization
        real(8):: fmaxtol_splsad !tolerance for maximum force component for splsad
        real(8):: fmaxtol_neb !tolerance for maximum force component for NEB
        character(20):: hybrid_splsad
        character(20):: doneb
        character(20):: docineb
        character(20):: typintpol !interpolation method for the maximization (cubic or quintic)
        character(20):: pickbestanchorpoints
        character(20):: runstat !new or restart (anchorposinp.xyz is needed if restart)
        character(10):: opt_method_saddle !SD or SDCG or SDDIIS or LBFGS or FIRE
        !-----------------------------------------------------------------------
        !parameters of [dynamics]
        logical:: avail_dynamics=.false.
        real(8):: dt_dynamics=-1.d0
        real(8):: temp_dynamics
        real(8):: init_temp_dynamics =0.d0
        real(8):: highest_frequency  =10.d0
        integer:: ntherm  = 2
        integer:: nmd_dynamics=0
        integer:: nfreq_dynamics=0
        character(20):: md_method_dynamics='unknown'
        logical:: print_force_dynamics=.false.
        logical:: restart_dynamics=.false.
        logical:: fix_cm_dynamics=.false.
        logical:: vflip_dynamics=.false.
        logical:: wall_repulsion_dynamics=.false.
        real(8):: time_dynamics = 0.d0
        !-----------------------------------------------------------------------
        !parameters of [genconf]
        logical:: avail_genconf=.false.
        character(50):: subtask_genconf='unknown'
        integer:: nat_add_genconf=0
        character(10):: cal_pot_genconf='no'
        character(10):: sat_genconf='unknown'
        real(8):: amargin_genconf=0.d0
        real(8):: dmin_genconf=-1.d0
        real(8):: dmax_genconf=-1.d0
        integer:: npoint_genconf=0
        real(8):: fbmin_genconf=-1.d0
        real(8):: fbmax_genconf=-1.d0
        logical :: variable_cell_genconf= .false.
        logical :: nonorthogonal_genconf= .false.
        !-----------------------------------------------------------------------
        !parameters of [conf_comp]
        logical:: avail_conf_comp=.false.
        real(8):: tol_conf_comp
        !-----------------------------------------------------------------------
        !parameters of [testforces]
        logical:: avail_testforces=.false.
        character(50):: testforces_method='unknown'
        !-----------------------------------------------------------------------
        !parameters of [single_point]
        logical:: avail_single_point=.false.
        logical:: print_force_single_point=.false.
        character(50):: frmt_single_point='unknown'
        !-----------------------------------------------------------------------
        !parameters of [ewald]
        integer:: nsp_ewald=-1
        real(8):: alpha_ewald=-1.d0
        real(8):: ecut_ewald=-1.d0
        real(8):: ecutz_ewald=-1.d0
        real(8):: rcut_ewald=-1.d0
        real(8):: rgcut_ewald=-1.d0
        real(8):: vu_ewald=0.d0
        real(8):: vl_ewald=0.d0
        real(8):: vu_ac_ewald=0.d0
        real(8):: frequency_ewald=0.d0
        real(8):: gnrmtol_eem=1.d-7
        real(8):: tolerance_ewald = 1.d-6
        real(8):: efield !external electric field
        real(8):: dielec_const
        real(8):: dielec_const1
        real(8):: dielec_const2
        logical :: ewald=.false.
        logical :: cell_ortho=.false.
        character(256):: bias_type='no'
        character(50):: psolver='unknown'
        logical:: cal_polar= .false.
        logical:: ecut_auto=.false.
        !-----------------------------------------------------------------------
        !parameters of [misc]
        logical:: avail_misc=.false.
        character(50):: subtask_misc='unknown'
        real(8):: gaussian_width=-1.d0
        character(50):: boundcond_misc='unknown'
        character(50):: posinp_misc='unknown'
        !-----------------------------------------------------------------------
        !Socket
        logical:: usesocket=.false.                 ! Use sockets to send the the results to
        integer:: inisock_inet, inisock_port        ! socket ID & address of the socket
        CHARACTER(LEN=1024) :: inisock_host
        !-----------------------------------------------------------------------
        !Minhocao global module
        integer:: verb                  !0: very little output, 1: normal output, 2: folders for geopt and md, 3: output stress and forces
        integer:: md_algo               !Algorithm for VCMD: 1=PR, 2=Cleveland, 3=Wentzcovitch
        real(8):: md_presscomp          !Pressure compensation during MD by substracting the kinetic energy pressure from the external pressure
        integer:: md_integrator         !Integrator for VCMD: 1=Verlet, 2=Velocity-Verlet, 3=Beeman
        logical:: auto_dtion_md         !If true, the timestep during MD will be adjusted during run
        logical:: energy_conservation   !Only used in fixed cell MD
        integer:: nit_per_min           !Target number of md steps per md minimum crossing
        real(8):: dtion_md              !Initial timestep for MD
        logical:: auto_mdmin            !If true, the mdmin parameter will be adjusted during run
        integer:: mdmin                 !Number of enthalpy minima crossed unit stop MD
        integer:: mdmin_min,mdmin_max   !min,max number of enthalpy minima crossed unit stop MD, only if automatically determined
        real(8):: bmass                 !Cell mass during MD and FIRE
        logical::  geopt_ext            !At the moment only used for siesta: if true, the external geometry optimizer is used
        real(8):: alphax_lat , alphax_at !Stepsize for BFGS of the atomic and lattice coordinates
        logical:: auto_soft             !If true, the softening stepsize will be adjusted during run
        logical:: mol_soften            !Switch on molecular softening
        real(8):: alpha_lat, alpha_at   !Stepsize for softening the atomic and lattice coordinates
        logical:: auto_kpt              !Currently a dummy variable
        integer:: ka,kb,kc              !The number of kpoints in each dimension
        real(8):: dkpt1,dkpt2           !Precisions of the kpt mesh if generated automatically
        !-----------------------------------------------------------------------
        integer:: bc                    !1: periodic, 2:free, 3:surface/slab
        real(8):: target_pressure_gpa !Target pressures
        real(8):: target_pressure_habohr  !Target pressures
        logical:: findsym               !If true, findsym will be used to get symmetry informations on the fly
        logical:: finddos               !If true, the DOS at the Fermi level will be evaluated at the end of every geometry optimization
        logical:: usewf_md,usewf_geopt,usewf_soften !Defines when the wavefunctions should be reused in the next step
        logical:: core_rep              !If or if not to add a purely repulsive force on top of the atoms
        integer:: correctalg            !Method to perform cell corrections
!        integer:: ntime_md              !Maximum number of iterations during MD
        !-----------------------------------------------------------------------
        logical:: use_confine           !if true, confinement is enable, otherwise disabled
        integer:: nconfine       !number of different confinements
        integer,allocatable:: conf_dim(:)    !1,2 or 3 for each of the 3 dimensions latvec(:,1), latvec(:,2), latvec(:,3)
        integer,allocatable:: conf_av(:)     !0: no confinement (should never occur)
                                             !1: confinement with respect to a fixed value along latvec(:,i)
                                             !2: confinement with respect to the average
        integer,allocatable:: conf_exp(:)    !The polynomial order for each confinement
        real(8),allocatable:: conf_cut(:)    !The cutoff distance from each confinement equilibrium
        real(8),allocatable:: conf_eq(:)     !Equlibrium position of confinement along the confinement direction, will be filled to average or fixed value
        integer,allocatable:: conf_nat(:)    !How many atoms per confinement
        integer,allocatable:: conf_list(:,:) !List of atoms per confinement
        real(8),allocatable:: conf_prefac(:) !The polynomial predactor for each confinement
        character(1),allocatable:: conf_cartred(:)!Cartesian or reduced coordinates, only if conf_eq is provided. d,D,r,R for reduced, C,c,K,k for cartesian
        !-----------------------------------------------------------------------
        character(20):: fp_method_ch
        real(8):: fp_rcut
        real(8):: fp_dbin
        real(8):: fp_sigma
        integer:: fp_nl
        integer:: fp_14_m
        real(8):: fp_14_w1
        real(8):: fp_14_w2
        integer:: fp_at_nmax
        integer:: fp_17_natx_sphere
        integer:: fp_17_lseg
        character(len=2) :: fp_17_orbital
        real(8):: fp_17_nex_cutoff
        real(8):: fp_17_width_cutoff
        integer :: fp_18_principleev = 6
        integer :: fp_18_lseg!=1
        integer :: fp_18_molecules=4
        integer :: fp_18_expaparameter = 4
        integer :: fp_18_molecules_sphere = 50
        real*8  :: fp_18_width_cutoff = 1.d0
        real*8  :: fp_18_width_overlap = 1.d0
        integer :: fp_18_nex_cutoff = 3
        character:: fp_18_orbital
        !-----------------------------------------------------------------------
        real(8):: sock_ecutwf(2)
        integer:: sock_inet, sock_port        ! socket ID & address of the socket
        character(len=1024):: sock_host
        !-----------------------------------------------------------------------
        integer::qbfgs_bfgs_ndim!=1
        real(8)::qbfgs_trust_radius_max!=0.5d0
        real(8)::qbfgs_trust_radius_min!=1.d-5
        real(8)::qbfgs_trust_radius_ini!=0.02D0
        real(8)::qbfgs_w_1!=0.05D0
        real(8)::qbfgs_w_2!=0.5D0
        !-----------------------------------------------------------------------
        integer:: nat                   !Number of atoms
        integer:: ntypat_global         !Number of atom types
        real(8),allocatable:: znucl(:)
        character(2),allocatable:: char_type(:) 
        logical:: voids                 !If or if not to use void creating LJ particles in the cell
        real(8),allocatable:: amu(:)
        real(8),allocatable:: rcov(:)
        logical,allocatable:: fixat(:)
        logical:: fixlat(7)
        integer,allocatable:: typat_global(:)
        integer,allocatable:: llist(:)
        integer,allocatable:: lhead(:)
        integer,allocatable:: fragarr(:)
        integer,allocatable:: fragsize(:)
        integer:: vasp_kpt_mode         !If 1, the kpoint mesh is defined by mesh length, else a monkhorst pack mesh is generated (only vasp)
        integer:: abinit_kpt_mode       !If 1, the kpoint mesh is defined by kptrlen length, else a monkhorst pack mesh is generated (only abinit)
        integer:: siesta_kpt_mode       !If 1, the kpoint mesh is defined by cutoff length, else a monkhorst pack mesh is generated (only siesta)
        !-----------------------------------------------------------------------
        !Block fit_elecpot
        integer:: lcn
        integer:: iat_plot
        real(8), allocatable:: qt(:,:)
        real(8), allocatable:: at(:,:)
        real(8):: alphax_q_fit_elecpot
        real(8):: alphax_a_fit_elecpot
        real(8):: alphax_r_fit_elecpot
        real(8):: pot_rmse_tol
        logical:: cutoff_fit_elecpot
        !-----------------------------------------------------------------------
        type(dictionary), pointer :: dict_user=>null()
        type(dictionary), pointer :: dict=>null()
        type(dictionary), pointer :: subdict=>null()
        type(dictionary), pointer :: subsubdict=>null()
    end type typ_parini
end module mod_parini
!*****************************************************************************************
