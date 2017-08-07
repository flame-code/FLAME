!*****************************************************************************************
module mod_parini
    use dictionaries
    use mod_opt, only: typ_paropt
    implicit none
    type typ_parini
        logical:: exists_yaml_file
        !-----------------------------------------------------------------------
        !parameters of [main]
        character(50):: task='unknown'
        character(256):: cwd='unknown'
        logical:: two_level_geopt=.false.
        integer:: iverbose=0
        integer:: iseed=-2
        character(50):: types_main='unknown'
        integer:: ntypat=-1
        integer:: ltypat(20)=-1
        integer:: iatomnum(20)=-1
        character(5):: stypat(20)='unknown'
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
        !parameters of [saddle_1s_opt]
        logical:: avail_saddle_1s_opt=.false.
        type(typ_paropt):: paropt_saddle_1s_opt
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
        character(50):: bias_potential='no'
        character(50):: bias_field='no'
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
        character(50):: psolver_ann='unknown'
        character(50):: symfunc_type_ann='behler'
        character(50):: symfunc='only_calculate'
        integer:: nstep_ekf=100
        integer:: nstep_cep=200
        integer:: nat_force=0
        real(8):: ampl_rand=1.d0
        real(8):: rgnrmtol=-1.d0
        real(8):: qgnrmtol=-1.d0
        real(8):: etol_ann !the tolerance difference of energies of two configuration
        real(8):: dtol_ann !distance between two FP
        logical:: normalization_ann=.false.
        logical:: prefit_ann=.false.
        logical:: read_forces_ann
        !-----------------------------------------------------------------------
        !parameters of [saddle_1s]
        logical:: avail_saddle_1s=.false.
        character(256):: str_moving_atoms_rand_saddle_1s
        real(8):: dimsep_saddle_1s=-1.d0
        real(8):: ampl_saddle_1s=-1.d0
        !-----------------------------------------------------------------------
        !parameters of [dynamics]
        logical:: avail_dynamics=.false.
        real(8):: dt_dynamics=-1.d0
        real(8):: temp_dynamics
        real(8):: init_temp_dynamics =0.d0
        integer:: nmd_dynamics=0
        character(20):: md_method_dynamics='unknown'
        logical:: print_force_dynamics=.false.
        logical:: restart_dynamics=.false.
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
        real(8):: hx_ewald=-1.d0
        real(8):: hy_ewald=-1.d0
        real(8):: hz_ewald=-1.d0
        real(8):: alpha_ewald=-1.d0
        real(8):: ecut_ewald=250.d0
        real(8):: rcut_ewald=-1.d0
        real(8):: rgcut_ewald=-1.d0
        integer:: nsp_ewald=-1
        real(8):: vu_ewald=0.d0
        real(8):: vl_ewald=0.d0
        real(8):: vu_ac_ewald=0.d0
        real(8):: frequency_ewald=0.d0
        real(8):: gnrmtol_eem=1.d-7
        logical :: ewald=.false.
        real(8):: tolerance_ewald = 1.d-6
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
        !-----------------------------------------------------------------------
        integer::qbfgs_bfgs_ndim!=1
        real(8)::qbfgs_trust_radius_max!=0.5d0
        real(8)::qbfgs_trust_radius_min!=1.d-5
        real(8)::qbfgs_trust_radius_ini!=0.02D0
        real(8)::qbfgs_w_1!=0.05D0
        real(8)::qbfgs_w_2!=0.5D0
        !-----------------------------------------------------------------------
        type(dictionary), pointer :: dict_user
        type(dictionary), pointer :: dict
        type(dictionary), pointer :: subdict
        type(dictionary), pointer :: subsubdict
    end type typ_parini
end module mod_parini
!*****************************************************************************************
