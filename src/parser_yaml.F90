!*****************************************************************************************
subroutine yaml_get_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use yaml_output
    implicit none
    type(typ_parini), intent(inout):: parini
    !nullify(parini%dict_user)
    !nullify(parini%dict)
    !nullify(parini%subdict)
    !nullify(parini%subsubdict)
    !-------------------------------------------------------
    call set_dict_parini_default(parini)
    !call yaml_comment('DEFAULT INPUT FILE ',hfill='~')
    !call yaml_dict_dump(parini%dict)
    !-------------------------------------------------------
    call set_dict_parini_user(parini)
    !call yaml_comment('USER INPUT FILE',hfill='~')
    !call yaml_dict_dump(parini%dict_user)
    !-------------------------------------------------------
    call dict_update(parini%dict,parini%dict_user)
    call yaml_comment('MERGED INPUT FILE',hfill='~')
    call yaml_dict_dump(parini%dict)
    call yaml_comment('',hfill='~')
    !-------------------------------------------------------
    parini%subdict=>dict_iter(parini%dict)
    do while(associated(parini%subdict))
        select case(trim(dict_key(parini%subdict)))
        case("main")
            call yaml_get_main_parameters(parini)
        case("minhopp")
            call yaml_get_minhopp_parameters(parini)
        case("geopt")
            call yaml_get_geopt_parameters(parini)
        case("geopt_prec")
            call yaml_get_geopt_prec_parameters(parini)
        case("saddle_opt")
            call yaml_get_saddle_opt_parameters(parini)
        case("saddle")
            call yaml_get_saddle_parameters(parini)
        case("potential")
            call yaml_get_potential_parameters(parini)
        case("ann")
            call yaml_get_ann_parameters(parini)
        case("dynamics")
            call yaml_get_dynamics_parameters(parini)
        case("bader")
            call yaml_get_bader_parameters(parini)
        case("genconf")
            call yaml_get_genconf_parameters(parini)
        case("conf_comp")
            call yaml_get_conf_comp_parameters(parini)
        case("testforces")
            call yaml_get_testforces_parameters(parini)
        case("single_point")
            call yaml_get_single_point_parameters(parini)
        case("fingerprint")
            call yaml_get_fingerprint_parameters(parini)
        case("misc")
            call yaml_get_misc_parameters(parini)
        case("fit_elecpot")
            call yaml_get_fit_elecpot_parameters(parini)
        end select
        parini%subdict=>dict_next(parini%subdict)
    end do
    nullify(parini%subdict)

    call check_nonoptional_parameters(parini)

    call dict_free(parini%dict)
    nullify(parini%dict)
    call dict_free(parini%dict_user)
    nullify(parini%dict_user)
end subroutine yaml_get_parameters
!*****************************************************************************************
subroutine yaml_get_main_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries, dict_set => set
    use defs_basis, only: HaBohr3_GPA
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    integer:: itype, verbosity_mh
    if(dict_size(parini%subdict)<1) stop 'ERROR: main block in flame_in.yaml is empty.'
    parini%task=parini%subdict//"task"
    parini%types_main=parini%subdict//"types"
    call set_atomc_types_info(parini)
    parini%two_level_geopt=parini%subdict//"two_level_geopt"
    parini%iverbose=parini%subdict//"verbosity"
    parini%rng_type=parini%subdict//"rng_type"
    parini%iseed=parini%subdict//"seed"
    parini%nrun_lammps=parini%subdict//"nrun_lammps"
    if(has_key(parini%subdict,"verbosity_mh")) then
        verbosity_mh=parini%subdict//"verbosity_mh"
        call dict_set(parini%subdict//"verbose",verbosity_mh)
    endif
    parini%verb=parini%subdict//"verbose"
    parini%params_new=parini%subdict//"params_new"
    parini%nat=parini%subdict//"nat"
    if(trim(parini%task)=='minhocao' .and. parini%nat<1) then
        write(*,*) 'ERROR: task=minhocao, did you set nat in input file?'
    endif
    parini%target_pressure_gpa=parini%subdict//"pressure"
    parini%target_pressure_habohr=parini%target_pressure_gpa/HaBohr3_GPA
    parini%ntypat_global=parini%ntypat
    if(trim(parini%task)=='minhocao') then
    if(.not.allocated(parini%znucl)) then ; allocate(parini%znucl(parini%ntypat_global),source=0.d0) ; endif
    if(.not.allocated(parini%amu)  ) then ; allocate(parini%amu(parini%ntypat_global),source=0.d0) ; endif
    if(.not.allocated(parini%rcov)  ) then ; allocate(parini%rcov(parini%ntypat_global),source=0.d0) ; endif
    if(.not.allocated(parini%char_type)) then; allocate(parini%char_type(parini%ntypat_global),source="  ") ; endif
    !if(.not.allocated(parini%char_type)) then
    !    allocate(parini%char_type(parini%ntypat_global))
    !    do itype=1,parini%ntypat_global
    !        parini%char_type(itype)="  "
    !    enddo
    !endif
    if(.not.allocated(parini%typat_global)) then; allocate(parini%typat_global(parini%nat),source=0) ; endif
    if(.not.allocated(parini%fixat)) then; allocate(parini%fixat(parini%nat),source=.false.) ; endif
    if(.not.allocated(parini%fragarr)) then; allocate(parini%fragarr(parini%nat),source=0) ; endif
    parini%str_typat_global=parini%subdict//"typat"
    read(parini%str_typat_global,*) parini%typat_global(:)
    !Get the correct atomic masses and atomic character
    do itype=1,parini%ntypat_global
        parini%char_type(itype)=trim(parini%stypat(itype))
        call symbol2znucl(parini%amu(itype),parini%rcov(itype),parini%char_type(itype),parini%znucl(itype))
    enddo
    if(has_key(parini%subdict,"znucl")) then
        parini%znucl=parini%subdict//"znucl"
    endif
    if(has_key(parini%subdict,"amass")) then
        parini%amu=parini%subdict//"amass"
    endif
    if(has_key(parini%subdict,"rcov")) then
        parini%rcov=parini%subdict//"rcov"
    endif
    parini%findsym=parini%subdict//"findsym"
    parini%finddos=parini%subdict//"finddos"
    endif !end of if on trim(parini%task)=='minhocao'
end subroutine yaml_get_main_parameters
!*****************************************************************************************
subroutine yaml_get_minhopp_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    !real(8):: alpha_lat_in, alpha_at_in
    if(dict_size(parini%subdict)<1) stop 'ERROR: minhopp block in flame_in.yaml is empty.'
    parini%nstep_minhopp=parini%subdict//"nstep"
    parini%nsoften_minhopp=parini%subdict//"nsoften"
    parini%mdmin_minhopp=parini%subdict//"mdmin"
    parini%minter_minhopp=parini%subdict//"minter"
    parini%nrandoff_minhopp=parini%subdict//"nrandoff"
    parini%npminx_minhopp=parini%subdict//"npminx"
    parini%etoler_minhopp=parini%subdict//"etoler"
    parini%eref_minhopp=parini%subdict//"eref"
    parini%ekinmax_minhopp=parini%subdict//"ekinmax"
    parini%alpha1_minhopp=parini%subdict//"alpha1"
    parini%alpha2_minhopp=parini%subdict//"alpha2"
    parini%beta1_minhopp=parini%subdict//"beta1"
    parini%beta2_minhopp=parini%subdict//"beta2"
    parini%beta3_minhopp=parini%subdict//"beta3"
    parini%trajectory_minhopp=parini%subdict//"trajectory"
    parini%print_force_minhopp=parini%subdict//"print_force"
    parini%auto_soft=parini%subdict//"auto_soft"
    !if(.not.parini%auto_soft) then
        parini%alpha_lat=parini%subdict//"alpha_lat"
    !endif
    !if(.not.parini%auto_soft) then
        parini%alpha_at=parini%subdict//"alpha_at"
    !endif
    parini%mol_soften=parini%subdict//"mol_soften"
end subroutine yaml_get_minhopp_parameters
!*****************************************************************************************
subroutine yaml_get_opt_parameters(parini,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use dictionaries
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt
    !local variales
    paropt%approach=parini%subdict//"method"
    paropt%fmaxtol=parini%subdict//"fmaxtol"
    paropt%alphax=parini%subdict//"alphax"
    paropt%condnum=parini%subdict//"condnum"
    paropt%precaution=parini%subdict//"precaution"
    paropt%lprint=parini%subdict//"lprint"
    paropt%dt_start=parini%subdict//"dt_start"
    paropt%nit=parini%subdict//"nit"
    paropt%dxmax=parini%subdict//"dxmax"
    paropt%anoise=parini%subdict//"anoise"
    paropt%nsatur=parini%subdict//"nsatur"
    paropt%cellrelax=parini%subdict//"cellrelax"
    paropt%funits=parini%subdict//"funits"
    paropt%print_force=parini%subdict//"print_force"
    paropt%trajectory=parini%subdict//"trajectory"
    paropt%nhist=parini%subdict//"nhist"
    paropt%dtmin=parini%subdict//"dt_min"
    paropt%dtmax=parini%subdict//"dt_max"
    paropt%strfact=parini%subdict//"strfact"
end subroutine yaml_get_opt_parameters
!*****************************************************************************************
subroutine yaml_get_geopt_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: geopt block in flame_in.yaml is empty.'
    parini%alphax_lat=parini%subdict//"hesslat"
    parini%alphax_at=parini%subdict//"hessat"
    parini%geopt_ext=parini%subdict//"geoext"
    parini%qbfgs_bfgs_ndim=parini%subdict//"qbfgsndim"
    parini%qbfgs_trust_radius_ini=parini%subdict//"qbfgstri"
    parini%qbfgs_trust_radius_min=parini%subdict//"qbfgstrmin"
    parini%qbfgs_trust_radius_max=parini%subdict//"qbfgstrmax"
    parini%qbfgs_w_1=parini%subdict//"qbfgsw1"
    parini%qbfgs_w_2=parini%subdict//"qbfgsw2"
    parini%paropt_geopt%maxrise=parini%subdict//"maxrise"
    parini%paropt_geopt%cutoffRatio=parini%subdict//"sqnmcutoff"
    parini%paropt_geopt%steepthresh=parini%subdict//"sqnmsteep"
    parini%paropt_geopt%trustr=parini%subdict//"sqnmtrustr"
    call yaml_get_opt_parameters(parini,parini%paropt_geopt)
end subroutine yaml_get_geopt_parameters
!*****************************************************************************************
subroutine yaml_get_geopt_prec_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: geopt_prec block in flame_in.yaml is empty.'
    call yaml_get_opt_parameters(parini,parini%paropt_geopt_prec)
end subroutine yaml_get_geopt_prec_parameters
!*****************************************************************************************
subroutine yaml_get_saddle_opt_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: saddle_opt block in flame_in.yaml is empty.'
    call yaml_get_opt_parameters(parini,parini%paropt_saddle_opt)
end subroutine yaml_get_saddle_opt_parameters
!*****************************************************************************************
subroutine yaml_get_saddle_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: saddle block in flame_in.yaml is empty.'
    parini%method_saddle=parini%subdict//"method"
    parini%str_moving_atoms_rand_saddle=parini%subdict//"list_random_displace"
    parini%dimsep_saddle=parini%subdict//"dimsep"
    parini%ampl_saddle=parini%subdict//"ampl"
    parini%np_splsad=parini%subdict//"np_splsad"
    parini%np_neb=parini%subdict//"np_neb"
    parini%ns2_splsad=parini%subdict//"ns2"
    parini%vdtol_splsad=parini%subdict//"vdtol"
    parini%dt_saddle=parini%subdict//"dt"
    parini%htol_splsad=parini%subdict//"htol"
    parini%alphax_saddle=parini%subdict//"alphax"
    parini%hybrid_splsad=parini%subdict//"hybrid"
    parini%typintpol=parini%subdict//"typintpol"
    parini%pickbestanchorpoints=parini%subdict//"pickbestanchorpoints"
    parini%runstat=parini%subdict//"runstat"
    parini%doneb=parini%subdict//"doneb"
    parini%docineb=parini%subdict//"docineb"
    parini%max_fcalls=parini%subdict//"fcalls_max"
    parini%fmaxtol_splsad=parini%subdict//"fmaxtol_splsad"
    parini%fmaxtol_neb=parini%subdict//"fmaxtol_neb"
    parini%opt_method_saddle=parini%subdict//"opt_method"
    parini%dbar=parini%subdict//"dbar"
    parini%alphax_bs=parini%subdict//"stepsize"
    parini%fnrmtol_coarse=parini%subdict//"fnrmtol_coarse"
    parini%nstep_bs=parini%subdict//"nstep"
    parini%bar_contract=parini%subdict//"bar_contract"
    parini%contr_dbar=parini%subdict//"dbar_contracted"
    parini%fnrmtol_contracted=parini%subdict//"fnrmtol_contracted"
    parini%nstep_contract=parini%subdict//"nstep_contract"
end subroutine yaml_get_saddle_parameters
!*****************************************************************************************
subroutine yaml_get_potential_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    integer:: kpt_abc(3)
    real(8):: dkpt_12(2)
    if(dict_size(parini%subdict)<1) stop 'ERROR: potential block in flame_in.yaml is empty.'
    parini%potential_potential=parini%subdict//"potential"
    parini%cal_charge=parini%subdict//"cal_charge"
    parini%potential_potential_sec=parini%subdict//"potential_sec"
    parini%potential_ann_boundcheck=parini%subdict//"ann_boundcheck"
    parini%component_ff=parini%subdict//"component_ff"
    parini%inisock_inet=parini%subdict//"sockinet"
    parini%inisock_port=parini%subdict//"sockport"
    parini%inisock_host=parini%subdict//"sockhost"
    parini%sock_inet=parini%subdict//"ipiinet"
    parini%sock_port=parini%subdict//"ipiport"
    parini%sock_host=parini%subdict//"ipihost"
    parini%sock_ecutwf=parini%subdict//"ipiecutwf"
    parini%sock_inet=parini%subdict//"sockinet"
    parini%sock_port=parini%subdict//"sockport"
    parini%sock_host=parini%subdict//"sockhost"
    parini%sock_ecutwf=parini%subdict//"sockecutwf"
    parini%bc=parini%subdict//"boundary"
    parini%drift_potential=parini%subdict//"drift"
    parini%add_repulsive=parini%subdict//"add_repulsive"
    parini%voids=parini%subdict//"voids"
    parini%core_rep=parini%subdict//"core_rep"
    parini%usewf_geopt=parini%subdict//"usewfgeo"
    parini%usewf_soften=parini%subdict//"usewfsoft"
    parini%usewf_md=parini%subdict//"usewfmd"
    parini%auto_kpt=parini%subdict//"auto_kpt"
    kpt_abc=parini%subdict//"kptmesh"
    if(has_key(parini%subdict,"kptmesh")) then
        parini%ka=kpt_abc(1)
        parini%kb=kpt_abc(2)
        parini%kc=kpt_abc(3)
    endif
    dkpt_12=parini%subdict//"kptden"
    if(has_key(parini%subdict,"kptden")) then
        parini%dkpt1=dkpt_12(1)
        parini%dkpt2=dkpt_12(2)
    endif
    if(has_key(parini%subdict,"ewald")) then
        parini%subsubdict => parini%subdict//"ewald"
        call yaml_get_ewald_parameters(parini)
        nullify(parini%subsubdict)
    endif
    if(has_key(parini%subdict,"confine")) then
        parini%subsubdict => parini%subdict//"confine"
        call yaml_get_confinement_parameters(parini)
        nullify(parini%subsubdict)
    endif
end subroutine yaml_get_potential_parameters
!*****************************************************************************************
subroutine yaml_get_ann_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: ann block in flame_in.yaml is empty.'
    parini%subtask_ann=parini%subdict//"subtask"
    parini%optimizer_ann=parini%subdict//"optimizer"
    parini%approach_ann=parini%subdict//"approach"
    parini%symfunc=parini%subdict//"symfunc"
    parini%nstep_opt_ann=parini%subdict//"nstep_opt"
    parini%nstep_cep=parini%subdict//"nstep_cep"
    parini%alphax_q=parini%subdict//"alphax_q"
    parini%alphax_r=parini%subdict//"alphax_r"
    parini%nconf_rmse=parini%subdict//"nconf_rmse"
    parini%ampl_rand=parini%subdict//"ampl_rand"
    parini%symfunc_type_ann=parini%subdict//"symfunc_type"
    parini%syslinsolver_ann=parini%subdict//"syslinsolver"
    parini%rgnrmtol=parini%subdict//"rgnrmtol"
    parini%qgnrmtol=parini%subdict//"qgnrmtol"
    parini%etol_ann=parini%subdict//"etol"
    parini%dtol_ann=parini%subdict//"dtol"
    parini%normalization_ann=parini%subdict//"normalization"
    parini%bondbased_ann=parini%subdict//"bondbased"
    parini%prefit_ann=parini%subdict//"prefit"
    parini%prefit_centt_ann=parini%subdict//"prefit_centt"
    parini%read_forces_ann=parini%subdict//"read_forces"
    parini%restart_param=parini%subdict//"restart_param"
    parini%restart_iter=parini%subdict//"restart_iter"
    parini%print_energy=parini%subdict//"print_energy"
    parini%fit_hoppint=parini%subdict//"fit_hoppint"
    parini%save_symfunc_force_ann=parini%subdict//"save_symfunc_force"
    parini%weight_hardness=parini%subdict//"weight_hardness"
    parini%save_symfunc_behnam=parini%subdict//"save_symfunc_behnam"
    parini%free_bc_direct=parini%subdict//"freeBC_direct"
    parini%ftol_ann=parini%subdict//"ftol"
end subroutine yaml_get_ann_parameters
!*****************************************************************************************
subroutine yaml_get_dynamics_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    integer:: mdmin_in
    if(dict_size(parini%subdict)<1) stop 'ERROR: dynamics block in flame_in.yaml is empty.'
    parini%nmd_dynamics=parini%subdict//"nmd"
    parini%dt_dynamics=parini%subdict//"dt"
    parini%temp_dynamics=parini%subdict//"temp"
    parini%init_temp_dynamics=parini%subdict//"init_temp"
    parini%highest_frequency=parini%subdict//"highest_freq"
    parini%ntherm=parini%subdict//"ntherm"
    parini%md_method_dynamics=parini%subdict//"md_method"
    parini%print_force_dynamics=parini%subdict//"print_force"
    parini%restart_dynamics=parini%subdict//"restart"
    parini%nfreq_dynamics=parini%subdict//"nfreq"
    parini%md_algo=parini%subdict//"algo"
    parini%md_integrator=parini%subdict//"integrator"
    parini%md_presscomp=parini%subdict//"presscomp"
    parini%bmass=parini%subdict//"cellmass"
    parini%auto_mdmin=parini%subdict//"auto_mdmin"
    parini%mdmin_min=parini%subdict//"mdmin_min"
    parini%mdmin_max=parini%subdict//"mdmin_max"
    mdmin_in=parini%subdict//"mdmin_init"
    if(.not.parini%auto_mdmin) then
        parini%mdmin=mdmin_in
    else
        parini%mdmin=max(mdmin_in,parini%mdmin_min)
    endif
    parini%auto_dtion_md=parini%subdict//"auto_mddt"
    !if(.not.parini%auto_dtion_md) then
        parini%dtion_md=parini%subdict//"dt_init"
    !endif
    parini%nit_per_min=parini%subdict//"nit_per_min"
    parini%energy_conservation=parini%subdict//"encon"
end subroutine yaml_get_dynamics_parameters
!*****************************************************************************************
subroutine yaml_get_bader_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: bader block in flame_in.yaml is empty.'
    parini%approach_bader=parini%subdict//"method"
    parini%filename_bader=parini%subdict//"filename"
    parini%vacuum_bader=parini%subdict//"vacuum"
end subroutine yaml_get_bader_parameters
!*****************************************************************************************
subroutine yaml_get_genconf_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: genconf block in flame_in.yaml is empty.'
    parini%subtask_genconf=parini%subdict//"subtask"
    parini%cal_pot_genconf=parini%subdict//"cal_pot"
    parini%nat_add_genconf=parini%subdict//"nat_add"
    parini%sat_genconf=parini%subdict//"sat"
    parini%amargin_genconf=parini%subdict//"amargin"
    parini%dmin_genconf=parini%subdict//"dmin"
    parini%dmax_genconf=parini%subdict//"dmax"
    parini%npoint_genconf=parini%subdict//"npoint"
    parini%fbmin_genconf=parini%subdict//"fbmin"
    parini%fbmax_genconf=parini%subdict//"fbmax"
    parini%variable_cell_genconf=parini%subdict//"variable_cell"
    parini%nonorthogonal_genconf=parini%subdict//"nonorthogonal"
end subroutine yaml_get_genconf_parameters
!*****************************************************************************************
subroutine yaml_get_conf_comp_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: conf_comp block in flame_in.yaml is empty.'
    parini%tol_conf_comp=parini%subdict//"tol"
end subroutine yaml_get_conf_comp_parameters
!*****************************************************************************************
subroutine yaml_get_testforces_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: testforces block in flame_in.yaml is empty.'
    parini%testforces_method=parini%subdict//"method"
end subroutine yaml_get_testforces_parameters
!*****************************************************************************************
subroutine yaml_get_single_point_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: single_point block in flame_in.yaml is empty.'
    parini%print_force_single_point=parini%subdict//"print_force"
    parini%frmt_single_point=parini%subdict//"format"
    parini%usesocket=parini%subdict//"usesocket"
    parini%inisock_inet=parini%subdict//"sockinet"
    parini%inisock_port=parini%subdict//"sockport"
    parini%inisock_host=parini%subdict//"sockhost"
end subroutine yaml_get_single_point_parameters
!*****************************************************************************************
subroutine yaml_get_ewald_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subsubdict)<1) stop 'ERROR: ewald block in flame_in.yaml is empty.'
    parini%alpha_ewald=parini%subsubdict//"alpha"
    parini%ecut_ewald=parini%subsubdict//"ecut"
    parini%ecutz_ewald=parini%subsubdict//"ecutz"
    parini%ecut_auto=parini%subsubdict//"ecut_auto"
    parini%rcut_ewald=parini%subsubdict//"rcut"
    parini%rgcut_ewald=parini%subsubdict//"rgcut"
    parini%nsp_ewald=parini%subsubdict//"nsp"
    parini%vu_ewald=parini%subsubdict//"plane_voltageu"
    parini%vl_ewald=parini%subsubdict//"plane_voltagel"
    parini%vu_ac_ewald=parini%subsubdict//"plane_voltageu_ac"
    parini%frequency_ewald=parini%subsubdict//"frequency"
    parini%gnrmtol_eem=parini%subsubdict//"gnrmtol"
    parini%ewald=parini%subsubdict//"ewald"
    parini%tolerance_ewald=parini%subsubdict//"ewald_tol"
    parini%efield=parini%subsubdict//"external_field"
    parini%bias_type=parini%subsubdict//"bias_type"
    parini%psolver=parini%subsubdict//"psolver"
    parini%cell_ortho=parini%subsubdict//"cell_ortho"
    parini%dielec_const=parini%subsubdict//"dielec_const"
    parini%dielec_const1=parini%subsubdict//"dielec_const1"
    parini%dielec_const2=parini%subsubdict//"dielec_const2"
    parini%cal_polar=parini%subsubdict//"cal_polar"
end subroutine yaml_get_ewald_parameters
!*****************************************************************************************
subroutine yaml_get_misc_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: misc block in flame_in.yaml is empty.'
    parini%subtask_misc=parini%subdict//"subtask"
    parini%gaussian_width=parini%subdict//"gaussian_width"
    parini%boundcond_misc=parini%subdict//"boundcond"
    parini%posinp_misc=parini%subdict//"posinp"
end subroutine yaml_get_misc_parameters
!*****************************************************************************************
subroutine yaml_get_confinement_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    integer:: iat, iconfine, nat_max, llen
    character(20):: strkey
    character(120):: strmess
    !return !TO_BE_CORRECTED
    if(dict_size(parini%subsubdict)<1) stop 'ERROR: confinement block in flame_in.yaml is empty.'
    parini%nconfine=parini%subsubdict//"nconfine"
    parini%use_confine=parini%subsubdict//"confinement"
    if(.not.allocated(parini%conf_dim) .and. parini%use_confine)    then; allocate(parini%conf_dim    (parini%nconfine),source=0)            ; endif
    if(.not.allocated(parini%conf_av) .and. parini%use_confine)     then; allocate(parini%conf_av     (parini%nconfine),source=0)            ; endif
    if(.not.allocated(parini%conf_exp) .and. parini%use_confine)    then; allocate(parini%conf_exp    (parini%nconfine),source=0)            ; endif
    if(.not.allocated(parini%conf_prefac) .and. parini%use_confine) then; allocate(parini%conf_prefac (parini%nconfine),source=0.d0)         ; endif
    if(.not.allocated(parini%conf_cut) .and. parini%use_confine)    then; allocate(parini%conf_cut    (parini%nconfine),source=0.d0)         ; endif
    if(.not.allocated(parini%conf_eq) .and. parini%use_confine)     then; allocate(parini%conf_eq     (parini%nconfine),source=0.d0)         ; endif
    if(.not.allocated(parini%conf_nat) .and. parini%use_confine)    then; allocate(parini%conf_nat    (parini%nconfine),source=0)            ; endif
    if(.not.allocated(parini%conf_cartred) .and. parini%use_confine)then; allocate(parini%conf_cartred(parini%nconfine),source="C")          ; endif
    strmess='ERROR: Provide confinement parameters either as a scalar ' &
        //'or an array of length nconfine: erroneous key= '
    strkey='cartred'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_cartred=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    strkey='dim'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_dim=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    strkey='exp'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_exp=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    strkey='prefac'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_prefac=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    strkey='cut'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_cut=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    strkey='av'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_av=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    strkey='eq'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_eq=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    strkey='nat'
    llen=dict_len(parini%subsubdict//strkey)
    if(llen==0 .or. llen==parini%nconfine) then
        if(parini%use_confine) then
        parini%conf_nat=parini%subsubdict//strkey
        endif
    else
        write(*,'(a,1x,a)') trim(strmess),trim(strkey)
        stop
    endif
    if(parini%use_confine) then
    nat_max=maxval(parini%conf_nat(1:parini%nconfine))
    if(.not.allocated(parini%conf_list))   then; allocate(parini%conf_list   (nat_max,parini%nconfine),source=0) ; endif
    if(trim(dict_value(parini%subsubdict//"atoms"))=='all') then
        !write(*,*) 'AAA ',dict_value(parini%subsubdict//"atoms")
        do iconfine=1,parini%nconfine
            do iat=1,parini%conf_nat(iconfine)
                parini%conf_list(iat,iconfine)=iat
            enddo
        enddo
    else
        do iconfine=1,parini%nconfine
            do iat=1,parini%conf_nat(iconfine)
                parini%conf_list(iat,iconfine)=parini%subsubdict//"atoms"//(iconfine-1)//(iat-1)
                !write(*,*) parini%conf_list(iat,iconfine)
            enddo
        enddo
    endif
    endif
end subroutine yaml_get_confinement_parameters
!*****************************************************************************************
subroutine yaml_get_fingerprint_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: fingerprint block in flame_in.yaml is empty.'
    parini%fp_method_ch=parini%subdict//"method"
    parini%fp_rcut=parini%subdict//"rcut"
    parini%fp_dbin=parini%subdict//"dbin"
    parini%fp_sigma=parini%subdict//"sigma"
    parini%fp_nl=parini%subdict//"nl"
    parini%fp_14_m=parini%subdict//"power"
    parini%fp_14_w1=parini%subdict//"gaussfac1"
    parini%fp_14_w2=parini%subdict//"gaussfac2"
    parini%fp_at_nmax=parini%subdict//"atnmax"
    parini%fp_17_natx_sphere=parini%subdict//"natx"
    parini%fp_17_orbital=parini%subdict//"orbital"
    parini%fp_18_orbital=parini%fp_17_orbital
    parini%fp_17_nex_cutoff=parini%subdict//"nexcut"
    parini%fp_18_nex_cutoff=int(parini%fp_17_nex_cutoff)
    parini%fp_18_principleev=parini%subdict//"principleev"
    parini%fp_18_molecules=parini%subdict//"molecules"
    parini%fp_18_expaparameter=parini%subdict//"expa"
    parini%fp_18_molecules_sphere=parini%subdict//"molsphere"
    parini%fp_18_width_cutoff=parini%subdict//"widthcut"
    parini%fp_18_width_overlap=parini%subdict//"widthover"
end subroutine yaml_get_fingerprint_parameters
!*****************************************************************************************
subroutine yaml_get_fit_elecpot_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    integer:: igto, ityp
    real(8):: ttq, tta
    if(dict_size(parini%subdict)<1) stop 'ERROR: fit_elecpot block in flame_in.yaml is empty.'
    parini%lcn=parini%subdict//"ngto"
    if(parini%lcn<1 .and. trim(parini%subtask_misc)=='fit_elecpot') then
        stop 'ERROR: parini%lcn<1 in flame_in.yaml'
    endif
    parini%iat_plot=parini%subdict//"iat_plot"
    parini%qt=f_malloc([1.to.parini%lcn,1.to.parini%ntypat],id='parini%q_per_type')
    parini%at=f_malloc([1.to.parini%lcn,1.to.parini%ntypat],id='parini%gwe_per_type')
    if(trim(parini%subtask_misc)=='fit_elecpot') then
        !write(*,*) dict_len(parini%subdict//"q_per_type")
        if(dict_len(parini%subdict//"q_per_type")/=parini%ntypat) then
            write(*,*) 'ERROR: incorrect length of list q_per_type in flame_in.yaml'
            stop
        endif
        if(dict_len(parini%subdict//"gwe_per_type")/=parini%ntypat) then
            write(*,*) 'ERROR: incorrect length of list gwe_per_type in flame_in.yaml'
            stop
        endif
        do ityp=0,parini%ntypat-1
            !write(*,*) dict_len(parini%subdict//"q_per_type"//ityp)
            if(dict_len(parini%subdict//"q_per_type"//ityp)/=parini%lcn) then
                write(*,*) 'ERROR: incorrect length of sublist q_per_type in flame_in.yaml'
                stop
            endif
            if(dict_len(parini%subdict//"gwe_per_type"//ityp)/=parini%lcn) then
                write(*,*) 'ERROR: incorrect length of sublist gwe_per_type in flame_in.yaml'
                stop
            endif
            do igto=0,parini%lcn-1
                parini%qt(igto+1,ityp+1)=parini%subdict//"q_per_type"//ityp//igto
                parini%at(igto+1,ityp+1)=parini%subdict//"gwe_per_type"//ityp//igto
            enddo
        enddo
    else
        !call yaml_dict_dump(parini%subdict)
        ttq=parini%subdict//"q_per_type"
        tta=parini%subdict//"gwe_per_type"
        parini%qt=ttq
        parini%at=tta
    endif
    parini%alphax_q_fit_elecpot=parini%subdict//"alphax_q"
    parini%alphax_a_fit_elecpot=parini%subdict//"alphax_a"
    parini%alphax_r_fit_elecpot=parini%subdict//"alphax_r"
    parini%pot_rmse_tol=parini%subdict//"pot_rmse_tol"
    parini%cutoff_fit_elecpot=parini%subdict//"cutoff_fit_elecpot"
end subroutine yaml_get_fit_elecpot_parameters
!*****************************************************************************************
subroutine set_dict_parini_default(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use yaml_parse
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    type(dictionary), pointer :: dict=>null()
    external :: get_input_variables_definition
    call yaml_parse_database(dict,get_input_variables_definition)
    call dict_copy(parini%dict,dict//0)
    call dict_free(dict)
    nullify(dict)
end subroutine set_dict_parini_default
!*****************************************************************************************
subroutine set_dict_parini_user(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use yaml_parse
    use dynamic_memory
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    type(dictionary), pointer :: dict=>null()
    character, dimension(:), allocatable :: fbuf
    character(len=*), parameter:: fname="flame_in.yaml"
    integer(kind = 8) :: cbuf, cbuf_len
    call getFileContent(cbuf,cbuf_len,fname,len_trim(fname))
    fbuf=f_malloc0_str(1,int(cbuf_len),id='fbuf')
    call copyCBuffer(fbuf,cbuf,cbuf_len)
    call freeCBuffer(cbuf)
    !then parse the user's input file
    call yaml_parse_from_char_array(dict,fbuf)
    call f_free_str(1,fbuf)
    call dict_copy(parini%dict_user,dict//0)
    call dict_free(dict)
    nullify(dict)
end subroutine set_dict_parini_user
!*****************************************************************************************
subroutine check_nonoptional_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(in):: parini
    !local variales
    type(dictionary), pointer :: subdict_parini=>null()
    !-------------------------------------------------------
    if(.not. has_key(parini%dict_user,"main")) then
        stop 'ERROR: flame_in.yaml must contain main block.'
    endif
    subdict_parini => parini%dict_user//"main"
    if(.not. has_key(subdict_parini,"task")) then
        stop 'ERROR: flame_in.yaml must contain task key in the main block.'
    endif
    if(.not. has_key(subdict_parini,"types")) then
        stop 'ERROR: flame_in.yaml must contain types key in the main block.'
    endif
    !-------------------------------------------------------
    nullify(subdict_parini)
end subroutine check_nonoptional_parameters
!*****************************************************************************************
