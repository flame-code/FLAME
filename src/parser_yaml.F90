!*****************************************************************************************
subroutine yaml_get_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use yaml_output
    implicit none
    type(typ_parini), intent(inout):: parini
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
            call set_atomc_types_info(parini)
        case("minhopp")
            call yaml_get_minhopp_parameters(parini)
        case("geopt")
            call yaml_get_geopt_parameters(parini)
        case("geopt_prec")
            call yaml_get_geopt_prec_parameters(parini)
        case("saddle_1s_opt")
            call yaml_get_saddle_1s_opt_parameters(parini)
        case("saddle_1s")
            call yaml_get_saddle_1s_parameters(parini)
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
        case("misc")
            call yaml_get_misc_parameters(parini)
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
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: main block in flame_in.yaml is empty.'
    parini%task=parini%subdict//"task"
    parini%types_main=parini%subdict//"types"
    parini%two_level_geopt=parini%subdict//"two_level_geopt"
    parini%iverbose=parini%subdict//"verbosity"
    parini%iseed=parini%subdict//"seed"
end subroutine yaml_get_main_parameters
!*****************************************************************************************
subroutine yaml_get_minhopp_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
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
end subroutine yaml_get_opt_parameters
!*****************************************************************************************
subroutine yaml_get_geopt_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: geopt block in flame_in.yaml is empty.'
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
subroutine yaml_get_saddle_1s_opt_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: saddle_1s_opt block in flame_in.yaml is empty.'
    call yaml_get_opt_parameters(parini,parini%paropt_saddle_1s_opt)
end subroutine yaml_get_saddle_1s_opt_parameters
!*****************************************************************************************
subroutine yaml_get_saddle_1s_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: saddle_1s block in flame_in.yaml is empty.'
    parini%str_moving_atoms_rand_saddle_1s=parini%subdict//"list_random_displace"
    parini%dimsep_saddle_1s=parini%subdict//"dimsep"
    parini%ampl_saddle_1s=parini%subdict//"ampl"
end subroutine yaml_get_saddle_1s_parameters
!*****************************************************************************************
subroutine yaml_get_potential_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: potential block in flame_in.yaml is empty.'
    parini%potential_potential=parini%subdict//"potential"
    parini%bias_potential=parini%subdict//"bias"
    parini%bias_field=parini%subdict//"bias_field"
    parini%cal_charge=parini%subdict//"cal_charge"
    parini%potential_potential_sec=parini%subdict//"potential_sec"
    parini%potential_ann_boundcheck=parini%subdict//"ann_boundcheck"
    parini%component_ff=parini%subdict//"component_ff"
    parini%inisock_inet=parini%subdict//"sockinet"
    parini%inisock_port=parini%subdict//"sockport"
    parini%inisock_host=parini%subdict//"sockhost"
    parini%drift_potential=parini%subdict//"drift"
    parini%add_repulsive=parini%subdict//"add_repulsive"
    if(has_key(parini%subdict,"ewald")) then
        parini%subsubdict => parini%subdict//"ewald"
        call yaml_get_ewald_parameters(parini)
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
    parini%nstep_ekf=parini%subdict//"nstep_ekf"
    parini%nstep_cep=parini%subdict//"nstep_cep"
    parini%nat_force=parini%subdict//"nat_force"
    parini%ampl_rand=parini%subdict//"ampl_rand"
    parini%symfunc_type_ann=parini%subdict//"symfunc_type"
    parini%syslinsolver_ann=parini%subdict//"syslinsolver"
    parini%psolver_ann=parini%subdict//"psolver"
    parini%rgnrmtol=parini%subdict//"rgnrmtol"
    parini%qgnrmtol=parini%subdict//"qgnrmtol"
    parini%etol_ann=parini%subdict//"etol"
    parini%dtol_ann=parini%subdict//"dtol"
    parini%normalization_ann=parini%subdict//"normalization"
    parini%bondbased_ann=parini%subdict//"bondbased"
    parini%prefit_ann=parini%subdict//"prefit"
end subroutine yaml_get_ann_parameters
!*****************************************************************************************
subroutine yaml_get_dynamics_parameters(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    if(dict_size(parini%subdict)<1) stop 'ERROR: dynamics block in flame_in.yaml is empty.'
    parini%nmd_dynamics=parini%subdict//"nmd"
    parini%dt_dynamics=parini%subdict//"dt"
    parini%temp_dynamics=parini%subdict//"temp"
    parini%init_temp_dynamics=parini%subdict//"init_temp"
    parini%md_method_dynamics=parini%subdict//"md_method"
    parini%print_force_dynamics=parini%subdict//"print_force"
    parini%restart_dynamics=parini%subdict//"restart"
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
    parini%hx_ewald=parini%subsubdict//"hx"
    parini%hy_ewald=parini%subsubdict//"hy"
    parini%hz_ewald=parini%subsubdict//"hz"
    parini%alpha_ewald=parini%subsubdict//"alpha"
    parini%ecut_ewald=parini%subsubdict//"ecut"
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
subroutine set_dict_parini_default(parini)
    use mod_parini, only: typ_parini
    use dictionaries
    use yaml_parse
    implicit none
    type(typ_parini), intent(inout):: parini
    !local variales
    type(dictionary), pointer :: dict
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
    type(dictionary), pointer :: dict
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
    type(dictionary), pointer :: subdict_parini
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
