!*****************************************************************************************
subroutine get_main_parameters(file_ini,parini)
    use mod_task, only: typ_file_ini
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[main]')
    if(file_ini%iline_header==0) then
        write(*,'(a)') 'ERROR: [main] block not available in input.ini,'
        write(*,'(a)') '       [main] block and task in that are mandatory inputs.'
        stop
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'task',char_var=parini%task)
        call get_one_param(file_ini,'two_level_geopt',log_var=parini%two_level_geopt)
        call get_one_param(file_ini,'verbosity',int_var=parini%iverbose)
        call get_one_param(file_ini,'seed',int_var=parini%iseed)
        call get_one_param(file_ini,'types',char_line_var=parini%types_main)
        !call get_one_param(file_ini,'types',char_line_var=parini%stypat_genconf)
    enddo
end subroutine get_main_parameters
!*****************************************************************************************
subroutine get_minhopp_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[minhopp]')
    if(file_ini%iline_header==0) then
        parini%avail_minhopp=.false.
        if(trim(parini%task)=='minhopp') then
            write(*,'(a)') 'WARNING: [minhopp] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_minhopp=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'nstep',int_var=parini%nstep_minhopp)
        call get_one_param(file_ini,'nsoften',int_var=parini%nsoften_minhopp)
        call get_one_param(file_ini,'mdmin',int_var=parini%mdmin_minhopp)
        call get_one_param(file_ini,'minter',int_var=parini%minter_minhopp)
        call get_one_param(file_ini,'nrandoff',int_var=parini%nrandoff_minhopp)
        call get_one_param(file_ini,'npminx',int_var=parini%npminx_minhopp)
        call get_one_param(file_ini,'etoler',real_var=parini%etoler_minhopp)
        call get_one_param(file_ini,'eref',real_var=parini%eref_minhopp)
        call get_one_param(file_ini,'ekinmax',real_var=parini%ekinmax_minhopp)
        call get_one_param(file_ini,'alpha1',real_var=parini%alpha1_minhopp)
        call get_one_param(file_ini,'alpha2',real_var=parini%alpha2_minhopp)
        call get_one_param(file_ini,'beta1',real_var=parini%beta1_minhopp)
        call get_one_param(file_ini,'beta2',real_var=parini%beta2_minhopp)
        call get_one_param(file_ini,'beta3',real_var=parini%beta3_minhopp)
        call get_one_param(file_ini,'trajectory',log_var=parini%trajectory_minhopp)
        call get_one_param(file_ini,'print_force',log_var=parini%print_force_minhopp)
    enddo
end subroutine get_minhopp_parameters
!*****************************************************************************************
subroutine get_opt_param(file_ini,paropt)
    use mod_task, only: typ_file_ini
    use mod_parser_ini, only: get_one_param, get_header_location
    use mod_opt, only: typ_paropt
    use mod_saddle, only: str_moving_atoms_rand, ampl, dimsep
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_paropt), intent(inout):: paropt
    !local variables
    call get_one_param(file_ini,'method',char_var=paropt%approach)
    call get_one_param(file_ini,'fmaxtol',real_var=paropt%fmaxtol)
    call get_one_param(file_ini,'alphax',real_var=paropt%alphax)
    call get_one_param(file_ini,'condnum',real_var=paropt%condnum)
    call get_one_param(file_ini,'precaution',char_var=paropt%precaution)
    call get_one_param(file_ini,'lprint',log_var=paropt%lprint)
    call get_one_param(file_ini,'dt_start',real_var=paropt%dt_start)
    call get_one_param(file_ini,'nit',int_var=paropt%nit)
    call get_one_param(file_ini,'dxmax',real_var=paropt%dxmax)
    call get_one_param(file_ini,'anoise',real_var=paropt%anoise)
    call get_one_param(file_ini,'nsatur',int_var=paropt%nsatur)
    call get_one_param(file_ini,'cellrelax',log_var=paropt%cellrelax)
    call get_one_param(file_ini,'funits',real_var=paropt%funits)
    call get_one_param(file_ini,'print_force',log_var=paropt%print_force)
    call get_one_param(file_ini,'trajectory',log_var=paropt%trajectory)
    call get_one_param(file_ini,'nhist',int_var=paropt%nhist)
end subroutine get_opt_param
!*****************************************************************************************
subroutine get_geopt_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[geopt]')
    if(file_ini%iline_header==0) then
        parini%avail_geopt=.false.
        if(trim(parini%task)=='minhopp' .or. trim(parini%task)=='geopt' .or. &
            trim(parini%task)=='saddle') then
            write(*,'(a)') 'WARNING: [geopt] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_geopt=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_opt_param(file_ini,parini%paropt_geopt)
    enddo
end subroutine get_geopt_parameters
!*****************************************************************************************
subroutine get_geopt_prec_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[geopt_prec]')
    if(file_ini%iline_header==0) then
        parini%avail_geopt_prec=.false.
        if(parini%two_level_geopt .and. (trim(parini%task)=='minhopp' .or. trim(parini%task)=='geopt' .or. &
            trim(parini%task)=='saddle')) then
            write(*,'(a)') 'WARNING: [geopt_prec] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_geopt_prec=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_opt_param(file_ini,parini%paropt_geopt_prec)
    enddo
end subroutine get_geopt_prec_parameters
!*****************************************************************************************
subroutine get_saddle_opt_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[saddle_opt]')
    if(file_ini%iline_header==0) then
        parini%avail_saddle_opt=.false.
        if(trim(parini%task)=='saddle') then
            write(*,'(a)') 'WARNING: [saddle_opt] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_saddle_opt=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_opt_param(file_ini,parini%paropt_saddle_opt)
    enddo
end subroutine get_saddle_opt_parameters
!*****************************************************************************************
subroutine get_saddle_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[saddle]')
    if(file_ini%iline_header==0) then
        parini%avail_saddle=.false.
        if(trim(parini%task)=='saddle') then
            write(*,'(a)') 'WARNING: [saddle] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_saddle=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'list_random_displace',char_line_var=parini%str_moving_atoms_rand_saddle)
        call get_one_param(file_ini,'dimsep',real_var=parini%dimsep_saddle)
        call get_one_param(file_ini,'ampl',real_var=parini%ampl_saddle)
    enddo
end subroutine get_saddle_parameters
!*****************************************************************************************
subroutine get_potential_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[potential]')
    if(file_ini%iline_header==0) then
        parini%avail_potential=.false.
        if(trim(parini%task)=='minhopp' .or. trim(parini%task)=='geopt' .or. &
            trim(parini%task)=='saddle' .or. trim(parini%task)=='dynamics' .or. &
            trim(parini%task)=='genconf') then
            write(*,'(a)') 'WARNING: [potential] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_potential=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'potential',char_var=parini%potential_potential)
        call get_one_param(file_ini,'bias_type',char_var=parini%bias_type)
        call get_one_param(file_ini,'cal_charge',log_var=parini%cal_charge)
        call get_one_param(file_ini,'potential_sec',char_var=parini%potential_potential_sec)
        call get_one_param(file_ini,'ann_boundcheck',char_var=parini%potential_ann_boundcheck)
        call get_one_param(file_ini,'component_ff',char_var=parini%component_ff)
        call get_one_param(file_ini,'sockinet',int_var=parini%inisock_inet)
        call get_one_param(file_ini,'sockport',int_var=parini%inisock_port)
        call get_one_param(file_ini,'sockhost',char_var=parini%inisock_host)
        call get_one_param(file_ini,'drift',log_var=parini%drift_potential)
        call get_one_param(file_ini,'add_repulsive',log_var=parini%add_repulsive)
    enddo
end subroutine get_potential_parameters
!*****************************************************************************************
subroutine get_ann_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[ann]')
    if(file_ini%iline_header==0) then
        parini%avail_ann=.false.
        if(trim(parini%task)=='ann') then
            write(*,'(a)') 'WARNING: [ann] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_ann=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'subtask',char_var=parini%subtask_ann)
        call get_one_param(file_ini,'optimizer',char_var=parini%optimizer_ann)
        call get_one_param(file_ini,'approach',char_var=parini%approach_ann)
        call get_one_param(file_ini,'symfunc',char_var=parini%symfunc)
        call get_one_param(file_ini,'nstep_opt',int_var=parini%nstep_opt_ann)
        call get_one_param(file_ini,'nstep_cep',int_var=parini%nstep_cep)
        call get_one_param(file_ini,'nconf_rmse',int_var=parini%nconf_rmse)
        call get_one_param(file_ini,'ampl_rand',real_var=parini%ampl_rand)
        call get_one_param(file_ini,'symfunc_type',char_var=parini%symfunc_type_ann)
        call get_one_param(file_ini,'syslinsolver',char_var=parini%syslinsolver_ann)
        call get_one_param(file_ini,'etol',real_var=parini%etol_ann)
        call get_one_param(file_ini,'dtol',real_var=parini%dtol_ann)
        call get_one_param(file_ini,'normalization',log_var=parini%normalization_ann)
        call get_one_param(file_ini,'bondbased',log_var=parini%bondbased_ann)
        call get_one_param(file_ini,'prefit',log_var=parini%prefit_ann)
        call get_one_param(file_ini,'restart_param',log_var=parini%restart_param)
        call get_one_param(file_ini,'restart_iter',int_var=parini%restart_iter)
    enddo
end subroutine get_ann_parameters
!*****************************************************************************************
subroutine get_dynamics_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    use mod_dynamics, only: dt, nmd,nfreq
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[dynamics]')
    if(file_ini%iline_header==0) then
        parini%avail_dynamics=.false.
        if(trim(parini%task)=='dynamics') then
            write(*,'(a)') 'WARNING: [dynamics] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_dynamics=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'nmd',int_var=parini%nmd_dynamics)
        call get_one_param(file_ini,'nfreq',int_var=parini%nfreq_dynamics)
        call get_one_param(file_ini,'dt',real_var=parini%dt_dynamics)
        call get_one_param(file_ini,'temp',real_var=parini%temp_dynamics)
        call get_one_param(file_ini,'init_temp',real_var=parini%init_temp_dynamics)
        call get_one_param(file_ini,'md_method',char_var=parini%md_method_dynamics)
        call get_one_param(file_ini,'print_force',log_var=parini%print_force_dynamics)
        call get_one_param(file_ini,'restart',log_var=parini%restart_dynamics)
        call get_one_param(file_ini,'fix_cm',log_var=parini%fix_cm_dynamics)
        call get_one_param(file_ini,'vflip',log_var=parini%vflip_dynamics)
        call get_one_param(file_ini,'wall_repulsion',log_var=parini%wall_repulsion_dynamics)
    enddo
end subroutine get_dynamics_parameters
!*****************************************************************************************
subroutine get_bader_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[bader]')
    if(file_ini%iline_header==0) then
        parini%avail_bader=.false.
        if(trim(parini%task)=='bader') then
            write(*,'(a)') 'WARNING: [bader] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_bader=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'method',char_var=parini%approach_bader)
        call get_one_param(file_ini,'filename',char_var=parini%filename_bader)
        call get_one_param(file_ini,'vacuum',char_var=parini%vacuum_bader)
    enddo
end subroutine get_bader_parameters
!*****************************************************************************************
subroutine get_genconf_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    use mod_genconf, only: typ_genconf
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[genconf]')
    if(file_ini%iline_header==0) then
        parini%avail_genconf=.false.
        if(trim(parini%task)=='genconf') then
            write(*,'(a)') 'WARNING: [genconf] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_genconf=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'subtask',char_var=parini%subtask_genconf)
        call get_one_param(file_ini,'cal_pot',char_var=parini%cal_pot_genconf)
        call get_one_param(file_ini,'nat_add',int_var=parini%nat_add_genconf)
        call get_one_param(file_ini,'sat',char_var=parini%sat_genconf)
        call get_one_param(file_ini,'amargin',real_var=parini%amargin_genconf)
        call get_one_param(file_ini,'dmin',real_var=parini%dmin_genconf)
        call get_one_param(file_ini,'dmax',real_var=parini%dmax_genconf)
        call get_one_param(file_ini,'npoint',int_var=parini%npoint_genconf)
        call get_one_param(file_ini,'fbmin',real_var=parini%fbmin_genconf)
        call get_one_param(file_ini,'fbmax',real_var=parini%fbmax_genconf)
        call get_one_param(file_ini,'variable_cell',log_var=parini%variable_cell_genconf)
        call get_one_param(file_ini,'nonorthogonal',log_var=parini%nonorthogonal_genconf)
    enddo
end subroutine get_genconf_parameters
!*****************************************************************************************
subroutine get_conf_comp_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[conf_comp]')
    if(file_ini%iline_header==0) then
        parini%avail_conf_comp=.false.
        if(trim(parini%task)=='conf_comp') then
            write(*,'(a)') 'WARNING: [conf_comp] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_conf_comp=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'tol',real_var=parini%tol_conf_comp)
    enddo
end subroutine get_conf_comp_parameters
!*****************************************************************************************
subroutine get_testforces_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[testforces]')
    if(file_ini%iline_header==0) then
        parini%avail_testforces=.false.
        if(trim(parini%task)=='testforces') then
            write(*,'(a)') 'WARNING: [testforces] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_testforces=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'method',char_var=parini%testforces_method)
    enddo
end subroutine get_testforces_parameters
!*****************************************************************************************
subroutine get_single_point_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[single_point]')
    if(file_ini%iline_header==0) then
        parini%avail_single_point=.false.
        if(trim(parini%task)=='single_point') then
            write(*,'(a)') 'WARNING: [single_point] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_single_point=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'print_force',log_var=parini%print_force_single_point)
        call get_one_param(file_ini,'format',char_var=parini%frmt_single_point)
        call get_one_param(file_ini,'usesocket',log_var=parini%usesocket)
        call get_one_param(file_ini,'sockinet',int_var=parini%inisock_inet)
        call get_one_param(file_ini,'sockport',int_var=parini%inisock_port)
        call get_one_param(file_ini,'sockhost',char_var=parini%inisock_host)
    enddo
end subroutine get_single_point_parameters
!*****************************************************************************************
subroutine get_ewald_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[ewald]')
    if(file_ini%iline_header==0) then
        parini%avail_single_point=.false.
        if(trim(parini%task)=='ewald') then
            write(*,'(a)') 'WARNING: [ewald] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_single_point=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'alpha',real_var=parini%alpha_ewald)
        call get_one_param(file_ini,'ecut',real_var=parini%ecut_ewald)
        call get_one_param(file_ini,'ecutz',real_var=parini%ecutz_ewald)
        call get_one_param(file_ini,'rcut',real_var=parini%rcut_ewald)
        call get_one_param(file_ini,'rgcut',real_var=parini%rgcut_ewald)
        call get_one_param(file_ini,'nsp',int_var=parini%nsp_ewald)
        call get_one_param(file_ini,'plane_voltageu',real_var=parini%vu_ewald)
        call get_one_param(file_ini,'plane_voltagel',real_var=parini%vl_ewald)
        call get_one_param(file_ini,'plane_voltageu_ac',real_var=parini%vu_ac_ewald)
        call get_one_param(file_ini,'frequency',real_var=parini%frequency_ewald)
        call get_one_param(file_ini,'gnrmtol',real_var=parini%gnrmtol_eem)
        call get_one_param(file_ini,'ewald',log_var=parini%ewald)
        call get_one_param(file_ini,'ewald_tol',real_var=parini%tolerance_ewald)
        call get_one_param(file_ini,'psolver',char_var=parini%psolver)
        call get_one_param(file_ini,'cell_ortho',log_var=parini%cell_ortho)
    enddo
end subroutine get_ewald_parameters
!*****************************************************************************************
subroutine get_misc_parameters(file_ini,parini)
    use mod_parini, only: typ_parini
    use mod_parser_ini, only: get_one_param, get_header_location, split_line
    use mod_task, only: typ_file_ini
    implicit none
    type(typ_file_ini), intent(inout):: file_ini
    type(typ_parini), intent(inout):: parini
    !local variables
    integer:: iline
    call get_header_location(file_ini,'[misc]')
    if(file_ini%iline_header==0) then
        parini%avail_misc=.false.
        if(trim(parini%task)=='misc') then
            write(*,'(a)') 'WARNING: [misc] block not available in input.ini, default values will be used.'
        endif
        return
    else
        parini%avail_misc=.true.
    endif
    do iline=file_ini%iline_header+1,file_ini%iline_next_header-1
        file_ini%iline=iline
        if(file_ini%stat_line_is_read(file_ini%iline)) cycle
        call split_line(file_ini)
        call get_one_param(file_ini,'subtask',char_var=parini%subtask_misc)
        call get_one_param(file_ini,'gaussian_width',real_var=parini%gaussian_width)
        call get_one_param(file_ini,'boundcond',char_var=parini%boundcond_misc)
        call get_one_param(file_ini,'posinp',char_var=parini%posinp_misc)
    enddo
end subroutine get_misc_parameters
!*****************************************************************************************
