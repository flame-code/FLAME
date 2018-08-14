!*****************************************************************************************
subroutine geopt(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_potential, only: fcalls, perfstatus, potential
    use mod_processors, only: iproc, nproc
    use mod_const, only: ang2bohr, bohr2ang, ev2ha
    use yaml_output
    implicit none
    type(typ_parini), intent(inout):: parini !poscar_getsystem must be called from parser
   ! type(typ_parini), intent(in):: parini
    !local variables
    integer:: iat, ierr, i
    real(8):: count_initial
    type(typ_paropt):: paropt, paropt_prec
    type(typ_atoms):: atoms
    type(typ_file_info):: file_info
    real(8):: rat_int(3,1000), fat_int(3,1000)
    logical:: acf_exists
    !character(5)::sat(natmax)
    inquire(file='posinp.acf',exist=acf_exists)
    if(acf_exists) then
        call acf_read(parini,'posinp.acf',1,atoms=atoms)
    else
        call read_poscar_for_single_point(parini,atoms)
    endif
    call set_ndof(atoms)
    potential=trim(parini%potential_potential)
    call init_geopt(parini,paropt,paropt_prec)
    paropt%trajectory=parini%paropt_geopt%trajectory
    paropt%filename='traj_mgo.acf'
    paropt%eps=1.d-8
    !write(*,*) 'POTENTIAL ',trim(potential)
    !stop
    !paropt%cellrelax=.true.
    if(paropt%cellrelax) then
        call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%rat,rat_int)
        !call vc_init_potential_forces(atoms)
        call init_potential_forces(parini,atoms)
        paropt%funits=250.d0
        call vc_minimize(parini,iproc,atoms,paropt)
        !do i=1,1
        !call cal_potential_forces_vc(iproc,atoms%nat,rat_int,atoms%cellvec,atoms%pressure,atoms%fat,atoms%celldv,atoms%stress,atoms%epot,atoms%enth)
        !call fxyz_cart2int(atoms%nat,atoms%fat,atoms%cellvec,fat_int)
        !write(*,'(a,f20.10)') 'EPOT ',atoms%epot
        !!write(*,'(a,3f20.10)') 'STRESS ',atoms%stress(1:3,1)
        !!write(*,'(a,3f20.10)') 'STRESS ',atoms%stress(1:3,2)
        !!write(*,'(a,3f20.10)') 'STRESS ',atoms%stress(1:3,3)
        !write(*,'(a,3f20.10)') 'STRESS ',atoms%celldv(1:3,1)
        !write(*,'(a,3f20.10)') 'STRESS ',atoms%celldv(1:3,2)
        !write(*,'(a,3f20.10)') 'STRESS ',atoms%celldv(1:3,3)
        !rat_int=rat_int+2.d-5*fat_int
        !atoms%cellvec=atoms%cellvec+1.d-5*atoms%celldv
        !enddo
        file_info%filename_positions='posout.acf'
        file_info%file_position='new'
        file_info%print_force=.true.
        call acf_write(file_info,atoms=atoms,strkey='posout')
        stop
    endif
    call init_potential_forces(parini,atoms)
    if(trim(potential)=='plato') then
        perfstatus='normal'
    endif
    !stop
    paropt%funits=5.d0
    if(parini%two_level_geopt) then
        paropt_prec%trajectory=parini%paropt_geopt_prec%trajectory
        paropt_prec%filename='traj_pgo.acf'
        !paropt_prec%nit=100
        paropt_prec%lprint=.true.
        !alphax not given in [geopt_prec] so we use the one in [geopt]
        if(paropt_prec%alphax<0.d0) then
            paropt_prec%alphax=paropt%alphax
        endif
        paropt_prec%condnum=paropt%condnum
        paropt_prec%precaution='high'
        paropt_prec%sdsaturation='yes'
        paropt_prec%funits=6.d0
        !paropt_prec%dt_start=paropt%dt_start
        call initminimize(paropt_prec)
    endif
    call initminimize(paropt)
    !if(iproc==0) write(*,*) 'RAT ',rat(1,2),rat(2,2),rat(3,2)
    !call cal_potential_forces(iproc,3*nat,rat,fat,epot)
    !rat(1,1)=rat(1,1)+1.d-2
    !call cal_potential_forces(iproc,3*nat,rat,fat,epot)
    count_initial=fcalls
    if(parini%two_level_geopt) then
        call minimize(parini,iproc,atoms,paropt_prec)
        call finalminimize(paropt_prec)
    endif
    !write(71,'(i5)') atoms%nat
    !do iat=1,atoms%nat
    !    write(71,'(3es23.13)') atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
    !call cal_potential_forces(iproc,3*atoms%nat,atoms%rat,atoms%fat,atoms%epot)
    call minimize(parini,iproc,atoms,paropt)
    file_info%filename_positions='posout.acf'
    file_info%file_position='new'
    file_info%print_force=.true.
    if(iproc==0) then
        call acf_write(file_info,atoms=atoms,strkey='posout')
    endif
    call finalminimize(paropt)
    if(iproc==0) then
        !write(*,'(a,l,i7)') 'converged ',paropt%converged,int(fcalls-count_initial)
        call yaml_mapping_open('task geopt')
        call yaml_map('converged',.true.)
        call yaml_map('total energy and force evaluations',int(fcalls-count_initial),fmt='(i5)')
        call yaml_mapping_close()
    endif
    call final_potential_forces(parini,atoms)
    call atom_deallocate_old(atoms)
end subroutine geopt
!*****************************************************************************************
subroutine init_geopt(parini,paropt,paropt_prec)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_paropt), intent(inout):: paropt, paropt_prec
    !local variables
    paropt=parini%paropt_geopt
    paropt_prec=parini%paropt_geopt_prec
end subroutine init_geopt
!*****************************************************************************************
