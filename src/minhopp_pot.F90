!*****************************************************************************************
subroutine setpot_init(parini,atoms_curr,paropt,paropt_prec)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_opt, only: typ_paropt
    use mod_processors, only: nproc, iproc
    use mod_potential, only: potcode, perfstatus
    use mod_potential, only: init_potential_forces
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms_curr
    type(typ_paropt), intent(inout):: paropt, paropt_prec
    !local variables
    call init_potential_forces(parini,atoms_curr)
    if(trim(potcode)=='lj') then
        paropt_prec%funits=200.d0
        paropt%funits=200.d0
    elseif(trim(potcode)=='bigdft') then
        stop 'ERROR: I do not remember why two_level_geopt so I commented the following line'
        !parini%two_level_geopt=.true.
    elseif(trim(potcode)=='vasp') then
        stop 'ERROR: I do not remember why two_level_geopt so I commented the following line'
        !parini%two_level_geopt=.true.
    elseif(trim(potcode)=='siesta') then
        stop 'ERROR: I do not remember why two_level_geopt so I commented the following line'
        !parini%two_level_geopt=.true.
        perfstatus='normal'
    endif
end subroutine setpot_init
!*****************************************************************************************
subroutine setpot_final(parini,atoms_curr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potential, only: fini_potential_forces
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms_curr
    !local variables
    call fini_potential_forces(parini,atoms_curr)
end subroutine setpot_final
!*****************************************************************************************
subroutine setpot_mdescape
    use mod_potential, only: potcode, perfstatus, single_point_calculation
    implicit none
    !local variables
    if(trim(potcode)=='bigdft') then
        perfstatus='fast'
    elseif(trim(potcode)=='vasp') then
        perfstatus='fast'
        single_point_calculation=.true.
    elseif(trim(potcode)=='siesta') then
        !call final_potential_forces
        perfstatus='fast'
        !call init_potential_forces(nat,rathopp,sat,cellvec,iproc,nproc)
    endif
end subroutine setpot_mdescape
!*****************************************************************************************
subroutine setpot_soften
    use mod_potential, only: potcode, perfstatus, single_point_calculation
    implicit none
    !local variables
    if(trim(potcode)=='bigdft') then
        perfstatus='accurate'
    elseif(trim(potcode)=='vasp') then
        perfstatus='soften'
        single_point_calculation=.true.
    endif
end subroutine setpot_soften
!*****************************************************************************************
subroutine setpot_geopt_prec
    use mod_potential, only: potcode, perfstatus, single_point_calculation
    implicit none
    !local variables
    if(trim(potcode)=='bigdft') then
        perfstatus='normal'
    elseif(trim(potcode)=='vasp') then
        perfstatus='normal'
        single_point_calculation=.false.
    elseif(trim(potcode)=='siesta') then
        !call final_potential_forces
        perfstatus='normal'
        !call init_potential_forces(atoms%nat,rat,sat,cellvec,iproc,nproc)
    endif
end subroutine setpot_geopt_prec
!*****************************************************************************************
subroutine setpot_geopt
    use mod_potential, only: potcode, perfstatus, single_point_calculation
    implicit none
    !local variables
    if(trim(potcode)=='bigdft') then
        perfstatus='accurate'
    elseif(trim(potcode)=='vasp') then
        single_point_calculation=.false.
        perfstatus='accurate'
    elseif(trim(potcode)=='siesta') then
        !call final_potential_forces
        perfstatus='accurate'
        !call init_potential_forces(atoms%nat,atoms%rat,sat,cellvec,iproc,nproc)
    endif
end subroutine setpot_geopt
!*****************************************************************************************
