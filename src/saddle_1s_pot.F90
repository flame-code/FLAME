!*****************************************************************************************
subroutine pot_initialize(parini,atoms,paropt,paropt_m)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_opt, only: typ_paropt
    use mod_processors, only: nproc, iproc
    use mod_potential, only: potential, perfstatus
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt, paropt_m
    !local variables
    call init_potential_forces(parini,atoms)
    if(trim(potential)=='lj') then
        paropt_m%funits=4.d0
        paropt%funits=4.d0
    elseif(trim(potential)=='bigdft') then
        stop 'ERROR: I do not remember why two_level_geopt so I commented the following line'
        !two_level_geopt=.true.
    elseif(trim(potential)=='vasp') then
        stop 'ERROR: I do not remember why two_level_geopt so I commented the following line'
        !two_level_geopt=.true.
    elseif(trim(potential)=='siesta') then
        stop 'ERROR: I do not remember why two_level_geopt so I commented the following line'
        !two_level_geopt=.true.
        perfstatus='normal'
    endif
!    use mod_atoms, only: typ_atoms
!#if defined(HAVE_VASP)
!    use potential_forces, only: init_potential_forces, fcalls, perfstatus, single_point_calculation
!#elif defined(HAVE_BigDFT)
!    use potential_forces, only: init_potential_forces, fcalls, perfstatus, final_potential_forces
!#elif defined(HAVE_LTB) || defined(HAVE_LJ)
!    use potential_forces, only: init_potential_forces, fcalls !, cal_potential_forces
!#elif defined(HAVE_QSC) || defined(HAVE_SIESTA)
!    use potential_forces, only: init_potential_forces, fcalls, perfstatus, final_potential_forces
!#else
!    stop 'ERROR: no potential is specified.'
!#endif
!    implicit none
!    type(typ_atoms), intent(in):: atoms
!    !local variables
!#if defined(HAVE_VASP)
!    perfstatus='normal'
!    call init_potential_forces(atoms%cellvec,atoms%nat,atoms%sat,atoms%atom_motion)
!#elif defined(HAVE_BigDFT)
!    call init_potential_forces(atoms%cellvec,atoms%nat,atoms%sat)
!#elif defined(HAVE_LTB) || defined(HAVE_QSC)
!    call init_potential_forces(atoms%cellvec,atoms%nat,atoms%sat)
!#elif defined(HAVE_LJ)
!    call init_potential_forces
!    !paropt%funits=4.d0
!    !paropt_target%funits=4.d0
!#elif defined(HAVE_SIESTA)
!    perfstatus='normal'
!    call init_potential_forces(atoms%nat,rat,atoms%sat,atoms%cellvec,iproc,nproc)
!#else
!    stop 'ERROR: no potential is specified.'
!#endif
end subroutine pot_initialize
!*****************************************************************************************
