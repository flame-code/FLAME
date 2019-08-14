!*****************************************************************************************
subroutine init_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potential, only: fcalls_sec, potential_sec
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    fcalls_sec=0.d0
    if(trim(potential_sec)=='lj') then
        call init_lennardjones
    elseif(trim(potential_sec)=='ltb') then
        call init_lenosky_tb(atoms)
    elseif(trim(potential_sec)=='mpmd') then
        call init_potential_forces_mpmd(atoms)
    elseif(trim(potential_sec)=='qsc') then
        call init_cal_potential_forces_qsc(atoms)
    elseif(trim(potential_sec)=='ann') then
        call init_potential_ann(parini,atoms)
    elseif(trim(potential_sec)=='bigdft') then
        call init_cal_potential_forces_bigdft(atoms)
    elseif(trim(potential_sec)=='vasp') then
        call init_cal_potential_forces_vasp(atoms)
    elseif(trim(potential_sec)=='siesta') then
#if defined(HAVE_SIESTA)
        call init_cal_potential_forces_siesta(atoms)
#else
        stop 'ERROR: Alborz is not linked with siesta during compilation.'
#endif
    else
        stop 'ERROR: potential_sec is unknown'
    endif
end subroutine init_potential_forces_sec
!*****************************************************************************************
subroutine cal_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potential, only: potential_sec, fcalls_sec
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    if(trim(potential_sec)=='lj') then
        call lennardjones(atoms)
    elseif(trim(potential_sec)=='ltb') then
        call lenosky_tb(parini,atoms)
    elseif(trim(potential_sec)=='mpmd') then
        call cal_potential_forces_mpmd(atoms)
    elseif(trim(potential_sec)=='qsc') then
        call cal_potential_forces_qsc(atoms)
    elseif(trim(potential_sec)=='ann') then
        call cal_potential_ann(parini,atoms)
    elseif(trim(potential_sec)=='bigdft') then
        call cal_potential_forces_bigdft(atoms)
    elseif(trim(potential_sec)=='vasp') then
        call cal_potential_forces_vasp(atoms)
    elseif(trim(potential_sec)=='siesta') then
#if defined(HAVE_SIESTA)
        call cal_potential_forces_siesta(atoms)
#else
        stop 'ERROR: Alborz is not linked with siesta during compilation.'
#endif
    else
        stop 'ERROR: potential_sec is unknown'
    endif
    fcalls_sec=fcalls_sec+1
end subroutine cal_potential_forces_sec
!*****************************************************************************************
subroutine final_potential_forces_sec(atoms)
    use mod_atoms, only: typ_atoms
    use mod_potential, only: potential_sec
    implicit none
    type(typ_atoms), intent(in):: atoms
    !local variables
    if(trim(potential_sec)=='lj') then
        call init_lennardjones
    elseif(trim(potential_sec)=='ltb') then
    elseif(trim(potential_sec)=='mpmd') then
    elseif(trim(potential_sec)=='qsc') then
        call final_potential_forces_qsc
    elseif(trim(potential_sec)=='ann') then
        call final_potential_ann
    elseif(trim(potential_sec)=='bigdft') then
        call final_potential_forces_bigdft
    elseif(trim(potential_sec)=='vasp') then
        call final_potential_forces_vasp
    elseif(trim(potential_sec)=='siesta') then
#if defined(HAVE_SIESTA)
        call final_potential_forces_siesta
#else
        stop 'ERROR: Alborz is not linked with siesta during compilation.'
#endif
    else
        stop 'ERROR: potential_sec is unknown'
    endif
end subroutine final_potential_forces_sec
!*****************************************************************************************
