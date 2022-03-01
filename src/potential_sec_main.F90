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
    select case(trim(potential_sec))
        case('netsock')
            call init_netsock(parini)
        case('lj')
            call init_lennardjones
        case('ltb')
            call init_lenosky_tb(atoms)
        case('mpmd')
            call init_potential_forces_mpmd(atoms)
        case('qsc')
            call init_potential_forces_qsc(atoms)
        case('ann')
            call init_potential_ann(parini,atoms)
        case('bigdft')
            call init_potential_forces_bigdft(atoms)
        case('vasp')
            call init_potential_forces_vasp(atoms)
        case('forcefield')
            call init_potential_forces_ff(parini,atoms)
        case('vcblj')
            call init_lennardjones_vc(atoms%nat,atoms%sat)
        case('dftb')
            call init_potential_forces_dftb(atoms)
        case('siesta')
#if defined(HAVE_SIESTA)
            call init_potential_forces_siesta(atoms)
#else
            stop 'ERROR: FLAME is not linked with siesta during compilation.'
#endif
        case default
            stop 'ERROR: potential_sec is unknown'
    end select
end subroutine init_potential_forces_sec
!*****************************************************************************************
subroutine cal_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, get_rat, update_ratp, set_rat, update_rat
    use mod_processors, only: iproc
    use mod_potential, only: potential_sec, fcalls_sec
    use dynamic_memory
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: iat
    real(8), allocatable:: ratred(:,:), rat_backup(:,:)
    allocate(rat_backup(3,atoms%nat))
    call get_rat(atoms,rat_backup)
    call update_ratp(atoms)
    if(trim(atoms%boundcond)=='bulk') then
        allocate(ratred(3,atoms%nat))
        call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,ratred)
        call backtocell_alborz(atoms%nat,atoms%cellvec,ratred)
        call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,ratred,atoms%ratp)
        call update_rat(atoms,upall=.true.)
        deallocate(ratred)
    elseif(trim(atoms%boundcond)=='slab') then
        allocate(ratred(3,atoms%nat))
        call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,ratred)
        do iat=1,atoms%nat
            ratred(1,iat)=modulo(modulo(ratred(1,iat),1.d0),1.d0)
            ratred(2,iat)=modulo(modulo(ratred(2,iat),1.d0),1.d0)
        enddo
        call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,ratred,atoms%ratp)
        call update_rat(atoms,upall=.true.)
        deallocate(ratred)
    endif
    do iat=1,atoms%nat
        atoms%fat(1,iat)=0.d0
        atoms%fat(2,iat)=0.d0
        atoms%fat(3,iat)=0.d0
    enddo
    select case(trim(potential_sec))
        case('netsock')
            call cal_potential_forces_netsock(atoms)
        case('lj')
            call lennardjones(atoms)
        case('ltb')
            call lenosky_tb(parini,atoms)
        case('mpmd')
            call cal_potential_forces_mpmd(atoms)
        case('qsc')
            call cal_potential_forces_qsc(atoms)
        case('ann')
            call cal_potential_ann(parini,atoms)
        case('bigdft')
            call cal_potential_forces_bigdft(atoms)
        case('vasp')
            call cal_potential_forces_vasp(atoms)
        case('forcefield')
            call cal_potential_forces_ff(parini,atoms)
        case('vcblj')
            allocate(ratred(3,atoms%nat))
            call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,ratred)
            call lennardjones_vc(iproc,atoms%nat,ratred,atoms%cellvec,atoms%pressure,atoms%fat,atoms%celldv,atoms%stress,atoms%epot,atoms%enth)
            deallocate(ratred)
        case('dftb')
            call cal_potential_forces_dftb(atoms)
        case('siesta')
#if defined(HAVE_SIESTA)
            call cal_potential_forces_siesta(atoms)
#else
            stop 'ERROR: FLAME is not linked with siesta during compilation.'
#endif
        case default
            stop 'ERROR: potential_sec is unknown'
    end select
    fcalls_sec=fcalls_sec+1
    call set_rat(atoms,rat_backup,setall=.true.)
    deallocate(rat_backup)
end subroutine cal_potential_forces_sec
!*****************************************************************************************
subroutine fini_potential_forces_sec(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potential, only: potential_sec
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    !local variables
    select case(trim(potential_sec))
        case('netsock')
            call final_netsock
        case('lj')
            call init_lennardjones
        case('ltb')
        case('mpmd')
        case('qsc')
            call final_potential_forces_qsc
        case('ann')
            call final_potential_ann
        case('bigdft')
            call final_potential_forces_bigdft
        case('vasp')
            call final_potential_forces_vasp
        case('forcefield')
            call final_potential_forces_ff(parini,atoms)
        case('dftb')
            call final_potential_forces_dftb
        case('siesta')
#if defined(HAVE_SIESTA)
            call final_potential_forces_siesta
#else
            stop 'ERROR: FLAME is not linked with siesta during compilation.'
#endif
        case default
            stop 'ERROR: potential_sec is unknown'
    end select
end subroutine fini_potential_forces_sec
!*****************************************************************************************
