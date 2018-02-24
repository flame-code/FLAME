!*****************************************************************************************
subroutine forcefield_init(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potential, only: poisson
    use mod_potential, only: shortrange
    use mod_electrostatics, only: typ_poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    type(typ_poisson):: poisson_rough
    integer:: ind
    !integer:: iat
    call set_typat(atoms)
    call set_qat(atoms)
    call set_interaction(atoms,shortrange)
    !-------------------------------------------------------
    ind=index(parini%component_ff,'buck')
    if(ind>0) then
        call set_buckingham(atoms,shortrange%tosifumi)
        poisson%spline%do_tosifumi=.true.

        deallocate (shortrange%interaction,shortrange%qq)
        allocate (shortrange%interaction(atoms%ntypat,atoms%ntypat),shortrange%qq(shortrange%ntypinter))
        shortrange%interaction=shortrange%tosifumi%interaction(1:shortrange%ntypinter,1:shortrange%ntypinter)
        shortrange%qq=shortrange%tosifumi%eee(1:shortrange%ntypinter)
    endif
    !-------------------------------------------------------
    ind=index(parini%component_ff,'coulomb')
    if(ind>0) then
        if(trim(atoms%boundcond)=='free') then
        elseif(trim(atoms%boundcond)=='slab') then
        call init_hartree(parini,atoms,poisson)
        shortrange%alpha=poisson%alpha
        poisson%spline%do_coulomb=.true.
        else
            stop 'ERROR: coulomb interaction not implemented for this type of BC'
        endif
    endif
    !-------------------------------------------------------
    !if(ind>0 .and. trim(atoms%boundcond)=='slab') then
    !call shortrange_init(atoms,shortrange,poisson%linked_lists,poisson%spline)
    !endif
    !-------------------------------------------------------
    !do iat=1,atoms%nat
    !    write(*,'(2i4,f6.1)') iat,atoms%itypat(iat),atoms%qat(iat)
    !enddo
    !stop
    ind=index(parini%component_ff,'tosifumi_old')
    if(ind>0) then
        call set_tosifumi(atoms,shortrange%tosifumi)
    endif
    !-------------------------------------------------------
    ind=index(parini%component_ff,'tosifumi')
    if(ind>0) then
        call set_tosifumi(atoms,shortrange%tosifumi)
        poisson%spline%do_tosifumi=.true.
    endif
    !-------------------------------------------------------
    if(poisson%spline%do_coulomb) then
    call shortrange_init(atoms,shortrange,poisson%linked_lists,poisson%spline)
    endif
end subroutine forcefield_init
!*****************************************************************************************
subroutine calculate_forces_energy_ff(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    !use mod_shortrange, only: typ_tosifumi
    use mod_potential, only: shortrange
    use mod_potential, only: poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: ind
    real(8):: epotshort
    ind=index(parini%component_ff,'coulomb')
    if(ind>0) then
        if(trim(atoms%boundcond)=='free') then
            call coulomb_free_direct(atoms)
        elseif(trim(atoms%boundcond)=='slab') then
        call calculate_forces_energy(parini,poisson,atoms)
        call cal_shortenergy(parini,shortrange,atoms,poisson%linked_lists,poisson%spline,poisson%alpha,poisson%cell,epotshort)
        write(*,'(a50,e32.15)') 'epotshort',epotshort
        atoms%epot=atoms%epot+epotshort
        endif
    endif
    ind=index(parini%component_ff,'tosifumi_old')
    if(ind>0) then
        call calenergyforces(atoms,shortrange%tosifumi)
    endif
end subroutine calculate_forces_energy_ff
!*****************************************************************************************
subroutine forcefield_final(parini,atoms)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    use mod_potential, only: poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    !local variables
    integer:: ind
    ind=index(parini%component_ff,'coulomb')
    if(ind>0) then
        if(trim(atoms%boundcond)=='free') then
        elseif(trim(atoms%boundcond)=='slab') then
        call destruct_poisson(parini,atoms,poisson)
        endif
    endif
    if(ind>0 .and. trim(atoms%boundcond)=='slab') then
    call shortrange_final(poisson%linked_lists,poisson%spline)
    endif
    ind=index(parini%component_ff,'tosifumi_old')
    if(ind>0) then
    endif
end subroutine forcefield_final
!*****************************************************************************************
