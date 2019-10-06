!*****************************************************************************************
subroutine init_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variabbles
    !integer:: nsp
    !real(8):: cell(3), hgxy_rough, hgz_rough, rgcut_rough, rcut
    !call readinputparam(nrowx,nrowy,nlayer,hgxy_rough,hgz_rough,g,rcut,rgcut_rough)
    call forcefield_init(parini,atoms)
end subroutine init_potential_forces_ff
!*****************************************************************************************
subroutine cal_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    !use mod_potential, only: poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(inout):: atoms
    !local variabbles
    integer:: iat
    !real(8):: cell(3)
    call calculate_forces_energy_ff(parini,atoms)
end subroutine cal_potential_forces_ff
!*****************************************************************************************
subroutine final_potential_forces_ff(parini,atoms)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms
    !use mod_potential, only: poisson
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    !call fini_hartree(poisson)
    call forcefield_final(parini,atoms)
end subroutine final_potential_forces_ff
!*****************************************************************************************
