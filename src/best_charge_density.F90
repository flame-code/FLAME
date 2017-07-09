!*****************************************************************************************
subroutine best_charge_density(parini)
    use mod_interface
    use mod_parini, only: typ_parini
    !use mod_electrostatics, only: typ_poisson_p3d, typ_ewald_p3d
    !use mod_atoms, only: typ_atoms
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    !type(typ_poisson_p3d):: poisson_p3d
    !type(typ_ewald_p3d):: ewald_p3d
    !type(typ_atoms):: atoms
    !call cube_read('rho.cube',atoms,poisson_p3d%typ_poisson)
end subroutine best_charge_density
!*****************************************************************************************
