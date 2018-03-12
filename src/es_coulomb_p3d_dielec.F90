subroutine dielec_potener_forces(parini,poisson,atoms,epot_dielec)
    use mod_interface
    use mod_electrostatics, only: typ_poisson
    use mod_atoms, only: typ_atoms
    use mod_potential, only: potential 
    use mod_parini, only: typ_parini
    use dynamic_memory
    implicit none
    type(typ_poisson), intent(inout):: poisson
    type(typ_atoms), intent(inout):: atoms
    type(typ_parini), intent(in):: parini
    real(8), intent(out):: epot_dielec
    !local variables
    !Do be added by Farhad.
end subroutine dielec_potener_forces
!*****************************************************************************************
