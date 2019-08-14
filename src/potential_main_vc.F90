!*****************************************************************************************
subroutine vc_init_potential_forces(atoms)
    use mod_atoms, only: typ_atoms
    use mod_potential, only: fcalls, potential
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variables
    fcalls=0
    if(trim(potential)=='vcblj') then
        call init_lennardjones_vc(atoms%nat,atoms%sat)
    else
        stop 'ERROR: potential is unknown'
    endif
end subroutine vc_init_potential_forces
!*****************************************************************************************
!(nat,latvec,xred,fxyz,celldv,stress,pressure,etot,enth,count1)
subroutine cal_potential_forces_vc(iproc,nat,rat,cellvec,pressure,fat,celldv,stress,epot,enth)
    use mod_potential, only: potential, fcalls
    implicit none
    integer, intent(in):: iproc, nat
    real(8), intent(in):: rat(3,nat), cellvec(3,3), pressure
    real(8), intent(inout):: fat(3,nat), celldv(3,3), stress(3,3), epot, enth
    !local variables
    !integer::
    if(trim(potential)=='vcblj') then
        call lennardjones_vc(iproc,nat,rat,cellvec,pressure,fat,celldv,stress,epot,enth)
    else
        stop 'ERROR: potential is unknown'
    endif
    fcalls=fcalls+1
end subroutine cal_potential_forces_vc
!*****************************************************************************************
