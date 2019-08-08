!*****************************************************************************************
subroutine init_lennardjones
    implicit none
end subroutine init_lennardjones
!*****************************************************************************************
!energy and forces for Lennard Jones potential
subroutine lennardjones(atoms)
    use mod_atoms, only: typ_atoms, update_ratp
    implicit none
    type(typ_atoms), intent(inout):: atoms
    !local variabbles
    integer:: iat, jat
    real(8):: xiat, yiat, ziat, dx, dy, dz, rsq, rinvsq, rinv4, rinv6, rinv12
    real(8):: tt, t1, t2, t3
    atoms%epot=0.d0
    atoms%fat(1:3,1:atoms%nat)=0.d0
    call update_ratp(atoms)
    do iat=1,atoms%nat
        xiat=atoms%ratp(1,iat)
        yiat=atoms%ratp(2,iat)
        ziat=atoms%ratp(3,iat)
        do jat=iat+1,atoms%nat
            dx=atoms%ratp(1,jat)-xiat
            dy=atoms%ratp(2,jat)-yiat
            dz=atoms%ratp(3,jat)-ziat
            rsq=dx**2+dy**2+dz**2
            rinvsq=1.d0/rsq
            rinv4=rinvsq**2;rinv6=rinvsq*rinv4;rinv12=rinv6**2
            atoms%epot=atoms%epot+4.d0*(rinv12-rinv6)
            tt=24.d0*rinvsq*(2.d0*rinv12-rinv6)
            t1=dx*tt;t2=dy*tt;t3=dz*tt
            atoms%fat(1,iat)=atoms%fat(1,iat)-t1  
            atoms%fat(2,iat)=atoms%fat(2,iat)-t2  
            atoms%fat(3,iat)=atoms%fat(3,iat)-t3  
            atoms%fat(1,jat)=atoms%fat(1,jat)+t1
            atoms%fat(2,jat)=atoms%fat(2,jat)+t2
            atoms%fat(3,jat)=atoms%fat(3,jat)+t3
        enddo
    enddo
end subroutine lennardjones
!*****************************************************************************************
