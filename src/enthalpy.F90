subroutine get_enthalpy(latvec,energy,pressure,enthalpy)
!This routine will compute the enthalpy within the given units
implicit none
integer:: natin,iat
character(40):: filename,units
real(8):: acell(3),v(3,3),ucvol,pressure,latvec(3,3),energy,enthalpy

!latvec(:,1)=acell(1)*rprim(:,1)
!latvec(:,2)=acell(2)*rprim(:,2)
!latvec(:,3)=acell(3)*rprim(:,3)

!Compute cell volume
 v=latvec
 ucvol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
        v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
enthalpy=energy+pressure*ucvol
end subroutine
