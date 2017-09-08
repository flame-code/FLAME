module interface_tinker
  use global
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    init_tinker,  &
    tinker, &
    count_tinker
  integer:: count_tinker = 0
 

contains
subroutine init_tinker(parini,nat,xred,latvec)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer, parameter:: maxbond=8
integer:: nat,iat,nat_tinker,tinker_typat(nat),tinker_bonded(maxbond,nat),j,n
real(8):: rxyz(3,nat),xred(3,nat),latvec(3,3),latvec_ang(3,3),latvec_rot(3,3),dist_ang(6)
character(2):: tinker_char_type(nat)
logical:: file_exists
character(300):: all_line


!tinker.type contains the info of the tinker atom types, including bonding of each atom to others!
  tinker_bonded=0
  INQUIRE(FILE="tinker.type", EXIST=file_exists)
     if(file_exists) then
     open(unit=26,file="tinker.type")
       do iat=1,nat
          read(26,'(a300)') all_line
          n=len(all_line)
          read(all_line(1:n),*,err=60,end=60) tinker_char_type(iat),tinker_typat(iat),(tinker_bonded(j,iat),j=1,8) 
          60 continue
       enddo 
     close(26)
     else
       do iat=1,nat
          tinker_char_type(iat)=char_type(parini%typat_global(iat))
          tinker_typat(iat)=parini%typat_global(iat)
       enddo
     endif
!Initiallize tinker
write(*,*) "Initiallizing tinker"
latvec_ang=latvec*Bohr_Ang
call latvec2dist_ang(dist_ang,latvec_ang,pi)
call dist_ang2latvec(dist_ang,latvec_rot,pi)
call rxyz_int2cart(latvec_rot,xred,rxyz,nat)
call mhm_tinker_init(nat,maxbond,parini%bc,ntypat,rxyz,dist_ang,tinker_char_type,tinker_typat,tinker_bonded)
end subroutine


!********************************************************************************
subroutine tinker(latvec,xred,fcart,energy,strten)
implicit none
integer:: iat
real(8):: latvec(3,3),xred(3,nat),fcart(3,nat),energy,strten(6)
real(8):: dist_ang(6),latvec_ang(3,3),rxyz(3,nat),virial(3,3),latvec_rot(3,3),latvec_rot_inv(3,3),rotmat(3,3),vol
latvec_ang=latvec*Bohr_Ang
call latvec2dist_ang(dist_ang,latvec_ang,pi)
!!!a=sqrt(dot_product(latvec_ang(:,1),latvec_ang(:,1)))
!!!b=sqrt(dot_product(latvec_ang(:,2),latvec_ang(:,2)))
!!!c=sqrt(dot_product(latvec_ang(:,3),latvec_ang(:,3)))
!!!alpha=dot_product(latvec_ang(:,2),latvec_ang(:,3))/b/c
!!!beta=dot_product(latvec_ang(:,3),latvec_ang(:,1))/c/a
!!!gamma=dot_product(latvec_ang(:,1),latvec_ang(:,2))/a/b
call dist_ang2latvec(dist_ang,latvec_rot,pi)
call rxyz_int2cart(latvec_rot,xred,rxyz,nat)
call mhm_tinker_update(rxyz,dist_ang(1),dist_ang(2),dist_ang(3),dist_ang(4),dist_ang(5),dist_ang(6))
call mhm_tinker_energyandforces(energy,fcart,virial)
!Rotate forces back to original cell
call rotmat_fcart_stress(latvec_ang,latvec_rot,rotmat)
do iat=1,nat
   fcart(:,iat)=matmul(rotmat,fcart(:,iat))
enddo
fcart=-fcart/Ha_eV*Bohr_Ang
energy=energy/Ha_eV
strten(1) = virial(1,1)
strten(2) = virial(2,2)
strten(3) = virial(3,3)
strten(6) = virial(2,1)
strten(5) = virial(3,1)
strten(4) = virial(3,2)
call getvol(latvec_rot,vol)
strten=strten/vol
strten=strten/Ha_eV*Bohr_Ang**3
call rotate_stresstensor(strten,rotmat)
end subroutine
!********************************************************************************
end module
