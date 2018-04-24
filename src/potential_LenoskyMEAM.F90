module interface_lenosky_meam
  use global
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    lenosky_meam
 

contains
!********************************************************************************
subroutine lenosky_meam(parini,latvec,xred0,iprec,ka,kb,kc,fcart,energy,strten)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: i,n_in(3),offset_in(3),do_kpt_in,cellsearch,tb_or_meam,n_silicon,nec1,nec2,nec3,iprec,ka,kb,kc,iat
real(8), parameter:: cut=9.d0
real(8):: xred(3,parini%nat),rxyz(3,parini%nat),xred0(3,parini%nat),fcart(3,parini%nat),latvec(3,3),latvec_ang(3,3),strten(6),stress(3,3),energy,count
real(8):: vol
!Only one type of atoms allowed for MEAM
if(parini%ntypat_global.gt.1) stop "Only one type of atoms allowed in MEAM"

!Number of kpoints
n_in(1)=1
n_in(2)=1
n_in(3)=1

!Disable k-points
   do_kpt_in=0
   offset_in(:)=0

!TB or MEAM
tb_or_meam=1 !1 for Si, 2 for Ti

latvec_ang=latvec*Bohr_Ang
xred=xred0
call n_rep_dim(latvec_ang,cut,nec1,nec2,nec3)
cellsearch=max(nec1,max(nec2,nec3))

call backtocell(parini%nat,latvec_ang,xred)
call rxyz_int2cart(latvec_ang,xred,rxyz,parini%nat)

n_silicon=parini%nat
call  lenoskytb(parini%nat,rxyz,fcart,energy,count,n_silicon,latvec_ang(:,1),&
      &latvec_ang(:,2),latvec_ang(:,3),stress(:,1),stress(:,2),stress(:,3),&
      &tb_or_meam,cellsearch,n_in(1),n_in(2),n_in(3),offset_in(1),offset_in(2),offset_in(3),do_kpt_in)
!  stress=-stress
  energy=energy/Ha_eV
  fcart=fcart/Ha_eV*Bohr_Ang
!Get full stress matrix
call getvol(latvec_ang,vol)
stress=matmul(stress,transpose(latvec_ang))/vol
!This matrix MUST be symmetric!!!
!!write(*,*) stress(:,1)           
!!write(*,*) stress(:,2)           
!!write(*,*) stress(:,3)           
  strten(1) = stress(1,1)
  strten(2) = stress(2,2)
  strten(3) = stress(3,3)
  strten(6) = stress(2,1)
  strten(5) = stress(3,1)
  strten(4) = stress(3,2)
!This is not very clear yet...
  strten=strten/Ha_eV*Bohr_Ang**3
end subroutine
!********************************************************************************
subroutine n_rep_dim(latvec,cut,nec1,nec2,nec3)
!This subroutine will return how many periodic expansions for each lattice vector direction are necessary for the periodic boundary conditions
!with for the given cut. nec1,nec2,nec3 for latvec(:,1),latvec(:,2),latvec(:,3)
implicit none
real*8 :: latvec(3,3),cut,nvec(3,3),point(3),point0(3),dist(3),eps,dd
integer:: i
integer:: nec1,nec2,nec3
! eps=1.d-6
nec1=0
nec2=0
nec3=0
call nveclatvec(latvec,nvec)
point0=(/0.d0,0.d0,0.d0/)
do i=1,3
call dist2plane(latvec(:,mod(i+1,3)+1),nvec(:,i),point0,dist(i))
! write(*,*) "cut",i,cut, dist
enddo
dist=abs(dist)
nec1=int(cut/dist(2))+1
nec2=int(cut/dist(3))+1
nec3=int(cut/dist(1))+1
end subroutine

!*********************************************************************************
end module
