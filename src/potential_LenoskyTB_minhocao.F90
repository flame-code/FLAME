module interface_lenosky_tb
  use global
  use defs_basis
  !use cell_utils

  implicit none

  private
  public :: &
    check_lenosky_tb,  &
    lenosky_tb
 

contains
subroutine check_lenosky_tb(parini,n_silicon)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
integer:: n_silicon,iat
logical:: in_h
!Check typat for consistentcy and provide n_silicon
in_h=.false.
n_silicon=0
do iat=1,parini%nat
   if(int(parini%znucl(parini%typat_global(iat))).ne.1.and.int(parini%znucl(parini%typat_global(iat))).ne.14.and.int(parini%znucl(parini%typat_global(iat))).lt.200) stop "Lenosky TB only implemented for Si and H"
   if(.not.parini%voids.and.int(parini%znucl(parini%typat_global(iat))).ge.200) stop "LJ particles only allowed when using voids"
   if(int(parini%znucl(parini%typat_global(iat)))==14) n_silicon=n_silicon+1
   if(int(parini%znucl(parini%typat_global(iat)))==14.and.in_h) stop "Lenosky TB: First Si, then H"
   if(int(parini%znucl(parini%typat_global(iat)))==1) in_h=.true.
enddo
end subroutine

!********************************************************************************
subroutine lenosky_tb(parini,latvec,xred0,iprec,ka,kb,kc,fcart,energy,strten)
use mod_parini, only: typ_parini
!All H-Atoms have to be at the end of the rxyz-array
implicit none
type(typ_parini), intent(in):: parini
integer:: i,n_in(3),offset_in(3),do_kpt_in,cellsearch,tb_or_meam,n_silicon,nec1,nec2,nec3,iprec,ka,kb,kc,iat
real(8), parameter:: cut=5.24d0
real(8):: xred(3,parini%nat),rxyz(3,parini%nat),xred0(3,parini%nat),fcart(3,parini%nat),latvec(3,3),latvec_ang(3,3),strten(6),stress(3,3),energy,count
real(8):: vol

!Check if the atomic types are allright
call check_lenosky_tb(parini,n_silicon)

!Number of kpoints
n_in(1)=ka
n_in(2)=kb
n_in(3)=kc

!Enable k-points
if(ka==1.and.kb==1.and.kc==1) then 
   do_kpt_in=0
!Offset
   offset_in(:)=0
else
   do_kpt_in=1
!Offset
   offset_in(:)=1
endif

!TB or MEAM
tb_or_meam=0

call check_lenosky_tb(parini,n_silicon)

latvec_ang=latvec*Bohr_Ang
xred=xred0
call n_rep_dim(latvec_ang,cut,nec1,nec2,nec3)
cellsearch=max(nec1,max(nec2,nec3))

call backtocell(parini%nat,latvec_ang,xred)
call rxyz_int2cart(latvec_ang,xred,rxyz,parini%nat)

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
