!!!module confinement
!!!implicit none
!!!!Confinement parameters
!!!integer:: nconfine       !number of different confinements
!!!integer,allocatable:: conf_dim(:)    !1,2 or 3 for each of the 3 dimensions latvec(:,1), latvec(:,2), latvec(:,3)
!!!integer,allocatable:: conf_av(:)     !0: no confinement (should never occur)
!!!                                     !1: confinement with respect to a fixed value along latvec(:,i)
!!!                                     !2: confinement with respect to the average
!!!integer,allocatable:: conf_exp(:)    !The polynomial order for each confinement
!!!real(8),allocatable:: conf_prefac(:) !The polynomial predactor for each confinement
!!!real(8),allocatable:: conf_cut(:)    !The cutoff distance from each confinement equilibrium
!!!real(8),allocatable:: conf_eq(:)     !Equlibrium position of confinement along the confinement direction, will be filled to average or fixed value
!!!integer,allocatable:: conf_nat(:)    !How many atoms per confinement
!!!integer,allocatable:: conf_list(:,:) !List of atoms per confinement
!!!end module
!!!
!!!program driver
!!!implicit none
!!!integer:: nat, iat
!!!real(8),allocatable:: xcart(:,:),xred(:,:),forces(:,:)
!!!real(8):: latvec(3,3),energy
!!!character(40):: filename
!!!!Read some stupid file
!!!open(unit=2,file="posconf")
!!!read(2,*) nat
!!!allocate(xcart(3,nat),xred(3,nat),forces(3,nat))
!!!read(2,*) latvec(:,1)
!!!read(2,*) latvec(:,2)
!!!read(2,*) latvec(:,3)
!!!do iat=1,nat
!!!  read(2,*) xcart(:,iat)
!!!enddo
!!!close(2)
!!!!Convert to xred
!!!call rxyz_cart2int(latvec,xred,xcart,nat)
!!!write(*,*) "xred"
!!!write(*,*) xred
!!!!Inint confinemant
!!!filename="confine.in"
!!!call init_confinement(nat,filename)
!!!!get confined energy forces
!!!call confinement_energy_forces(nat,xred,latvec,energy,forces)
!!!write(*,*) "Energy",energy
!!!do iat=1,nat
!!!write(*,*) forces(:,iat)
!!!enddo
!!!
!!!
!!!end program

!subroutine init_confinement(nat,filename)
!use confinement, only: nconfine,conf_dim,conf_av,conf_exp,conf_prefac,conf_cut,conf_eq,conf_nat,conf_list,conf_cartred
!use global, only: units
!use defs_basis, only: Bohr_Ang, Ha_eV
!implicit none
!integer:: i,io,iconf,nat
!character(40):: filename,my_fmt,allatoms
!character(200):: line
!open(unit=67,file=trim(filename))
!read(67,*,iostat=io) nconfine
!if(io.lt.0) stop "Could not read nconfine"
!allocate(conf_dim(nconfine),conf_av(nconfine),conf_exp(nconfine),conf_prefac(nconfine),&
!        &conf_cut(nconfine),conf_eq(nconfine),conf_list(nat,nconfine),conf_nat(nconfine),conf_cartred(nconfine))
!conf_cartred="C"
!!Read each block of confinements
!conf_eq=0.d0
!write(*,'(a,a)') " # Constraints in the internal, atomic units. Input interpreted as ",trim(units) 
!do iconf=1,nconfine
!  read(67,'(a200)',end=99) line
!  !Dimension
!  read(line,*,iostat=io) conf_dim(iconf),conf_exp(iconf),conf_prefac(iconf),conf_cut(iconf),conf_av(iconf)
!  if(conf_av(iconf)==1) then
!    read(line,*,iostat=io) &
!    &conf_dim(iconf),conf_exp(iconf),conf_prefac(iconf),conf_cut(iconf),conf_av(iconf),conf_eq(iconf),conf_cartred(iconf)
!    if(conf_cartred(iconf).ne."C".and.conf_cartred(iconf).ne."c".and.conf_cartred(iconf).ne."k".and.conf_cartred(iconf).ne."K"&
!    &.and.conf_cartred(iconf).ne."r".and.conf_cartred(iconf).ne."R".and.conf_cartred(iconf).ne."D".and.conf_cartred(iconf).ne."d")&
!    stop "Provide cartesian or reduced indicator for the equilirium position"
!  endif
!
!  if(io.lt.0) then
!    write(*,*) "Could not read the confinement specifications for No. ",iconf
!    stop
!  endif
!  !Convert here the units if necessary
!  if(trim(units)=="angstroem") then
!   conf_prefac(iconf)=conf_prefac(iconf)/Ha_eV
!   if(conf_cartred(iconf).eq."C".or.conf_cartred(iconf).eq."c".or.conf_cartred(iconf).eq."k".or.conf_cartred(iconf).eq."K") then   
!      conf_eq(iconf)=conf_eq(iconf)/Bohr_Ang
!   endif
!   conf_cut(iconf)=conf_cut(iconf)/Bohr_Ang
!  endif
!  !Read nat
!  !Its possible to select all atoms simply by writing "All" or "all" instead of nat
!  read(67,'(a200)',end=99) line
!  read(line,*) allatoms
!  if(trim(allatoms)=="All".or.trim(allatoms)=="all") then
!    conf_nat(iconf)=nat
!    do i=1,nat
!      conf_list(i,iconf)=i
!    enddo
!  else
!    read(line,*,iostat=io) conf_nat(iconf);if(io.lt.0) stop "Could not read conf_nat"
!    read(67,*,iostat=io) conf_list(1:conf_nat(iconf),iconf);if(io.lt.0) stop "Could not read conf_list"
!  endif
!  !Read list
!  write(my_fmt,'(a,i5,a)') '(a,',conf_nat(iconf),'(i4))'
!  write(*, '(a,i3,a,i3,i3,es15.7,es15.7,i3,es15.7)') ' # Constraint No. ',iconf,' with dim, exp, coeff, cutoff, av, eq: ',&
!  &conf_dim(iconf),conf_exp(iconf),conf_prefac(iconf),conf_cut(iconf),conf_av(iconf),conf_eq(iconf)
!  write(*,my_fmt) ' # Containing atoms: ', (conf_list(i,iconf), i = 1, conf_nat(iconf))
!enddo
!99 continue
!close(67)
!!Check some input stuff
!if(minval(conf_dim(:)).lt.1.or.maxval(conf_dim(:)).gt.3) stop "Wrong dimension confined"
!if(minval(conf_av(:)).lt.1.or.maxval(conf_av(:)).gt.2) stop "Wrong average in confinement"
!end subroutine

subroutine init_confinement_parser(parini)
use mod_parini, only: typ_parini
use global, only: units
use defs_basis, only: Bohr_Ang, Ha_eV
implicit none
type(typ_parini), intent(inout):: parini
integer:: i,io,iconf
character(40):: filename,my_fmt,allatoms
character(200):: line
!Read each block of confinements
!write(*,'(a)') " # Initiallizing confinements"
write(*,'(a,a)') " # Constraints in the internal, atomic units. Input interpreted as ",trim(units) 
do iconf=1,parini%nconfine
  !Dimension
!  read(line,*,iostat=io) conf_dim(iconf),conf_exp(iconf),conf_prefac(iconf),conf_cut(iconf),conf_av(iconf)
  if(parini%conf_av(iconf)==1) then
    if(parini%conf_cartred(iconf).ne."C".and.parini%conf_cartred(iconf).ne."c".and.parini%conf_cartred(iconf).ne."k".and.parini%conf_cartred(iconf).ne."K"&
    &.and.parini%conf_cartred(iconf).ne."r".and.parini%conf_cartred(iconf).ne."R".and.parini%conf_cartred(iconf).ne."D".and.parini%conf_cartred(iconf).ne."d")&
    stop "Provide cartesian or reduced indicator for the equilirium position"
  endif
  !Convert here the units if necessary
  if(trim(units)=="angstroem") then
   parini%conf_prefac(iconf)=parini%conf_prefac(iconf)/Ha_eV
   if(parini%conf_cartred(iconf).eq."C".or.parini%conf_cartred(iconf).eq."c".or.parini%conf_cartred(iconf).eq."k".or.parini%conf_cartred(iconf).eq."K") then   
      parini%conf_eq(iconf)=parini%conf_eq(iconf)/Bohr_Ang
   endif
   parini%conf_cut(iconf)=parini%conf_cut(iconf)/Bohr_Ang
  endif
  !Read nat
  !Its possible to select all atoms simply by writing "All" or "all" instead of nat
  !Read list
  write(my_fmt,'(a,i5,a)') '(a,',parini%conf_nat(iconf),'(i4))'
  write(*, '(a,i3,a,i3,i3,es15.7,es15.7,i3,es15.7)') ' # Constraint in atomic units, No. ',iconf,' with dim, exp, coeff, cutoff, av, eq: ',&
  &parini%conf_dim(iconf),parini%conf_exp(iconf),parini%conf_prefac(iconf),parini%conf_cut(iconf),parini%conf_av(iconf),parini%conf_eq(iconf)
!  write(*,my_fmt) ' # Containing atoms: ', (conf_list(i,iconf), i = 1, conf_nat(iconf))
enddo
!Check some input stuff
if(minval(parini%conf_dim(:)).lt.1.or.maxval(parini%conf_dim(:)).gt.3) stop "Wrong dimension confined"
if(minval(parini%conf_av(:)).lt.1.or.maxval(parini%conf_av(:)).gt.2) stop "Wrong average in confinement"
end subroutine

subroutine confinement_energy_forces(parini,nat,xred,latvec,energy,forces,strten)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: nat,iconf,iat
real(8):: xred(3,nat),latvec(3,3),energy,forces(3,nat),dist,dist_av,nvec(3,3),point0(3),point(3)
real(8):: xcart(3,nat),tt,flat(3,3),xred_ppoint(3),str(3,3),strten(6),vol,fcart_all(3),ft(3)
real(8),allocatable:: conf_eq(:)
!Initiallize
strten=0.d0
energy=0.d0
forces=0.d0
flat=0.d0
!Create cartesian coordinates
do iat=1,nat
  point=xred(:,iat)
!!  point(1)=modulo(modulo(point(1),1.d0),1.d0) 
!!  point(2)=modulo(modulo(point(2),1.d0),1.d0) 
!!  point(3)=modulo(modulo(point(3),1.d0),1.d0) 
  xcart(:,iat)=matmul(latvec,point(:))
enddo
allocate(conf_eq      (parini%nconfine))
conf_eq(1:parini%nconfine)=parini%conf_eq(1:parini%nconfine)
!First compute the average if not provided by input
call nveclatvec(latvec,nvec)
point0=0.d0
do iconf=1,parini%nconfine
   if(parini%conf_av(iconf)==2) then
     dist_av=0.d0
     do iat=1,parini%conf_nat(iconf)
      call dist2plane(xcart(:,parini%conf_list(iat,iconf)),nvec(:,mod(parini%conf_dim(iconf),3)+1),point0,dist)
      dist_av=dist_av+dist
     enddo
     dist_av=dist_av/real(parini%conf_nat(iconf),8)
     conf_eq(iconf)=dist_av
   endif
enddo

!Run over all confinements and compute the forces and energies
do iconf=1,parini%nconfine
! write(*,*) "EQ",conf_eq(iconf)
 if(parini%conf_cartred(iconf).eq."R".or.parini%conf_cartred(iconf).eq."r".or.parini%conf_cartred(iconf).eq."D".or.parini%conf_cartred(iconf).eq."d") then
  point0=0.d0
  call dist2plane(latvec(:,parini%conf_dim(iconf)),nvec(:,mod(parini%conf_dim(iconf),3)+1),point0,dist)
  point0(:)=nvec(:,mod(parini%conf_dim(iconf),3)+1)*conf_eq(iconf)*dist
 else 
  point0(:)=nvec(:,mod(parini%conf_dim(iconf),3)+1)*conf_eq(iconf)
 endif
! write(*,*) "POINT0",point0(:)
 fcart_all=0.d0
 do iat=1,parini%conf_nat(iconf)
!Compute the distance of the atom to the equilibrium plane
    call dist2plane(xcart(:,parini%conf_list(iat,iconf)),nvec(:,mod(parini%conf_dim(iconf),3)+1),point0,dist)
    call rxyz_cart2int(latvec,xred_ppoint,point0,1)
!    write(*,*)"dist",dist
!    write(*,*)"xcart",conf_list(iat,iconf),xcart(:,conf_list(iat,iconf))
!Compute energy and forces and stress
    if(abs(dist).gt.parini%conf_cut(iconf)) then
      energy=energy+parini%conf_prefac(iconf)*(abs(dist)-parini%conf_cut(iconf))**parini%conf_exp(iconf)  
      tt=parini%conf_prefac(iconf)*(abs(dist)-parini%conf_cut(iconf))**(parini%conf_exp(iconf)-1)*parini%conf_exp(iconf)
      ft(1)=-tt*dist/abs(dist)*nvec(1,mod(parini%conf_dim(iconf),3)+1)
      ft(2)=-tt*dist/abs(dist)*nvec(2,mod(parini%conf_dim(iconf),3)+1)
      ft(3)=-tt*dist/abs(dist)*nvec(3,mod(parini%conf_dim(iconf),3)+1)
      forces(1,parini%conf_list(iat,iconf))=forces(1,parini%conf_list(iat,iconf))+ft(1)!-tt*dist/abs(dist)*nvec(1,mod(conf_dim(iconf),3)+1)
      forces(2,parini%conf_list(iat,iconf))=forces(2,parini%conf_list(iat,iconf))+ft(2)!-tt*dist/abs(dist)*nvec(2,mod(conf_dim(iconf),3)+1)
      forces(3,parini%conf_list(iat,iconf))=forces(3,parini%conf_list(iat,iconf))+ft(3)!-tt*dist/abs(dist)*nvec(3,mod(conf_dim(iconf),3)+1)
      write(*,'(a,i5,a,i5,a,3(es15.7))') " # Confinement No: ",iconf,", atom : ,",iat,&
            &", force along confinement :",dot_product(ft,nvec(:,mod(parini%conf_dim(iconf),3)+1)) 
      fcart_all=fcart_all+forces(:,parini%conf_list(iat,iconf))
      call conf_latforce(latvec,parini%conf_dim(iconf),xred(:,parini%conf_list(iat,iconf)),xred_ppoint,str)
      flat=flat+str*tt
    endif     
!    write(*,*)"fcart",conf_list(iat,iconf),forces(:,conf_list(iat,iconf))
 enddo
!Project out translation if conf_av(iconf)==2
 if(parini%conf_av(iconf)==2) then
 fcart_all=fcart_all/real(parini%conf_nat(iconf),8)
 do iat=1,parini%conf_nat(iconf)
    forces(:,parini%conf_list(iat,iconf))=forces(:,parini%conf_list(iat,iconf))-fcart_all
 enddo
 endif
 

enddo
call getvol(latvec,vol)
flat=matmul(flat,transpose(latvec))/vol
!This matrix MUST be symmetric!!!
!write(*,*) flat(:,1)           
!write(*,*) flat(:,2)           
!write(*,*) flat(:,3)           
  strten(1) = flat(1,1)
  strten(2) = flat(2,2)
  strten(3) = flat(3,3)
  strten(6) = flat(2,1)
  strten(5) = flat(3,1)
  strten(4) = flat(3,2)
  deallocate(conf_eq)
end subroutine
!!
!!!-**************************************
!! subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
!! !This subrouine will convert the internal coordinates into cartesian coordinates
!! implicit none
!! real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
!! integer:: nat,iat
!! call invertmat(latvec,latvecinv,3)
!! do iat=1,nat
!!  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
!! enddo
!! end subroutine rxyz_cart2int
!!
!!!************************************************************************************
!!
!! subroutine invertmat(mat,matinv,n)
!! implicit none
!! real(8),intent(in) :: mat(n,n)
!! integer               :: n
!! real(8),allocatable   :: WORK(:)
!! real(8)               :: matinv(n,n),det(3),a(n,n),div
!! integer               :: IPIV(n), INFO
!! integer               :: LDWORK
!! !Here only for a 3*3 matrix
!! if (n==3) then
!! a=mat
!! div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+&
!! &a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
!! div=1.d0/div
!!      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
!!      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
!!      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
!!      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
!!      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
!!      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
!!      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
!!      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
!!      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
!! else
!! !General n*n matrix 
!! matinv=mat
!! allocate(WORK(n))
!! call  DGETRF( n, n, matinv, n, IPIV, INFO )
!! if (info.ne.0) stop "Error in DGETRF"
!! LDWORK=-1
!! call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
!! LDWORK=WORK(1)
!! deallocate(WORK)
!! allocate(WORK(LDWORK))
!! call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
!! if (info.ne.0) stop "Error in DGETRI"
!! endif
!! end
!!!**********************************************************************************************
!!
!! subroutine dist2plane(point,nvec,ppoint,dist)
!! !This subroutine will calculate  the distance between a plane and a point in space
!! !The point is 'point', the normalized normal vector of the plane is 'nvec', 'ppoint' is an arbitrary point on the plane
!! !and the output is the distance 'dist'  
!! real(8), intent(in) :: point(3),nvec(3),ppoint(3)
!! real(8), intent(out):: dist
!! integer             :: i
!! real(8)             :: p,nvectmp(3)
!! nvectmp(:)=nvec(:)!/sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
!! p=DOT_PRODUCT(nvectmp,ppoint) 
!! p=-p
!! dist=DOT_PRODUCT(nvectmp,point)+p
!! end subroutine
!!
!!!**********************************************************************************************
!! subroutine nveclatvec(latvec,nvec)
!! !Will calculate the normalized normal vector to the 3 planes of the cell
!! implicit none
!! real*8, intent(in) :: latvec(3,3)
!! real*8, intent(out):: nvec(3,3)
!! real*8             :: a(3),b(3),crossp(3),norm
!! integer:: i
!! do i=1,3
!! a=latvec(:,i)
!! b=latvec(:,mod(i,3)+1)
!! call cross_product(a,b,crossp)
!! norm=dsqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
!! nvec(:,i)=crossp(:)/norm
!! enddo
!! end subroutine
!! subroutine cross_product(a,b,crossp)
!! !a very simple implementation of the cross product
!! implicit none
!! real(8)::a(3),b(3)
!! real(8)::crossp(3)
!! crossp(1)=a(2)*b(3)-a(3)*b(2)
!! crossp(2)=a(3)*b(1)-a(1)*b(3)
!! crossp(3)=a(1)*b(2)-a(2)*b(1)
!! return
!! end subroutine

subroutine conf_latforce(latvec,conf_dim,xred_point,xred_ppoint,str)
implicit none
real(8):: latvec(3,3),xred_point(3),xred_ppoint(3),str(3,3)
integer:: conf_dim
real(8):: sd1,sd2,sd3,h11,h21,h31,h12,h22,h32,h13,h23,h33

sd1=xred_point(1) - xred_ppoint(1)
sd2=xred_point(2) - xred_ppoint(2)
sd3=xred_point(3) - xred_ppoint(3)
h11=latvec(1,1)
h12=latvec(1,2)
h13=latvec(1,3)
h21=latvec(2,1)
h22=latvec(2,2)
h23=latvec(2,3)
h31=latvec(3,1)
h32=latvec(3,2)
h33=latvec(3,3)



select case(conf_dim)
case(1)
str(1,1) = &
&((sd1)*(h22*h33 - h23*h32)*(((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 +            &
&(h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/     &
&((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) +    &
&h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0)))/             &
&(((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 +            &
&(h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 +     &
&(h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/     &
&((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))**2)**(0.5d0)*((h12*h23 - h13*h22)**2 +& 
&(h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))
 
str(2,1) = &
&-((sd1)*(h12*h33 - h13*h32)*(((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 +         &
&(h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) +                &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 -            &
&h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0)))/(((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33  &
&- h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))**2)**(0.5d0)*((h12*h23  &
&- h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))
 
str(3,1) =&
&((sd1)*(h12*h23 - h13*h22)*(((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 &
&- h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0)))/(((((h22*h33 -         &
&h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 -   &
&h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 -   &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))**2)**(0.5d0)*((h12*h23 - h13*h22)**2 + (h12*h33 &
&- h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))
 
str(1,2) =&
&-((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 &
&- h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - &
&h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 -   &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))*((h33*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - (h23*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) -            &
&((sd2)*(h22*h33 - h23*h32))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +        &
&((2*h23*(h12*h23 - h13*h22) + 2*h33*(h12*h33 - h13*h32))*(h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) - ((2*h23*(h12*h23 - &
&h13*h22) + 2*h33*(h12*h33 - h13*h32))*(h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) + ((2*h23*(h12*h23 - h13*h22) + 2*h33*(h12*h33  &
&- h13*h32))*(h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 &
&+ (h22*h33 - h23*h32)**2)**(1.5d0))))/((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 -            &
&h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0))**2)**(0.5d0)
 
str(2,2) =&
&((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33  &
&- h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - &
&h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 -   &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))*((h33*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - (h13*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) -            &
&((sd2)*(h12*h33 - h13*h32))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +        &
&((2*h13*(h12*h23 - h13*h22) - 2*h33*(h22*h33 - h23*h32))*(h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) - ((2*h13*(h12*h23 - &
&h13*h22) - 2*h33*(h22*h33 - h23*h32))*(h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) + ((2*h13*(h12*h23 - h13*h22) - 2*h33*(h22*h33  &
&- h23*h32))*(h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 &
&+ (h22*h33 - h23*h32)**2)**(1.5d0))))/((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 -            &
&h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0))**2)**(0.5d0)
 
str(3,2) =&
&((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33  &
&- h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - &
&h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 -   &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))*((h13*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - (h23*(h11*(sd1) +      &
&h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +            &
&((sd2)*(h12*h23 - h13*h22))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +        &
&((2*h13*(h12*h33 - h13*h32) + 2*h23*(h22*h33 - h23*h32))*(h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) - ((2*h13*(h12*h33 - &
&h13*h32) + 2*h23*(h22*h33 - h23*h32))*(h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) + ((2*h13*(h12*h33 - h13*h32) + 2*h23*(h22*h33  &
&- h23*h32))*(h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 &
&+ (h22*h33 - h23*h32)**2)**(1.5d0))))/((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 -            &
&h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0))**2)**(0.5d0)
 
str(1,3) =&
&((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33  &
&- h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - &
&h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 -   &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))*((h32*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - (h22*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +            &
&((sd3)*(h22*h33 - h23*h32))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +        &
&((2*h22*(h12*h23 - h13*h22) + 2*h32*(h12*h33 - h13*h32))*(h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) - ((2*h22*(h12*h23 - &
&h13*h22) + 2*h32*(h12*h33 - h13*h32))*(h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) + ((2*h22*(h12*h23 - h13*h22) + 2*h32*(h12*h33  &
&- h13*h32))*(h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 &
&+ (h22*h33 - h23*h32)**2)**(1.5d0))))/((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 -            &
&h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0))**2)**(0.5d0)
 
str(2,3) =&
&-((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 &
&- h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - &
&h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 -   &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))*((h32*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - (h12*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +            &
&((sd3)*(h12*h33 - h13*h32))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +        &
&((2*h12*(h12*h23 - h13*h22) - 2*h32*(h22*h33 - h23*h32))*(h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) - ((2*h12*(h12*h23 - &
&h13*h22) - 2*h32*(h22*h33 - h23*h32))*(h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) + ((2*h12*(h12*h23 - h13*h22) - 2*h32*(h22*h33  &
&- h23*h32))*(h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 &
&+ (h22*h33 - h23*h32)**2)**(1.5d0))))/((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 -            &
&h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0))**2)**(0.5d0)
 
str(3,3) =&
&((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33  & 
&- h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - &
&h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 -   &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0))*((h22*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - (h12*(h21*(sd1) +      &
&h22*(sd2) + h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) +            &
&((sd3)*(h12*h23 - h13*h22))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) -        &
&((2*h12*(h12*h33 - h13*h32) + 2*h22*(h22*h33 - h23*h32))*(h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) + ((2*h12*(h12*h33 - &
&h13*h32) + 2*h22*(h22*h33 - h23*h32))*(h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(1.5d0)) - ((2*h12*(h12*h33 - h13*h32) + 2*h22*(h22*h33  &
&- h23*h32))*(h12*h23 - h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 &
&+ (h22*h33 - h23*h32)**2)**(1.5d0))))/((((h22*h33 - h23*h32)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h12*h23 -            &
&h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) - ((h12*h33 - h13*h32)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 - h23*h32)**2)**(0.5d0) + ((h12*h23 -            &
&h13*h22)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h12*h23 - h13*h22)**2 + (h12*h33 - h13*h32)**2 + (h22*h33 -              &
&h23*h32)**2)**(0.5d0))**2)**(0.5d0)
 
case(2)
str(1,1) = &
&-((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 &
&- h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))*((h33*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - (h23*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) -            &
&((sd1)*(h21*h33 - h23*h31))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +        &
&((2*h23*(h11*h23 - h13*h21) + 2*h33*(h11*h33 - h13*h31))*(h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) - ((2*h23*(h11*h23 - &
&h13*h21) + 2*h33*(h11*h33 - h13*h31))*(h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) + ((2*h23*(h11*h23 - h13*h21) + 2*h33*(h11*h33  &
&- h13*h31))*(h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 &
&+ (h21*h33 - h23*h31)**2)**(1.5d0))))/((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 -            &
&h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(2,1) = &
&((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33  &
&- h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))*((h33*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - (h13*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) -            &
&((sd1)*(h11*h33 - h13*h31))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +        &
&((2*h13*(h11*h23 - h13*h21) - 2*h33*(h21*h33 - h23*h31))*(h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) - ((2*h13*(h11*h23 - &
&h13*h21) - 2*h33*(h21*h33 - h23*h31))*(h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) + ((2*h13*(h11*h23 - h13*h21) - 2*h33*(h21*h33  &
&- h23*h31))*(h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 &
&+ (h21*h33 - h23*h31)**2)**(1.5d0))))/((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 -            &
&h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(3,1) = &
&((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33  &
&- h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))*((h13*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - (h23*(h11*(sd1) +      &
&h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +            &
&((sd1)*(h11*h23 - h13*h21))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +        &
&((2*h13*(h11*h33 - h13*h31) + 2*h23*(h21*h33 - h23*h31))*(h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) - ((2*h13*(h11*h33 - &
&h13*h31) + 2*h23*(h21*h33 - h23*h31))*(h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) + ((2*h13*(h11*h33 - h13*h31) + 2*h23*(h21*h33  &
&- h23*h31))*(h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 &
&+ (h21*h33 - h23*h31)**2)**(1.5d0))))/((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 -            &
&h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(1,2) = &
&((sd2)*(h21*h33 - h23*h31)*(((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 &  
&- h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0)))/(((((h21*h33 -         &
&h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 -   &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))**2)**(0.5d0)*((h11*h23 - h13*h21)**2 + (h11*h33 &
&- h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))
 
str(2,2) = &
&-((sd2)*(h11*h33 - h13*h31)*(((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 +         &
&(h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) +                &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 -            &
&h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0)))/(((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33  &
&- h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))**2)**(0.5d0)*((h11*h23  &
&- h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))
 
str(3,2) = &
&((sd2)*(h11*h23 - h13*h21)*(((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 & 
&- h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0)))/(((((h21*h33 -         &
&h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 -   &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))**2)**(0.5d0)*((h11*h23 - h13*h21)**2 + (h11*h33 &
&- h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))
 
str(1,3) = &
&((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33  &
&- h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))*((h31*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - (h21*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +            &
&((sd3)*(h21*h33 - h23*h31))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +        &
&((2*h21*(h11*h23 - h13*h21) + 2*h31*(h11*h33 - h13*h31))*(h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) - ((2*h21*(h11*h23 - &
&h13*h21) + 2*h31*(h11*h33 - h13*h31))*(h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) + ((2*h21*(h11*h23 - h13*h21) + 2*h31*(h11*h33  &
&- h13*h31))*(h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 &
&+ (h21*h33 - h23*h31)**2)**(1.5d0))))/((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 -            &
&h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(2,3) = &
&-((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 & 
&- h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))*((h31*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - (h11*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +            &
&((sd3)*(h11*h33 - h13*h31))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +        &
&((2*h11*(h11*h23 - h13*h21) - 2*h31*(h21*h33 - h23*h31))*(h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) - ((2*h11*(h11*h23 - &
&h13*h21) - 2*h31*(h21*h33 - h23*h31))*(h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) + ((2*h11*(h11*h23 - h13*h21) - 2*h31*(h21*h33  &
&- h23*h31))*(h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 &
&+ (h21*h33 - h23*h31)**2)**(1.5d0))))/((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 -            &
&h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(3,3) =  &
&((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33  &
&- h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - &
&h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 -   &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0))*((h21*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - (h11*(h21*(sd1) +      &
&h22*(sd2) + h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) +            &
&((sd3)*(h11*h23 - h13*h21))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) -        &
&((2*h11*(h11*h33 - h13*h31) + 2*h21*(h21*h33 - h23*h31))*(h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) + ((2*h11*(h11*h33 - &
&h13*h31) + 2*h21*(h21*h33 - h23*h31))*(h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(1.5d0)) - ((2*h11*(h11*h33 - h13*h31) + 2*h21*(h21*h33  &
&- h23*h31))*(h11*h23 - h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 &
&+ (h21*h33 - h23*h31)**2)**(1.5d0))))/((((h21*h33 - h23*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h23 -            &
&h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) - ((h11*h33 - h13*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 - h23*h31)**2)**(0.5d0) + ((h11*h23 -            &
&h13*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h23 - h13*h21)**2 + (h11*h33 - h13*h31)**2 + (h21*h33 -              &
&h23*h31)**2)**(0.5d0))**2)**(0.5d0)

case(3) 
str(1,1) = &
&-((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 & 
&- h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))*((h32*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - (h22*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) -            &
&((sd1)*(h21*h32 - h22*h31))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +        &
&((2*h22*(h11*h22 - h12*h21) + 2*h32*(h11*h32 - h12*h31))*(h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) - ((2*h22*(h11*h22 - &
&h12*h21) + 2*h32*(h11*h32 - h12*h31))*(h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) + ((2*h22*(h11*h22 - h12*h21) + 2*h32*(h11*h32  &
&- h12*h31))*(h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 &
&+ (h21*h32 - h22*h31)**2)**(1.5d0))))/((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 -            &
&h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(2,1) = &
&((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32  &
&- h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))*((h32*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - (h12*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) -            &
&((sd1)*(h11*h32 - h12*h31))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +        &
&((2*h12*(h11*h22 - h12*h21) - 2*h32*(h21*h32 - h22*h31))*(h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) - ((2*h12*(h11*h22 - &
&h12*h21) - 2*h32*(h21*h32 - h22*h31))*(h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) + ((2*h12*(h11*h22 - h12*h21) - 2*h32*(h21*h32  &
&- h22*h31))*(h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 &
&+ (h21*h32 - h22*h31)**2)**(1.5d0))))/((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 -            &
&h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(3,1) = &
&((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32  &
&- h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))*((h12*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - (h22*(h11*(sd1) +      &
&h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +            &
&((sd1)*(h11*h22 - h12*h21))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +        &
&((2*h12*(h11*h32 - h12*h31) + 2*h22*(h21*h32 - h22*h31))*(h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) - ((2*h12*(h11*h32 - &
&h12*h31) + 2*h22*(h21*h32 - h22*h31))*(h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) + ((2*h12*(h11*h32 - h12*h31) + 2*h22*(h21*h32  &
&- h22*h31))*(h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 &
&+ (h21*h32 - h22*h31)**2)**(1.5d0))))/((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 -            &
&h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(1,2) = &
&((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32  &
&- h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))*((h31*(h21*(sd1) + h22*(sd2) +                  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - (h21*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +            &
&((sd2)*(h21*h32 - h22*h31))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +        &
&((2*h21*(h11*h22 - h12*h21) + 2*h31*(h11*h32 - h12*h31))*(h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) - ((2*h21*(h11*h22 - &
&h12*h21) + 2*h31*(h11*h32 - h12*h31))*(h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) + ((2*h21*(h11*h22 - h12*h21) + 2*h31*(h11*h32  &
&- h12*h31))*(h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 &
&+ (h21*h32 - h22*h31)**2)**(1.5d0))))/((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 -            &
&h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(2,2) = &
&-((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 & 
&- h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))*((h31*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - (h11*(h31*(sd1) +      &
&h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +            &
&((sd2)*(h11*h32 - h12*h31))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +        &
&((2*h11*(h11*h22 - h12*h21) - 2*h31*(h21*h32 - h22*h31))*(h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) - ((2*h11*(h11*h22 - &
&h12*h21) - 2*h31*(h21*h32 - h22*h31))*(h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) + ((2*h11*(h11*h22 - h12*h21) - 2*h31*(h21*h32  &
&- h22*h31))*(h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 &
&+ (h21*h32 - h22*h31)**2)**(1.5d0))))/((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 -            &
&h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(3,2) = &
&((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32  &
&- h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))*((h21*(h11*(sd1) + h12*(sd2) +                  &
&h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - (h11*(h21*(sd1) +      &
&h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) +            &
&((sd2)*(h11*h22 - h12*h21))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) -        &
&((2*h11*(h11*h32 - h12*h31) + 2*h21*(h21*h32 - h22*h31))*(h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) +                    &
&h13*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) + ((2*h11*(h11*h32 - &
&h12*h31) + 2*h21*(h21*h32 - h22*h31))*(h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/(2*((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(1.5d0)) - ((2*h11*(h11*h32 - h12*h31) + 2*h21*(h21*h32  &
&- h22*h31))*(h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/(2*((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 &
&+ (h21*h32 - h22*h31)**2)**(1.5d0))))/((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 -            &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) +  &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 -            &
&h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0))**2)**(0.5d0)
 
str(1,3) = &
&((sd3)*(h21*h32 - h22*h31)*(((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 &
&- h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0)))/(((((h21*h32 -         &
&h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 -   &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))**2)**(0.5d0)*((h11*h22 - h12*h21)**2 + (h11*h32 &
&- h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))
 
str(2,3) = &
&-((sd3)*(h11*h32 - h12*h31)*(((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 +         & 
&(h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) +                &
&h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 -            &
&h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0)))/(((((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32  &
&- h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))**2)**(0.5d0)*((h11*h22  &
&- h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))
 
str(3,3) = &
&((sd3)*(h11*h22 - h12*h21)*(((h21*h32 - h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 &  
&- h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) +  &
&h33*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0)))/(((((h21*h32 -         &
&h22*h31)*(h11*(sd1) + h12*(sd2) + h13*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 -              &
&h22*h31)**2)**(0.5d0) - ((h11*h32 - h12*h31)*(h21*(sd1) + h22*(sd2) + h23*(sd3)))/((h11*h22 - h12*h21)**2 + (h11*h32 -   &
&h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0) + ((h11*h22 - h12*h21)*(h31*(sd1) + h32*(sd2) + h33*(sd3)))/((h11*h22 -   &
&h12*h21)**2 + (h11*h32 - h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))**2)**(0.5d0)*((h11*h22 - h12*h21)**2 + (h11*h32 &
&- h12*h31)**2 + (h21*h32 - h22*h31)**2)**(0.5d0))
end select
end subroutine

 
