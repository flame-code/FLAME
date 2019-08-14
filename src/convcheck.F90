subroutine convcheck(parini,nat,latvec_in,fcart_in,strten_in,target_pressure_habohr,strfact,fmax,fmax_at,fmax_lat,tolmxf,iexit)
use mod_parini, only: typ_parini
use defs_basis
implicit none
type(typ_parini), intent(in):: parini
integer:: nat, iexit,iat,istr,i
real(8):: latvec_in(3,3),fcart_in(3,nat),strten_in(6),target_pressure_habohr,fmax,dstr(6)
real(8):: tolmxf,strtarget(6),strfact,fmax_at,fmax_lat
real(8):: gradtarget(3,3),strmattarget(3,3),latvect(3,3),latvectinv(3,3),flattarget(3,3)
real(8):: strmat(3,3),flat(3,3),dflat(3,3),dstrmat(3,3),dist_ang(6)
!Compute maximal component of forces, EXCLUDING any fixed components
fmax_at=0.0d0
fmax_lat=0.0d0
if(.not.any(parini%fixat)) then
 do iat=1,nat
   do i=1,3
       if(abs(fcart_in(i,iat)) >= fmax_at ) fmax_at=abs(fcart_in(i,iat))
   end do
 end do
else
 do iat=1,nat
   do i=1,3
     if (.not.parini%fixat(iat)) then
       if( abs(fcart_in(i,iat)) >= fmax_at ) fmax_at=abs(fcart_in(i,iat))
     end if
   end do
 end do
endif


strtarget=0.d0
strtarget(1:3)=-target_pressure_habohr
if(.not.any(parini%fixlat)) then
 dstr(:)=strten_in(:)-strtarget(:)
!Evaluate the convergence
 do istr=1,6
     if(abs(dstr(istr))*strfact >= fmax_lat ) fmax_lat=abs(dstr(istr))*strfact
 end do
elseif(parini%fixlat(7)) then
 dstr(:)=strten_in(:)-strtarget(:)
 fmax_lat=strfact*abs(((dstr(1)+dstr(2)+dstr(3))/3.d0))
! write(*,'(a,5(es25.15))') "Pressure", HaBohr3_GPa*(strten_in(1)+strten_in(2)+strten_in(3))/3.d0,HaBohr3_GPa*strten_in(1),HaBohr3_GPa*strten_in(2),HaBohr3_GPa*strten_in(3),abs(((dstr(1)+dstr(2)+dstr(3))/3.d0))*HaBohr3_GPa
! call dist_latvec2ang(dist_ang,latvec_in,pi)
! write(*,'(a,6(es25.15))') "angdeg", dist_ang
elseif(all(parini%fixlat(1:6))) then
 fmax_lat=0.d0
else
!Convert sigma target to the gradient target
!Get full stress matrix of target
 strmattarget(1,1)=strtarget(1)
 strmattarget(2,2)=strtarget(2)
 strmattarget(3,3)=strtarget(3)
 strmattarget(1,2)=strtarget(6)
 strmattarget(2,1)=strtarget(6)
 strmattarget(1,3)=strtarget(5)
 strmattarget(3,1)=strtarget(5)
 strmattarget(2,3)=strtarget(4)
 strmattarget(3,2)=strtarget(4)
!Convert the target stress to the target gradient, assuming cell volume to be 1
 latvect(:,1)=latvec_in(1,:)
 latvect(:,2)=latvec_in(2,:)
 latvect(:,3)=latvec_in(3,:)
 call invertmat(latvect,latvectinv,3)
 flattarget=matmul(strmattarget,latvectinv)
!Get full stress matrix
 strmat(1,1)=strten_in(1)
 strmat(2,2)=strten_in(2)
 strmat(3,3)=strten_in(3)
 strmat(1,2)=strten_in(6)
 strmat(2,1)=strten_in(6)
 strmat(1,3)=strten_in(5)
 strmat(3,1)=strten_in(5)
 strmat(2,3)=strten_in(4)
 strmat(3,2)=strten_in(4)
!Convert the target stress to the target gradient, assuming cell volume to be 1
 flat=matmul(strmat,latvectinv)
!Eliminate fixed lattice components from the lattice gradient difference
 dflat=flat-flattarget
 call elim_fixed_lat(parini,latvec_in,dflat)
!Transform back
 dstrmat=-matmul(dflat,latvect)
!Extract components
 dstr(1)=dstrmat(1,1)
 dstr(2)=dstrmat(2,2)
 dstr(3)=dstrmat(3,3)
 dstr(4)=dstrmat(2,3)
 dstr(5)=dstrmat(1,3)
 dstr(6)=dstrmat(1,2)
!Evaluate the convergence
 do istr=1,6
     if(abs(dstr(istr))*strfact >= fmax_lat ) fmax_lat=abs(dstr(istr))*strfact
 end do
endif
fmax=max(fmax_at,fmax_lat)
 iexit=0
 if(fmax.lt.tolmxf) iexit=1
!write(*,*) "FLAT;FAT", fmax_lat,fmax_at
end subroutine
