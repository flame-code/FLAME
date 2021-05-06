subroutine propagate(parini,nat,xred,latvec0,dxred,dlatvec,xredout,latvecout)
use mod_parini, only: typ_parini
!This subroutine will propagate the coordinates of the atoms and the cells according to the
!value of fixed or free degrees of freedom according to
!xred=xred+dxred,latvec=latvec+dlatvec
implicit none
type(typ_parini), intent(in):: parini
integer::nat,i,iat,j
real(8):: xred(3,nat),latvec(3,3),dxred(3,nat),dlatvec(3,3),xredout(3,nat),latvecout(3,3),len1,len2
real(8):: orig_angle(3),new_angle(3),axis(3),rotmat(3,3),center(3),latvec0(3,3)

latvec=latvec0
!write(*,*) "In propagate", fixlat
!Eliminate components not to be changed
if(any(parini%fixlat)) call elim_fixed_lat(parini,latvec,dlatvec)
if(any(parini%fixat))  call elim_fixed_at(parini,nat,dxred)

!Propagate
xredout=xred+dxred
latvecout=latvec+dlatvec
!if(fixlat(7)) write(*,*) "We have fixed cell shape"


!We need to fix the angles since we propagated along a linear direction, not rotational
if(any(parini%fixlat(4:6)).and.(.not.all(parini%fixlat(1:6)))) then
do j=1,10 !Do iteratively
!  write(*,*) "We have fixed cell angles",fixlat(4:6)
  do i=1,3
  if(parini%fixlat(3+i)) then
!First get the original angle
  if(j==1)   orig_angle(i)=acos(dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i+1,3)+1))/&
             &sqrt(dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))*&
             &dot_product(latvec(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))))
             new_angle(i) =acos(dot_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i+1,3)+1))/&
             &sqrt(dot_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i,3)+1))*&
             &dot_product(latvecout(:,modulo(i+1,3)+1),latvecout(:,modulo(i+1,3)+1))))
   write(*,'(a,2(es25.15),1x,i5)') " # Old and New angle: ",orig_angle(i),new_angle(i),j
   call cross_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i+1,3)+1),axis)
   axis=axis/sqrt(dot_product(axis,axis))

   len1=sqrt(dot_product(latvecout(:,modulo(i,3)+1),latvecout(:,modulo(i,3)+1)))
   len2=sqrt(dot_product(latvecout(:,modulo(i+1,3)+1),latvecout(:,modulo(i+1,3)+1)))

   center=0.5d0*(latvecout(:,modulo(i,3)+1)/len1+latvecout(:,modulo(i+1,3)+1)/len2)
   center=center/sqrt(dot_product(center,center))

   call rotation(rotmat,-0.5d0*orig_angle(i),axis)
   latvecout(:,modulo(i,3)+1)=matmul(rotmat,center)*len1
   call rotation(rotmat,0.5d0*orig_angle(i),axis)
   latvecout(:,modulo(i+1,3)+1)=matmul(rotmat,center)*len2
   endif
   enddo
   latvec=latvecout
enddo
endif

!If the cell length were to be fixed, we have to rescale again
if(any(parini%fixlat(1:3))) then
!  write(*,*) "We have fixed cell length",fixlat(1:3)
  do i=1,3
  if(parini%fixlat(i)) then
   len1=dot_product(latvec0(:,i),latvec0(:,i))
   len2=dot_product(latvecout(:,i),latvecout(:,i))
   latvecout(:,i)=latvecout(:,i)*sqrt(len1/len2)
  endif
  enddo
endif
end subroutine propagate
