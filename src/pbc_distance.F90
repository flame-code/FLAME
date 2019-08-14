subroutine pbc_distance0(latvec,xred_1,xred_2,distance2,dxyz)
!This routine computes the distance of the reduced coordinates xred_1 and xred_2 and applies 
!periodic boundary conditions to them and computes the squared distance
implicit none
integer:: i
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3),dxyz(3)
diff=xred_2-xred_1
do i=1,3
    if(diff(i).le.-0.5d0) then
        diff(i)=diff(i)+1.d0
    elseif(diff(i).gt.0.5d0) then
        diff(i)=diff(i)-1.d0
    endif
enddo
dxyz=matmul(latvec,diff)
distance2=dot_product(dxyz,dxyz)
end subroutine

!************************************************************************************

subroutine pbc_distance1(latvec,xred_1,xred_2,distance2)
!This routine computes the distance of the reduced coordinates xred_1 and xred_2 and applies 
!periodic boundary conditions to them and computes the squared distance. Minimal image convention
implicit none
integer:: i
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3)
diff=xred_2-xred_1
do i=1,3
    if(diff(i).le.-0.5d0) then
        diff(i)=diff(i)+1.d0
    elseif(diff(i).gt.0.5d0) then
        diff(i)=diff(i)-1.d0
    endif
enddo
distance2=dot_product(matmul(latvec,diff),matmul(latvec,diff))
end subroutine

!************************************************************************************

subroutine pbc_distance2(latvec,xred_1,xcart_1,xred_2,xcart_2,distance2)
!This routine computes the distance of the reduced coordinates xred_1 and xred_2 and applies 
!periodic boundary conditions to them and computes the squared distance
implicit none
integer:: i,k,l,m
real(8):: xred_1(3),xred_2(3),diff(3),distance2,latvec(3,3),xcart_1(3),xcart_2(3),xcart_tmp(3),xcart_20(3),xcart_10(3)
real(8):: xred_10(3),xred_20(3),distance2_tmp
xred_10(1)=modulo(modulo(xred_1(1),1.d0),1.d0)
xred_10(2)=modulo(modulo(xred_1(2),1.d0),1.d0)
xred_10(3)=modulo(modulo(xred_1(3),1.d0),1.d0)
xcart_10=matmul(latvec,xred_10)
xred_20(1)=modulo(modulo(xred_2(1),1.d0),1.d0)
xred_20(2)=modulo(modulo(xred_2(2),1.d0),1.d0)
xred_20(3)=modulo(modulo(xred_2(3),1.d0),1.d0)
xcart_20=matmul(latvec,xred_20)
distance2=1.d5
do k=-1,1
do l=-1,1
do m=-1,1
  xcart_tmp=xcart_20+real(k,8)*latvec(:,1)+real(l,8)*latvec(:,2)+real(m,8)*latvec(:,3)
  diff=xcart_10-xcart_tmp
  distance2_tmp=dot_product(diff,diff)
  if(distance2_tmp.lt.distance2) then
    distance2=distance2_tmp
    xcart_2=xcart_tmp+xcart_1-xcart_10
  endif
enddo
enddo
enddo
end subroutine
