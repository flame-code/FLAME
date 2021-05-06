subroutine inertia_tensor(nat,xcart,cmass,amass,intens)
!This routine computes the inertia tensor with respect to the center of mass of a system with nat atoms
implicit none
integer:: nat,iat,i,j
real(8):: xcart(3,nat),amass(nat),intens(3,3),cmass(3),xtmp(3),dist2
real(8):: delta_kronecker

intens=0.d0
do iat=1,nat
xtmp=xcart(:,iat)-cmass(:)
dist2=dot_product(xtmp,xtmp)
  do i=1,3
  do j=1,i
     intens(i,j)=intens(i,j)+amass(iat)*(dist2*delta_kronecker(i,j)-xtmp(i)*xtmp(j))
  enddo
  enddo
enddo

intens(1,2)=intens(2,1)
intens(1,3)=intens(3,1)
intens(2,3)=intens(3,2)
end subroutine
