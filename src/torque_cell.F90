subroutine torque_cell(latvec0,vlat,torquenrm)
implicit none
real(8), intent(in)    :: latvec0(3,3)
real(8), intent(inout) :: vlat(3,3),torquenrm
real(8) :: torque(3),crossp(3),sx,sy,sz,tmax,cx,cy,cz
integer :: i,ii,it,itmax
real(8) :: unitmat(3,3),rotmat(3,3),rotmatall(3,3),xaxis(3),axis(3),latvec(3,3),tnorm,angle,axisnorm
       latvec=latvec0
!Initialize

       torque=0.d0
       do i=1,3
       call cross_product(latvec(:,i),vlat(:,i),crossp)
       torque=torque+crossp
       enddo
       torquenrm=sqrt(torque(1)**2+torque(2)**2+torque(3)**2)
end subroutine torque_cell
