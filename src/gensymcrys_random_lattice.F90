subroutine random_lattice(LATSGP,CRYSSYS,BRAV,LATVEC,LATVEC_PRIM,TARGET_VOL)
!This subroutine will generate the conventional and primitive lattice matrix,
!with correct volume scaling 
implicit none
INTEGER:: LATSGP, CRYSSYS, NGUESS, I, BRAV
real(8):: LATVEC(3,3),LATVEC_PRIM(3,3),prim2conv_mat(3,3),conv_vol,vol,v(3,3)
real(8):: dist_ang(6),TARGET_VOL
real(8):: tolerance_dist,tolerance_ang,pi,ang_min,ang_max,minlat

pi=acos(-1.d0)
ang_min=60.d0
ang_max=120.d0
tolerance_dist=1.d-2
tolerance_ang=1.d0
minlat=0.3d0
select case(LATSGP)

case (  1:  2)!Triclinic crystal system
  do  
    call random_number(dist_ang(:))
    dist_ang(1:3)=dist_ang(1:3)*(1.d0-minlat)+minlat
    if(abs(dist_ang(1)-dist_ang(2)).lt.tolerance_dist) cycle 
    if(abs(dist_ang(2)-dist_ang(3)).lt.tolerance_dist) cycle 
    if(abs(dist_ang(3)-dist_ang(1)).lt.tolerance_dist) cycle 
    dist_ang(4)=ang_min+dist_ang(4)*(ang_max-ang_min) 
    dist_ang(5)=ang_min+dist_ang(5)*(ang_max-ang_min)
    dist_ang(6)=ang_min+dist_ang(6)*(ang_max-ang_min) 
    if(abs(dist_ang(4)-dist_ang(5)).lt.tolerance_ang) cycle
    if(abs(dist_ang(5)-dist_ang(6)).lt.tolerance_ang) cycle
    if(abs(dist_ang(6)-dist_ang(4)).lt.tolerance_ang) cycle
    exit
  enddo
case (  3: 15)!Monoclinic crystal system, b-unique
  do  
    call random_number(dist_ang(:))
    dist_ang(1:3)=dist_ang(1:3)*(1.d0-minlat)+minlat
    if(abs(dist_ang(1)-dist_ang(2)).lt.tolerance_dist) cycle 
    if(abs(dist_ang(2)-dist_ang(3)).lt.tolerance_dist) cycle 
    if(abs(dist_ang(3)-dist_ang(1)).lt.tolerance_dist) cycle 
    dist_ang(4)=90.d0 
    dist_ang(5)=ang_min+dist_ang(2)*(ang_max-ang_min)
    dist_ang(6)=90.d0
    exit
  enddo
case ( 16: 74)!Orthorhombic crystal system
  do  
    call random_number(dist_ang(:))
    dist_ang(1:3)=dist_ang(1:3)*(1.d0-minlat)+minlat
    if(abs(dist_ang(1)-dist_ang(2)).lt.tolerance_dist) cycle 
    if(abs(dist_ang(2)-dist_ang(3)).lt.tolerance_dist) cycle 
    if(abs(dist_ang(3)-dist_ang(1)).lt.tolerance_dist) cycle 
    dist_ang(4)=90.d0 
    dist_ang(5)=90.d0
    dist_ang(6)=90.d0
    exit
  enddo
case ( 75:142)!Tetragonal crystal system
  do  
    call random_number(dist_ang(:))
    dist_ang(1:3)=dist_ang(1:3)*(1.d0-minlat)+minlat
    dist_ang(2)=dist_ang(1) 
    if(abs(dist_ang(3)-dist_ang(1)).lt.tolerance_dist) cycle 
    dist_ang(4)=90.d0 
    dist_ang(5)=90.d0
    dist_ang(6)=90.d0
    exit
  enddo
case (143:167)!Trigonal crystal system, also hexagonal cells with rombohedral unit cells (146,148,155,160,161,166,167)
  do  
    call random_number(dist_ang(:))
    dist_ang(1:3)=dist_ang(1:3)*(1.d0-minlat)+minlat
    dist_ang(2)=dist_ang(1) 
    dist_ang(4)=90.d0 
    dist_ang(5)=90.d0
    dist_ang(6)=120.d0
    exit
  enddo
case (168:194)!Hexagonal crystal system
  do  
    call random_number(dist_ang(:))
    dist_ang(1:3)=dist_ang(1:3)*(1.d0-minlat)+minlat
!    if(dist_ang(1).lt.dist_ang(3).or.dist_ang(2).lt.dist_ang(3)) cycle  !make sure that the c-axis is the largest
    dist_ang(2)=dist_ang(1) 
    dist_ang(4)=90.d0 
    dist_ang(5)=90.d0
    dist_ang(6)=120.d0
    exit
  enddo
case (195:230)!Cubic crystal system
  do  
    call random_number(dist_ang(:))
    dist_ang(1:3)=dist_ang(1:3)*(1.d0-minlat)+minlat
    dist_ang(2)=dist_ang(1) 
    dist_ang(3)=dist_ang(1) 
    dist_ang(4)=90.d0 
    dist_ang(5)=90.d0
    dist_ang(6)=90.d0
    exit
  enddo
case default
stop "Wrong space group number"
end select
call dist_ang2latvec_gensymcrys(dist_ang,latvec,pi)
call make_cell_prim(latvec,latvec_prim,BRAV)
call scale_vol(TARGET_VOL,latvec_prim)
call prim2conv(BRAV,prim2conv_mat)

 v=prim2conv_mat
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
conv_vol=TARGET_VOL/vol
call scale_vol(conv_vol,latvec)
end subroutine

subroutine dist_ang2latvec_gensymcrys(dist_ang,latvec,pi)
!This subroutine will generate the lattice vector representation of the cell
!from the length/angle representation
implicit none
real(8):: dist_ang(6),latvec(3,3),pi,convang
convang=pi/180.d0

latvec(1,1)=dist_ang(1)
latvec(2,1)=0.d0
latvec(3,1)=0.d0
latvec(1,2)=dist_ang(2)*cos(convang*dist_ang(6))
latvec(2,2)=dist_ang(2)*sin(convang*dist_ang(6))
latvec(3,2)=0.d0
latvec(1,3)=dist_ang(3)*cos(convang*dist_ang(5))
latvec(2,3)=(dist_ang(2)*dist_ang(3)*cos(convang*dist_ang(4))-latvec(1,2)*latvec(1,3))/latvec(2,2)
latvec(3,3)=sqrt(dist_ang(3)**2-latvec(1,3)**2-latvec(2,3)**2)
end subroutine

subroutine scale_vol(TARGET_VOL,latvec)
!This subroutine will scale the cell in such a way that the volume of the cell will be the same as
!as the target volume TARGET_vol
implicit none
real(8):: latvec(3,3),TARGET_VOL,v(3,3),vol,ratio
 v=latvec
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
 ratio=TARGET_VOL/vol
 latvec(:,:)=latvec(:,:)*(ratio**(1.d0/3.d0))
end subroutine

subroutine make_cell_prim(latvec,prim_lat,BRAV)
!This subroutine will convert the given cell latvec (initially conventional cell) into a new, primitive cell prim_lat,
!according to the transformation proposed in international table of crystallography
!latvec_prim=latvec*P, where P=Mat(P-->A,B,C,F,I,H)
!For your info:
!xp=P^-1*xc, where xp and xc contain the reduced coordinates in the primitive and conventional cell, respectively.
implicit none
real(8):: latvec(3,3), prim_lat(3,3), conv2prim_mat(3,3), prim2conv_mat(3,3)
integer:: BRAV
call prim2conv(BRAV,prim2conv_mat)
prim_lat=matmul(latvec,prim2conv_mat)
end subroutine
