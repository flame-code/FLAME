 subroutine dist2line(point,ppoint1,ppoint2,dist)
 !This subroutine will calculate a the distance between a plane and a line in space
 !The point is 'point', 'ppoint1' and 'ppoint2' are arbitrary points on the line
 !and the output is the distance 'dist'  
 implicit none
 real(8), intent(in) :: point(3),ppoint1(3),ppoint2(3)
 real(8), intent(out):: dist
 integer             :: i
 real(8)             :: v1(3),v2(3),vnrm,crossp(3)
 v1(:)=ppoint2(:)-ppoint1(:)
 vnrm=sqrt(v1(1)*v1(1)+v1(2)*v1(2)+v1(3)*v1(3))
 v1(:)=v1(:)/vnrm
 v2(:)=point(:)-ppoint1(:)
 call cross_product(v1,v2,crossp)
 dist=sqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 end subroutine

!**********************************************************************************************

 subroutine dist2plane(point,nvec,ppoint,dist)
 !This subroutine will calculate  the distance between a plane and a point in space
 !The point is 'point', the normalized normal vector of the plane is 'nvec', 'ppoint' is an arbitrary point on the plane
 !and the output is the distance 'dist'  
 implicit none
 real(8), intent(in) :: point(3),nvec(3),ppoint(3)
 real(8), intent(out):: dist
 integer             :: i
 real(8)             :: p,nvectmp(3)
 nvectmp(:)=nvec(:)!/sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
 p=DOT_PRODUCT(nvectmp,ppoint) 
 p=-p
 dist=DOT_PRODUCT(nvectmp,point)+p
 end subroutine

!**********************************************************************************************

subroutine dist_ang2latvec(dist_ang,latvec,pi)
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
subroutine dist_latvec2ang(dist_ang,latvec,pi)
!This subroutine will generate the angdeg represenation of the cell from the lattice vectors
implicit none
real(8):: dist_ang(6),latvec(3,3),pi,convang
convang=180.d0/pi
dist_ang(1)=sqrt(dot_product(latvec(:,1),latvec(:,1)))
dist_ang(2)=sqrt(dot_product(latvec(:,2),latvec(:,2)))
dist_ang(3)=sqrt(dot_product(latvec(:,3),latvec(:,3)))
dist_ang(4)=acos(dot_product(latvec(:,2),latvec(:,3))/(dist_ang(2)*dist_ang(3)))*convang
dist_ang(5)=acos(dot_product(latvec(:,3),latvec(:,1))/(dist_ang(3)*dist_ang(1)))*convang
dist_ang(6)=acos(dot_product(latvec(:,1),latvec(:,2))/(dist_ang(1)*dist_ang(2)))*convang
end subroutine

!**********************************************************************************************

 subroutine dproj2latvec(dproj,latvec)
 !This subroutine will convert the distance and projective representation of 
 !a periodic cell (dxx,dyx,dyy,dzx,dzy,dzz) into a 
 !lattice vektor format (vec1(:,1),vec2(:,2),vec3(:,3)) with dxx oriented into x direction
 implicit none
 real*8:: dproj(6),latvec(3,3)

 latvec(:,:)=0.d0
 latvec(1,1)=dproj(1)
 latvec(1,2)=dproj(2)
 latvec(2,2)=dproj(3)
 latvec(1,3)=dproj(4)
 latvec(2,3)=dproj(5)
 latvec(3,3)=dproj(6)
 return
 end subroutine

!**********************************************************************************************

 subroutine latvec2dproj(dproj,latvec,rotmat,rxyz,nat)
 !This subroutine will convert the lattice vector representation of thei
 !periodic cell (vec1,vec2,vec3) into the projective representation (dxx,dyx,dyy,dzx,dzy,dzz)
 !The cell will thus be rotated. The rotational matrix is stored in rotmat as an operator rotmat
 !and the atomic position rxyz are transformed into the new coordination sizstem as well
 implicit none
 integer,intent(in)  :: nat
 integer             :: iat,i
 real*8,intent(inout):: dproj(6),latvec(3,3),rotmat(3,3),rxyz(3,nat)
 real*8  :: tempvec(3),rotmat1(3,3),rotmat2(3,3),crossp(3),alpha,latvect(3,3)
 real*8  :: eps,axe(3),vnrm,rxyzt(3),rotmatt(3,3)
 eps=1.d-6
 !Calculating dxx
 dproj(1)=sqrt(latvec(1,1)*latvec(1,1)+latvec(2,1)*latvec(2,1)+latvec(3,1)*latvec(3,1))

 !Calculate the first rotation to align the first axis in x direction
 rotmat1(:,:)=0.d0
 do i=1,3
    rotmat1(i,i)=1.d0
 enddo
 tempvec(:)=0.d0
 tempvec(1)=1.d0
 !tempvec is the x-unit vector
 call cross_product(latvec(:,1),tempvec,crossp)
! if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. abs(crossp(3)).lt.eps*1.d-1) goto 1001 !no rotation needed
 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. abs(crossp(3)).lt.eps*1.d-1) then
!Check if the latvec points along the positive x-rection
    alpha=dot_product(tempvec(:),latvec(:,1))/sqrt(dot_product(latvec(:,1),latvec(:,1)))
    if(alpha.gt.0.d0) then
        goto 1001 !no rotation needed
    else
!Rotate along the y-axis
    axe(:)=0.d0
    axe(2)=1.d0
    alpha=dacos(-1.d0)
    call rotation(rotmat1,alpha,axe)
    latvec(:,:)=matmul(rotmat1(:,:),latvec(:,:))
        goto 1001
    endif
 endif
 vnrm=sqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 axe(:)=crossp(:)/vnrm
 alpha=dacos(dot_product(tempvec,latvec(:,1))/dproj(1))
 call rotation(rotmat1,alpha,axe)
 latvec(:,:)=matmul(rotmat1(:,:),latvec(:,:))
! call DGEMM('N','N',3,3,3,1.d0,rotmat,3,latvec,3,0.d0,latvect,3)
! latvec=latvect

 1001 continue
 if(latvec(2,1).gt.eps .or. latvec(3,1).gt.eps) then
 write(*,*) "Error in 1. rotation",latvec(2,1),latvec(3,1)
 stop
 endif
 !Calculate the second rotation to align the second axis in xy plane 
 rotmat2(:,:)=0.d0
 do i=1,3
    rotmat2(i,i)=1.d0
 enddo
! axe(:)=0.d0
! axe(1)=1.d0
 tempvec(:)=latvec(:,2)
 tempvec(1)=0.d0
 call cross_product(tempvec,(/0.d0,1.d0,0.d0/),axe)

! if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1 .and. tempvec(2).gt.0.d0) goto 1002 !no rotation needed
 if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1) then
 if (tempvec(2).gt.0.d0) then  !     goto 1002 !no rotation needed
!Check if the latvec points along the positive x-rection
    goto 1002 !no rotation needed
 else
!Rotate along the y-axis
    axe(:)=0.d0
    axe(1)=1.d0
    alpha=dacos(-1.d0)
    call rotation(rotmat2,alpha,axe)
    latvec(:,:)=matmul(rotmat2(:,:),latvec(:,:))
        goto 1002
 endif
 endif

 vnrm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
 axe=axe/vnrm
 call cross_product(axe,latvec(:,2),crossp)

 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002
 vnrm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
 axe=axe/vnrm
 call cross_product(axe,latvec(:,2),crossp)
! if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002 !no rotation needed
 if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1 .and. tempvec(2).gt.0.d0) goto 1002 !no rotation needed
 vnrm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
 axe=axe/vnrm
 call cross_product(axe,latvec(:,2),crossp)
! if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002 !no rotation needed

 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) then
!Check if the latvec points along the positive y-rection
    alpha=dot_product(axe(:),latvec(:,2))/sqrt(dot_product(latvec(:,2),latvec(:,2)))
    if(alpha.gt.0.d0) then
        goto 1002 !no rotation needed
    else
!Rotate along the y-axis
    axe(:)=0.d0
    axe(1)=1.d0
    alpha=dacos(-1.d0)
    call rotation(rotmat2,alpha,axe)
    latvec(:,:)=matmul(rotmat2(:,:),latvec(:,:))
        goto 1002
    endif
 endif
 vnrm=sqrt(tempvec(2)*tempvec(2)+tempvec(3)*tempvec(3))
 alpha=dot_product((/0.d0,1.d0,0.d0/),tempvec(:))/vnrm
 alpha=dacos(alpha)
 call rotation(rotmat2,alpha,axe)
 latvec(:,:)=matmul(rotmat2(:,:),latvec(:,:))
! call DGEMM('N','N',3,3,3,1.d0,rotmat,3,latvec,3,0.d0,latvect,3)
! latvec=latvect
 1002 continue
 if(latvec(3,2).gt.eps) then! stop "Error in 2. rotation" 
! write(*,*) latvec(:,1)
! write(*,*) latvec(:,2)
 write(*,*) "Error in 2. rotation"
 stop
 endif
 if(latvec(3,3).lt.0.d0) stop "Error in orientation of the cell"

 !The total rotational matrix:
 rotmat=matmul(rotmat2(:,:),rotmat1(:,:))
! call DGEMM('N','N',3,3,3,1.d0,rotmat2,3,rotmat1,3,0.d0,rotmatt,3)
! rotmat=rotmatt 

 !Apply rotation on all atoms
 do iat=1,nat
     rxyz(:,iat)=matmul(rotmat,rxyz(:,iat))
!     call DGEMM('N','N',3,3,3,1.d0,rotmat,3,latvec,3,0.d0,latvect,3)
!     rxyz(:,iat)=rxyzt(:)
 enddo

 !Calculate all other elements of dproj
 dproj(2)=latvec(1,2)
 dproj(3)=latvec(2,2)
 dproj(4)=latvec(1,3)
 dproj(5)=latvec(2,3)
 dproj(6)=latvec(3,3)
 end subroutine

!**********************************************************************************************

 subroutine latvec2acell_rprim(latvec,acell,rprim)
!This routine will split up the latvec into two parts, 
!acell and rprim, where rprim is normalized to 1
 implicit none
 real(8):: latvec(3,3), rprim(3,3), acell(3)
 integer:: i
 do i=1,3
   acell(i)=sqrt(latvec(1,i)**2+latvec(2,i)**2+latvec(3,i)**2)
   rprim(:,i)=latvec(:,i)/acell(i)
 enddo
 end subroutine latvec2acell_rprim

!**********************************************************************************************

subroutine getvol(latvec,vol)
implicit none
real(8):: latvec(3,3),v(3,3),vol
 v=latvec
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
end subroutine

!**********************************************************************************************

 subroutine acell_rprim2latvec(latvec,acell,rprim)
!This routine will combine 
!acell and rprim to latvec
 implicit none
 real(8):: latvec(3,3), rprim(3,3), acell(3)
 integer:: i
 do i=1,3
   latvec(:,i)=acell(i)*rprim(:,i)
 enddo
 end subroutine acell_rprim2latvec

!**********************************************************************************************

 subroutine nveclatvec(latvec,nvec)
 !Will calculate the normalized normal vector to the 3 planes of the cell
 implicit none
 real*8, intent(in) :: latvec(3,3)
 real*8, intent(out):: nvec(3,3)
 real*8             :: a(3),b(3),crossp(3),norm
 integer:: i
 do i=1,3
 a=latvec(:,i)
 b=latvec(:,mod(i,3)+1)
 call cross_product(a,b,crossp)
 norm=dsqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 nvec(:,i)=crossp(:)/norm
 enddo
 end subroutine

!**********************************************************************************************

 subroutine rxyz_cart2int(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
 call invertmat(latvec,latvecinv,3)
 do iat=1,nat
  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
 enddo
 end subroutine rxyz_cart2int

!**********************************************************************************************

 subroutine rxyz_int2cart(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3)
 integer:: nat,iat
 do iat=1,nat
  rxyzcart(:,iat)=matmul(latvec,rxyzint(:,iat))
 enddo
 end subroutine rxyz_int2cart 

!**********************************************************************************************

 subroutine fxyz_cart2int(nat,fxyz_cart,fxyz_int,latvec)
 !This subrtouine will transform theforces initially in the cartesian system into the internal coordinates with respect to the
 !cell vectors provided in the latvec
 implicit none
 real(8):: fxyz_cart(3,nat),fxyz_int(3,nat),latvec(3,3),transmat(3,3)
 integer:: nat,iat
 transmat(1,:)=latvec(:,1)
 transmat(2,:)=latvec(:,2)
 transmat(3,:)=latvec(:,3)
 do iat=1,nat
 fxyz_int(:,iat)=matmul(transmat,fxyz_cart(:,iat))
 enddo
 end subroutine fxyz_cart2int 

!**********************************************************************************************

subroutine k_expansion(parini,latvec,xred,ka,kb,kc,k_latvec,k_xcart)
use mod_parini, only: typ_parini

!This routine expands the cell defined by latevec to a supercell
!of dimension ka,kb,kc. The atomic positions in real space
!will be returned in k_xcart. k_nat will then be the number of
!all atoms in the supercell and is ka*kb*kc*nat
implicit none
type(typ_parini), intent(in):: parini
real(8):: latvec(3,3),k_latvec(3,3),k_xcart(3,parini%nat,ka,kb,kc),xred(3,parini%nat) 
integer:: iat,k,l,m,ka,kb,kc
do k=1,ka
do l=1,kb
do m=1,kc
do iat=1,parini%nat
   k_xcart(:,iat,k,l,m)=matmul(latvec,xred(:,iat))+&
   &real(k-1,8)*latvec(:,1)+real(l-1,8)*latvec(:,2)+real(m-1,8)*latvec(:,3)
enddo
enddo
enddo
enddo
k_latvec(:,1)=real(ka,8)*latvec(:,1)
k_latvec(:,2)=real(kb,8)*latvec(:,2)
k_latvec(:,3)=real(kc,8)*latvec(:,3)
end subroutine

!**********************************************************************************************

 subroutine strten2flat(strten,flat,latvec,press)
 !flat is the force on the lettice vector per unit cell volume
 implicit none
 real(8):: strten(6),flat(3,3),latvec(3,3),press,pressmat(3,3),str_matrix(3,3),latvect(3,3),latvectinv(3,3),vol
 pressmat=0.d0
 pressmat(1,1)=press
 pressmat(2,2)=press
 pressmat(3,3)=press
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
           flat=str_matrix+pressmat
           flat=-flat
 latvect(:,1)=latvec(1,:)
 latvect(:,2)=latvec(2,:)
 latvect(:,3)=latvec(3,:)
 call invertmat(latvect,latvectinv,3)
 flat=matmul(flat,latvectinv)
! I guess i forgot to divide through volume
 call getvol(latvec,vol)
 flat=flat*vol
 end subroutine

!**********************************************************************************************

subroutine rotate_stresstensor(strten,rotmat)
!This subroutine will rotate the stress tensor by rotmat according to rotmat*stress*rotmat^T
implicit none
real(8):: strten(6),rotmat(3,3),stress(3,3)
        stress(1,1) =  strten(1) 
        stress(2,2) =  strten(2) 
        stress(3,3) =  strten(3) 
        stress(2,1) =  strten(6) 
        stress(3,1) =  strten(5) 
        stress(3,2) =  strten(4) 
        stress(1,2) =  stress(2,1)
        stress(1,3) =  stress(3,1)
        stress(2,3) =  stress(3,2)
           stress=matmul(rotmat,matmul(stress,transpose(rotmat)))
        strten(1) =  stress(1,1)
        strten(2) =  stress(2,2)
        strten(3) =  stress(3,3)
        strten(6) =  stress(2,1)
        strten(5) =  stress(3,1)
        strten(4) =  stress(3,2)
end subroutine

!**********************************************************************************************

 subroutine find_kpt(k1, k2, k3, lat, gridden)
! This code will define the KPT mesh based on the desired grid density
   implicit none
   integer, intent(out) :: k1,k2,k3
   real(8), intent(in)  :: lat(3,3), gridden
   
   integer :: i, j
   real(8) :: lat1, lat2, lat3
   real(8) :: angles(3), cos_arr(3)
   real(8) :: glat(3,3), crossp(3), vol, glen(3), a(3,3)
   real(8) :: pi
   
   pi = acos(-1.d0)
   
   lat1 = sqrt(lat(1,1)**2 + lat(2,1)**2 + lat(3,1)**2)
   lat2 = sqrt(lat(1,2)**2 + lat(2,2)**2 + lat(3,2)**2)
   lat3 = sqrt(lat(1,3)**2 + lat(2,3)**2 + lat(3,3)**2)
   
   cos_arr(1) = (lat(1,2)*lat(1,3) + lat(2,2)*lat(2,3) + lat(3,2)*lat(3,3))/(lat2*lat3)
   cos_arr(2) = (lat(1,1)*lat(1,3) + lat(2,1)*lat(2,3) + lat(3,1)*lat(3,3))/(lat1*lat3)
   cos_arr(3) = (lat(1,1)*lat(1,2) + lat(2,1)*lat(2,2) + lat(3,1)*lat(3,2))/(lat1*lat2)
   
   angles(1) = (acos(cos_arr(1))/pi)*180.d0
   angles(2) = (acos(cos_arr(2))/pi)*180.d0
   angles(3) = (acos(cos_arr(3))/pi)*180.d0
   
   call getvol(lat,vol)
   
   call cross_product(lat(:,2), lat(:,3), crossp(:))
   glat(:,1) = 2.d0*pi*crossp(:)/vol
   call cross_product(lat(:,3), lat(:,1), crossp(:))
   glat(:,2) = 2.d0*pi*crossp(:)/vol
   call cross_product(lat(:,1), lat(:,2), crossp(:))
   glat(:,3) = 2.d0*pi*crossp(:)/vol
   
   !Compute the correct kpts
   glen(1) = sqrt(glat(1,1)**2 + glat(2,1)**2 + glat(3,1)**2)
   glen(2) = sqrt(glat(1,2)**2 + glat(2,2)**2 + glat(3,2)**2)
   glen(3) = sqrt(glat(1,3)**2 + glat(2,3)**2 + glat(3,3)**2)
   
   call track_kpt(gridden, glen(1), k1)
   call track_kpt(gridden, glen(2), k2)
   call track_kpt(gridden, glen(3), k3)
   
 end subroutine find_kpt

!**********************************************************************************************
   
   subroutine track_kpt(gridden, glen, kpt)
     implicit none
     real(8), intent(in) :: gridden, glen
     integer :: kpt,j
     real(8) :: d_test
   real(8) :: pi
   
   pi = acos(-1.d0)
     
     kpt = int(glen/(gridden*2.d0*pi))
     if (kpt == 0) kpt = 1
     d_test=glen/(kpt*2.d0*pi)
     if (d_test.ge.gridden) then
       do j = 1, 25
         kpt = kpt + j
         d_test = glen/(kpt*2.d0*pi)
         if (d_test.le.gridden) exit
       enddo
     endif
   end subroutine track_kpt

!**********************************************************************************************

subroutine rotmat_fcart_stress(latvec_init,latvec_trans,rotmat)
!This subroutine will compute a rotation matrix, which transforms
!fcart_trans into the original orientation forces fcart by fcart=matmul(rotmat,fcart_trans)
!stress_trans into the original orientation stress by stress=rotmat*stress_trans*rotnat^T
implicit none
real(8):: latvec_init(3,3),latvec_trans(3,3),latvec_trans_inv(3,3),rotmat(3,3)
call invertmat(latvec_trans,latvec_trans_inv,3)
rotmat=matmul(latvec_init,latvec_trans_inv)
end subroutine

!**********************************************************************************************

 subroutine updaterxyz(latvecold,latvecnew,rxyz,nat)
 !This subroutine will update the atomic positions in the cartesian coordinates after the cell shape has been changed according
 !to the change in the lattice vectors thus keeping the relative coordinates of all atoms constant
 implicit none
 real(8), intent(in)   :: latvecold(3,3), latvecnew(3,3)
 real(8), intent(inout):: rxyz(3,nat)
 integer, intent(in)   :: nat
 real(8)               :: latvecold_inv(3,3),trafo(3,3)
 integer               :: iat
 call invertmat(latvecold,latvecold_inv,3)
 trafo=matmul(latvecnew,latvecold_inv) !Transform atomic positions according to cell
 do iat=1,nat
 rxyz(:,iat)=matmul(trafo,rxyz(:,iat))
 enddo
 end subroutine
