!*************************************************************************
!                                                                        *
! TOOLS FOR PERIODIC CELLS                                               *                   
!                                                                        *
!*************************************************************************

 subroutine backtocell_gensymcrys(nat,latvec,rxyz)
 !This subroutine will transform back all atoms into the periodic cell
 !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
 implicit none
 integer:: nat,i,iat,j
 real(8) :: latvec(3,3), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 logical:: neccesary
 !To really be on the safe side, the translation vector can be shortened by  a factor eps in order
 !to get the atom into the cell. 
!  eps=1.d-10
  eps=1.d-15
 ! eps=1.d0-eps
 count=0.d0
 neccesary=.true.
 do while(neccesary)
 neccesary=.false.
 count=count+1.d0
 !generate 3 normal vectors of the 3 planes
 call nveclatvec_gensymcrys(latvec,nvec)
 do iat=1,nat
 !3 planes through origin (xy,yz,zx)
 do i=1,3
 dist(i)=DOT_PRODUCT(rxyz(:,iat),nvec(:,i))
! if(dist(i).lt.0.d0) then
 if(dist(i).lt.-abs(dist(i))*eps) then
! write(*,*) "unten 1",i,iat
 rxyz(:,iat)=rxyz(:,iat)+latvec(:,mod(i+1,3)+1)!*eps
 neccesary=.true.
 endif

 !3 planes on top/side/back (xy,yz,zx)
 dist(i+3)=DOT_PRODUCT(rxyz(:,iat)-latvec(:,mod(i+1,3)+1),nvec(:,i))
! if(dist(i+3).gt.0.d0) then
 if(dist(i+3).gt.abs(dist(i+3))*eps) then
! write(*,*) "unten 1",i,iat
! if(dist(i+3).gt.eps) then
! write(*,*) "oben 2",i,iat
 rxyz(:,iat)=rxyz(:,iat)-latvec(:,mod(i+1,3)+1)!*eps
 neccesary=.true.
 endif
 enddo
 enddo
 if(count.gt.1.d6) stop "Too many iterations in back-to-cell"
 enddo
 end subroutine


!************************************************************************************

 subroutine cross_product_gensymcrys(a,b,crossp)
 !a very simple implementation of the cross product
 implicit none
 real(8)::a(3),b(3)
 real(8)::crossp(3)
 crossp(1)=a(2)*b(3)-a(3)*b(2)
 crossp(2)=a(3)*b(1)-a(1)*b(3)
 crossp(3)=a(1)*b(2)-a(2)*b(1)
 return
 end subroutine

!************************************************************************************

 subroutine latvec2dproj_gensymcrys(dproj,latvec,rotmat,rxyz,nat)
 !This subroutine will convert the lattice vector representation of thei
 !periodic cell (vec1,vec2,vec3) into the projective representation (dxx,dyx,dyy,dzx,dzy,dzz)
 !The cell will thus be rotated. The rotational matrix is stored in rotmat as an operator rotmat
 !and the atomic position rxyz are transformed into the new coordination sizstem as well
 implicit none
 integer,intent(in)  :: nat
 integer             :: iat,i
 real*8,intent(inout):: dproj(6),latvec(3,3),rotmat(3,3),rxyz(3,nat)
 real*8  :: tempvec(3),rotmat1(3,3),rotmat2(3,3),crossp(3),alpha,latvect(3,3)
 real*8  :: eps,axe(3),norm,rxyzt(3),rotmatt(3,3)
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
 call cross_product_gensymcrys(latvec(:,1),tempvec,crossp)
 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. abs(crossp(3)).lt.eps*1.d-1) goto 1001 !no rotation needed
 norm=sqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 axe(:)=crossp(:)/norm
 alpha=dacos(dot_product(tempvec,latvec(:,1))/dproj(1))
 call rotation_gensymcrys(rotmat1,alpha,axe)
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
 call cross_product_gensymcrys(tempvec,(/0.d0,1.d0,0.d0/),axe)
 if (abs(axe(1)).lt.eps*1.d-1 .and. abs(axe(2)).lt.eps*1.d-1 .and. abs(axe(3)).lt.eps*1.d-1 .and. tempvec(2).gt.0.d0) goto 1002 !no rotation needed
 norm=sqrt(axe(1)*axe(1)+axe(2)*axe(2)+axe(3)*axe(3))
 axe=axe/norm
 call cross_product_gensymcrys(axe,latvec(:,2),crossp) 
 if (abs(crossp(1)).lt.eps*1.d-1 .and. abs(crossp(2)).lt.eps*1.d-1 .and. crossp(3).gt.0.d0) goto 1002 !no rotation needed
 norm=sqrt(tempvec(2)*tempvec(2)+tempvec(3)*tempvec(3))
 alpha=dot_product((/0.d0,1.d0,0.d0/),tempvec(:))/norm
 alpha=dacos(alpha)
 call rotation_gensymcrys(rotmat2,alpha,axe)
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

!************************************************************************************

 subroutine rotation_gensymcrys(rotmat,angle,axe)
 !This subroutine will calculate the rotational matrix rotmat for a
 !3-dim vector around an axis 'axe' by the angle 'angle'.
 real(8), intent(in) :: angle
 real(8), intent(in) :: axe(3)
 real(8):: rotator(3,3)
 real(8):: rotmat(3,3)
 
 !Define Rotation Matrix
 rotator(1,1)=dcos(angle)+(axe(1)**2)*(1.d0-dcos(angle))
 rotator(1,2)=axe(1)*axe(2)*(1.d0-dcos(angle))-axe(3)*dsin(angle)
 rotator(1,3)=axe(1)*axe(3)*(1.d0-dcos(angle))+axe(2)*dsin(angle)
 
 rotator(2,1)=axe(2)*axe(1)*(1.d0-dcos(angle))+axe(3)*dsin(angle)
 rotator(2,2)=dcos(angle)+(axe(2)**2)*(1.d0-dcos(angle))
 rotator(2,3)=axe(2)*axe(3)*(1.d0-dcos(angle))-axe(1)*dsin(angle)
 
 rotator(3,1)=axe(3)*axe(1)*(1.d0-dcos(angle))-axe(2)*dsin(angle)
 rotator(3,2)=axe(3)*axe(2)*(1.d0-dcos(angle))+axe(1)*dsin(angle)
 rotator(3,3)=dcos(angle)+(axe(3)**2)*(1.d0-dcos(angle))
 rotmat(:,:)=rotator(:,:)
 
 !do i=1,3
 !   vector2(i)=rotator(i,1)*vector(1)+rotator(i,2)*vector(2)+rotator(i,3)*vector(3)
 !enddo
 !vector(:)=vector2(:)
 end subroutine rotation_gensymcrys

!************************************************************************************

 subroutine dist2plane_gensymcrys(point,nvec,ppoint,dist)
 !This subroutine will calculate  the distance between a plane and a point in space
 !The point is 'point', the normalized normal vector of the plane is 'nvec', 'ppoint' is an arbitrary point on the plane
 !and the output is the distance 'dist'  
 real(8), intent(in) :: point(3),nvec(3),ppoint(3)
 real(8), intent(out):: dist
 integer             :: i
 real(8)             :: p,nvectmp(3)
 nvectmp(:)=nvec(:)!/sqrt(nvec(1)*nvec(1)+nvec(2)*nvec(2)+nvec(3)*nvec(3))
 p=DOT_PRODUCT(nvectmp,ppoint) 
 p=-p
 dist=abs(DOT_PRODUCT(nvectmp,point)+p)
 end subroutine

!************************************************************************************

 subroutine nveclatvec_gensymcrys(latvec,nvec)
 !Will calculate the normalized normal vector to the 3 planes of the cell
 implicit none
 real*8, intent(in) :: latvec(3,3)
 real*8, intent(out):: nvec(3,3)
 real*8             :: a(3),b(3),crossp(3),norm
 integer:: i
 do i=1,3
 a=latvec(:,i)
 b=latvec(:,mod(i,3)+1)
 call cross_product_gensymcrys(a,b,crossp)
 norm=dsqrt(crossp(1)*crossp(1)+crossp(2)*crossp(2)+crossp(3)*crossp(3))
 nvec(:,i)=crossp(:)/norm
 enddo
 end

!************************************************************************************

 subroutine expand_gensymcrys(rxyz,rxyzout,transvecall,latvec,nat)
 !This subroutine will expand the unit cell into 26 periodic cells and store them in rxyzout
 implicit none
 real*8, intent(in)  :: rxyz(3,nat),latvec(3,3)
 integer, intent(in) :: nat
 real*8, intent(out) :: rxyzout(3,nat,3,3,3) !26 periodic images plus the main cell
 integer             :: iat,iplane,icorner,iedge,m,k,l
 real*8,intent(inout):: transvecall(3,3,3,3)!,(transvecp(3,6),transvecc(3,8),transvece(3,12)

 do m=-1,1
    do k=-1,1
       do l=-1,1
       transvecall(:,l+2,k+2,m+2)=real(l,8)*latvec(:,1)+real(k,8)*latvec(:,2)+real(m,8)*latvec(:,3)
       enddo
    enddo
 enddo

 do m=1,3
    do k=1,3
       do l=1,3
       do iat=1,nat
       rxyzout(:,iat,l,k,m)=rxyz(:,iat)+transvecall(:,l,k,m)
       enddo
       enddo
    enddo
 enddo
 end

!************************************************************************************
 
 subroutine invertmat_gensymcrys(mat,matinv,n)
 implicit none
 real(8),intent(in) :: mat(n,n)
 integer               :: n
 real(8),allocatable   :: WORK(:)
 real(8)               :: matinv(n,n),det(3),a(n,n),div
 integer               :: IPIV(n), INFO 
 integer               :: LDWORK
 !Here only for a 3*3 matrix
! a=mat
! div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)) 
! div=1.d0/div
!      matinv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
!      matinv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
!      matinv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
!      matinv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
!      matinv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
!      matinv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
!      matinv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
!      matinv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
!      matinv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
 !General n*n matrix 
 matinv=mat 
 allocate(WORK(n))
 call  DGETRF( n, n, matinv, n, IPIV, INFO )
 if (info.ne.0) stop "Error in DGETRF"
 LDWORK=-1
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 LDWORK=WORK(1)
 deallocate(WORK)
 allocate(WORK(LDWORK))
 call  DGETRI( n, matinv, n, IPIV, WORK,LDWORK , INFO )
 if (info.ne.0) stop "Error in DGETRI"
 end subroutine 

!************************************************************************************

 subroutine latvec2metric(latvec,metric)
 !This subroutine will calculate the scalar products between every lattice vector and store them
 !in the symmetric matrix metric. Diagonal elements represent the length of the lattice vectors,
 !All others represent the angles
 implicit none
 real(8),intent(in)    :: latvec(3,3)
 real(8),intent(out)   :: metric(3,3)
 metric(1,1)=dot_product(latvec(:,1),latvec(:,1))
 metric(2,2)=dot_product(latvec(:,2),latvec(:,2))
 metric(3,3)=dot_product(latvec(:,3),latvec(:,3))
 metric(1,2)=dot_product(latvec(:,1),latvec(:,2))
 metric(1,3)=dot_product(latvec(:,1),latvec(:,3))
 metric(2,3)=dot_product(latvec(:,2),latvec(:,3))
 metric(2,1)=metric(1,2)
 metric(3,1)=metric(1,3)
 metric(3,2)=metric(2,3)
 end subroutine

!************************************************************************************

 subroutine checkmetric(metric1,metric2,alpha,info)
 !This subroutine will compare the metric matrices metric1 and metric2 and check if the change in any part
 !of the metrices exceeds more than alpha (for example alpha=10%). The difference is calculated as follows:
 !err(i,j)=abs(metric1(i,j)-metric2(i,j))/Vol(metric1)
 !info=0 if err(ij).le.alpha for all ij
 !info=1 if err(ij).gt.alpha for i=j
 !info=2 if err(ij).gt.alpha for i.ne.j
 !info=3 if err(ij).gt.alpha for all ij
 !info=info+5 if (vol1-vol2)*volinv.gt.alpha
 implicit none
 real(8), intent(in):: metric1(3,3),metric2(3,3),alpha
 integer            :: info,i,j,infoii,infoij
 real(8)            :: err(3,3), det, volinv,vol1,vol2, a(3,3)
 info=0
 infoii=0
 infoij=0
 a=metric1
 det=a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
     a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
 vol1=sqrt(det)
 a=metric2
 det=a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
     a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
 vol2=sqrt(det)

 volinv=1.d0/vol1

 do i=1,3
   do j=1,3
   err(i,j)=abs(metric1(i,j)-metric2(i,j))*volinv
   if(i==j.and.infoii==0) then
     if(err(i,j).gt.alpha) infoii=1
   elseif(i.ne.j.and.infoij==0) then
     if(err(i,j).gt.alpha) infoij=2
   endif
   if(infoij.ne.0.and.infoii.ne.0) goto 1001
   enddo
 enddo

 1001 continue
 info=infoij+infoii
 if((vol1-vol2)*volinv.gt.alpha) info=info+5
 end subroutine

!************************************************************************************

 subroutine rescale_cellforces(stress,latvec,alphapar,alphashear)
 !This subroutine will split the forces on the cell into a part parallel to the cell vectors and scale this part
 !with alphapar and a part in within the plane spanned by the two other cell vectors scaled with alphashear. 
 !The output will be written back to "stress"
 implicit none
 real(8), intent(inout):: stress(3,3)
 real(8), intent(in)   :: latvec(3,3),alphapar,alphashear
 !Has to be implemented first
 end subroutine

!************************************************************************************

 subroutine updaterxyz_gensymcrys(latvecold,latvecnew,rxyz,nat)
 !This subroutine will update the atomic positions in the cartesian coordinates after the cell shape has been changed according
 !to the change in the lattice vectors thus keeping the relative coordinates of all atoms constant
 implicit none
 real(8), intent(in)   :: latvecold(3,3), latvecnew(3,3)
 real(8), intent(inout):: rxyz(3,nat)
 integer, intent(in)   :: nat
 real(8)               :: latvecold_inv(3,3),trafo(3,3)
 integer               :: iat 
 call invertmat_gensymcrys(latvecold,latvecold_inv,3)
 trafo=matmul(latvecnew,latvecold_inv) !Transform atomic positions according to cell
 do iat=1,nat
 rxyz(:,iat)=matmul(trafo,rxyz(:,iat))
 enddo
 end subroutine
 
!************************************************************************************

 subroutine r_symm(rxyz,latvec,nat)
 !This subroutine will symmetrize the lattice vectors. It will also transform the atomic position into the rotated system
 implicit none
 real(8), intent(inout):: rxyz(3,nat), latvec(3,3)
 integer, intent(in)   :: nat
 real(8)               :: metric(3,3),evals(3),dmat(3,3),metrict(3,3),latvecnew(3,3)
 real(8),allocatable   :: work(:)
 integer               :: lwork, info
 !Transform cell to latve*latvec_transposed=metric
 call latvec2metric(latvec,metric)
 !Calculate eigenvalues and eigenvectors of s2
      lwork=100
      allocate(work(lwork))
      call DSYEV( 'V', 'L', 3, metric, 3, evals, WORK, LWORK, INFO )
      if(info.ne.0) stop 'Error in DSYEV'
 !Take sqrt of the diagonal elements
      evals(1)=sqrt(evals(1));evals(2)=sqrt(evals(2));evals(3)=sqrt(evals(3))
      dmat=0.d0
 !Generate diagonal matrix
      dmat(1,1)=evals(1)
      dmat(2,2)=evals(2)
      dmat(3,3)=evals(3)
 !Transform back by multiplying the diagonal matrix with the unitary matrices and its transposed from lef and right
      metrict(:,1)=metric(1,:)
      metrict(:,2)=metric(2,:)
      metrict(:,3)=metric(3,:)
      dmat=matmul(metric,dmat)
 !Update lattice vector and the positions     
      latvecnew=matmul(dmat,metrict)
      call updaterxyz_gensymcrys(latvec,latvecnew,rxyz,nat)
      latvec=latvecnew
 end subroutine

!************************************************************************************

 subroutine sigma_sigmat(latvec,sigma,sigmat)
 !This subroutine will calulate the sigma tensor which consists of the three unnormalized normal vectors sigma=[v2xv3,v3xv1,v1xv2]
 !sigmat=sigma(transposed)*sigma
 !This code is useful for the Wentzcovitch MD
 implicit none
 real(8), intent(in) :: latvec(3,3)
 real(8), intent(out):: sigma(3,3),sigmat(3,3)
 real(8)             :: a(3),b(3),crossp(3)
 integer:: i
 do i=1,3
 a=latvec(:,i)
 b=latvec(:,mod(i,3)+1)
 call cross_product_gensymcrys(a,b,crossp)
 sigma(:,i)=crossp(:)
 enddo
 sigmat(:,1)=sigma(1,:)
 sigmat(:,2)=sigma(2,:)
 sigmat(:,3)=sigma(3,:)
 end subroutine
 
!************************************************************************************
 
 subroutine rxyz_int2cart_gensymcrys(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3)
 integer:: nat,iat
 do iat=1,nat
  rxyzcart(:,iat)=matmul(latvec,rxyzint(:,iat))
 enddo
 end subroutine rxyz_int2cart_gensymcrys

!************************************************************************************

 subroutine rxyz_cart2int_gensymcrys(latvec,rxyzint,rxyzcart,nat)
 !This subrouine will convert the internal coordinates into cartesian coordinates
 implicit none
 real(8):: rxyzint(3,nat), rxyzcart(3,nat),latvec(3,3),latvecinv(3,3)
 integer:: nat,iat
 call invertmat_gensymcrys(latvec,latvecinv,3)
 do iat=1,nat
  rxyzint(:,iat)=matmul(latvecinv,rxyzcart(:,iat))
 enddo
 end subroutine rxyz_cart2int_gensymcrys

!************************************************************************************

subroutine getvol_gensymcrys(cellvec,vol)
    implicit none
    real(8), intent(in):: cellvec(3,3)
    real(8), intent(out):: vol
    vol=cellvec(1,1)*cellvec(2,2)*cellvec(3,3)-cellvec(1,1)*cellvec(2,3)*cellvec(3,2)- &
        cellvec(1,2)*cellvec(2,1)*cellvec(3,3)+cellvec(1,2)*cellvec(2,3)*cellvec(3,1)+ &
        cellvec(1,3)*cellvec(2,1)*cellvec(3,2)-cellvec(1,3)*cellvec(2,2)*cellvec(3,1)
    if(vol<0.d0) stop 'ERROR: Negative volume!'
end subroutine getvol_gensymcrys


