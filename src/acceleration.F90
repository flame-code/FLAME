!**********************************************************************************************

subroutine acceleration(pressure,accpos,acclat,accvol,vpos,vlat,vvol,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
implicit none
integer:: iat,i,j,md_type,nat
real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
real(8):: accvol,vvol,vol_1_3
!Get volume
           a=latvec
           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!md_type=1
!Convert cartesian forces to reduced forces on atoms
        call invertmat(latvec,tmplat,3)
        do iat=1,nat
          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
        enddo
!Get full stress matrix
           str_matrix(1,1)=strten(1)
           str_matrix(2,2)=strten(2)
           str_matrix(3,3)=strten(3)
           str_matrix(1,2)=strten(6)
           str_matrix(2,1)=strten(6)
           str_matrix(1,3)=strten(5)
           str_matrix(3,1)=strten(5)
           str_matrix(2,3)=strten(4)
           str_matrix(3,2)=strten(4)
!Update velocity part of stress tensor
           velmat=0.d0
           do iat=1,nat
              vpostmp=matmul(latvec,vpos(:,iat))
              do i=1,3
                 do j=1,3
                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
                 enddo
              enddo
           enddo
if(md_type.ge.1.and.md_type.le.3) then
!Compute sigma
        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!Compute the atomic acceleration
!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
          do i=1,3
            lattrans(:,i)=latvec(i,:)
            latdottrans(:,i)=vlat(i,:)
          enddo
          gdot=matmul(latdottrans,latvec)+matmul(lattrans,vlat)
          g=matmul(lattrans,latvec)
          call invertmat(g,ginv,3)
          gtot=matmul(ginv,gdot)
!Total acceleration
          do iat=1,nat
          accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat) - matmul(gtot,vpos(:,iat))
!          if(fixlat(7))  accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)
          enddo
elseif(md_type==4) then
!Andersen
!Be careful: the positions are not the real positions but scaled with vol**(1/3)
          vol_1_3=vol**(1.d0/3.d0)
          do iat=1,nat
              accpos(:,iat)=(fcart(:,iat)/amass(iat)-2.d0*vpos(:,iat)*vvol/(vol*3.d0))/vol_1_3
          enddo
!Update velocity part of stress tensor
           velmat=0.d0
           do iat=1,nat
              vpostmp=vol_1_3*vpos(:,iat)
              do i=1,3
                 do j=1,3
                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)*amass(iat) !Added mass on Oct 23 2015
                 enddo
              enddo
           enddo
endif

if(md_type==1) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)

elseif(md_type==2) then
!Fist term, same as in PR******************************
!Compute the acceleration of the cell
        term1=velmat/vol-str_matrix
!Here the pressure is applied
        term1=term1-pressure
!Scale with lattice volume and mass
        term1=term1/vol/latmass
!Combine it with the cell
        term1=matmul(term1,latvec)

!Second term*********************************************
        sigmatrans(:,1)=sigma(1,:)
        sigmatrans(:,2)=sigma(2,:)
        sigmatrans(:,3)=sigma(3,:)
        term2=matmul(sigmatrans,vlat)
        volvel=term2(1,1)+term2(2,2)+term2(3,3)
        term2=-2.d0*volvel/vol*vlat

!Third term**********************************************
!This term was taken from the moldyn manual. It is slightly different in the cleveland paper!!!
!        term3=matmul(sigmatrans,sigma)
!        term3=matmul(term3,latdottrans)
!        term3=matmul(term3,vlat)
!From cleveland
        term3=matmul(vlat,sigmatrans)
        term3=matmul(term3,sigma)
        term3=matmul(term3,latdottrans)
        trace3=term3(1,1)+term3(2,2)+term3(3,3)
        term3=1.d0/vol**2*trace3*latvec 

!Fourth term*********************************************
        term4=matmul(vlat,sigmatrans)
        term4=matmul(term4,vlat)
        term4=1.d0/vol*term4
 
!Fifth term**********************************************
        term5_1=matmul(vlat,sigmatrans)
        term5_1=matmul(term5_1,sigma)
        term5_1=matmul(term5_1,latdottrans)
        term5_2=matmul(sigma,latdottrans)
        term5_2=matmul(term5_2,vlat)
        term5_2=matmul(term5_2,sigmatrans)
        term5=1.d0/vol**2*matmul((term5_1-term5_2),latvec)
!Sum
        acclat=term1+term2+term3+term4+term5
elseif(md_type==3) then
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass*vol**(4.d0/3.d0)
!Multiply with sigma from left
        acclat=matmul(acclat,sigma)
        acclat=matmul(acclat,f0inv)
elseif(md_type==4) then
!Andersen MD
!Compute the acceleration of the cell
        acclat=velmat/vol-str_matrix
!Here the pressure is applied
        acclat=acclat-pressure
!Scale with lattice mass
        acclat=acclat/latmass
!Compute the hydrostatic pressure stuff
        accvol=(acclat(1,1)+acclat(2,2)+acclat(3,3))/3.d0
else
stop "Wrong option in MD"
endif
end subroutine

!**********************************************************************************************

!subroutine acceleration(pressure,accpos,acclat,vpos,vlat,strten,fcart,latvec,amass,latmass,f0inv,md_type,nat) 
!use global, only: fixlat
!implicit none
!integer:: iat,i,j,md_type
!real(8),dimension(3,nat):: accpos,vpos,fcart,fpos
!real(8),dimension(3,3)  :: acclat,vlat,latvec,tmplat,pressure,a,velmat,sigma,lattrans,latdottrans,gdot,g,ginv,gtot,str_matrix
!real(8),dimension(3,3)  :: term1,term2,term3,term4,term5,term5_1,term5_2,sigmatrans,f0inv
!real(8):: amass(nat),latmass,crossp(3),strten(6),vol,vpostmp(3),volvel,trace3
!!md_type=1
!!Convert cartesian forces to reduced forces on atoms
!        call invertmat(latvec,tmplat,3)
!        do iat=1,nat
!          fpos(:,iat)=matmul(tmplat,fcart(:,iat))
!        enddo
!!Get full stress matrix
!           str_matrix(1,1)=strten(1)
!           str_matrix(2,2)=strten(2)
!           str_matrix(3,3)=strten(3)
!           str_matrix(1,2)=strten(6)
!           str_matrix(2,1)=strten(6)
!           str_matrix(1,3)=strten(5)
!           str_matrix(3,1)=strten(5)
!           str_matrix(2,3)=strten(4)
!           str_matrix(3,2)=strten(4)
!!Get volume
!           a=latvec
!           vol= a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+&
!                a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1)
!!Update velocity part of stress tensor
!           velmat=0.d0
!           do iat=1,nat
!              vpostmp=matmul(latvec,vpos(:,iat))
!              do i=1,3
!                 do j=1,3
!                 velmat(i,j)=velmat(i,j)+vpostmp(i)*vpostmp(j)
!                 enddo
!              enddo
!           enddo
!!Compute sigma
!        call cross_product(latvec(:,2),latvec(:,3),crossp); sigma(:,1)=crossp
!        call cross_product(latvec(:,3),latvec(:,1),crossp); sigma(:,2)=crossp
!        call cross_product(latvec(:,1),latvec(:,2),crossp); sigma(:,3)=crossp
!!Compute the atomic acceleration
!!For this we first need to calculate the matrix g^-1*dg/dt=ginv*gdot=gtot
!          do i=1,3
!            lattrans(:,i)=latvec(i,:)
!            latdottrans(:,i)=vlat(i,:)
!          enddo
!          gdot=matmul(latdottrans,latvec)+matmul(lattrans,vlat)
!          g=matmul(lattrans,latvec)
!          call invertmat(g,ginv,3)
!          gtot=matmul(ginv,gdot)
!!Total acceleration
!          do iat=1,nat
!          accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat) - matmul(gtot,vpos(:,iat))
!!          if(fixlat(7))  accpos(:,iat)=1.d0/amass(iat)*fpos(:,iat)
!          enddo
!if(md_type==1) then
!!Compute the acceleration of the cell
!        acclat=velmat/vol-str_matrix
!!Here the pressure is applied
!        acclat=acclat-pressure
!!Scale with lattice mass
!        acclat=acclat/latmass
!!Multiply with sigma from left
!        acclat=matmul(acclat,sigma)
!
!elseif(md_type==2) then
!!Fist term, same as in PR******************************
!!Compute the acceleration of the cell
!        term1=velmat/vol-str_matrix
!!Here the pressure is applied
!        term1=term1-pressure
!!Scale with lattice volume and mass
!        term1=term1/vol/latmass
!!Combine it with the cell
!        term1=matmul(term1,latvec)
!
!!Second term*********************************************
!        sigmatrans(:,1)=sigma(1,:)
!        sigmatrans(:,2)=sigma(2,:)
!        sigmatrans(:,3)=sigma(3,:)
!        term2=matmul(sigmatrans,vlat)
!        volvel=term2(1,1)+term2(2,2)+term2(3,3)
!        term2=-2.d0*volvel/vol*vlat
!
!!Third term**********************************************
!!This term was taken from the moldyn manual. It is slightly different in the cleveland paper!!!
!!        term3=matmul(sigmatrans,sigma)
!!        term3=matmul(term3,latdottrans)
!!        term3=matmul(term3,vlat)
!!From cleveland
!        term3=matmul(vlat,sigmatrans)
!        term3=matmul(term3,sigma)
!        term3=matmul(term3,latdottrans)
!        trace3=term3(1,1)+term3(2,2)+term3(3,3)
!        term3=1.d0/vol**2*trace3*latvec 
!
!!Fourth term*********************************************
!        term4=matmul(vlat,sigmatrans)
!        term4=matmul(term4,vlat)
!        term4=1.d0/vol*term4
! 
!!Fifth term**********************************************
!        term5_1=matmul(vlat,sigmatrans)
!        term5_1=matmul(term5_1,sigma)
!        term5_1=matmul(term5_1,latdottrans)
!        term5_2=matmul(sigma,latdottrans)
!        term5_2=matmul(term5_2,vlat)
!        term5_2=matmul(term5_2,sigmatrans)
!        term5=1.d0/vol**2*matmul((term5_1-term5_2),latvec)
!!Sum
!        acclat=term1+term2+term3+term4+term5
!elseif(md_type==3) then
!!Compute the acceleration of the cell
!        acclat=velmat/vol-str_matrix
!!Here the pressure is applied
!        acclat=acclat-pressure
!!Scale with lattice mass
!        acclat=acclat/latmass*vol**(4.d0/3.d0)
!!Multiply with sigma from left
!        acclat=matmul(acclat,sigma)
!        acclat=matmul(acclat,f0inv)
!else
!stop "Wrong option in MD"
!endif
!end subroutine
