 subroutine correct_latvec_oganov(latvec,pos_red,nat,iproc)
 !use cell_utils
 !This subroutine will use the algorithm proposed by oganov and glass (J.Phys,Cond.Mat 20,2008) to perform a transformation of the lattice vectors into an equivalent
 !system where the length of all cell vectors are similar (no nasty angles).
 implicit none
 real(8)              :: latvec(3,3),rxyz(3,nat),pos_red(3,nat)  
 logical              :: correct 
 real(8)              :: val_inter,norm_half,vnrm !The val_... are the absolute value of the projections
 real(8)              :: a(3),b(3),c(3) ! the three latticevectors
 real(8)              :: tempvec(3),v(3,3),vol,sign_inter
 integer              :: i,nat,counter,iproc

 call backtocell(nat,latvec,pos_red)
 call rxyz_int2cart(latvec,pos_red,rxyz,nat)
 a=latvec(:,1);b=latvec(:,2);c=latvec(:,3)
 !check volume
 v=latvec
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! write(*,*) "initial vol",vol
 !Check the 6 criteria ab,ba,ac,ca,bc,cb
 correct=.true.
 counter=0
 do while(correct) 
 counter=counter+1
 correct=.false.
 !ab
 val_inter=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 sign_inter=sign(1.d0,val_inter)
 vnrm=b(1)*b(1)+b(2)*b(2)+b(3)*b(3) 
 vnrm=dsqrt(vnrm)
 norm_half=vnrm*0.5d0
 val_inter=abs(val_inter/vnrm)
 if(val_inter.gt.norm_half) then
! write(*,*) "Transform ab"
 correct=.true.
 val_inter=val_inter/vnrm
 a(:)=a(:)-ceiling(val_inter)*sign_inter*b(:)
 latvec(:,1)=a(:)
 endif


 !ba
 val_inter=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
 sign_inter=sign(1.d0,val_inter)
 vnrm=a(1)*a(1)+a(2)*a(2)+a(3)*a(3) 
 vnrm=dsqrt(vnrm)
 norm_half=vnrm*0.5d0
 val_inter=abs(val_inter/vnrm)
 if(val_inter.gt.norm_half) then
! write(*,*) "Transform ba"
 correct=.true.
 val_inter=val_inter/vnrm
 b(:)=b(:)-ceiling(val_inter)*sign_inter*a(:)
 latvec(:,2)=b(:)
 endif


 !ac
 val_inter=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 vnrm=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)
 vnrm=dsqrt(vnrm)
 norm_half=vnrm*0.5d0
 val_inter=abs(val_inter/vnrm)
 if(val_inter.gt.norm_half) then
! write(*,*) "Transform ac"
 correct=.true.
 val_inter=val_inter/vnrm
 a(:)=a(:)-ceiling(val_inter)*sign_inter*c(:)
 latvec(:,1)=a(:)
 endif


 !ca
 val_inter=a(1)*c(1)+a(2)*c(2)+a(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 vnrm=a(1)*a(1)+a(2)*a(2)+a(3)*a(3)
 vnrm=dsqrt(vnrm)
 norm_half=vnrm*0.5d0
 val_inter=abs(val_inter/vnrm)
 if(val_inter.gt.norm_half)then
! write(*,*) "Transform ca"
 correct=.true.
 val_inter=val_inter/vnrm
 c(:)=c(:)-ceiling(val_inter)*sign_inter*a(:)
 latvec(:,3)=c(:)
 endif


 !bc
 val_inter=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 vnrm=c(1)*c(1)+c(2)*c(2)+c(3)*c(3)
 vnrm=dsqrt(vnrm)
 norm_half=vnrm*0.5d0
 val_inter=abs(val_inter/vnrm)
 if(val_inter.gt.norm_half)then
! write(*,*) "Transform bc"
 correct=.true.
 val_inter=val_inter/vnrm
 b(:)=b(:)-ceiling(val_inter)*sign_inter*c(:)
 latvec(:,2)=b(:)
 endif


 !cb
 val_inter=b(1)*c(1)+b(2)*c(2)+b(3)*c(3)
 sign_inter=sign(1.d0,val_inter)
 vnrm=b(1)*b(1)+b(2)*b(2)+b(3)*b(3)
 vnrm=dsqrt(vnrm)
 norm_half=vnrm*0.5d0
 val_inter=abs(val_inter/vnrm)
 if(val_inter.gt.norm_half)then
! write(*,*) "Transform cb"
 correct=.true.
 val_inter=val_inter/vnrm
 c(:)=c(:)-ceiling(val_inter)*sign_inter*b(:)
 latvec(:,3)=c(:)
 endif
 
 enddo

 if(iproc==0) write(*,*)"# Number of correct cycles" , counter
 v=latvec
 vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! write(100,*) "final vol",vol
 call backtocell_cart(nat,latvec,rxyz)
 call rxyz_cart2int(latvec,pos_red,rxyz,nat)
 end subroutine
