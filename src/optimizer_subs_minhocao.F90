!******************************************
!                                         *
!GEOMETRY OPTIMIZERS                      *                   
!                                         *
!******************************************

module interfaces
!Module used for the geometry optimizer BFGS with ABINIT LINESEARCH
interface
   subroutine unit_matrix(mat)
     implicit none
     real(8),DIMENSION(:,:), INTENT(INOUT) :: mat
   end subroutine unit_matrix

   function vabs(v) result(res)
     implicit none
     real(8),dimension(:):: v
     real(8):: res
   end function vabs

   function outerprod(a,b)
     real(8),dimension(:),intent(in)::a,b
     real(8),dimension(size(a),size(b))::outerprod
   end function outerprod
end interface
end module interfaces
!>   Geometry optimization, parametrisation routine.
subroutine geopt_init()
  use minpar
  implicit none

  parmin_bfgs%approach  = 'unknown'
  parmin_bfgs%iter      = 0
  parmin_bfgs%iflag     = 0
  parmin_bfgs%verbosity = 1
  parmin_bfgs%MSAVE=7
  parmin_bfgs%MP=16
  parmin_bfgs%LP=16
  parmin_bfgs%MAXFEV=10
  parmin_bfgs%GTOL=9.d-1
  parmin_bfgs%XTOL=1.d-15
  parmin_bfgs%FTOL=1.d-6
  parmin_bfgs%STPMIN=1.d-20
  parmin_bfgs%STPMAX=20.d0
  parmin_bfgs%DIAGCO=.FALSE.
  parmin_bfgs%IWRITE=.FALSE.

END SUBROUTINE geopt_init


!************************************************************************************

subroutine GEOPT_RBFGS_MHM(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter)
!subroutine bfgs_driver_atoms(latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,counter,fmax_tol)
 use global, only: ntypat,znucl,typat,char_type,units
 use defs_basis
 use minpar

!SUBROUTINE dfpmin_pos(nat,latvec,rxyz,fxyz,stress,pressure,etot,fnrmtol,iter,count)
use interfaces
use mod_parini, only: typ_parini
IMPLICIT NONE
type(typ_parini), intent(in):: parini
!REAL(8), INTENT(IN) :: fnrmtol!gtol 
REAL(8) :: fret, counter
REAL(8), INTENT(INOUT) :: xred_in(3*parini%nat),latvec_in(9),fcart_in(3*parini%nat),strten_in(6),etot_in
INTEGER, PARAMETER :: ITMAX=4000
REAL(8), PARAMETER :: STPMX=1.0d0,EPS=epsilon(xred_in),TOLX=4.0d0*EPS
!   Given a starting point p that is a vector of length N , the Broyden-Fletcher-Goldfarb-Shanno
!   variant of Davidon-Fletcher-Powell minimization is performed on a function func, using its
!   gradient as calculated by a routine dfunc. The convergence requirement on zeroing the
!   gradient is input as gtol. Returned quantities are p (the location of the minimum), iter
!   (the number of iterations that were performed), and fret (the minimum value of the
!   function). The routine lnsrch is called to perform approximate line minimizations.
!   Parameters: ITMAX is the maximum allowed number of iterations; STPMX is the scaled
!   maximum step length allowed in line searches; EPS is the machine precision; TOLX is the
!   convergence criterion on x values.
INTEGER :: its,assert_eq,i
LOGICAL :: check
REAL(8) :: den,fac,fad,fae,fp,stpmax,sumdg,sumxi
REAL(8):: dg(3*parini%nat+9),g(3*parini%nat+9),hdg(3*parini%nat+9),pnew(3*parini%nat+9),xi(3*parini%nat+9),p(3*parini%nat+9)
REAL(8):: tp(3*parini%nat+9),tg(3*parini%nat+9),dvin(3*parini%nat+9),vout(3*parini%nat+9),vout_prev(3*parini%nat+9)
REAL(8):: vin_min(3*parini%nat+9),vin(3*parini%nat+9)
REAL(8):: vout_min(3*parini%nat+9),dedv_min(3*parini%nat+9)
REAL(8), DIMENSION(3*parini%nat+9,3*parini%nat+9) :: hessin,hessin0
REAL(8) :: alpha_pl
REAL(8) :: gamma0,gammax,lambda_1,lambda_2,tfp,dedv_1,dedv_2,etotal_1,etotal_2,dedv_predict
REAL(8) :: d2edv2_1,d2edv2_2,d2edv2_predict,etotal_predict,lambda_predict
INTEGER :: choice,status,sumstatus,iprec,iexit
REAL(8) :: latvec0(9),rxyz0(3*parini%nat),eval(3*parini%nat+9),fmax,fmax_at,fmax_lat,pressure
LOGICAL :: getwfk
REAL(8) :: ent_pos_0,enthalpy,en0000
character(40)::filename
character(4) ::fn4
logical:: multiprec
real(8):: tolmxf_switch 
!multiprec is hardcoded and, if true, starts a geopt with iprec==2, and then switches 
!to iprec==1 when the fmax==tolmxf_switch. The switch only occurs once
 multiprec=.true.
 tolmxf_switch=10.d0*parini%paropt_geopt%fmaxtol

 counter=0.d0
write(*,'(a,es15.7,es15.7)') " # BFGS BETAX, BETAX_LAT: ", parmin_bfgs%betax, parmin_bfgs%betax_lat
 
pressure=parini%target_pressure_habohr

open(unit=16,file="geopt.mon")
alpha_pl=1.d-0
gammax=1.d0
fret=0.d0
check=.true.
!call rxyz_cart2int(latvec,p(1:3*nat),rxyz,nat)
p(1:3*parini%nat)=xred_in(:)
p(3*parini%nat+1:3*parini%nat+9)=latvec_in(:)
!Calculate starting function value and gradient.
getwfk=.false.
iprec=1
call get_BFGS_forces_max(parini,p,g,fp,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
!INITIAL STEP, STILL THE SAME STRUCTURE AS INPUT
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)') 0
       filename="posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS AC ",0,fp,fmax,fmax_lat,fmax_at,0.d0,iprec
!*********************************************************************
   iexit=0
   if(fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
   if(iexit==1) then
   write(*,'(a)') " # BFGS converged before entering optimization"
!   call wtpos_inter(nat,rxyz,latvec,555)
   RETURN 
   endif

!Initial iprec after running the first force call
 if(multiprec) iprec=2

!call energyandforces(nat,latvec,rxyz,fxyz,stress,pressure,fp,count)
write(16,*) "Initial energy",fp
!call wtpos_inter(nat,rxyz,latvec,500)
!call fxyz_cart2int(nat,fxyz,g(1:3*nat),latvec)
g(3*parini%nat+1:3*parini%nat+9)=g(3*parini%nat+1:3*parini%nat+9)*alpha_pl
g=-g
call unit_matrix(hessin,3*parini%nat+9) !Initialize inverse Hessian to the unit matrix.

!Initialize Hessian diagonal elements
hessin=hessin*parmin_bfgs%betax
!Initialize Hessian diagonal elements for the cell variables
do i=3*parini%nat+1,3*parini%nat+9
hessin(i,i)=parmin_bfgs%betax_lat
enddo

xi=-matmul(hessin,g)
!Main loop over the iterations.
! call wtpos_inter(nat,rxyz,latvec,500)
 sumstatus=0

hessin0=hessin
do its=1,ITMAX

 gamma0=1.d0*gammax
1001 continue
 do i=1,3*parini%nat+9
    tp(i)=p(i)+gamma0*xi(i)
 enddo
 write(16,*) "Length of movement along xi", vabs(gamma0*xi)
 if (vabs(gamma0*xi).gt.7.d-1) then
 gamma0=gamma0*0.5d0
 goto 1001
 endif

 if(parini%usewf_geopt) then
     getwfk=.true.
 else
     getwfk=.false.
 endif
 if(iprec==2.and.multiprec.and.fmax.lt.tolmxf_switch) then
     getwfk=.false.
     iprec=1
 endif
 counter=counter+1.d0
 call get_BFGS_forces_max(parini,tp,tg,tfp,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
 call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)') (its)*2-1 
       filename="posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call fcart_in,strten_in,&
            &write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS LS ",(its)*2-1,fp,fmax,fmax_lat,fmax_at,gamma0,iprec
!*********************************************************************
! call rxyz_int2cart(tp(3*nat+1:3*nat+9),tp(1:3*nat),rxyz,nat)
! call energyandforces(nat,tp(3*nat+1:3*nat+9),rxyz,fxyz,stress,pressure,tfp,count)
! call wtpos(nat,tp(3*nat+1:3*nat+9),rxyz,real(its,8))
! call fxyz_cart2int(nat,fxyz,tg,tp(3*nat+1:3*nat+9))
! tg(3*nat+1:3*nat+9)=stress(:)*alpha_pl


 tg=-tg 
 if (fmax.lt.1.d-3.or.sumstatus.gt.10000)then
 choice=1
 else
 choice=4
 endif 
 1000 continue
 lambda_1=1.0d0       ; lambda_2=0.0d0
 etotal_1=tfp      ; etotal_2=fp
! dvin(:)=vin(:)-vin_prev(:)
 dvin(:)=tp(:)-p(:)
 vout=tg
 vout_prev=g
! vout_prev=xi
 dedv_1=dot_product(vout,dvin)
 dedv_2=dot_product(vout_prev,dvin)
 call findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)

if(status==2) sumstatus=sumstatus+status
if(status==3) then
choice=1
call findmin(choice,dedv_1,dedv_2,dedv_predict,&
& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)
endif

!NEW: limit the step to go to 1.d0
lambda_predict=min(lambda_predict,1.5d0)
lambda_predict=max(lambda_predict,-1.d0)


!Generates vin at the minimum, and an interpolated vout, modified
!to have the right value of dedv_predict, from findmin.
! vin_min(:)=vin_prev(:)+lambda_predict*dvin(:)
 vin_min(:)=p(:)+lambda_predict*dvin(:)
 vout_min(:)=vout_prev(:)+lambda_predict*(vout(:)-vout_prev(:))
 dedv_min=dedv_2+lambda_predict*(dedv_1-dedv_2)
!Modify vout_min in order to impose dedv_predict
 vout_min(:)=vout_min(:)+dvin(:)*(dedv_predict-dedv_min)/dot_product(dvin,dvin)
 vin(:)=vin_min(:)
 pnew(:)=vin(:)

!   Update the line direction, and the current point.
   xi=pnew-p
   p=pnew
   dg=g       !Save the old gradient,

   if(parini%usewf_geopt) then
       getwfk=.true.
   else
       getwfk=.false.
   endif
   counter=counter+1.d0
   call get_BFGS_forces_max(parini,p,g,fp,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
   call get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
!MHM: Write output to file in every step***********************************
       write(*,*) "Pressure, Energy",pressure,etot_in
       ent_pos_0=fp
       en0000=fp-ent_pos_0
       write(fn4,'(i4.4)')  its*2
       filename="posgeopt."//fn4//".ascii"
       units=units
       write(*,*) "# Writing the positions in :",filename
       call write_atomic_file_ascii(parini,filename,parini%nat,units,xred_in,latvec_in,fcart_in,strten_in,&
            &char_type(1:ntypat),ntypat,typat,etot_in,pressure,fp,en0000)
       write(*,'(a,i4,4(1x,es17.8),1x,es9.2,1x,i4)') " # GEOPT BFGS AC ",its*2,fp,fmax,fmax_lat,fmax_at,lambda_predict,iprec
!*********************************************************************
!   call rxyz_int2cart(p(3*nat+1:3*nat+9),p(1:3*nat),rxyz,nat)
!   call energyandforces(nat,p(3*nat+1:3*nat+9),rxyz,fxyz,stress,pressure,etot,count)
!   call fxyz_cart2int(nat,fxyz,g,p(3*nat+1:3*nat+9))
!   g(3*nat+1:3*nat+9)=stress(:)*alpha_pl

!   fp=etot  
   g=-g
   den=max(fret,1.0d0)
   write(16,'(a,1x,I5,1x,1pe21.14,1x,1pe12.5,1x,1pe12.5,1x,1pe12.5)') "  &
   &  BFGS_all",its,fp,fmax,fmax_at,fmax_lat
   iexit=0
   if(fmax.lt.parini%paropt_geopt%fmaxtol) iexit=1
!   if (fnrm < fnrmtol) then  !Test for convergence on zero gradient.
!   latvec=p(3*nat+1:3*nat+9)
!   call backtocell(nat,latvec,rxyz)
   if(iexit==1) then
   write(*,'(a,i5,a)') " # BFGS converged in ",its*2," iterations"
!   call wtpos_inter(nat,rxyz,latvec,555)
   RETURN 
   endif
   if(int(counter).gt.parini%paropt_geopt%nit) then
   write(*,'(a,i5)') " # BFGS did not converg in steps: ", int(counter)
   RETURN
   endif

   dg=g-dg                !Compute difference of gradients,
   hdg=matmul(hessin,dg)  !and difference times current matrix.
   fac=dot_product(dg,xi) ! Calculate dot products for the denominators.
   fac=dot_product(dg,xi)
   fae=dot_product(dg,hdg)
   sumdg=dot_product(dg,dg)
   sumxi=dot_product(xi,xi)
    if (fac > sqrt(EPS*sumdg*sumxi)) then !Skip update if fac not sufficiently positive.
        fac=1.0d0/fac
        fad=1.0d0/fae
        dg=fac*xi-fad*hdg                 !Vector that makes BFGS different from DFP.
        hessin=hessin+fac*outerprod(xi,xi)
        hessin=hessin-fad*outerprod(hdg,hdg)
        hessin=hessin+fae*outerprod(dg,dg)
    end if
!Now calculate the next direction to go
    xi=-matmul(hessin,g)
!and go back for another iteration.
end do
write(*,*) "Too many iterations"

!stop "Too many iterations"
END SUBROUTINE

!************************************************************************************

real(8) FUNCTION vabs(v)
implicit none
real(8),dimension(:):: v
real(8)::sumv
integer::i
sumv=0.d0
do i=1,size(v)
sumv=sumv+v(i)*v(i)
enddo
vabs=sqrt(sumv)
!read(*,*) i
END FUNCTION vabs

!************************************************************************************

function outerprod(a,b)
implicit none
real(8),dimension(:),intent(in)::a,b
real(8),dimension(size(a),size(b))::outerprod
integer:: i,j
!outerprod=spread(a,dim=2,ncopies=size(b))*spread(b,dim=1,ncopies=size(a))
do i=1,size(a)
do j=1,size(b)
outerprod(i,j)=a(i)*b(j)
enddo
enddo
end function

!************************************************************************************

SUBROUTINE unit_matrix(mat)
implicit none
real(8),DIMENSION(:,:), INTENT(INOUT) :: mat
INTEGER :: i,n
!Action:
!Sets the diagonal components of mat to unity, all other components to zero.
!When mat is square, this will be the unit matrix; otherwise, a unit matrix
!with appended rows or columns of zeros.
mat(:,:)=0.0d0
n=min(size(mat,1),size(mat,2))
do i=1,n
mat(i,i)=1.0d0
end do
END SUBROUTINE unit_matrix

!************************************************************************************

FUNCTION assert_eq(n1,n2,n3,n4,string)
implicit none
CHARACTER(LEN=*), INTENT(IN) :: string
INTEGER, INTENT(IN) :: n1,n2,n3,n4
INTEGER :: assert_eq
!Action:
!Embedding program dies gracefully with an error message if any of the
!integer arguments are not equal to the first. Otherwise, return the value of
!the first argument. Typical use is for enforcing equality on the sizes of arrays
!passed to a subprogram. nrutil implements and overloads forms with 1, 2,
!3, and 4 integer arguments.
if (n1==n2.and.n2==n3.and.n3==n4) then
assert_eq=n1
else
write (*,*) "error: an assert_eq failed with this tag:", string
STOP "program terminated by assert_eq"
end if
END FUNCTION assert_eq

!************************************************************************************

subroutine findmin(choice,dedv_1,dedv_2,dedv_predict,&
!{\src2tex{textfont=tt}}
!!****f* ABINIT/findmin
!!
!! NAME
!! findmin
!!
!! FUNCTION
!! Compute the minimum of a function whose value
!! and derivative are known at two points,
!! using different algorithms.
!! Also deduce different quantities at this predicted
!! point, and at the two other points
!!
!! COPYRIGHT
!! Copyright (C) 1998-2009 ABINIT group (XG, GMR)
!! This file is distributed under the terms of the
!! GNU General Public License, see ~abinit/COPYING
!! or http://www.gnu.org/copyleft/gpl.txt .
!! For the initials of contributors, see ~abinit/doc/developers/contributors.txt .
!!

!! INPUTS
!! choice=1,uses a linear interpolation of the derivatives
!!       =2,uses a quadratic interpolation based on the
!!        values of the function, and the second derivative at mid-point
!!       =3,uses a cubic interpolation
!!       =4,uses a quartic interpolation, with the supplementary
!!          condition that the second derivative vanishes at one and
!!          only one point (See Schlegel, J. Comp. Chem. 3, 214 (1982).
!!          For this option, lambda_1 must be 1 (new point),
!!          and lambda_2 must be 0 (old point).
!!          Also, if the derivative at the new point is more negative
!!          than the derivative at the old point, the predicted
!!          point cannot correspond to a minimum, but will be lambda=2.5d0,
!!          if the energy of the second point is lower than the energy
!!          of the first point.
!! etotal_1=first value of the function
!! etotal_2=second value of the function
!! dedv_1=first value of the derivative
!! dedv_2=second value of the derivative
!! lambda_1=first value of the argument
!! lambda_2=second value of the argument
!!
!! OUTPUT
!! dedv_predict=predicted value of the derivative (usually zero,
!!  except if choice=4, if it happens that a minimum cannot be located,
!!  and a trial step is taken)
!! d2edv2_predict=predicted value of the second derivative (not if choice=4)
!! d2edv2_1=first value of the second derivative (not if choice=4)
!! d2edv2_2=second value of the second derivative (not if choice=4)
!! etotal_predict=predicted value of the function
!! lambda_predict=predicted value of the argument
!! status= 0 if everything went normally ;
!!         1 if negative second derivative
!!         2 if some other problem
!!!!My part:: if choice=4 and ee=0 then status=3
!!
!! PARENTS
!!      brdene,scfcge
!!
!! CHILDREN
!!      leave_new,wrtout
!!
!! SOURCE

& d2edv2_1,d2edv2_2,d2edv2_predict,&
& etotal_1,etotal_2,etotal_predict,&
& lambda_1,lambda_2,lambda_predict,status)


!This section has been created automatically by the script Abilint (TD).
!Do not modify the following lines by hand.
!End of the abilint section

 implicit none

!Arguments ------------------------------------
!scalars
 integer,intent(in) :: choice
 integer,intent(out) :: status
 real(8),intent(in) :: dedv_1,dedv_2,etotal_1,etotal_2,lambda_1,lambda_2
 real(8),intent(out) :: d2edv2_1,d2edv2_2,d2edv2_predict,dedv_predict
 real(8),intent(out) :: etotal_predict,lambda_predict

!Local variables-------------------------------
!scalars
 integer :: printvol
 real(8) :: aa,bb,bbp,cc,ccp,d2edv2_mid,d_lambda,dd,dedv_2bis,dedv_mid1
 real(8) :: dedv_mid2,discr,ee,eep,etotal_2bis,lambda_shift,sum1,sum2,sum3,uu
 real(8) :: uu3,vv,vv3
 real(8) :: tol12
 character(len=500) ::message 

! *************************************************************************
 tol12=0.000000000001d0
 open(unit=16,file="geopt.mon")
!DEBUG
!write(6,*)' findmin : enter'
!write(6,*)' choice,lambda_1,lambda_2=',choice,lambda_1,lambda_2
!ENDDEBUG

 status=0
 d_lambda=lambda_1-lambda_2

!DEBUG
!do choice=3,1,-1
!ENDDEBUG

 if(choice==3)then

! Evaluate cubic interpolation
! etotal = aa + bb * lambda + cc * lambda**2 + dd * lambda**3
  dedv_mid1=(dedv_1+dedv_2)/2.0d0
  dedv_mid2=(etotal_1-etotal_2)/d_lambda
  d2edv2_mid=(dedv_1-dedv_2)/d_lambda
  dd=2.0d0 * ( dedv_mid1 - dedv_mid2 ) / (d_lambda**2)
  cc=0.5d0 * ( d2edv2_mid - 3.0d0*dd*(lambda_1+lambda_2) )
  bb=dedv_2 - 2*cc*lambda_2 - 3*dd*lambda_2**2
  aa=etotal_2 - bb*lambda_2 - cc*lambda_2**2 - dd*lambda_2**3

! Find the lambda at the minimum
  discr=cc*cc-3*bb*dd
  if(discr<0.0d0)then
   write(16, '(a,a,/,a)' )"  ",&
&   ' findmin : BUG -',&
&   '  The 2nd degree equation has no root (choice=3).'
!   call wrtout(06,*,'COLL')
!   call leave_new('COLL')
  end if
  discr=sqrt(discr)
! The root that gives a minimum corresponds to  +discr
  lambda_predict=(-cc+discr)/(3.0d0*dd)

! Predict etotal at that lambda
  etotal_predict=aa+lambda_predict*(bb&
&  +lambda_predict*(cc&
&  +lambda_predict* dd  ))
  dedv_predict=bb+2.0d0*cc*lambda_predict&
&  +3.0d0*dd*lambda_predict**2
  d2edv2_1=2*cc+6*dd*lambda_1
  d2edv2_2=2*cc+6*dd*lambda_2
  d2edv2_predict=2*cc+6*dd*lambda_predict

 else if(choice==4)then

  if(abs(lambda_1-1.0d0)>tol12 .or. abs(lambda_2)>tol12) then
   write(16, '(a,a,/,a)' )"  ",&
&   ' findmin : BUG -',&
&   '  For choice=4, lambda_1 must be 1 and lambda_2 must be 0.'
!   call wrtout(06,*,'COLL')
!   call leave_new('COLL')
  end if

! Evaluate quartic interpolation
! etotal = aa + bb * lambda + cc * lambda**2 + dd * lambda**3 + ee * lambda**4
! Impose positive second derivative everywhere, with
! one point where it vanishes :  3*dd**2=8*cc*ee
  aa=etotal_2
  bb=dedv_2
  sum1=etotal_1-aa-bb
  sum2=dedv_1-bb
  sum3=sum2-2.0d0*sum1

! Build the discriminant of the associated 2nd degree equation
  discr=sum2**2-3.0d0*sum3**2
  if(discr<0.0d0 .or. sum2<0.0d0)then

!  Even if there is a problem, try to keep going ...
   write(16, '(a,a,/,a)' )"  ",&
&   ' findmin : WARNING -',&
&   '  The 2nd degree equation has no positive root (choice=4).'
   status=2
!   call wrtout(06,*,'COLL')
   if(etotal_1<etotal_2)then
    write(16, '(a,a,/,a,/,a)' )"  ",&
&    ' findmin : COMMENT -',&
&    '  Will continue, since the new total energy is lower',&
&    '  than the old. Take a larger step in the same direction.'
!    call wrtout(06,*,'COLL')
    lambda_predict=1.5d0
   else
    write(16, '(a,a,/,a,/,a,/,a)' )"  ",&
&    ' findmin : COMMENT -',&
&    '  There is a problem, since the new total energy is larger',&
&    '  than the old (choice=4).',&
&    '  I take a point between the old and new, close to the old .'
!    call wrtout(06,*,'COLL')
    lambda_predict=0.25d0
   end if
!  Mimick a zero-gradient lambda, in order to avoid spurious
!  action of the inverse hessian (the next line would be a realistic estimation)
   dedv_predict=0.0d0
!  dedv_predict=dedv_2+lambda_predict*(dedv_1-dedv_2)
!  Uses the energies, and the gradient at lambda_2
   etotal_predict=etotal_2+dedv_2*lambda_predict&
&   +(etotal_1-etotal_2-dedv_2)*lambda_predict**2

  else

!  Here, there is an acceptable solution to the 2nd degree equation
   discr=sqrt(discr)
!  The root that gives the smallest ee corresponds to  -discr
!  This is the one to be used: one aims at modelling the
!  behaviour of the function as much as possible with the
!  lowest orders of the polynomial, not the quartic term.
   ee=(sum2-discr)*0.5d0
   dd=sum3-2.0d0*ee
   cc=sum1-dd-ee

!  My additional part
   if(ee==0.d0) then
   status=3
   return
   endif
!  END My additional part   



!  DEBUG
!  write(6,*)'aa,bb,cc,dd,ee',aa,bb,cc,dd,ee
!  ENDDEBUG

!  Now, must find the unique root of
!  $0 = bb + 2*cc * lambda + 3*dd * lambda^2 + 4*ee * lambda^3$
!  This root is unique because it was imposed that the second derivative
!  of the quartic polynomial is everywhere positive.
!  First, remove the quadratic term, by a shift of lambda
!  lambdap=lambda-lambda_shift
!  $0 = bbp + ccp * lambdap + eep * lambdap^3$
   eep=4.0d0*ee
   lambda_shift=-dd/(4.0d0*ee)
   ccp=2.0d0*cc-12.0d0*ee*lambda_shift**2
   bbp=bb+ccp*lambda_shift+eep*lambda_shift**3

!  DEBUG
!  write(6,*)'bbp,ccp,eep,lambda_shift',bbp,ccp,eep,lambda_shift
!  ENDDEBUG

!  The solution of a cubic polynomial equation is as follows :
   discr=(bbp/eep)**2+(4.0d0/27.0d0)*(ccp/eep)**3
   write(16,*)"Discr",discr,(4.0d0/27.0d0)*(ccp/eep)**3,(bbp/eep)**2
   if(discr.lt.0.d0) then
   status=3
   return
   endif
!  In the present case, discr will always be positive
   discr=sqrt(discr)
   uu3=0.5d0*(-bbp/eep+discr) ; uu=sign((abs(uu3))**(1.0d0/3.0d0),uu3)
   vv3=0.5d0*(-bbp/eep-discr) ; vv=sign((abs(vv3))**(1.0d0/3.0d0),vv3)
   lambda_predict=uu+vv
   write(16,*) "Shift",lambda_shift,lambda_predict,bbp,eep,discr
!  Restore the shift
   lambda_predict=lambda_predict+lambda_shift
   etotal_predict=aa+bb*lambda_predict+cc*lambda_predict**2+&
&   dd*lambda_predict**3+ee*lambda_predict**4
   dedv_predict=bb+2.0d0*cc*lambda_predict+3.0d0*dd*lambda_predict**2+&
&   4.0d0*ee*lambda_predict**3
   d2edv2_1=2*cc+6*dd*lambda_1+12*ee*lambda_1**2
   d2edv2_2=2*cc+6*dd*lambda_2+12*ee*lambda_2**2
   d2edv2_predict=2*cc+6*dd*lambda_predict+12*ee*lambda_predict**2

  end if

 else if(choice==1) then
 write(16, '(a,i3)' )'   line minimization, algorithm ',choice

! Use the derivative information to predict lambda
  d2edv2_mid=(dedv_1-dedv_2)/d_lambda
  lambda_predict=lambda_2-dedv_2/d2edv2_mid
  dedv_predict=dedv_2+(lambda_predict-lambda_2)*d2edv2_mid
  d2edv2_1=d2edv2_mid
  d2edv2_2=d2edv2_mid
  d2edv2_predict=d2edv2_mid
! also use the first energy to predict new energy
  etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&  +0.5d0*d2edv2_1*(lambda_predict-lambda_1)**2
  etotal_2bis=etotal_1+dedv_1*(lambda_2-lambda_1)&
&  +0.5d0*d2edv2_1*(lambda_2-lambda_1)**2

  if(d2edv2_mid<0.0d0)then
   write(16, '(a,a,/,a,es18.10,a)' ) "  ",&
&   ' findmin : WARNING -',&
&   '  (scfcge) The second derivative is negative, equal to',d2edv2_mid        ,'.'
!   call wrtout(6,*,'COLL')
   status=1
  end if

 else if(choice==2) then

! Use energies and first derivative information
! etotal = aa + bb * lambda + cc * lambda**2
  dedv_mid2=(etotal_1-etotal_2)/d_lambda
  cc=(dedv_1-dedv_mid2)/d_lambda
  lambda_predict=lambda_1-0.5d0*dedv_1/cc
  d2edv2_1=2*cc
  d2edv2_2=d2edv2_1
  d2edv2_predict=d2edv2_1
  if(d2edv2_predict<0.0d0)then
   write(16, '(a,a,/,a,es18.10,a,/,a)' ) "  ",&
&   ' findmin : WARNING -',&
&   '  (scfcge) The second derivative is negative, equal to',d2edv2_predict,'.',&
&   '  (scfcge) => Pivoting                     '
!  call wrtout(6,*,'COLL')
   status=1
   if(etotal_2 < etotal_1)then
    lambda_predict=lambda_2-0.5d0*(lambda_1-lambda_2)
   else
    lambda_predict=lambda_1-0.5d0*(lambda_2-lambda_1)
   end if
  end if
  dedv_predict=dedv_1+(lambda_predict-lambda_1)*d2edv2_1
  dedv_2bis=dedv_1+(lambda_2-lambda_1)*d2edv2_1
  etotal_predict=etotal_1+dedv_1*(lambda_predict-lambda_1)&
&  +0.5d0*d2edv2_1*(lambda_predict-lambda_1)**2

 end if
 printvol=1
 if(choice==4)printvol=2
 if(printvol==2)then
  write(16, '(a,i3)' )'   line minimization, algorithm ',choice
! call wrtout(6,*,'COLL')
  write(16, '(a,a)' )'                        lambda      etotal ',&
&  '           dedv        d2edv2    '
! call wrtout(6,*,'COLL')
  write(16, '(a,es12.4,es18.10,2es12.4)' )&
&  '   old point         :',lambda_2,etotal_2,dedv_2,d2edv2_2
! call wrtout(6,*,'COLL')
  write(16, '(a,es12.4,es18.10,2es12.4)' )&
&  '   new point         :',lambda_1,etotal_1,dedv_1,d2edv2_1
! call wrtout(6,*,'COLL')
  write(16, '(a,es12.4,es18.10,2es12.4)' )&
&  '   predicted point   :',&
&  lambda_predict,etotal_predict,dedv_predict,d2edv2_predict
! call wrtout(6,*,'COLL')
  if(choice==1) then
   write(16, '(a,es10.4)' ) &
&   ' consistency check :    etotal_2 =',etotal_2bis
!  call wrtout(6,*,'COLL')
  end if
  if(choice==2) then
   write(16, '(a,es10.4)' ) &
&   ' consistency check :    dedv_2 =',dedv_2bis
!  call wrtout(6,*,'COLL')
  end if
  write(16, '(a)' ) ' '
!  call wrtout(6,*,'COLL')
 else if(printvol==1)then
  write(16, '(a,es12.4,a,es18.10)' ) &
&  '   findmin : lambda_predict ',lambda_predict,&
&  '   etotal_predict ',etotal_predict
!  call wrtout(6,*,'COLL')
 end if

end subroutine findmin


subroutine get_BFGS_forces_max(parini,pos_all,force_all,enthalpy,getwfk,iprec,latvec_in,xred_in,etot_in,fcart_in,strten_in)
use mod_parini, only: typ_parini
!This routine hides away all cumbersome conversion of arrays in lattice and positions and forces and stresses
!such that they can be directly passed on to bfgs. It also outputs the enthalpy instead of the energy
implicit none
type(typ_parini), intent(in):: parini
integer:: iprec,iat
real(8):: pos_all(3*parini%nat+9)
real(8):: force_all(3*parini%nat+9)
real(8):: enthalpy,pressure,vol
real(8):: xred_in(3,parini%nat),fcart_in(3,parini%nat),strten_in(6),etot_in,latvec_in(3,3),transformed(3,3),transformed_inv(3,3)
real(8):: str_matrix(3,3),flat(3,3),pressure_mat(3,3),tmplat(3,3),sigma(3,3),crossp(3),stressvol(3,3)
logical:: getwfk

!Here we perform the force call
!****************************************************************************************************************        
!****************************************************************************************************************        
!Get the current acell_in and rprim_in, and also the atomic positions
!Generate a set of variables containing all degrees of freedome
    do iat=1,parini%nat
      xred_in(:,iat)=pos_all((iat-1)*3+1:iat*3)
    enddo
    do iat=1,3
      latvec_in(:,iat)=pos_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)
    enddo
       call get_energyandforces_single(parini,latvec_in,xred_in,fcart_in,strten_in,etot_in,iprec,getwfk)
!****************************************************************************************************************   
!Conversion of forces is more complicate: 
!start with atomic forces
        call fxyz_cart2int(parini%nat,fcart_in,force_all(1:3*parini%nat),latvec_in)
!now the stresses   
!Convert strten into "stress", which is actually the forces on the cell vectors
       str_matrix(1,1)=strten_in(1)
       str_matrix(2,2)=strten_in(2)
       str_matrix(3,3)=strten_in(3)
       str_matrix(1,2)=strten_in(6)
       str_matrix(2,1)=strten_in(6)
       str_matrix(1,3)=strten_in(5)
       str_matrix(3,1)=strten_in(5)
       str_matrix(2,3)=strten_in(4)
       str_matrix(3,2)=strten_in(4)
       call getvol(latvec_in,vol)
       transformed(:,1)=latvec_in(1,:)
       transformed(:,2)=latvec_in(2,:)
       transformed(:,3)=latvec_in(3,:)
       call invertmat(transformed,transformed_inv,3)
       flat=(-vol*matmul(str_matrix,transformed_inv))
       pressure=parini%target_pressure_habohr
       call stress_volume(latvec_in,vol,pressure,stressvol)
       flat=flat+stressvol
!Finally, write those values into fxyz
        do iat=1,3
          force_all(3*parini%nat+(iat-1)*3+1:3*parini%nat+iat*3)=flat(:,iat)
        enddo
!Get the enthalpy
        call get_enthalpy(latvec_in,etot_in,pressure,enthalpy)
end subroutine


subroutine get_fmax(parini,fcart_in,strten_in,fmax,fmax_at,fmax_lat)
use mod_parini, only: typ_parini
type(typ_parini), intent(in):: parini
implicit none
integer:: iat,i,istr
real(8):: fcart_in(3,parini%nat),strten_in(6),fmax,fmax_at,fmax_lat
real(8):: dstr(6), strtarget(6)

!!Compute maximal component of forces, EXCLUDING any fixed components
 fmax=0.d0
 fmax_at=0.d0
 fmax_lat=0.d0
 do iat=1,parini%nat
   do i=1,3
!     if (dtsets(1)%iatfix(i,iat) /= 1) then
       if( abs(fcart_in(i,iat)) >= fmax_at ) fmax_at=abs(fcart_in(i,iat))
!     end if
   end do
 end do
 strtarget=0.d0
 strtarget(1:3)=-parini%target_pressure_habohr
 dstr(:)=strten_in(:)-strtarget(:)
!Eventually take into account the stress
 do istr=1,6
     if(abs(dstr(istr))*parini%paropt_geopt%strfact >= fmax_lat ) fmax_lat=abs(dstr(istr))*parini%paropt_geopt%strfact
 end do
 fmax=max(fmax_at,fmax_lat)
end subroutine

        subroutine stress_volume(latvec,vol,pressure,stressvol)
        !This subroutine will compute the additional components to the negative
        !derivative of the enthalpy with respect to the cell variables
        !For the derivatives with respect to hij, the lattive vector components of the
        !lattece-matrix h, the relation on http://en.wikipedia.org/wiki/Determinant is used:
        !\frac{\partial \det(h)}{\partial h_{ij}}= \det(h)(h^{-1})_{ji}. 
        implicit none
        integer:: i,j
        real(8):: latvec(3,3),vol,stressvol(3,3),inv_latvec(3,3),pressure
        stressvol=0.d0
        call invertmat(latvec,inv_latvec,3)
        do i=1,3
           stressvol(i,1)=-vol*inv_latvec(1,i)*pressure
           stressvol(i,2)=-vol*inv_latvec(2,i)*pressure
           stressvol(i,3)=-vol*inv_latvec(3,i)*pressure
        enddo
        end subroutine
