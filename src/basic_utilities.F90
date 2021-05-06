!*****************************************************************************************
module mod_rng
    implicit none
    private
    integer:: iseed=11
interface random_number_generator_simple
    module procedure random_number_generator_simple_single
    module procedure random_number_generator_simple_array
end interface random_number_generator_simple
    public:: random_number_generator_simple
contains
subroutine random_number_generator_simple_single(rnd)
    implicit none
    real(8), intent(out):: rnd
    !local variables
    integer:: irnd1, irnd2, nlarge
    real(8):: rlarge
    nlarge=huge(1)
    rlarge=real(nlarge,kind=8)
    irnd1=modulo((57*iseed+1),nlarge)
    irnd2=modulo((57*irnd1+1),nlarge)
    iseed=irnd2
    rnd=real(irnd2,kind=8)/rlarge
end subroutine random_number_generator_simple_single
subroutine random_number_generator_simple_array(n,rnd)
    implicit none
    integer, intent(in):: n
    real(8), intent(out):: rnd(n)
    !local variables
    integer:: i, irnd1, irnd2, nlarge
    real(8):: rlarge
    nlarge=huge(1)
    rlarge=real(nlarge,kind=8)
    do i=1,n
        irnd1=modulo((57*iseed+1),nlarge)
        irnd2=modulo((57*irnd1+1),nlarge)
        iseed=irnd2
        rnd(i)=real(irnd2,kind=8)/rlarge
    enddo
end subroutine random_number_generator_simple_array
end module mod_rng
!*****************************************************************************************
module mod_utils
    use mod_rng
    implicit none
end module mod_utils
!*****************************************************************************************
subroutine elim_moment_alborz(nat,atomic_vector)
  implicit none
  integer, intent(in):: nat
  real(8), intent(inout):: atomic_vector(3,nat)
  !local variables
  integer:: iat
  real(8):: sx, sy, sz
  sx=0.d0 ; sy=0.d0 ; sz=0.d0
  do iat=1,nat
     sx=sx+atomic_vector(1,iat)
     sy=sy+atomic_vector(2,iat)
     sz=sz+atomic_vector(3,iat)
  enddo
  sx=sx/real(nat,8) ; sy=sy/real(nat,8) ; sz=sz/real(nat,8)
  do iat=1,nat
     atomic_vector(1,iat)=atomic_vector(1,iat)-sx
     atomic_vector(2,iat)=atomic_vector(2,iat)-sy
     atomic_vector(3,iat)=atomic_vector(3,iat)-sz
  enddo
end subroutine elim_moment_alborz
!*****************************************************************************************
subroutine elim_moment_mass(nat,atomic_vector,atomic_mass)
  implicit none
  integer, intent(in):: nat
  real(8), intent(inout):: atomic_vector(3,nat)
  real(8), intent(in):: atomic_mass(nat)
  !local variables
  integer:: iat
  real(8):: sx, sy, sz, mass_tot, t1 
  sx=0.d0 ; sy=0.d0 ; sz=0.d0 
  mass_tot=0.d0
  do iat=1,nat
     t1=atomic_mass(iat)
     mass_tot=mass_tot+t1
     sx=sx+t1*atomic_vector(1,iat)
     sy=sy+t1*atomic_vector(2,iat)
     sz=sz+t1*atomic_vector(3,iat)
  enddo
  sx=sx/mass_tot ; sy=sy/mass_tot ; sz=sz/mass_tot
  do iat=1,nat
     atomic_vector(1,iat)=atomic_vector(1,iat)-sx
     atomic_vector(2,iat)=atomic_vector(2,iat)-sy
     atomic_vector(3,iat)=atomic_vector(3,iat)-sz
  enddo
end subroutine elim_moment_mass
!*****************************************************************************************
subroutine calc_rotation_eigenvectors(nat,rat0,vrot)
  implicit none
  integer, intent(in) :: nat
  real(8), dimension(3*nat), intent(in) :: rat0
  !real(8), dimension(3*nat), intent(inout) :: fat
  real(8), dimension(3*nat,3), intent(out) :: vrot
  !local variables
  character(len=*), parameter :: subname='elim_torque_reza_alborz'
  integer :: i,iat,i_all,i_stat
  real(8) :: vrotnrm,cmx,cmy,cmz,alpha,totmass, DNRM2
  !this is an automatic array but it should be allocatable
  real(8), dimension(3) :: evaleria
  real(8), dimension(3,3) :: teneria
  real(8), dimension(3*nat) :: rat
  real(8), dimension(:), allocatable :: amass
  
  allocate(amass(nat),stat=i_stat)

  rat=rat0
  amass(1:nat)=1.d0
  !project out rotations
  totmass=0.d0
  cmx=0.d0 
  cmy=0.d0
  cmz=0.d0
  do i=1,3*nat-2,3
     iat=(i+2)/3
     cmx=cmx+amass(iat)*rat(i+0)
     cmy=cmy+amass(iat)*rat(i+1)
     cmz=cmz+amass(iat)*rat(i+2)
     totmass=totmass+amass(iat)
  enddo
  cmx=cmx/totmass 
  cmy=cmy/totmass 
  cmz=cmz/totmass
  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)-cmx
     rat(i+1)=rat(i+1)-cmy
     rat(i+2)=rat(i+2)-cmz
  enddo

  call moment_of_inertia_alborz(nat,rat,teneria,evaleria)
  do iat=1,nat
     i=iat*3-2
     call mycross(teneria(1,1),rat(i),vrot(i,1))
     call mycross(teneria(1,2),rat(i),vrot(i,2))
     call mycross(teneria(1,3),rat(i),vrot(i,3))
  enddo
  call normalizevector_alborz(3*nat,vrot(1,1))
  call normalizevector_alborz(3*nat,vrot(1,2))
  call normalizevector_alborz(3*nat,vrot(1,3))
  
  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)+cmx
     rat(i+1)=rat(i+1)+cmy
     rat(i+2)=rat(i+2)+cmz
  enddo

  vrotnrm=DNRM2(3*nat,vrot(1,1),1)
  if (vrotnrm /= 0.d0) vrot(1:3*nat,1)=vrot(1:3*nat,1)/vrotnrm
  vrotnrm=DNRM2(3*nat,vrot(1,2),1)
  if (vrotnrm /= 0.d0) vrot(1:3*nat,2)=vrot(1:3*nat,2)/vrotnrm
  vrotnrm=DNRM2(3*nat,vrot(1,3),1)
  if (vrotnrm /= 0.d0) vrot(1:3*nat,3)=vrot(1:3*nat,3)/vrotnrm
  !do i=1,3
  !   alpha=0.d0  
  !   if(abs(evaleria(i))>1.d-10) then
  !      alpha=dot_product(vrot(:,i),fat(:))
  !      fat(:)=fat(:)-alpha*vrot(:,i) 
  !   endif
  !enddo
  !i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
end subroutine calc_rotation_eigenvectors
!*****************************************************************************************
!> Eliminate the translational forces before calling this subroutine!!!
!! Main subroutine: Input is nat (number of atoms), rat0 (atomic positions) and fat (forces on atoms)
!! The atomic positions will be returned untouched
!! In fat, the rotational forces will be eliminated with respect to the center of mass. 
!! All atoms are treated equally (same atomic mass) 
subroutine elim_torque_reza_alborz(nat,rat0,fat)
  implicit none
  integer, intent(in) :: nat
  real(8), dimension(3*nat), intent(in) :: rat0
  real(8), dimension(3*nat), intent(inout) :: fat
  !local variables
  character(len=*), parameter :: subname='elim_torque_reza_alborz'
  integer :: i,iat,i_all,i_stat
  real(8) :: vrotnrm,cmx,cmy,cmz,alpha,totmass, DNRM2
  !this is an automatic array but it should be allocatable
  real(8), dimension(3) :: evaleria
  real(8), dimension(3,3) :: teneria
  real(8), dimension(3*nat) :: rat
  real(8), dimension(3*nat,3) :: vrot
  real(8), dimension(:), allocatable :: amass
  
  allocate(amass(nat),stat=i_stat)

  rat=rat0
  amass(1:nat)=1.d0
  !project out rotations
  totmass=0.d0
  cmx=0.d0 
  cmy=0.d0
  cmz=0.d0
  do i=1,3*nat-2,3
     iat=(i+2)/3
     cmx=cmx+amass(iat)*rat(i+0)
     cmy=cmy+amass(iat)*rat(i+1)
     cmz=cmz+amass(iat)*rat(i+2)
     totmass=totmass+amass(iat)
  enddo
  cmx=cmx/totmass 
  cmy=cmy/totmass 
  cmz=cmz/totmass
  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)-cmx
     rat(i+1)=rat(i+1)-cmy
     rat(i+2)=rat(i+2)-cmz
  enddo

  call moment_of_inertia_alborz(nat,rat,teneria,evaleria)
  do iat=1,nat
     i=iat*3-2
     call mycross(teneria(1,1),rat(i),vrot(i,1))
     call mycross(teneria(1,2),rat(i),vrot(i,2))
     call mycross(teneria(1,3),rat(i),vrot(i,3))
  enddo
  call normalizevector_alborz(3*nat,vrot(1,1))
  call normalizevector_alborz(3*nat,vrot(1,2))
  call normalizevector_alborz(3*nat,vrot(1,3))
  
  do i=1,3*nat-2,3
     rat(i+0)=rat(i+0)+cmx
     rat(i+1)=rat(i+1)+cmy
     rat(i+2)=rat(i+2)+cmz
  enddo

  vrotnrm=DNRM2(3*nat,vrot(1,1),1)
  if (vrotnrm /= 0.d0) vrot(1:3*nat,1)=vrot(1:3*nat,1)/vrotnrm
  vrotnrm=DNRM2(3*nat,vrot(1,2),1)
  if (vrotnrm /= 0.d0) vrot(1:3*nat,2)=vrot(1:3*nat,2)/vrotnrm
  vrotnrm=DNRM2(3*nat,vrot(1,3),1)
  if (vrotnrm /= 0.d0) vrot(1:3*nat,3)=vrot(1:3*nat,3)/vrotnrm
  do i=1,3
     alpha=0.d0  
     if(abs(evaleria(i))>1.d-10) then
        alpha=dot_product(vrot(:,i),fat(:))
        fat(:)=fat(:)-alpha*vrot(:,i) 
     endif
  enddo
  !i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
end subroutine elim_torque_reza_alborz
!*****************************************************************************************
subroutine mycross(a,b,c)
  implicit none
  real(8), dimension(3), intent(in):: a,b
  real(8), dimension(3), intent(out):: c
  c(1)=a(2)*b(3)-b(2)*a(3)
  c(2)=a(3)*b(1)-b(3)*a(1)
  c(3)=a(1)*b(2)-b(1)*a(2)
end subroutine mycross
!*****************************************************************************************
subroutine moment_of_inertia_alborz(nat,rat,teneria,evaleria)
  implicit none
  integer, intent(in) :: nat
  real(8), dimension(3,nat), intent(in) :: rat
  real(8), dimension(3), intent(out) :: evaleria
  real(8), dimension(3,3), intent(out) :: teneria
  !local variables
  character(len=*), parameter :: subname='moment_of_inertia_alborz'
  integer, parameter::lwork=100
  integer :: iat,info,i_all,i_stat
  real(8) :: tt
  real(8), dimension(lwork) :: work
  real(8), dimension(:), allocatable :: amass
  allocate(amass(nat),stat=i_stat)
  !positions relative to center of geometry
  amass(1:nat)=1.d0
  !calculate inertia tensor
  teneria(1:3,1:3)=0.d0
  do iat=1,nat
     tt=amass(iat)
     teneria(1,1)=teneria(1,1)+tt*(rat(2,iat)*rat(2,iat)+rat(3,iat)*rat(3,iat))
     teneria(2,2)=teneria(2,2)+tt*(rat(1,iat)*rat(1,iat)+rat(3,iat)*rat(3,iat))
     teneria(3,3)=teneria(3,3)+tt*(rat(1,iat)*rat(1,iat)+rat(2,iat)*rat(2,iat))
     teneria(1,2)=teneria(1,2)-tt*(rat(1,iat)*rat(2,iat))
     teneria(1,3)=teneria(1,3)-tt*(rat(1,iat)*rat(3,iat))
     teneria(2,3)=teneria(2,3)-tt*(rat(2,iat)*rat(3,iat))
     teneria(2,1)=teneria(1,2)
     teneria(3,1)=teneria(1,3)
     teneria(3,2)=teneria(2,3)
  enddo
  !diagonalize inertia tensor
  call DSYEV('V','L',3,teneria,3,evaleria,work,lwork,info)
  i_all=-product(shape(amass))*kind(amass)
  deallocate(amass,stat=i_stat)
end subroutine moment_of_inertia_alborz
!*****************************************************************************************
subroutine normalizevector_alborz(n,v)
    implicit none
    integer, intent(in):: n
    real(8), intent(inout):: v(n)
    !local variables
    integer:: i
    real(8):: vnrm
    !integer, save:: icall=0
    !icall=icall+1
    vnrm=0.d0
    do i=1,n
       vnrm=vnrm+v(i)**2
    enddo
    vnrm=sqrt(vnrm)
    !write(21,'(i5,es24.15)') icall,vnrm
    if(vnrm/=0.d0) v(1:n)=v(1:n)/vnrm
end subroutine normalizevector_alborz
!*****************************************************************************************
subroutine calnorm(n,v,vnrm)
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: v(n)
    real(8), intent(out):: vnrm
    !local variables
    integer:: i
    vnrm=0.d0
    do i=1,n
        vnrm=vnrm+v(i)**2
    enddo
    vnrm=sqrt(vnrm)
end subroutine calnorm
!*****************************************************************************************
subroutine calmaxforcecomponent(n,v,vmax)
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: v(n)
    real(8), intent(out):: vmax
    !local variables
    integer:: i
    vmax=0.d0
    do i=1,n
        vmax=max(vmax,abs(v(i)))
    enddo
end subroutine calmaxforcecomponent
!*****************************************************************************************
subroutine rxyz_cart2int_alborz(nat,latvec,rxyzcart,rxyzint)
    !This subrouine will convert the internal coordinates into cartesian coordinates
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: rxyzcart(3,nat), latvec(3,3)
    real(8), intent(out):: rxyzint(3,nat)
    !local variables
    integer:: iat
    real(8):: latvecinv(3,3)
    real(16):: latvec_q(3,3), latvecinv_q(3,3)
    latvec_q(1:3,1:3)=real(latvec(1:3,1:3),kind=16)
    call invertmat_alborz_qp(latvec_q,latvecinv_q)
    latvecinv(1:3,1:3)=real(latvecinv_q(1:3,1:3),kind=8)
    do iat=1,nat
        rxyzint(1,iat)=latvecinv(1,1)*rxyzcart(1,iat)+latvecinv(1,2)*rxyzcart(2,iat)+latvecinv(1,3)*rxyzcart(3,iat)
        rxyzint(2,iat)=latvecinv(2,1)*rxyzcart(1,iat)+latvecinv(2,2)*rxyzcart(2,iat)+latvecinv(2,3)*rxyzcart(3,iat)
        rxyzint(3,iat)=latvecinv(3,1)*rxyzcart(1,iat)+latvecinv(3,2)*rxyzcart(2,iat)+latvecinv(3,3)*rxyzcart(3,iat)
    enddo
end subroutine rxyz_cart2int_alborz
!*****************************************************************************************
subroutine rxyz_int2cart_alborz(nat,cellvec,rat_int,rat_cart)
    implicit none
    integer, intent(in):: nat
    real(8), intent(in):: cellvec(3,3), rat_int(3,nat)
    real(8), intent(inout):: rat_cart(3,nat)
    !local variables
    integer:: iat
    do iat=1,nat
        rat_cart(1:3,iat)=matmul(cellvec,rat_int(1:3,iat))
    enddo
end subroutine rxyz_int2cart_alborz
!*****************************************************************************************
subroutine invertmat_alborz(a,ainv)
    implicit none
    real(8),intent(in):: a(3,3)
    real(8),intent(out):: ainv(3,3)
    !local variables
    real(8):: div
    !integer:: ipiv(3), info, ldwork
    !real(8), allocatable:: WORK(:)
    div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+ &
         a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1))
    div=1.d0/div
    ainv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
    ainv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
    ainv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
    ainv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
    ainv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
    ainv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
    ainv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
    ainv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
    ainv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
    !General n*n matrix 
    !ainv=mat
    !allocate(WORK(n))
    !call  DGETRF( n, n, ainv, n, IPIV, INFO )
    !if (info.ne.0) stop "Error in DGETRF"
    !LDWORK=-1
    !call  DGETRI( n, ainv, n, IPIV, WORK,LDWORK , INFO )
    !LDWORK=WORK(1)
    !deallocate(WORK)
    !allocate(WORK(LDWORK))
    !call  DGETRI( n, ainv, n, IPIV, WORK,LDWORK , INFO )
    !if (info.ne.0) stop "Error in DGETRI"
end subroutine invertmat_alborz
!*****************************************************************************************
subroutine invertmat_alborz_qp(a,ainv)
    use mod_defs
    implicit none
    real(kind=fqp),intent(in):: a(3,3)
    real(kind=fqp),intent(out):: ainv(3,3)
    !local variables
    real(16):: div
    !integer:: ipiv(3), info, ldwork
    !real(8), allocatable:: WORK(:)
    div=(a(1,1)*a(2,2)*a(3,3)-a(1,1)*a(2,3)*a(3,2)-a(1,2)*a(2,1)*a(3,3)+ &
         a(1,2)*a(2,3)*a(3,1)+a(1,3)*a(2,1)*a(3,2)-a(1,3)*a(2,2)*a(3,1))
    div=1.0_fqp/div
    ainv(1,1) = (a(2,2)*a(3,3)-a(2,3)*a(3,2))*div
    ainv(1,2) =-(a(1,2)*a(3,3)-a(1,3)*a(3,2))*div
    ainv(1,3) = (a(1,2)*a(2,3)-a(1,3)*a(2,2))*div
    ainv(2,1) =-(a(2,1)*a(3,3)-a(2,3)*a(3,1))*div
    ainv(2,2) = (a(1,1)*a(3,3)-a(1,3)*a(3,1))*div
    ainv(2,3) =-(a(1,1)*a(2,3)-a(1,3)*a(2,1))*div
    ainv(3,1) = (a(2,1)*a(3,2)-a(2,2)*a(3,1))*div
    ainv(3,2) =-(a(1,1)*a(3,2)-a(1,2)*a(3,1))*div
    ainv(3,3) = (a(1,1)*a(2,2)-a(1,2)*a(2,1))*div
end subroutine invertmat_alborz_qp
!*****************************************************************************************
subroutine convertupper(str)
    character(*), intent(inout):: str
    integer:: i
    do i=1,len(str)
        select case(str(i:i))
            case("a":"z")
            str(i:i)=achar(iachar(str(i:i))-32)
        end select
    end do
end subroutine convertupper
!*****************************************************************************************
subroutine convertlower(str)
    character(*), intent(inout) :: str
    integer:: i
    do i=1,len(str)
        select case(str(i:i))
            case("A":"Z")
            str(i:i)=achar(iachar(str(i:i))+32)
        end select
    end do
end subroutine convertlower
!*****************************************************************************************
subroutine check_whether_time_exceeded
    use mod_task, only: time_start, time_exceeded
    use mod_processors, only: iproc
    implicit none
    !local variables
    real(8):: cpulimit, time_now
    integer:: k
    open(unit=55,file='CPUlimit',status='unknown')
    read(55,*,iostat=k) cpulimit
    if(k==0) then !k=0 means there was no error, nor was EOF encountered.
        call cpu_time(time_now)
        if(time_now-time_start>cpulimit) then
            write(*,'(a,i4)') 'CPU time exceeded: iproc',iproc
            time_exceeded=.true.
        endif
    endif
    close(55)
end subroutine check_whether_time_exceeded
!*****************************************************************************************
subroutine expdist(n,x)
    !generates n random numbers distributed according to  exp(-x)
    implicit none
    integer, intent(in):: n
    real(8), intent(out):: x(n)
    !local variables
    integer:: i
    real(8):: tt, ss
    !on Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps.
    real(8), parameter::eps=1.d-8
    do i=1,n
        call random_number(ss)
        tt=eps+(1.d0-2.d0*eps)*real(ss,8)
        x(i)=log(tt)
    enddo
end subroutine expdist
!*****************************************************************************************
subroutine gausdist_alborz(n,x)
    !generates n random numbers distributed according to  exp(-.5*x**2)
    implicit none
    integer, intent(in) ::n
    real(8), intent(out) :: x(n)
    !local variables
    integer:: i
    real(4):: s1, s2, t1, t2, tt, twopi
    !On Intel the random_number can take on the values 0. and 1.. To prevent overflow introduce eps.
    real(8), parameter:: eps=1.d-8
    twopi=8.d0*atan(1.d0)
    do i=1,n-1,2
        call random_number(s1)
        t1=eps+(1.d0-2.d0*eps)*real(s1,8)
        call random_number(s2)
        t2=real(s2,8)
        tt=sqrt(-2.d0*log(t1))
        x(i)=tt*cos(twopi*t2)
        x(i+1)=tt*sin(twopi*t2)
    enddo
    call random_number(s1)
    t1=eps+(1.d0-2.d0*eps)*real(s1,8)
    call random_number(s2)
    t2=real(s2,8)
    tt=sqrt(-2.d0*log(t1))
    x(n)=tt*cos(twopi*t2)
end subroutine gausdist_alborz
!*****************************************************************************************
subroutine randdist(a,n,x)
    !create a uniform random numbers in interval [-a,a]
    implicit none
    real(8), intent(in) :: a
    integer, intent(in):: n
    real(8), intent(out) :: x(n)
    !local variables
    integer:: i
    real(8):: tt1, tt2
    tt1=a*2.d0
    do i=1,n
        call random_number(tt2)
        x(i)=(tt2-0.5d0)*tt1
    enddo
end subroutine randdist
!*****************************************************************************************
subroutine randdist_simple(a,n,x)
    !create a uniform random numbers in interval [-a,a]
    use mod_utils
    implicit none
    real(8), intent(in) :: a
    integer, intent(in):: n
    real(8), intent(out) :: x(n)
    !local variables
    integer:: i
    real(8):: tt1, tt2
    tt1=a*2.d0
    do i=1,n
        call random_number_generator_simple(tt2)
        x(i)=(tt2-0.5d0)*tt1
    enddo
end subroutine randdist_simple
!*****************************************************************************************
subroutine hpsort(n,ra)
    implicit real*8 (a-h,o-z)
    real*8 ::ra(n)
    if (n.lt.2) return
    l=n/2+1
    ir=n
    do
        if(l.gt.1) then
            l=l-1
            rra=ra(l)
        else
            rra=ra(ir)
            ra(ir)=ra(1)
            ir=ir-1
            if(ir.eq.1) then
                ra(1)=rra
                return
            endif
        endif
        i=l
        j=l+l
        do
            if(j.le.ir) then
                if(j.lt.ir) then
                    if(ra(j).lt.ra(j+1))  j=j+1
                endif
                if(rra.lt.ra(j)) then
                    ra(i)=ra(j)
                    i=j
                    j=j+j
                else
                    j=ir+1
                endif
            else
                exit
            endif
        enddo
      ra(i)=rra
    enddo
end subroutine hpsort
!*****************************************************************************************
function flm_index(str1,str2) result(ind)
    implicit none
    character(*), intent(in):: str1, str2
    !local variables
    integer:: ind, indp, len_str
    character(1000):: str
    len_str=len(str1)
    if(len_str>999) then
        write(*,*) 'ERROR: flm_index works only for string whose lengths'
        write(*,*) '       are less than 256 character'
        stop
    endif
    str=trim(str1)//','
    ind=0
    do
        indp=index(str,str2)
        ind=ind+indp
        if(indp<1) then
            ind=0
            return
        endif
        str=str(indp:)
        indp=scan(str,',')-1
        if(indp<1) then
            ind=0
            return
        endif
        !write(*,*) 'BBB ',str(1:indp),'   ',trim(str2)
        if(str(1:indp)==trim(str2)) return
        str(1:indp+1)=''
    enddo
end function flm_index
!*****************************************************************************************
!subroutine projtransout(n,v)
!    implicit none
!    integer::n,istat,i
!    real(8)::v(n),t1,t2,t3,DDOT
!    real(8), allocatable::dirx(:),diry(:),dirz(:)
!    allocate(dirx(n),stat=istat);if(istat/=0) stop 'ERROR: failure deallocating dirx'
!    allocate(diry(n),stat=istat);if(istat/=0) stop 'ERROR: failure deallocating diry'
!    allocate(dirz(n),stat=istat);if(istat/=0) stop 'ERROR: failure deallocating dirz'
!    dirx=0.d0
!    do i=1,n,3
!        dirx(i)=1.d0
!    enddo
!    call normalizevector(n,dirx)
!    !call atom_normalizevector(n/3,dirx)
!    diry=0.d0
!    do i=1,n,3
!        diry(i+1)=1.d0
!    enddo
!    call normalizevector(n,diry)
!    !call atom_normalizevector(n/3,diry)
!    dirz=0.d0
!    do i=1,n,3
!        dirz(i+2)=1.d0
!    enddo
!    call normalizevector(n,dirz)
!    !call atom_normalizevector(n/3,dirz)
!    !-------------------------------------------------------
!    t1=dot_product(v,dirx)
!    t2=dot_product(v,diry)
!    t3=dot_product(v,dirz)
!    !write(61,*) t1,t2,t3
!    !-------------------------------------------------------
!    t1=dot_product(v,dirx)
!    !write(*,*) 't1',t1
!    v(1:n)=v(1:n)-t1*dirx(1:n)
!    t1=dot_product(v,diry)
!    !write(*,*) 't1',t1
!    v(1:n)=v(1:n)-t1*diry(1:n)
!    t1=dot_product(v,dirz)
!    !write(*,*) 't1',t1
!    v(1:n)=v(1:n)-t1*dirz(1:n)
!    deallocate(dirx,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating dirx'
!    deallocate(diry,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating diry'
!    deallocate(dirz,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating dirz'
!end subroutine projtransout
!*****************************************************************************************
!subroutine projrotout(n,nr,x,v)
!    implicit none
!    integer::n,nr,istat,i
!    real(8)::x(n),v(nr),t1,cmx,cmy,cmz,t2,t3,DDOT
!    real(8), allocatable::pxy(:),pyz(:),pxz(:)
!    !write(*,*) 'n,nr',n,nr
!    allocate(pxy(nr),stat=istat);if(istat/=0) stop 'ERROR: failure allocating pxy'
!    allocate(pyz(nr),stat=istat);if(istat/=0) stop 'ERROR: failure allocating pyz'
!    allocate(pxz(nr),stat=istat);if(istat/=0) stop 'ERROR: failure allocating pxz'
!    cmx=0.d0;cmy=0.d0;cmz=0.d0
!    do i=1,n,3
!        cmx=cmx+x(i+0)
!        cmy=cmy+x(i+1)
!        cmz=cmz+x(i+2)
!    enddo
!    cmx=3.d0*cmx/n ; cmy=3.d0*cmy/n ; cmz=3.d0*cmz/n
!    pxy=0.d0
!    do i=1,nr,3
!        pxy(i+0)=-(x(i+1)-cmy)
!        pxy(i+1)= (x(i+0)-cmx)
!    enddo
!    call normalizevector(nr,pxy)
!    !call atom_normalizevector(nr/3,pxy)
!    pxz=0.d0
!    do i=1,nr,3
!        pxz(i+0)=-(x(i+2)-cmz)
!        pxz(i+2)= (x(i+0)-cmx)
!    enddo
!    call normalizevector(nr,pxz)
!    !call atom_normalizevector(nr/3,pxz)
!    pyz=0.d0
!    do i=1,nr,3
!        pyz(i+1)=-(x(i+2)-cmz)
!        pyz(i+2)= (x(i+1)-cmy)
!    enddo
!    call normalizevector(nr,pyz)
!    !call atom_normalizevector(nr/3,pyz)
!    !-------------------------------------------------------
!    t1=DDOT(nr,v,1,pxy,1)
!    t2=DDOT(nr,v,1,pyz,1)
!    t3=DDOT(nr,v,1,pxz,1)
!    !write(51,*) t1,t2,t3
!    !-------------------------------------------------------
!    !t1=dot_product(v,pxy)
!    t1=DDOT(nr,v,1,pxy,1)
!    v(1:nr)=v(1:nr)-t1*pxy(1:nr)
!    !t1=dot_product(v,pyz)
!    t1=DDOT(nr,v,1,pyz,1)
!    v(1:nr)=v(1:nr)-t1*pyz(1:nr)
!    !t1=dot_product(v,pxz)
!    t1=DDOT(nr,v,1,pxz,1)
!    v(1:nr)=v(1:nr)-t1*pxz(1:nr)
!    !-------------------------------------------------------
!    deallocate(pxy,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating pxy'
!    deallocate(pyz,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating pyz'
!    deallocate(pxz,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating pxz'
!end subroutine projrotout
!*****************************************************************************************
subroutine elim_fixed_at(parini,nat,x)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: iat,nat
real(8):: x(3,nat)
!write(*,*) "# Eliminiate"
do iat=1,nat
  if(parini%fixat(iat)) then
     x(:,iat)=0.d0

  endif
!  write(*,*) x(:,iat)
enddo
end subroutine

!************************************************************************************

subroutine elim_fixed_lat(parini,latvec,x)
use mod_parini, only: typ_parini
implicit none
type(typ_parini), intent(in):: parini
integer:: i,k
real(8):: x(3,3),latvec(3,3),lenlat,tmpvec(3)
real(8):: len1,len2,len3,tmp1(3),tmp2(3),tmp3(3),tmp1len,tmp2len,tmp3len
!The fixlat has 7 components:
!a,b,c,alha,beta,gamma, and for fixed cell shape (volume fluctuation)

!We should perform the projection out of the forces self-consistently
do k=1,10
if(parini%fixlat(7)) then
!Here we implement cell fluctuation
!There are only forces along the cell vectors
   call diagcomp(latvec,x)
   return
elseif(all(parini%fixlat(1:6))) then
!When the whole cell is fixed
   x=0.d0
   return
elseif(parini%bc==2) then
!If the boundary condition is free
   x=0.d0
   return
else
!Treat the case where a, b, or c are fixed
do i=1,3
   if(parini%fixlat(i)) then
!Project out the component of x onto latvec
   lenlat=dot_product(latvec(:,i),latvec(:,i))
   tmpvec=dot_product(x(:,i),latvec(:,i))*latvec(:,i)/lenlat
   x(:,i)=x(:,i)-tmpvec(:)
   endif   
enddo

!Now eliminate the compoments of x which would change the angle between the lattice vectors
!The correct way:
do i=1,3
   if(parini%fixlat(i+3)) then
   len1=dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))
   len2=dot_product(latvec(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))
   call cross_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i+1,3)+1),tmpvec)
   call cross_product(tmpvec,latvec(:,modulo(i,3)+1),tmp1)
   call cross_product(tmpvec,latvec(:,modulo(i+1,3)+1),tmp2)
   tmp1len=dot_product(tmp1,tmp1)  
   tmp2len=dot_product(tmp2,tmp2)  
!Project out these components
   tmpvec=dot_product(x(:,modulo(i,3)+1),tmp1)*tmp1/tmp1len
   x(:,modulo(i,3)+1)=x(:,modulo(i,3)+1)-tmpvec
   tmpvec=dot_product(x(:,modulo(i+1,3)+1),tmp2)*tmp2/tmp2len
   x(:,modulo(i+1,3)+1)=x(:,modulo(i+1,3)+1)-tmpvec
   endif
enddo
endif
!write(*,*) "elim_fixed_lat norm", sqrt(sum(x*x))
enddo
!!The simple way:
!!This means only keep the components along the lattice vectors...
!do i=1,3
!   if(fixlat(i+3)) then
!!   write(*,*) i+3,modulo(i,3)+1,modulo(i+1,3)+1
!   lenlat=dot_product(latvec(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))
!   tmpvec=dot_product(x(:,modulo(i,3)+1),latvec(:,modulo(i,3)+1))*latvec(:,modulo(i,3)+1)/lenlat
!!   write(*,*) lenlat
!   x(:,modulo(i,3)+1)=tmpvec(:)
!   lenlat=dot_product(latvec(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))
!   tmpvec=dot_product(x(:,modulo(i+1,3)+1),latvec(:,modulo(i+1,3)+1))*latvec(:,modulo(i+1,3)+1)/lenlat
!!   write(*,*) lenlat
!   x(:,modulo(i+1,3)+1)=tmpvec(:)
!   endif   
!enddo

!!write(*,*) x(:,1)
!!write(*,*) x(:,2)
!!write(*,*) x(:,3)
end subroutine

!************************************************************************************

        subroutine elim_moment(nat,vxyz,atmass)
        implicit none
        real(8):: vxyz(3,nat),sx,sz,sy,atmass(nat)
        integer:: iat,nat       
 
        sx=0.d0 ; sy=0.d0 ; sz=0.d0
        do iat=1,nat
        sx=sx+vxyz(1,iat)*atmass(iat)
        sy=sy+vxyz(2,iat)*atmass(iat)
        sz=sz+vxyz(3,iat)*atmass(iat)
        enddo
        sx=sx/nat ; sy=sy/nat ; sz=sz/nat
        do iat=1,nat
        vxyz(1,iat)=vxyz(1,iat)-sx/atmass(iat)
        vxyz(2,iat)=vxyz(2,iat)-sy/atmass(iat)
        vxyz(3,iat)=vxyz(3,iat)-sz/atmass(iat)
        enddo
        return
        end subroutine

!************************************************************************************

subroutine elim_torque_cell(latvec0,vlat)
implicit none
real(8), intent(in)    :: latvec0(3,3)
real(8), intent(inout) :: vlat(3,3)
real(8) :: torque(3),crossp(3),sx,sy,sz,tmax,cx,cy,cz
integer :: i,ii,it,itmax
real(8) :: unitmat(3,3),rotmat(3,3),rotmatall(3,3),xaxis(3),axis(3),latvec(3,3),tnorm,angle,axisnorm
       latvec=latvec0
       itmax=5000
!Initialize
       unitmat=0.d0
       do i=1,3
       unitmat(i,i)=1.d0
       enddo
       rotmatall=unitmat
       xaxis=(/1.d0,0.d0,0.d0/)

it=0
do
it=it+1
       torque=0.d0
       do i=1,3
       call cross_product(latvec(:,i),vlat(:,i),crossp)
       torque=torque+crossp
       enddo
       !write(*,'(3(e11.3),i6)')torque,loop
       if (it.ge.itmax) goto 1001
       if (torque(1)**2+torque(2)**2+torque(3)**2.lt.1.d-22) goto 1000


        tnorm=sqrt(torque(1)**2+torque(2)**2+torque(3)**2)
        angle=-acos(dot_product(torque,xaxis)/tnorm)
        call cross_product(xaxis,torque,axis)
        axisnorm=sqrt(axis(1)**2+axis(2)**2+axis(3)**2)
        rotmat=unitmat
        if(axisnorm.gt.0.d0) then
        axis=axis/axisnorm
        call rotation(rotmat,angle,axis)
        endif
        rotmatall=matmul(rotmat,rotmatall)
        sy=0.d0 ; sz=0.d0
        do i=1,3
          latvec(:,i)=matmul(rotmat,latvec(:,i))
          vlat(:,i)=matmul(rotmat,vlat(:,i))
          sy=sy+latvec(2,i)**2
          sz=sz+latvec(3,i)**2
        enddo
         cx=torque(1)/(sz+sy)*0.9d0
         do i=1,3
         vlat(2,i)=vlat(2,i)+cx*latvec(3,i)
         vlat(3,i)=vlat(3,i)-cx*latvec(2,i)
         enddo
enddo
1001  write(100,'(a,3(e11.3))') 'WARNING REMAINING TORQUE',torque

1000  call invertmat(rotmatall,rotmat,3)
      do i=1,3
      vlat(:,i)=matmul(rotmat,vlat(:,i))
      enddo
      return
end subroutine elim_torque_cell

!************************************************************************************

subroutine diagcomp(latvec,x)
implicit none
real(8):: latvec(3,3),x(3,3),xnrm,latvect(3,3),latvecinv(3,3),sigma(3,3)
real(8):: len1,len2,len3,tmp1,tmp2,tmp3,tmpvec(3),vol
!There are only forces along the cell vectors
!Convert the target stress to the target gradient, assuming cell volume to be 1
   len1=sqrt(dot_product(latvec(:,1),latvec(:,1)))
   len2=sqrt(dot_product(latvec(:,2),latvec(:,2)))
   len3=sqrt(dot_product(latvec(:,3),latvec(:,3)))
   tmp1=dot_product(x(:,1),latvec(:,1))/len1
   tmp2=dot_product(x(:,2),latvec(:,2))/len2
   tmp3=dot_product(x(:,3),latvec(:,3))/len3
!!!!Convert the target stress to the target gradient, assuming cell volume to be 1
   latvect(:,1)=latvec(1,:)
   latvect(:,2)=latvec(2,:)
   latvect(:,3)=latvec(3,:)
   sigma=matmul(x,latvect)
   xnrm=sigma(1,1)+sigma(2,2)+sigma(3,3)
   xnrm=xnrm/3.d0
   sigma=0.d0
   sigma(1,1)=xnrm
   sigma(2,2)=xnrm
   sigma(3,3)=xnrm
   call invertmat(latvect,latvecinv,3)
   x=matmul(sigma,latvecinv)

   call getvol(latvec,vol)
   vol=vol**(1.d0/3.d0)
!   xnrm=xnrm/vol
!   xnrm=(tmp1/len1+tmp2/len2+tmp3/len3)/3.d0*vol
   xnrm=(x(1,1)+x(2,2)+x(3,3))/3.d0
   vol=1.d0/vol
   x(:,1)=latvec(:,1)*vol*xnrm
   x(:,2)=latvec(:,2)*vol*xnrm
   x(:,3)=latvec(:,3)*vol*xnrm
end subroutine

 subroutine backtocell(nat,latvec,rxyz_red)
 !This subroutine will transform back all atoms into the periodic cell
 !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
 implicit none
 integer:: nat,i,iat,j
 real(8) :: latvec(3,3), rxyz_red(3,nat), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 real(8) :: v(3,3),vol
 logical:: neccesary
!First check if the volume is positive
 call getvol(latvec,vol)
 if(vol.le.0.d0) stop "Negative volume during backtocell"
 do iat=1,nat
        rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.d0),1.d0)
        rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.d0),1.d0)
        rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.d0),1.d0)
 enddo
! v=latvec
! vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
!      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! if(vol.le.0.d0) stop "Negative volume during backtocell"
! call rxyz_int2cart(latvec,rxyz_red,rxyz,nat)
!
! !To really be on the safe side, the translation vector can be shortened by  a factor eps in order
! !to get the atom into the cell. 
!!  eps=1.d-10
!  eps=1.d-15
! ! eps=1.d0-eps
! count=0.d0
! neccesary=.true.
! do while(neccesary)
! neccesary=.false.
! count=count+1.d0
! !generate 3 normal vectors of the 3 planes
! call nveclatvec(latvec,nvec)
! do iat=1,nat
! !3 planes through origin (xy,yz,zx)
! do i=1,3
! dist(i)=DOT_PRODUCT(rxyz(:,iat),nvec(:,i))
!! if(dist(i).lt.0.d0) then
! if(dist(i).lt.-abs(dist(i))*eps) then
!! write(*,*) "unten 1",i,iat
! rxyz(:,iat)=rxyz(:,iat)+latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
!
! !3 planes on top/side/back (xy,yz,zx)
! dist(i+3)=DOT_PRODUCT(rxyz(:,iat)-latvec(:,mod(i+1,3)+1),nvec(:,i))
!! if(dist(i+3).gt.0.d0) then
! if(dist(i+3).gt.abs(dist(i+3))*eps) then
!! write(*,*) "unten 1",i,iat
!! if(dist(i+3).gt.eps) then
!! write(*,*) "oben 2",i,iat
! rxyz(:,iat)=rxyz(:,iat)-latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
! enddo
! enddo
! if(count.gt.1.d6) stop "Too many iterations in back-to-cell"
! enddo
! 
! call rxyz_cart2int(latvec,rxyz_red,rxyz,nat)
 end subroutine

!************************************************************************************

 subroutine backtocell_cart(nat,latvec,rxyz)
 !This subroutine will transform back all atoms into the periodic cell
 !defined by the 3 lattice vectors in latvec=[v1.v2.v3]
 implicit none
 integer:: nat,i,iat,j
 real(8) :: latvec(3,3), rxyz(3,nat), crossp(3),a(3),b(3), nvec(3,3), dist(6),eps,count
 real(8) :: v(3,3),vol,rxyz_red(3,nat)
 logical:: neccesary
 call getvol(latvec,vol)
   !First check if the volume is positive
   if(vol.le.0.d0) stop "Negative volume during backtocell"
   call rxyz_cart2int(latvec,rxyz_red,rxyz,nat)
   do iat=1,nat
       rxyz_red(1,iat)=modulo(modulo(rxyz_red(1,iat),1.d0),1.d0)
       rxyz_red(2,iat)=modulo(modulo(rxyz_red(2,iat),1.d0),1.d0)
       rxyz_red(3,iat)=modulo(modulo(rxyz_red(3,iat),1.d0),1.d0)
   enddo
   call rxyz_int2cart(latvec,rxyz_red,rxyz,nat)
!!First check if the volume is positive
! v=latvec
! vol= v(1,1)*v(2,2)*v(3,3)-v(1,1)*v(2,3)*v(3,2)-v(1,2)*v(2,1)*v(3,3)+&
!      v(1,2)*v(2,3)*v(3,1)+v(1,3)*v(2,1)*v(3,2)-v(1,3)*v(2,2)*v(3,1)
! if(vol.le.0.d0) stop "Negative volume during backtocell"
!
! !To really be on the safe side, the translation vector can be shortened by  a factor eps in order
! !to get the atom into the cell. 
!!  eps=1.d-10
!  eps=1.d-15
! ! eps=1.d0-eps
! count=0.d0
! neccesary=.true.
! do while(neccesary)
! neccesary=.false.
! count=count+1.d0
! !generate 3 normal vectors of the 3 planes
! call nveclatvec(latvec,nvec)
! do iat=1,nat
! !3 planes through origin (xy,yz,zx)
! do i=1,3
! dist(i)=DOT_PRODUCT(rxyz(:,iat),nvec(:,i))
!! if(dist(i).lt.0.d0) then
! if(dist(i).lt.-abs(dist(i))*eps) then
!! write(*,*) "unten 1",i,iat
! rxyz(:,iat)=rxyz(:,iat)+latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
!
! !3 planes on top/side/back (xy,yz,zx)
! dist(i+3)=DOT_PRODUCT(rxyz(:,iat)-latvec(:,mod(i+1,3)+1),nvec(:,i))
!! if(dist(i+3).gt.0.d0) then
! if(dist(i+3).gt.abs(dist(i+3))*eps) then
!! write(*,*) "unten 1",i,iat
!! if(dist(i+3).gt.eps) then
!! write(*,*) "oben 2",i,iat
! rxyz(:,iat)=rxyz(:,iat)-latvec(:,mod(i+1,3)+1)!*eps
! neccesary=.true.
! endif
! enddo
! enddo
! if(count.gt.1.d6) stop "Too many iterations in back-to-cell"
! enddo
 end subroutine
