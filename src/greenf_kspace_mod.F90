!*****************************************************************************************
module mod_greenf_kspace
    use mod_spline_new, only: typ_spline_new
    implicit none
    private
    public:: typ_greenf_kspace

    type:: typ_greenf_kspace
        type(typ_spline_new):: spline_new
        integer:: npow
        real(8):: h, pi
        real(8):: c(0:10,2:6) !2<=npow<=6
        character(10), private:: str_method
        logical, private:: constructed=.false.
        contains
        procedure, public, pass(self):: init_greenf_kspace
        procedure, public, pass(self):: fini_greenf_kspace
        procedure, public, pass(self):: get_greenf_kspace_single
    end type typ_greenf_kspace
contains
!*****************************************************************************************
subroutine init_greenf_kspace(self,npow,str_method)
    implicit none
    class(typ_greenf_kspace), intent(inout):: self
    integer, intent(in):: npow
    character(*), intent(in):: str_method
    !local variables
    include 'fftw3.f'
    integer:: ip, np, i, n
    real(8), allocatable:: x(:), f(:)
    real(8), allocatable:: c(:)
    real(8):: pi, box, w, tt1
    integer(8):: planb
    if(npow<2) stop 'ERROR: I do not know if it works for npow<1'
    if(npow>6) stop 'ERROR: are you sure? this module not tested for npow>6'
    self%npow=npow
    pi=4.d0*atan(1.d0)
    self%pi=pi
    box=400*100.d0
    n=2000000
    self%h=2.d-2
    allocate(f(n))
    call dfftw_plan_dft_c2r_1d(planb,n,f(1),f(1),fftw_patient)
    do i=1,n
        w=2.d0*pi*((i-1)/2)/(real(n,kind=8)*self%h)
        !write(41,'(es14.5)') w
        if(i==1 .or. i==2) then
            f(i)=0.d0
        elseif(mod(i,2)==1) then
            tt1=w**npow
            if(tt1>100.d0) then
                f(i)=0.d0
            else
                f(i)=pi*w**(npow-1)*exp(-tt1)
            endif
        else
            f(i)=0.d0
        endif
    enddo
    call dfftw_execute(planb)
    do i=1,n
        f(i)=f(i)/box
    enddo
    np=n/50
    allocate(x(0:np))
    do ip=0,np
        x(ip)=real(ip,kind=8)*self%h
    enddo
    call self%spline_new%init_spline_new(np,x,'cubic')
    call self%spline_new%interpolate(f)
    self%str_method=str_method
    deallocate(f)
    deallocate(x)
    call dfftw_destroy_plan(planb)
    !-----------------------------------
    self%c( 0,2)=1.d0
    self%c( 1,2)=0.16666666666666667d0
    self%c( 2,2)=0.016666666666666667d0
    self%c( 3,2)=0.0011904761904761905d0
    self%c( 4,2)=0.000066137566137566138d0
    self%c( 5,2)=3.0062530062530063d-6
    self%c( 6,2)=1.1562511562511563d-7
    self%c( 7,2)=3.8541705208371875d-9
    self%c( 8,2)=1.133579564952114d-10
    self%c( 9,2)=2.9831041182950368d-12
    self%c(10,2)=7.1026288530834209d-14
    !-----------------------------------
    !-----------------------------------
    self%c( 0,4)=1.772453850905516d0
    self%c( 1,4)=0.16666666666666667d0
    self%c( 2,4)=0.0073852243787729834d0
    self%c( 3,4)=0.00019841269841269841d0
    self%c( 4,4)=3.6633057434389799d-6
    self%c( 5,4)=5.0104216770883438d-8
    self%c( 6,4)=5.3369838919565558d-10
    self%c( 7,4)=4.5882982390918899d-12
    self%c( 8,4)=3.2702107181106347d-14
    self%c( 9,4)=1.9729524591898391d-16
    self%c(10,4)=1.0245021046712515d-18
    !-----------------------------------
    !-----------------------------------
    self%c( 0,6)=2.6789385347077476d0
    self%c( 1,6)=0.2256863232377334d0
    self%c( 2,6)=0.0083333333333333333d0
    self%c( 3,6)=0.00017717847451770818d0
    self%c( 4,6)=2.4877240215799537d-6
    self%c( 5,6)=2.5052108385441719d-8
    self%c( 6,6)=1.9120529495565503d-10
    self%c( 7,6)=1.1505735105542392d-12
    self%c( 8,6)=5.6229145086910415d-15
    self%c( 9,6)=2.2838227524994151d-17
    self%c(10,6)=7.8530579618588898d-20
    !-----------------------------------
    self%constructed=.true.
end subroutine init_greenf_kspace
!*****************************************************************************************
subroutine fini_greenf_kspace(self)
    implicit none
    class(typ_greenf_kspace), intent(inout):: self
    !local variables
    call self%spline_new%fini_spline_new()
    self%constructed=.false.
end subroutine fini_greenf_kspace
!*****************************************************************************************
subroutine get_greenf_kspace_single(self,ak,sf,res)
    implicit none
    class(typ_greenf_kspace), intent(inout):: self
    real(8), intent(in):: ak, sf
    real(8), intent(out):: res
    !local variables
    integer:: mp
    real(8):: f, fd, fdd
    real(8):: t, t0, t1, t2, t3, t4, t5, t6, t7, t8, t9
    if(.not. self%constructed) stop 'ERROR: greenf_kspace not constructed'
    if(.not. (sf>0.d0)) stop 'ERROR: sf not greater than zero'
    if(ak<0.d0) stop 'ERROR: k<0.d0'
    t=ak/sf
    if(t>800.d0) stop 'ERROR: greenf_kspace does not work for k/sf>800'
    if(t>1.d0) then
        mp=floor(t/self%h)+1
        call self%spline_new%get_spline_new_single(mp,t,f,fd,fdd)
        res=4.d0*self%pi*(1.0-real(self%npow,kind=8)*f)/ak**2
    else !if(t>0.d0) then
        t0=t**2
        t1=t0*t0
        t2=t0*t1
        t3=t0*t2
        t4=t0*t3
        t5=t0*t4
        t6=t0*t5
        t7=t0*t6
        t8=t0*t7
        t9=t0*t8
        res= self%c( 0,self%npow)   &
            -self%c( 1,self%npow)*t0 &
            +self%c( 2,self%npow)*t1 &
            -self%c( 3,self%npow)*t2 &
            +self%c( 4,self%npow)*t3 &
            -self%c( 5,self%npow)*t4 &
            +self%c( 6,self%npow)*t5 &
            -self%c( 7,self%npow)*t6 &
            +self%c( 8,self%npow)*t7 &
            -self%c( 9,self%npow)*t8 &
            +self%c(10,self%npow)*t9
        res=res*4.d0*self%pi/(real(self%npow,kind=8)*sf**2)
    !else
    !    res=4.d0*self%pi*gamma(2.d0/real(self%npow,kind=8))/(real(self%npow,kind=8)*sf**2)
    endif
end subroutine get_greenf_kspace_single
!*****************************************************************************************
end module mod_greenf_kspace
!*****************************************************************************************
