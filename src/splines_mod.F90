!*****************************************************************************************
module mod_spline_new
    implicit none
    private
    public:: typ_spline_new

    type:: typ_spline_new
        integer:: np
        character(10), private:: str_method
        character(30), private:: str_stat='unknown'
        logical, private:: constructed=.false.
        real(8), private, allocatable:: s(:)
        real(8), private, allocatable:: h(:)
        real(8), private, allocatable:: e1(:)
        real(8), private, allocatable:: e2(:)
        real(8), private, allocatable:: c(:)
        real(8), private, allocatable:: y(:)
        contains
        procedure, public, pass(self):: init_spline_new
        procedure, public, pass(self):: fini_spline_new
        procedure, public, pass(self):: interpolate
        procedure, public, pass(self):: get_spline_new
        procedure, public, pass(self):: get_spline_new_single
        procedure, public, pass(self):: get_coeff
        procedure, private, nopass:: factor_cubic_second
        procedure, private, nopass:: inter_cubic_second
        procedure, private, nopass:: ffdfdd_cubic
    end type typ_spline_new

    !interface typ_spline_new
    !    module procedure:: init_spline_new
    !end interface typ_spline_new
contains
!*****************************************************************************************
subroutine init_spline_new(self,np,s,str_method)
    implicit none
    class(typ_spline_new), intent(inout):: self
    integer, intent(in):: np
    real(8), intent(in):: s(0:np)
    character(*), intent(in):: str_method
    !local variables
    integer:: ip
    self%np=np
    self%str_method=str_method
    !if(trim(self%str_method)=='old') then
    allocate(self%h(np),self%e1(np-1),self%e2(np-2),self%c(0:np))
    allocate(self%s(0:np),self%y(0:np))
    self%s(0:np)=s(0:np)
    do ip=1,np
        self%h(ip)=s(ip)-s(ip-1)
    enddo
    call factor_cubic_second(np,self%h,self%e1,self%e2)
    self%constructed=.true.
end subroutine init_spline_new
!*****************************************************************************************
subroutine fini_spline_new(self)
    implicit none
    class(typ_spline_new), intent(inout):: self
    !local variables
    deallocate(self%h,self%e1,self%e2,self%c)
    deallocate(self%s,self%y)
    self%str_stat='unknown'
    self%constructed=.false.
end subroutine fini_spline_new
!*****************************************************************************************
subroutine interpolate(self,y)
    implicit none
    class(typ_spline_new), intent(inout):: self
    real(8), intent(in):: y(0:self%np)
    !local variables
    if(.not. self%constructed) stop 'ERROR: spline_new not initialized'
    self%y(0:self%np)=y(0:self%np)
    call inter_cubic_second(self%np,y,self%h,self%e1,self%e2,self%c)
    self%str_stat='interpolated'
end subroutine interpolate
!*****************************************************************************************
subroutine get_spline_new(self)
    implicit none
    class(typ_spline_new), intent(inout):: self
    !local variables
    stop 'ERROR: get_spline_new not ready yet!'
end subroutine get_spline_new
!*****************************************************************************************
subroutine get_spline_new_single(self,mp,t,f,fd,fdd)
    implicit none
    class(typ_spline_new), intent(inout):: self
    integer, intent(in):: mp
    real(8), intent(in):: t
    real(8), intent(out):: f, fd, fdd
    !local variables
    if(trim(self%str_stat)/='interpolated') stop 'ERROR: spline_new not interpolated'
    if(mp<1) stop 'ERROR: mp<1 in get_spline_new_single'
    if(mp>self%np) stop 'ERROR: mp>self%np in get_spline_new_single'
    call ffdfdd_cubic(self%np,self%y,self%s,mp,self%h(mp),t,self%c,f,fd,fdd)
end subroutine get_spline_new_single
!*****************************************************************************************
subroutine get_coeff(self,c)
    implicit none
    class(typ_spline_new), intent(inout):: self
    real(8), intent(out):: c(0:self%np)
    !local variables
    c(0:self%np)=self%c(0:self%np)
end subroutine get_coeff
!*****************************************************************************************
subroutine factor_cubic_second(np,h,e1,e2)
    implicit none
    integer::np
    integer::ip,info
    real(8)::h(np),e1(np-1),e2(np-2)
    do ip=1,np-2;e1(ip)=2.d0*(h(ip+1)+h(ip));e2(ip)=h(ip+1);enddo
    e1(np-1)=2.d0*(h(np)+h(np-1))+h(np)
    e1(1)=e1(1)+h(1)
    call dpttrf(np-1,e1,e2,info)
    if(info/=0) write(*,*) 'ERROR: factorization failed: info,np',info,np
end subroutine factor_cubic_second
!*****************************************************************************************
subroutine inter_cubic_second(np,y,h,e1,e2,c)
    implicit none
    integer::np
    integer::ip,info
    real(8)::y(0:np),h(np),e1(np-1),e2(np-2),c(0:np)
    !real(8)::tt,dt,b(4),ipiv(4),hi,bt0,btn,p0,p1,p2,p3
    do ip=1,np-1;c(ip)=(y(ip+1)-y(ip))/h(ip+1)-(y(ip)-y(ip-1))/h(ip);enddo
    !write(*,*) c(1:np-1)
    call dpttrs(np-1,1,e1,e2,c(1),np-1,info)
    !write(*,*) e1(1:np-1)
    !stop
    if(info/=0) write(*,*) 'ERROR: solution of dpttrs failed: info',info
    c(0)=c(1);c(np)=c(np-1)
end subroutine inter_cubic_second
!*****************************************************************************************
subroutine ffdfdd_cubic(np,y,s,mp,hmp,t,c,f,fd,fdd)
    implicit none
    integer::np,mp
    real(8)::y(0:np),s(0:np),hmp,c(0:np),t,p0,p1,p2,p3,f,fd,fdd
    if(mp<1 .or. mp>np) stop 'ERROR: invalid mp in cubic evaluation'
    !hmp=s(mp)-s(mp-1)
    p3=(c(mp)-c(mp-1))/hmp
    p2=3.d0*(c(mp-1)*s(mp)-c(mp)*s(mp-1))/hmp
    p1=(3.d0*(c(mp)*s(mp-1)**2-c(mp-1)*s(mp)**2)+y(mp)-y(mp-1))/hmp+hmp*(c(mp-1)-c(mp))
    p0=(c(mp-1)*s(mp)**3-c(mp)*s(mp-1)**3-y(mp)*s(mp-1)+y(mp-1)*s(mp))/hmp+ &
        hmp*(c(mp)*s(mp-1)-c(mp-1)*s(mp))
    f=((p3*t+p2)*t+p1)*t+p0
    fd=(3.d0*p3*t+2.d0*p2)*t+p1
    fdd=6.d0*p3*t+2.d0*p2
    !write(*,'(4es14.5)') t,p3,p2,fdd
end subroutine ffdfdd_cubic
!*****************************************************************************************
end module mod_spline_new
!*****************************************************************************************
