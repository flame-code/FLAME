!*****************************************************************************************
!program erf_rinv_spline
!    implicit none
!    integer::nsp !number of spline nodes.
!    integer::npt !number of points for testing the spline.
!    real(8)::rcut !first and second cutoff for the function
!    real(8)::a !for the arguement of error function erf(r/(a*sqrt(2)))/r
!    !contains the coefficients for function
!    real(8), allocatable::fsp(:,:)
!    !contains the coefficients for derivative of function
!    real(8), allocatable::fdsp(:,:)
!    real(8)::r,hsp,hspinv,t,rhspinv,f,fd,spf,spfd
!    real(8)::errf,errfd,errfsp,errfspd,errfsprel,errfspdrel,r_withmaxerr
!    integer::ipt,isp
!    real(16)::fsd_tmp,hspq,aq,rq,fq,fdq
!    nsp=10**3
!    rcut=10.d0
!    a=2.d0
!    allocate(fdsp(0:3,0:nsp),fsp(0:4,0:nsp-1))
!    fdsp=0.d0;fsp=0.d0
!    call build_shortrange_spline(fdsp,fsp,nsp,rcut,a,hsp)
!    hspinv=1.d0/hsp
!    errfsp=0.d0;errfspd=0.d0;errf=0.d0;errfd=0.d0
!    npt=100*nsp
!    errfsprel=-1.d0
!    errfspdrel=-1.d0
!    r_withmaxerr=-1.d0
!    do ipt=1,npt
!        call random_number(r)
!        r=r*rcut
!        !if(r<hsp) cycle
!        !r=1.d-6
!        !write(*,*) 'r',r
!        !-----------------------------------------------------------------------
!        !following 5 command lines show how to use the spline.
!        rhspinv=r*hspinv
!        isp=floor(rhspinv) 
!        t=rhspinv-isp !t is canonical coordinates in [0,1)
!        spf=fsp(0,isp)+(fsp(1,isp)+(fsp(2,isp)+(fsp(3,isp)+fsp(4,isp)*t)*t)*t)*t
!        spfd=fdsp(0,isp)+(fdsp(1,isp)+(fdsp(2,isp)+fdsp(3,isp)*t)*t)*t
!        !-----------------------------------------------------------------------
!        !write(*,'(5e)') fsp(0,isp),fsp(1,isp),fsp(2,isp),fsp(3,isp),fsp(4,isp)
!        !-----------------------------------------------------------------------
!        rq=real(r,16)
!        aq=real(a,16)
!        hspq=real(hsp,16)
!        call func_funcder_funcsecder(rq,aq,hspq,fq,fdq,fsd_tmp)
!        f=real(fq,8)
!        fd=real(fdq,8)
!        !errfsp=errfsp+(f-spf)**2
!        !errf=errf+f**2
!        !errfspd=errfspd+(fd-spfd)**2
!        !errfd=errfd+fd**2
!        !errfsprel=max(errfsprel,abs(1.d0-spf/f))
!        !errfspdrel=max(errfspdrel,abs(1.d0-spfd/fd))
!        
!        if(abs(fd-spfd)>errfspdrel) then
!        !if(abs(f-spf)>errfsprel) then
!            errfsprel=abs(f-spf)
!            errfspdrel=abs(fd-spfd)
!            r_withmaxerr=r
!        endif
!        !errfsprel=max(errfsprel,abs(f-spf))
!        !errfspdrel=max(errfspdrel,abs(fd-spfd))
!    enddo
!    !errfsprel=sqrt(errfsp)/sqrt(errf)
!    !errfspdrel=sqrt(errfspd)/sqrt(errfd)
!    write(*,*) 'errfsprel=',errfsprel
!    write(*,*) 'errfspdrel=',errfspdrel
!    write(*,*) 'r_withmaxerr=',r_withmaxerr
!    deallocate(fdsp,fsp)
!end program erf_rinv_spline
!*****************************************************************************************
subroutine build_shortrange_spline(shortrange,spline,rcut,a)
    use mod_shortrange, only: typ_shortrange
    use mod_spline, only: typ_spline
    use mod_defs
    use yaml_output
    implicit none
    type(typ_shortrange), intent(in):: shortrange
    type(typ_spline), intent(inout):: spline
    real(8), intent(in):: rcut !first and second cutoff for the function
    real(8), intent(in):: a !prefacor in exponent of exponential function
    !local variables
    !following variables are quadruple precisions in order 
    !to have better accuracy.
    real(16):: fdspq(0:3,0:spline%nsp)
    real(16):: fspq(0:4,0:spline%nsp-1)
    real(16):: fdspq_1(0:3,0:spline%nsp)
    real(16):: fspq_1(0:4,0:spline%nsp-1)
    real(16):: fdspq_2(0:3,0:spline%nsp)
    real(16):: fspq_2(0:4,0:spline%nsp-1)
    real(16):: fdspq_3(0:3,0:spline%nsp)
    real(16):: fspq_3(0:4,0:spline%nsp-1)
    real(16):: fdspq_4(0:3,0:spline%nsp)
    real(16):: fspq_4(0:4,0:spline%nsp-1)
    real(16):: hspq 
    real(16):: aq, rcutq, qq, tt
    real(16):: a2, a3, a4, a5
    integer:: itypinter, isp
    external erf_over_r, one_over_r6, one_over_r8, exp_ar
    associate(nsp=>spline%nsp)
    rcutq=real(rcut,16)
    aq=real(a,16)
    hspq=rcutq/real(nsp,16)
    call build_spline(erf_over_r,rcutq,hspq,aq,nsp,fspq_1,fdspq_1)
    if(spline%do_tosifumi) then
        tt=1.0_fqp !will not be used
        call build_spline(one_over_r6,rcutq,hspq,tt,nsp,fspq_2,fdspq_2)
        call build_spline(one_over_r8,rcutq,hspq,tt,nsp,fspq_3,fdspq_3)
    endif
    do itypinter=1,shortrange%ntypinter
        qq=real(shortrange%qq(itypinter),16)
        fdspq(0:3,0:nsp)=qq*fdspq_1(0:3,0:nsp)
        fspq(0:4,0:nsp-1)=qq*fspq_1(0:4,0:nsp-1)
        if(spline%do_tosifumi) then
            a2=shortrange%tosifumi%ccc(itypinter)
            a3=shortrange%tosifumi%ddd(itypinter)
            a4=shortrange%tosifumi%aaa(itypinter)
            a5=shortrange%tosifumi%bbb(itypinter)
            fdspq_4=0.0_fqp
            fspq_4=0.0_fqp
            call build_spline(exp_ar,rcutq,hspq,a5,nsp,fspq_4,fdspq_4)
            fdspq(0:3,0:nsp)=fdspq(0:3,0:nsp)+a2*fdspq_2(0:3,0:nsp)+a3*fdspq_3(0:3,0:nsp)-a4*fdspq_4(0:3,0:nsp)
            fspq(0:4,0:nsp-1)=fspq(0:4,0:nsp-1)+a2*fspq_2(0:4,0:nsp-1)+a3*fspq_3(0:4,0:nsp-1)-a4*fspq_4(0:4,0:nsp-1)
            !fdspq(0:3,0:nsp)=fdspq(0:3,0:nsp) -a4*fdspq_4(0:3,0:nsp)
            !fspq(0:4,0:nsp-1)=fspq(0:4,0:nsp-1) -a4*fspq_4(0:4,0:nsp-1)
        endif
        spline%fdsp(0:3,0:nsp,itypinter)=real(fdspq(0:3,0:nsp),8)
        spline%fsp(0:4,0:nsp-1,itypinter)=real(fspq(0:4,0:nsp-1),8)
        !write(401,'(i5,2f15.10)') itypinter,a4,a5 
        !do isp=0,nsp-1
        !    write(300+itypinter,'(i5,5f15.10)') isp,spline%fsp(0,isp,itypinter),spline%fsp(1,isp,itypinter), &
        !        spline%fsp(2,isp,itypinter),spline%fsp(3,isp,itypinter),spline%fsp(4,isp,itypinter)
        !enddo
    enddo
    spline%hsp=real(hspq,8)
    call yaml_map('hsp',spline%hsp,fmt='(es22.14)')
    !write(*,*)'hsp=',spline%hsp
    end associate
end subroutine build_shortrange_spline 
            !fdspq(0:3,0:nsp)=fdspq(0:3,0:nsp)-a3*fdspq_3(0:3,0:nsp)
            !fspq(0:4,0:nsp-1)=fspq(0:4,0:nsp-1)-a3*fspq_3(0:4,0:nsp-1)
            !fdspq(0:3,0:nsp)=fdspq(0:3,0:nsp)+a4*fdspq_4(0:3,0:nsp)
            !fspq(0:4,0:nsp-1)=fspq(0:4,0:nsp-1)+a4*fspq_4(0:4,0:nsp-1)
!*****************************************************************************************
subroutine build_spline(cal_f_fd_fdd,rcutq,hspq,aq,nsp,fspq,fdspq)
    use mod_spline, only: typ_spline
    use mod_defs
    use yaml_output
    implicit none
    external:: cal_f_fd_fdd
    real(16), intent(in):: rcutq !first and second cutoff for the function
    real(16), intent(in):: hspq
    real(16), intent(in):: aq !prefacor in exponent of exponential function
    integer, intent(in):: nsp
    real(16), intent(out):: fdspq(0:3,0:nsp)
    real(16), intent(out):: fspq(0:4,0:nsp-1)
    !local variables
    !following variables are quadruple precisions in order 
    !to have better accuracy.
    real(16):: pi, onethird, fsp_int
    real(16):: del !shift to make short range part to zero, it is not the difference
    !between spline and the value with original function
    real(16):: r, fd0, fsd0, fd1, fsd1
    real(16):: f_tmp, fd_tmp, fsd_tmp, f_rcut
    integer:: isp
    pi=4.0_fqp*atan(1.0_fqp)
    onethird=1.0_fqp/3.0_fqp
    r=0.0_fqp
    call cal_f_fd_fdd(r,aq,hspq,f_tmp,fd0,fsd0)
    do isp=0,nsp
        r=hspq*real((isp+1),16)
        call cal_f_fd_fdd(r,aq,hspq,f_tmp,fd1,fsd1)
        !write(501,'(3es25.10)') r,fd1,fsd1
        fdspq(0,isp)=fd0
        fdspq(1,isp)=fsd0
        fdspq(2,isp)=3.0_fqp*(-fd0+fd1)-2.0_fqp*fsd0-fsd1
        fdspq(3,isp)=2.0_fqp*( fd0-fd1)+     fsd0+fsd1
        fd0=fd1
        fsd0=fsd1
    enddo
    !for the time being assume that the potential starts from this value
    fsp_int=1.0_fqp 
    do isp=0,nsp-1
        !the zeroth coefficient is the value of the integral
        fspq(0,isp)=fsp_int
        !the coefficients of the energy spline are obtained by integration
        fspq(1,isp)=hspq*fdspq(0,isp)
        fspq(2,isp)=hspq*fdspq(1,isp)*.5q0
        fspq(3,isp)=hspq*fdspq(2,isp)*onethird
        fspq(4,isp)=hspq*fdspq(3,isp)*.25q0
        !calculate the value of integral over the subinterval
        fsp_int=fsp_int+fspq(1,isp)+fspq(2,isp)+fspq(3,isp)+fspq(4,isp)
    enddo
    !correct the additive constant:
    !the value of the potential at rcut should be correct
    call cal_f_fd_fdd(rcutq,aq,hspq,f_rcut,fd_tmp,fsd_tmp)
    del=f_rcut-fsp_int
    call yaml_map('delta in build_spline',real(del,kind=8),fmt='(es22.14)')
    !write(*,*)'delta=',del
    do isp=0,nsp-1
        fspq(0,isp)=fspq(0,isp)+del
    enddo
end subroutine build_spline 
!*****************************************************************************************
subroutine erf_over_r(r,a,hsp,func,funcder,funcsecder)
    use mod_defs
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
    !local variables
    real(16):: pi, twosqrtinv
    !----------------------------------------------------------------
    !evaluating function at r
    !func=exp(-a*r)
    !evaluating first derivative of function at r
    !funcder=-a*func
    !evaluating second derivative of function at r
    !funcsecder=-a*funcder*hsp
    !----------------------------------------------------------------
    twosqrtinv=1.0_fqp/sqrt(2.0_fqp)
    pi=4.0_fqp*atan(1.0_fqp)
    if(r/a<1.e-10_fqp) then
        !evaluating function at r=0
        func=sqrt(2.0_fqp/pi)/a
        !evaluating first derivative of function at r=0
        funcder=0.0_fqp
        !evaluating second derivative of function at r=0
        funcsecder=-func/(3.0_fqp*a**2)*hsp
    else
        !evaluating function at r/=0
        func=erf(r*twosqrtinv/a)/r
        !evaluating first derivative of function at r/=0
        funcder=(sqrt(2.0_fqp/pi)*exp(-r**2/a**2*0.5q0)/a-func)/r
        !evaluating second derivative of function at r/=0
        funcsecder=(-sqrt(2.0_fqp/pi)*exp(-r**2/a**2*0.5q0)*(2.0_fqp*a**2+r*r)/a**3+2.0_fqp*func)/(r*r)*hsp
    endif
end subroutine erf_over_r
!*****************************************************************************************
!funcder=12.0_fqp*exp(-r12)/(sqrt(pi)*r)-(6.0_fqp*erf(r6))/r**7
!funcsecder=-144.0_fqp*exp(-r12)/(sqrt(pi)*r2)+exp(-r12)*(5.0_fqp*r4-144.0_fqp*r**16)/(sqrt(pi)*r6)+(42.0_fqp*erf(r6))/r8
subroutine one_over_r6(r,a,hsp,func,funcder,funcsecder)
    use mod_defs
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
    !local variables
    real(16):: pi, r2, r4, r6, r8, r12
    !----------------------------------------------------------------
    pi=4.0_fqp*atan(1.0_fqp)
    if(r<5.e-1_fqp) then
        !evaluating function at r=0
        func=2.0_fqp/sqrt(pi)
        !evaluating first derivative of function at r=0
        funcder=0.0_fqp
        !evaluating second derivative of function at r=0
        funcsecder=0.0_fqp
    else
        !evaluating function at r/=0
        r2=r**2
        r4=r2**2
        r6=r4*r2
        r8=r4**2
        r12=r6**2
        func=erf(r6)/r6
        !evaluating first derivative of function at r/=0
        funcder=12.0_fqp*exp(-r12)/(sqrt(pi)*r)-(6.0_fqp*erf(r6))/r**7
        !evaluating second derivative of function at r/=0
        funcsecder=(-144.0_fqp*exp(-r12)/(sqrt(pi)*r2)+ &
                   exp(-r12)*(5.0_fqp*r4-144.0_fqp*r**16)/(sqrt(pi)*r6)+(42.0_fqp*erf(r6))/r8)*hsp
    endif
end subroutine one_over_r6
!*****************************************************************************************
subroutine one_over_r8(r,a,hsp,func,funcder,funcsecder)
    use mod_defs
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
    !local variables
    real(16):: pi, r2, r4, r8, r9, r10, r16
    !----------------------------------------------------------------
    pi=4.0_fqp*atan(1.0_fqp)
    if(r<5.e-1_fqp) then
        !evaluating function at r=0
        func=2.0_fqp/sqrt(pi)
        !evaluating first derivative of function at r=0
        funcder=0.0_fqp
        !evaluating second derivative of function at r=0
        funcsecder=0.0_fqp
    else
        !evaluating function at r/=0
        r2=r**2
        r4=r2**2
        r8=r4**2
        r9=r8*r
        r10=r8*r2
        r16=r8**2
        func=erf(r8)/r8
        !evaluating first derivative of function at r/=0
        funcder=16.0_fqp*exp(-r16)/(sqrt(pi)*r)-(8.0_fqp*erf(r8))/r9
        !evaluating second derivative of function at r/=0
        funcsecder=((-144.0_fqp*exp(-r16)/r2-256.0_fqp/r16)/sqrt(pi)+(72.0_fqp*erf(r8))/r10)*hsp
    endif
end subroutine one_over_r8
!16.0_fqp*exp(-r16)/(sqrt(pi)*r)-(8.0_fqp*erf(r8))/r9
!(-144.0_fqp*exp(-r16)/r2-256.0_fqp/r16)/sqrt(pi)+(72*erf(r8))/r10
!*****************************************************************************************
subroutine exp_ar(r,a,hsp,func,funcder,funcsecder)
    use mod_defs
    implicit none 
    real(16), intent(in):: r
    real(16), intent(in):: a
    real(16), intent(in):: hsp
    real(16), intent(out):: func
    real(16), intent(out):: funcder
    real(16), intent(out):: funcsecder
    !local variables
    !----------------------------------------------------------------
    if(r<1.e-5_fqp) then
        !evaluating function at r=0
        func=1.0_fqp
        !evaluating first derivative of function at r=0
        funcder=-a
        !evaluating second derivative of function at r=0
        funcsecder=a**2*hsp
    else
        !evaluating function at r/=0
        func=exp(-a*r)
        !evaluating first derivative of function at r/=0
        funcder=-a*func
        !evaluating second derivative of function at r/=0
        funcsecder=a**2*func*hsp
    endif
        !write(601,'(3es25.15)') r,funcder,funcsecder
end subroutine exp_ar
!*****************************************************************************************
!function phi_fun(x)
!    ! computes the long range potential
!    implicit none
!    real(8)::phi_fun
!    phi_fun=erf(x*sqrt(.5d0))/x
!end function phi_fun
!!*****************************************************************************************
!function f(r)
!    !computes the exact force
!    implicit none
!    parameter(pi=3.14159265358979323846264338328d0)
!    rm=r*sqrt(.5d0)
!    f=( sqrt(2.d0/pi)*exp(-rm*rm)*r-erf(rm) )/(r*r)
!end function f
!!*****************************************************************************************
!function f1(r,hgrid)
!    ! retutns the exact force derivative
!    implicit none
!    parameter(pi=3.14159265358979323846264338328d0)
!    fac=sqrt(2.d0/pi)
!    rm=r*sqrt(.5d0)
!    f1=hgrid*(-fac*exp(-rm*rm)*r*(2.d0+r*r)+2.d0*erf(rm)   )/(r*r*r)
!end function f1
!*****************************************************************************************

