!*****************************************************************************************
subroutine mydfp(nr,x,epot,f,nwork,work,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer::nr,nwork,mf,my,ms,nrsq,iw1,iw2,iw3,info,i,j,l,mx
    real(8)::x(nr),f(nr),epot,work(nwork)
    integer, allocatable::ipiv(:)
    type(typ_paropt)::paropt
    real(8)::DDOT,DNRM2,tt1,tt2,de,fnrm,fmax
    real(8), save::epotold,alpha,alphamax
    logical, save::reset
    if(nwork/=nr*nr+3*nr+2*nr*nr+2*nr) then
        stop 'ERROR: size of work array is insufficient.'
    endif
    mf=nr*nr+1 !for force of previous iteration in wiki notation
    my=nr*nr+nr+1 !for y_k in wiki notation
    ms=nr*nr+2*nr+1 !for s_k in wiki notation
    iw1=nr*nr+3*nr+1 !work array to keep the hessian untouched
    iw2=nr*nr+3*nr+nr*nr+1 !for work array of DSYTRF
    iw3=nr*nr+3*nr+2*nr*nr+1 !for p_k in wiki notation
    mx=nr*nr+3*nr+2*nr*nr+nr+1 !for position of previous iteration
    nrsq=nr*nr
    if(paropt%iflag==0) then
        paropt%iflag=1
        paropt%iter=0
        epotold=epot
        alpha=10.d-1
        reset=.false.
        alphamax=1.d0
    else
        paropt%iter=paropt%iter+1
    endif
    de=epot-epotold
    call calnorm(nr,f,fnrm)
    call calmaxforcecomponent(nr,f,fmax)
    write(*,'(a10,i4,es23.15,es11.3,2es12.5,1es12.4)') &
        'BFGSMIN   ',paropt%iter,epot,de,fnrm,fmax,alpha
    !if(paropt%iter==1714) then
    !    do i=1,nr/3
    !        write(31,*) x(i*3-2),x(i*3-1),x(i*3-0)
    !    enddo
    !    stop
    !endif
    if(fmax<paropt%fmaxtol) then
        paropt%iflag=0
        write(*,'(a,i4,es23.15,2es12.5)') &
            'BFGS FINISHED: itfire,epot,fnrm,fmax ',paropt%iter,epot,fnrm,fmax
        return
    endif

    if(de>0.d0) then
        epot=epotold
        x(1:nr)=work(mx:mx-1+nr)
        f(1:nr)=work(mf:mf-1+nr)
        reset=.true.
        alpha=max(alpha*0.5d0/1.1d0,1.d-2)
    endif
    if(paropt%iter==0 .or. reset) then
        reset=.false.
        work(1:nr*nr)=0.d0
        do i=1,nr
            work(i+(i-1)*nr)=1.d0/paropt%alphax
        enddo
        work(iw3:iw3-1+nr)=paropt%alphax*f(1:nr)
    else
        work(ms:ms-1+nr)=x(1:nr)-work(mx:mx-1+nr)
        work(my:my-1+nr)=work(mf:mf-1+nr)-f(1:nr)
        write(21,*) paropt%iter,DNRM2(nr,work(my),1)
        tt1=DDOT(nr,work(my),1,work(ms),1)
        do i=1,nr
            tt2=0.d0
            do j=1,nr
                tt2=tt2+work(i+(j-1)*nr)*work(ms-1+j)
            enddo
            work(iw2-1+i)=tt2
        enddo
        tt2=DDOT(nr,work(ms),1,work(iw2),1)
        do i=1,nr
            do j=i,nr
                l=i+(j-1)*nr
                work(l)=work(l)+work(my-1+i)*work(my-1+j)/tt1-work(iw2-1+i)*work(iw2-1+j)/tt2
                work(j+(i-1)*nr)=work(l)
            enddo
        enddo
        allocate(ipiv(nr))
        work(iw1:iw1-1+nr*nr)=work(1:nr*nr)
        work(iw3:iw3-1+nr)=f(1:nr)
        !http://alcinoe.net/fortran/optim/optim.f90.html
        !http://www.netlib.no/netlib/lapack/double/dsytrf.f
        !http://www.netlib.no/netlib/lapack/double/dsytrs.f
        call DSYTRF('L',nr,work(iw1),nr,ipiv,work(iw2),nrsq,info)
        if(info/=0) then;write(*,*) 'ERROR: DSYTRF failed: info',info;stop;endif
        call DSYTRS('L',nr,1,work(iw1),nr,ipiv,work(iw3),nr,info)
        if(info/=0) then;write(*,*) 'ERROR: DSYTRS failed: info',info;stop;endif
        deallocate(ipiv)
    endif
    epotold=epot
    work(mf:mf-1+nr)=f(1:nr)
    work(mx:mx-1+nr)=x(1:nr)
    alpha=min(alphamax,alpha*1.1d0)
    x(1:nr)=x(1:nr)+alpha*work(iw3:iw3-1+nr)
end subroutine mydfp
!*****************************************************************************************
