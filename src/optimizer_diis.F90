!*****************************************************************************************
subroutine diisminimum(n,nr,x,epot,f,paropt,nwork,work)
    use mod_opt, only: typ_paropt
    implicit none
    integer:: n, nr, nwork, info, id, jd
    real(8):: x(n), f(n), epot, work(nwork), fnrm, dnrm2, ddot, fmax
    type(typ_paropt):: paropt
    real(8), save:: emin, fnrmlowest, epotold
    integer, save:: itdiis, ld, nd
    !local variables
    integer:: i
    real(8):: dx
    if(paropt%iflag==0) then
        paropt%iflag=1;itdiis=0;epotold=epot
        emin=1.d100;fnrmlowest=1.d100;ld=0;nd=0
    endif
    fnrm=dnrm2(nr,f,1)


    if(epot>emin+1.d-2*abs(emin) .or. fnrm>2.d0*fnrmlowest) then 
        write(*,'(a37,i8,2E10.2)') 'DIVERGENCE in DIIS, switch back to SD', &
            itdiis,(epot-emin)/abs(emin),fnrm/fnrmlowest
        call dcopy(nr,work((3*paropt%idsx+2)*nr+1),1,x,1)
        !call sdminimum(0,n,n,x,fnrmtol,f,epot,sdconverged)
        !emin=1.d100;fnrmlowest=1.d100;ld=0;nd=0;epotold=epot
        paropt%sdsaturated=.false.
        return
    endif



    nd=min(nd+1,paropt%idsx)
    ld=mod(ld,paropt%idsx)+1
    !call dcopy(n,x,1,xh(1,ld),1)
    call dcopy(nr,x,1,work(ld*nr+1),1)
    if(epot<emin) then 
        emin=epot;fnrmlowest=fnrm
        call dcopy(nr,x,1,work((3*paropt%idsx+2)*nr+1),1)
    endif
    call dcopy(nr,f,1,work((paropt%idsx+1)*nr+ld*nr+1),1)
    !call calmaxforcecomponent(nr,f,fmax)
    fmax=paropt%fmax
    write(*,'(a10,i4,e23.15,e11.3,2e12.5)') 'DIISMIN   ',itdiis,epot,epot-epotold,fnrm,fmax
    !if(fmax<paropt%fmaxtol) then
    if(paropt%converged) then
        write(*,'(a,i4,e23.15,2e12.5)') 'DIIS finished: ',itdiis,epot,fnrm,fmax
        paropt%iflag=0
        return
    endif
    call dcopy(nr,f,1,work((2*paropt%idsx+2)*nr+(ld-1)*nr+1),1)
    !set up DIIS matrix (upper triangle)
    if(itdiis>paropt%idsx-1) then !shift left up matrix
        do i=1,paropt%idsx-1;paropt%a(1:i,i,1)=paropt%a(2:i+1,i+1,1);enddo
    endif
    !calculate new line, use b as work array for summation
    do id=1,nd
        jd=mod(ld+id-1,nd)+1
        paropt%a(id,nd,1)=ddot(nr,work((2*paropt%idsx+2)*nr+(ld-1)*nr+1),1, &
        work((2*paropt%idsx+2)*nr+(jd-1)*nr+1),1)
    enddo
    do i=1,nd;paropt%a(i,i:nd,2)=paropt%a(i,i:nd,1);enddo !copy to work array
    paropt%a(1:nd,nd+1,2)=1.d0;paropt%a(nd+1,nd+1,2)=0.d0 !prepare boundary elements
    paropt%b(1:nd)=0.d0;paropt%b(nd+1)=1.d0 !prepare right hand side
    if(itdiis>0) then !solve linear system:(LAPACK)
        call DSYSV('U',nd+1,1,paropt%a(1,1,2),paropt%idsx+1,paropt%ipiv,paropt%b,paropt%idsx+1,paropt%a(1,1,3),(paropt%idsx+1)**2,info)
        if(info/=0) then;write(*,*) 'ERROR: DSYSV failed: info',info;stop;endif
    else
        paropt%b(1)=1.d0
    endif
    write(100,'(a,11(1pe11.2))')'DIIS weights',paropt%b(1:nd+1)
    !xh(1:nr,0)=0.d0;fh(1:nr,0)=0.d0 !new guess
    work(1:nr)=0.d0;work((paropt%idsx+1)*nr+1:(paropt%idsx+1)*nr+nr)=0.d0 !new guess
    do id=1,nd
        jd=mod(ld+id-1,nd)+1
        !xh(1:nr,0)=xh(1:nr,0)+b(id)*xh(1:nr,jd)
        work(1:nr)=work(1:nr)+paropt%b(id)*work(jd*nr+1:jd*nr+nr)  
        !fh(1:nr,0)=fh(1:nr,0)+b(id)*fh(1:nr,jd)
        work((paropt%idsx+1)*nr+1:(paropt%idsx+1)*nr+nr)=work((paropt%idsx+1)*nr+1:(paropt%idsx+1)*nr+nr)+&
            paropt%b(id)*work((paropt%idsx+1)*nr+jd*nr+1:(paropt%idsx+1)*nr+jd*nr+nr)
    enddo
    x(1:nr)=work(1:nr)+work((paropt%idsx+1)*nr+1:(paropt%idsx+1)*nr+nr)*paropt%alphax*2.d0
    !do i=1,nr
    !    dx=work((paropt%idsx+1)*nr+i)*paropt%alphax*1.d0
    !    dx=sign(min(5.d-2,abs(dx)),dx)
    !    x(i)=work(i)+dx
    !enddo
    epotold=epot
    itdiis=itdiis+1
end subroutine diisminimum
!*****************************************************************************************
