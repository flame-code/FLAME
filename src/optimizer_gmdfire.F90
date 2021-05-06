!*****************************************************************************************
subroutine gmdfire(nr,x,epot,f,work,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer:: nr
    real(8):: x(nr), epot, f(nr), de, DDOT, fnrm, fmax, vnrm, dt, p
    real(8):: tt, vnrmmax
    real(8):: work(5*nr) !1:nr velocities, nr+1:2*nr previous force
    real(8), save:: epotold, alpha
    integer, save:: ndown
    type(typ_paropt)::paropt
    if(paropt%iflag==0) then
        if(paropt%dt<-0.d0) then
            write(*,*) 'ERROR: time step in FIRE method must be set by user'
            return
        endif
        if(paropt%dtmax<0.d0) paropt%dtmax=20.d0*paropt%dt
        if(paropt%finc<0.d0) paropt%finc=1.1d0  !1.15d0
        if(paropt%fdec<0.d0) paropt%fdec=0.5d0  !0.2d0
        if(paropt%falpha<0.d0) paropt%falpha=0.99d0
        if(paropt%alphastart<0.d0) paropt%alphastart=0.1d0  !0.3d0
        if(paropt%alphax<0.d0) then
            write(*,*) 'ERROR: alphax in FIRE method must be set by user'
            return
        endif
        if(paropt%ndowntol<0) paropt%ndowntol=3
        !work(1:nr)=0.1d0*paropt%alphax*f(1:nr)
        work(1:nr)=0.d0
        work(nr+1:2*nr)=-f(1:nr)
        work(4*nr+1:5*nr)=-work(3*nr+1:4*nr)
        epotold=epot
        alpha=paropt%alphastart
        ndown=0
        paropt%iflag=1
    endif
    dt=paropt%dt
    if(ndown==0) work(nr+1:2*nr)=-f(1:nr)
    if(ndown==0) work(4*nr+1:5*nr)=-work(3*nr+1:4*nr)
    !work(1:nr)=work(1:nr)+0.5d0*dt*(f(1:nr)+work(nr+1:2*nr))
    work(1:nr)=work(1:nr)+0.5d0*dt*(work(3*nr+1:4*nr)+work(4*nr+1:5*nr))
    p=DDOT(nr,f,1,work,1)
    call calnorm(nr,work,vnrm)
    call calnorm(nr,f,fnrm)
    call calmaxforcecomponent(nr,f,fmax)
    de=epot-epotold
    write(*,'(a10,i5,es23.15,es11.3,2es12.5,3es12.4,i4,1es12.4)') &
        'FIREMIN   ',paropt%itfire,epot,de,fnrm,fmax,vnrm,dt,alpha,ndown,p
    write(21,'(a10,i5,es23.15,es11.3,2es12.5,3es12.4,i4,1es12.4)') &
        'FIREMIN   ',paropt%itfire,epot,de,fnrm,fmax,vnrm,dt,alpha,ndown,p
    if(fmax<paropt%fmaxtol) then
        paropt%iflag=0
        write(*,'(a,i4,es23.15,2es12.5)') &
            'FIRE FINISHED: itfire,epot,fnrm,fmax ',paropt%itfire,epot,fnrm,fmax
        return
    endif
    paropt%itfire=paropt%itfire+1
    epotold=epot
    work(2*nr+1:3*nr)=x(1:nr)
    !x(1:nr)=x(1:nr)+dt*work(1:nr)+0.5d0*dt**2*f(1:nr)
    x(1:nr)=x(1:nr)+dt*work(1:nr)+0.5d0*dt**2*work(3*nr+1:4*nr)
    !f(1:nr)=f(1:nr)/fnrm
    tt=min(alpha*vnrm/fnrm,5.d-2*paropt%alphax)
    work(1:nr)=(1.d0-alpha)*work(1:nr)+tt*f(1:nr) !*min(alpha*vnrm,5.d0*paropt%alphax*fnrm)
    !--------------------------------------
    !call calnorm(nr,work,tt)
    !!if(iproc==0) write(*,'(a,2es19.10)') 'fort56 ',tt,fnrm
    !vnrmmax=1.d-1*fnrm
    !if(tt>vnrmmax) then
    !    work(1:nr)=work(1:nr)*vnrmmax/tt
    !endif
    !--------------------------------------
    !if(p>0.d0) then
    if(.not. p<0.d0) then
        if(ndown>paropt%ndowntol) then
            paropt%dt=min(paropt%finc*paropt%dt,paropt%dtmax)
            !alpha=max(paropt%falpha*alpha,5.d-2) !10.d-2)
            alpha=max(paropt%falpha*alpha,1.d-10) !10.d-2)
            !alpha=paropt%falpha*alpha
        endif
        ndown=ndown+1
    else
        paropt%dt=paropt%fdec*paropt%dt
        x(1:nr)=work(2*nr+1:3*nr)
        f(1:nr)=work(nr+1:2*nr)
        work(3*nr+1:4*nr)=work(4*nr+1:5*nr)
        !work(1:nr)=0.1d0*paropt%alphax*f(1:nr)  !1.d-5*f(1:nr)
        work(1:nr)=0.d0
        alpha=paropt%alphastart
        ndown=0
        !work(nr+1:2*nr)=-f(1:nr)
    endif
    work(nr+1:2*nr)=f(1:nr)
    work(4*nr+1:5*nr)=work(3*nr+1:4*nr)
end subroutine gmdfire
!*****************************************************************************************
