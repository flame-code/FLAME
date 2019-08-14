!*****************************************************************************************
subroutine fire(parini,iproc,n,x,epot,f,work,paropt)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt, frmt_base
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: n, iproc
    real(8), intent(inout):: x(n), epot, f(n)
    real(8), intent(inout):: work(3*n) !1:n velocities, n+1:2*n previous force
    type(typ_paropt), intent(inout):: paropt
    !local variables
    character(59), parameter:: frmt='('//frmt_base//',3es12.4,i4,1es12.4,a)'
    real(8):: de, DDOT, fnrm, fmax, vnrm, dt, p
    real(8):: tt, vnrmmax
    character(14):: filename
    character(56):: comment
    if(paropt%iflag==0) call init_fire(n,f,epot,work,paropt)
    dt=paropt%dt
    if(paropt%ndown==0) work(n+1:2*n)=-f(1:n)
    work(1:n)=work(1:n)+0.5d0*dt*(f(1:n)+work(n+1:2*n))
    p=DDOT(n,f,1,work,1)
    call calnorm(n,work,vnrm)
    call calnorm(n,f,fnrm)
    !call calmaxforcecomponent(n,f,fmax)
    fmax=paropt%fmax
    de=epot-paropt%epotold
    !write(*,'(a10,i5,es23.15,es11.3,2es12.5,3es12.4,i4,1es12.4)') &
    !    'FIREMIN   ',paropt%itfire,epot,de,fnrm,fmax,vnrm,dt,paropt%alpha,paropt%ndown,p
    if(paropt%lprint .and. parini%iverbose>=0) then
        !write(*,'(a4,i3.3,1x,i5,es23.15,es11.3,2es12.5,3es12.4,i4,1es12.4,a)') &
        !write(*,frmt) &
        !    'MIN:',iproc,paropt%itfire,epot,de,fnrm,fmax,vnrm,dt,paropt%alpha,paropt%ndown,p,' FIRE'
        call yaml_sequence(advance='no')
        call yaml_mapping_open('FIRE',flow=.true.)
        call yaml_map('iter',paropt%itfire,fmt='(i5)')
        call yaml_map('epot',epot,fmt='(es20.12)')
        call yaml_map('de',de,fmt='(es9.1)')
        call yaml_map('fmax',fmax,fmt='(es10.3)')
        call yaml_map('fnrm',fnrm,fmt='(es10.3)')
        call yaml_map('vnrm',vnrm,fmt='(es10.3)')
        call yaml_map('alpha',paropt%alpha,fmt='(e12.4)')
        call yaml_map('dt',dt,fmt='(e12.4)')
        call yaml_map('ndown',paropt%ndown,fmt='(i5)')
        call yaml_map('power',p,fmt='(e12.4)')
        call yaml_mapping_close()
    endif
    !write(21,'(a10,i5,es23.15,es11.3,2es12.5,3es12.4,i4,1es12.4)') &
    !    'FIREMIN   ',paropt%itfire,epot,de,fnrm,fmax,vnrm,dt,paropt%alpha,paropt%ndown,p
    !if(fmax<paropt%fmaxtol) then
    if(paropt%converged) then
        paropt%iflag=0
        !write(*,'(a,i4,es23.15,2es12.5)') &
        !    'FIRE FINISHED: itfire,epot,fnrm,fmax ',paropt%itfire,epot,fnrm,fmax
        call yaml_sequence(advance='no')
        call yaml_mapping_open('FIRE FINISHED') !,label='id001')
        call yaml_map('success',.true.)
        call yaml_map('iter',paropt%itfire,fmt='(i5)')
        call yaml_map('epot',epot,fmt='(es20.12)')
        call yaml_map('fnrm',fnrm,fmt='(es12.5)')
        call yaml_map('fmax',fmax,fmt='(es12.5)')
        call yaml_mapping_close()
        return
    endif
    paropt%itfire=paropt%itfire+1
    if(paropt%itfire==paropt%nit) then
        paropt%iflag=0
        write(*,'(a)') 'FIRE: NO CONVERGENCE '
        return
    endif
    paropt%epotold=epot
    work(2*n+1:3*n)=x(1:n)
    x(1:n)=x(1:n)+dt*work(1:n)+0.5d0*dt**2*f(1:n)
    !f(1:n)=f(1:n)/fnrm
    tt=min(paropt%alpha*vnrm/(paropt%funits*fnrm),2.d-1*paropt%alphax)
    work(1:n)=(1.d0-paropt%alpha)*work(1:n)+tt*f(1:n) !*min(paropt%alpha*vnrm,5.d0*paropt%alphax*fnrm)
    !--------------------------------------
    !call calnorm(n,work,tt)
    !!if(iproc==0) write(*,'(a,2es19.10)') 'fort56 ',tt,fnrm
    !vnrmmax=1.d-1*fnrm
    !if(tt>vnrmmax) then
    !    work(1:n)=work(1:n)*vnrmmax/tt
    !endif
    !--------------------------------------
    !if(p>0.d0) then
    if(.not. p<0.d0) then
        if(paropt%ndown>paropt%ndowntol) then
            paropt%dt=min(paropt%finc*paropt%dt,paropt%dtmax)
            !paropt%alpha=max(paropt%falpha*paropt%alpha,5.d-2) !10.d-2)
            paropt%alpha=max(paropt%falpha*paropt%alpha,1.d-10) !10.d-2)
            !paropt%alpha=paropt%falpha*paropt%alpha
        endif
        paropt%ndown=paropt%ndown+1
    else
        paropt%dt=paropt%fdec*paropt%dt
        x(1:n)=work(2*n+1:3*n)
        f(1:n)=work(n+1:2*n)
        !work(1:n)=0.1d0*paropt%alphax*f(1:n)  !1.d-5*f(1:n)
        work(1:n)=0.d0
        paropt%alpha=paropt%alphastart
        paropt%ndown=0
        !work(n+1:2*n)=-f(1:n)
    endif
    work(n+1:2*n)=f(1:n)
end subroutine fire
!*****************************************************************************************
subroutine init_fire(n,f,epot,work,paropt)
    use mod_opt, only: typ_paropt
    use yaml_output
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: epot, f(n)
    real(8), intent(inout):: work(3*n)
    type(typ_paropt), intent(inout):: paropt
    !local variables
    if(paropt%dt_start<-0.d0) then
        paropt%dt_start=5.d-3
        !write(*,*) 'ERROR: time step in FIRE method must be set by user'
        !return
    endif
    paropt%dt=paropt%dt_start
    if(paropt%dtmax<0.d0) paropt%dtmax=25.d0*paropt%dt_start
    if(paropt%finc<0.d0) paropt%finc=1.1d0  !1.15d0
    if(paropt%fdec<0.d0) paropt%fdec=0.5d0  !0.2d0
    if(paropt%falpha<0.d0) paropt%falpha=0.99d0
    if(paropt%alphastart<0.d0) paropt%alphastart=0.1d0  !0.3d0
    if(paropt%funits<0.d0) paropt%funits=1.d0
    if(paropt%alphax<0.d0) then
        write(*,*) 'ERROR: alphax in FIRE method must be set by user'
        return
    endif
    if(paropt%ndowntol<0) paropt%ndowntol=3
    !work(1:n)=0.1d0*paropt%alphax*f(1:n)
    work(1:n)=0.d0
    work(n+1:2*n)=-f(1:n)
    paropt%epotold=epot
    paropt%alpha=paropt%alphastart
    paropt%ndown=0
    paropt%iflag=1
    paropt%itfire=0
    if(paropt%nit<0) paropt%nit=1000
    if(paropt%lprint .and. .not. paropt%param_reported) call report_param(paropt)
    call yaml_sequence_open('FIRE optimization iterations')
end subroutine init_fire
!*****************************************************************************************
