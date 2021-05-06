!*****************************************************************************************
subroutine cgminimum(iproc,n,nr,x,f,epot,paropt,nwork,work)
    use mod_opt, only: typ_paropt, frmt_base
    use yaml_output
    implicit none
    type(typ_paropt):: paropt
    integer, intent(in):: iproc, n, nr, nwork
    real(8):: x(n), f(n), epot, work(nwork)
    !local variables
    character(46), parameter:: frmt='('//frmt_base//',e12.4,a)'
    real(8):: fnrm, tt, y0, y1, DDOT, rlambda, t1, t2, t3, fmax, dx
    integer:: i
    t1=DDOT(nr,f,1,f,1)
    fnrm=sqrt(t1)
    if(paropt%iflag==0) call init_cgminimum(paropt,n,nr,f,nwork,work,epot,fnrm)
    if(paropt%dolinesearch) then
        !if((epot-paropt%epotitm1)>paropt%anoise) then
        if(fnrm>5.d0*paropt%fnrmitm1) then
            call yaml_sequence(advance='no')
            call yaml_mapping_open('back to SD')
            call yaml_map('iter',paropt%itcg,fmt='(i5)')
            call yaml_map('epot',epot,fmt='(es20.12)')
            call yaml_map('de',epot-paropt%epotitm1,fmt='(es20.12)')
            call yaml_mapping_close()
            !write(*,'(a35,i5,2e25.15)') 'back to SD:itcg,epot,epot-epotold', &
            !    paropt%itcg,epot,epot-paropt%epotitm1
            x(1:nr)=work(1:nr)
            paropt%sdsaturated=.false.
            paropt%iflag=0
            return
        endif
        t2=DDOT(nr,work(nr+1),1,f,1)
        t3=DDOT(nr,work(nr+1),1,work(nr+1),1)
        rlambda=(t1-t2)/t3 
        work(2*nr+1:2*nr+nr)=f(1:nr)+rlambda*work(2*nr+1:2*nr+nr)
        !call calmaxforcecomponent(nr,f,fmax)
        fmax=paropt%fmax
        call yaml_sequence(advance='no')
        call yaml_mapping_open('CG',flow=.true.)
        !call yaml_map('iproc',iproc,fmt='(i3.3)')
        call yaml_map('iter',paropt%itcg,fmt='(i5)')
        call yaml_map('epot',epot,fmt='(es20.12)')
        call yaml_map('de',epot-paropt%epotitm1,fmt='(es9.1)')
        call yaml_map('fmax',fmax,fmt='(es10.3)')
        call yaml_map('fnrm',fnrm,fmt='(es10.3)')
        call yaml_map('alpha/alphax',paropt%alpha/paropt%alphax,fmt='(e12.4)')
        !call yaml_map('method','SD',fmt='(a)')
        call yaml_mapping_close()
        !write(*,frmt) 'MIN:',iproc,paropt%itcg,epot, &
        !    epot-paropt%epotitm1,fnrm,fmax,paropt%alpha/paropt%alphax,' CG'
        !if(fmax<paropt%fmaxtol) then
        if(paropt%converged) then
            call yaml_sequence(advance='no')
            call yaml_mapping_open('CG FINISHED')
            call yaml_map('success',.true.)
            call yaml_map('iter',paropt%itcg,fmt='(i5)')
            call yaml_map('epot',epot,fmt='(es20.12)')
            call yaml_map('fnrm',fnrm,fmt='(es12.5)')
            call yaml_map('fmax',fmax,fmt='(es12.5)')
            call yaml_mapping_close()
            paropt%iflag=0
            return
        endif
        paropt%epotitm1=epot
        paropt%fnrmitm1=fnrm
        if(paropt%itcg==0) then
            paropt%alpha0=2.d0*paropt%alphax
        else
            paropt%alpha0=max(min(2.d0*paropt%alpha,4.d0*paropt%alphax),-paropt%alphax)
        endif
        if(paropt%itcg==paropt%nit) then
            !stop 'WARNING: I should do something here'
            !write(*,*) 'NO conv in CG after 500 its: switching back to SD',paropt%itcg,fnrm,epot
            x(1:nr)=work(1:nr)
            f(1:nr)=work(nr+1:nr+nr)
            !paropt%ifail=paropt%ifail+1
            !if(paropt%ifail==paropt%nfail) stop 'too many failures of cgminimum'
            paropt%iflag=0
            paropt%sdsaturated=.false.
            return
        endif
        work(1:nr)=x(1:nr)
        x(1:nr)=x(1:nr)+paropt%alpha0*work(2*nr+1:2*nr+nr)
        work(nr+1:nr+nr)=f(1:nr)
        paropt%dolinesearch=.false.
        paropt%itcg=paropt%itcg+1
        !return
    else
        y0=DDOT(nr,work(nr+1),1,work(2*nr+1),1)
        y1=DDOT(nr,f,1,work(2*nr+1),1)
        tt=y0/(y0-y1)
        !write(*,'(a,2(1x,e10.3),2x,e12.5)')  'y0,y1,y0/(y0-y1)',y0,y1,tt
        paropt%alpha=paropt%alpha0*max(min(tt,4.d0),-0.5d0)
        !paropt%alpha=paropt%alpha0*max(min(tt,2.d0),-0.25d0)
        !paropt%alpha=paropt%alpha0*max(min(tt,10.d0),-1.d0)
        !x(1:nr)=work(1:nr)+paropt%alpha*work(2*nr+1:2*nr+nr)
        do i=1,nr
            dx=paropt%alpha*work(2*nr+i)
            dx=sign(min(paropt%dxmax,abs(dx)),dx)
            x(i)=work(i)+dx
        enddo
        paropt%avgalpha=paropt%avgalpha+paropt%alpha/paropt%alphax
        paropt%avgnum=paropt%avgnum+1.d0
        paropt%dolinesearch=.true.
    endif !end of if for dolinesearch
    !write(*,*) 'average CG stepsize in terms of alphax',paropt%avgalpha/paropt%avgnum
end subroutine cgminimum
!*****************************************************************************************
subroutine init_cgminimum(paropt,n,nr,f,nwork,work,epot,fnrm)
    use mod_opt, only: typ_paropt
    use yaml_output
    implicit none
    type(typ_paropt), intent(inout):: paropt
    integer, intent(in):: n, nr, nwork
    real(8), intent(in):: f(n), epot, fnrm
    real(8), intent(inout):: work(nwork)
    paropt%iflag=1
    paropt%itcg=0
    paropt%dolinesearch=.true.
    paropt%alpha=0.d0
    paropt%epotitm1=epot
    paropt%fnrmitm1=fnrm
    !paropt%ifail=0
    paropt%avgalpha=0.d0
    paropt%avgnum=0.d0
    if(paropt%fmaxtol<0.d0) stop 'ERROR: fmaxtol<0, it is not set.'
    if(paropt%alphax<0.d0) stop 'ERROR: alphax<0, it is not set.'
    if(paropt%dxmax<0.d0) paropt%dxmax=0.1d0
    if(paropt%nit<0) paropt%nit=1000
    !if(paropt%nfail<0) paropt%nfail=10
    if(paropt%alpha0<0.d0) paropt%alpha0=2.d0*paropt%alphax
    work(nr+1:nr+nr)=f(1:nr)
    !rlambda=0.d0
    if(paropt%lprint .and. .not. paropt%param_reported) call report_param(paropt)
    call yaml_sequence_open('CG optimization iterations')
end subroutine init_cgminimum
!*****************************************************************************************
