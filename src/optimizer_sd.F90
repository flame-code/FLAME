!*****************************************************************************************
subroutine sdminimum(parini,iproc,nr,x,f,epot,paropt,nwork,work)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt, frmt_base
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, nr, nwork
    real(8), intent(inout):: x(nr), f(nr), work(nwork)
    real(8), intent(in):: epot
    type(typ_paropt), intent(inout):: paropt
    !local variables
    character(52), parameter:: frmt='('//frmt_base//',e12.4,i5,l2,a)'
    real(8):: fmax, fnrm, de1, de2, df1, df2, dx, tt1
    logical:: feedbackcondition
    integer:: i
    if(paropt%iflag==0) call init_sdminimum(paropt,nr,x,nwork,work)
    call calnorm(nr,f,fnrm)
    !call calmaxforcecomponent(nr,f,fmax)
    fmax=paropt%fmax
    de1=epot-paropt%epotitm1;de2=epot-2.d0*paropt%epotitm1+paropt%epotitm2
    df1=fnrm-paropt%fnrmitm1;df2=fnrm-2.d0*paropt%fnrmitm1+paropt%fnrmitm2
    if(paropt%itsd==0) de1=0.d0
    paropt%xmoved=.true.
    call what_is_condition_of_feedback(paropt,de1,df1,feedbackcondition)
    !if(care .and. de1>paropt%anoise) then
    if(paropt%care .and. feedbackcondition) then
        if(paropt%alpha<paropt%alphamin .and. paropt%itsd/=0) then
            write(*,'(a)') 'alpha getting too small, do not care anymore if energy goes up'
            paropt%care=.false.
         else
            x(1:nr)=work(1:nr)
            f(1:nr)=work(nr+1:nr+nr)
            paropt%xmoved=.false.
        endif
    endif
    !write(*,'(5(1x,e11.4),1x,i3)') fnrm/fnrmitm1, de1,de2,df1,df2,isatur
    !if(care .and. paropt%itsd>5 .and. alpha==alphax .and. fnrm/fnrmitm1>0.8d0 &
    call test_saturation(paropt,de1,de2,df2,fnrm)
    if(paropt%lprint .and. parini%iverbose>=0) then
        tt1=0.5d0*(sign(1.d0,real(paropt%itsd-1,8))+1.d0)*paropt%alpha/paropt%alphax
        !write(*,'(a4,i3.3,1x,i5,es23.15,e11.3,2e12.5,e12.4,i5,l2,a)') 'MIN:',iproc,paropt%itsd,epot,de1, &
        !a4,i3.3,1x,i5,es20.12,es11.3,2es12.5 //,e12.4,i5,l2,a
        call yaml_sequence(advance='no')
        call yaml_mapping_open('SD',flow=.true.)
        !call yaml_map('iproc',iproc,fmt='(i3.3)')
        call yaml_map('iter',paropt%itsd,fmt='(i5)')
        call yaml_map('epot',epot,fmt='(es20.12)')
        call yaml_map('de',de1,fmt='(es9.1)')
        call yaml_map('fmax',fmax,fmt='(es10.3)')
        call yaml_map('fnrm',fnrm,fmt='(es10.3)')
        call yaml_map('alpha/alphax',tt1,fmt='(e12.4)')
        call yaml_map('isatur',paropt%isatur,fmt='(i5)')
        call yaml_map('xmoved',paropt%xmoved,fmt='(l2)')
        !call yaml_map('method','SD',fmt='(a)')
        call yaml_mapping_close()
        !write(*,frmt) 'MIN:',iproc,paropt%itsd,epot,de1, &
        !fnrm,fmax,tt1,paropt%isatur,paropt%xmoved,' SD'
    endif
    !if(paropt%itsd==126) then
    !close(74)
    !write(74,'(i5)') nr/3
    !do i=1,nr,3
    !    write(74,'(3es23.13)') x(i),x(i+1),x(i+2)
    !enddo
    !write(*,'(3es23.13)') x(1),x(2),x(3)
    !endif
    !if(fmax<paropt%fmaxtol) then
    if(paropt%converged) then
        !paropt%iflag=0
        !paropt%alpha=-1.d0
        !write(*,'(a,i6,e23.15,2e12.5)') &
        !    'SD FINISHED: itsd,epot,fnrm,fmax',paropt%itsd,epot,fnrm,fmax
        call yaml_sequence(advance='no')
        call yaml_mapping_open('SD FINISHED') !,label='id001')
        call yaml_map('success',.true.)
        call yaml_map('iter',paropt%itsd,fmt='(i5)')
        call yaml_map('epot',epot,fmt='(es20.12)')
        call yaml_map('fnrm',fnrm,fmt='(es12.5)')
        call yaml_map('fmax',fmax,fmt='(es12.5)')
        call yaml_mapping_close()
        call final_sdminimum(paropt)
        return
    endif
    if(paropt%isatur>paropt%nsatur .and. .not. paropt%sdsaturated &
        .and. trim(paropt%sdsaturation)=='yes') then
        call yaml_sequence(advance='no')
        call yaml_mapping_open('SD SATURATED') !,label='id001')
        call yaml_map('iter',paropt%itsd,fmt='(i5)')
        call yaml_map('epot',epot,fmt='(es20.12)')
        call yaml_map('fnrm',fnrm,fmt='(es12.5)')
        call yaml_map('fmax',fmax,fmt='(es12.5)')
        call yaml_mapping_close()
        call final_sdminimum(paropt)
        !write(*,'(a,i6,e23.15,2e12.5)') &
        !    'SD SATURATED: itsd,epot,fnrm,fmax',paropt%itsd,epot,fnrm,fmax
        paropt%sdsaturated=.true.
        !paropt%alpha=-1.d0
        !paropt%itsd=0
        if(trim(paropt%approach)/='SD') then
            paropt%iflag=0
            return
        endif
        call final_sdminimum(paropt)
        return
    endif
    !call what_is_condition_of_feedback(paropt,de1,df1,feedbackcondition)
    !if(care .and. de1>paropt%anoise) then
    if(paropt%care .and. feedbackcondition) then
        paropt%alpha=5.d-1*paropt%alpha
    endif
    if(paropt%xmoved) then
        paropt%epotitm2=paropt%epotitm1;paropt%epotitm1=epot
        paropt%fnrmitm2=paropt%fnrmitm1;paropt%fnrmitm1=fnrm
        paropt%alpha=min(1.2d0*paropt%alpha,paropt%alphamax)
        work(1:nr)=x(1:nr)
        work(nr+1:nr+nr)=f(1:nr)
    endif
    !x(1:nr)=x(1:nr)+paropt%alpha*f(1:nr)
    do i=1,nr
        dx=paropt%alpha*f(i)
        dx=sign(min(paropt%dxmax,abs(dx)),dx)
        x(i)=x(i)+dx
    enddo
    if(.not. paropt%care .and. paropt%alpha>2.d0*paropt%alphamin) then
        paropt%care=.true.
        write(*,'(a)') 'sdminimum starts to care whether energy goes up' 
    endif
    paropt%itsd=paropt%itsd+1
    if(paropt%itsd>=paropt%nit) then 
        paropt%iflag=-1
        x(1:nr)=work(1:nr)
        f(1:nr)=work(nr+1:nr+nr)
        write(*,'(a,e23.15,e12.5)') 'SD: NO CONVERGENCE: fnrm,epot ',epot,fmax
        call final_sdminimum(paropt)
    endif
end subroutine sdminimum
!*****************************************************************************************
subroutine init_sdminimum(paropt,nr,x,nwork,work)
    use mod_opt, only: typ_paropt
    use yaml_output
    implicit none
    type(typ_paropt), intent(inout):: paropt
    integer, intent(in):: nr, nwork
    real(8), intent(in):: x(nr)
    real(8), intent(inout):: work(nwork)
    !local variables
    !write(*,*) len(trim(paropt%approach))
    if(paropt%fmaxtol<0.d0) stop 'ERROR: fmaxtol<0, it is not set.'
    if(paropt%alphax<0.d0) stop 'ERROR: alphax<0, it is not set.'
    paropt%iflag=1
    paropt%itsd=0
    if(paropt%anoise<0.d0) then
        paropt%anoise=epsilon(paropt%anoise)
    endif
    if(trim(paropt%approach)=='SD') then
        if(trim(paropt%sdsaturation)=='unknown') paropt%sdsaturation='no'
    else
        paropt%sdsaturation='yes'
    endif
    if(paropt%nit<0) paropt%nit=1000
    if(paropt%funits<0.d0) paropt%funits=1.d0
    if(paropt%fnrmtolsatur<0.d0) paropt%fnrmtolsatur=4.d-2*paropt%funits !paropt%fmaxtol**0.1d0
    if(paropt%nsatur<0) paropt%nsatur=5
    if(paropt%feedback<0) paropt%feedback=1
    if(paropt%maxforcecall<0) paropt%maxforcecall=100
    if(paropt%alpha<0.d0) paropt%alpha=0.5d0*paropt%alphax
    !if(paropt%alpha<0.d0) paropt%alpha=1.d-1*paropt%alphax
    if(paropt%alphamin<0.d0) paropt%alphamin=1.d-1*paropt%alphax
    if(paropt%alphamax<0.d0) paropt%alphamax=2.0d0*paropt%alphax
    if(paropt%dxmax<0.d0) paropt%dxmax=0.1d0
    paropt%optional_control_on_saturation=.true.
    paropt%sdsaturated=.false.
    paropt%epotitm2=1.d50
    paropt%fnrmitm2=1.d50
    paropt%epotitm1=1.d50
    paropt%fnrmitm1=1.d50
    paropt%care=.true.
    paropt%isatur=0
    work(1:nr)=x(1:nr)
    work(nr+1:2*nr)=0.d0
    !write(*,'(a,i4,2es15.7)') 'nsatur,fnrmtolsatur,alphax ', &
    !    paropt%nsatur,paropt%fnrmtolsatur,paropt%alphax
    if(paropt%lprint .and. .not. paropt%param_reported) call report_param(paropt)
    call yaml_sequence_open('SD optimization iterations')
end subroutine init_sdminimum
!*****************************************************************************************
subroutine what_is_condition_of_feedback(paropt,de1,df1,feedbackcondition)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(in):: paropt
    real(8), intent(in):: de1, df1
    logical, intent(out):: feedbackcondition
    if(paropt%feedback==1) then
        if(de1>paropt%anoise) then
            feedbackcondition=.true.
        else
            feedbackcondition=.false.
        endif
    elseif(paropt%feedback==2) then
        if(df1>paropt%anoise) then
            feedbackcondition=.true.
        else
            feedbackcondition=.false.
        endif
    else
        feedbackcondition=.true.
    endif
end subroutine what_is_condition_of_feedback
!*****************************************************************************************
subroutine test_saturation(paropt,de1,de2,df2,fnrm)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    !logical, intent(in):: care
    real(8), intent(in):: de1, de2, df2, fnrm !, fnrmitm1
    !integer, intent(inout):: isatur
    !local variables
    logical:: c1, c2, c3, c4, c5, c6, c7, c8, c9
    c1=.false. ; c2=.false. ; c3=.false. ; c4=.false.
    c5=.false. ; c6=.false. ; c7=.false. ; c8=.false.
    if(paropt%care) c1=.true.
    if(paropt%itsd>5) c2=.true.
    if(paropt%feedback==1) then
        if(de1>-3.5d-3) c3=.true.
        if(de1<paropt%anoise) c4=.true.
        if(de2>-2.d0*paropt%anoise) c5=.true.
        !write(81,*) de2,paropt%anoise,c5
    else
        c3=.true.
        c4=.true.
        c5=.true.
    endif
    if(fnrm/paropt%fnrmitm1>0.9d0) c6=.true.
    if(fnrm<paropt%fnrmtolsatur) c7=.true.
    if(df2>-2.d0*paropt%anoise) c8=.true.
    c9=paropt%optional_control_on_saturation
    !write(*,'(a,i6,9l3,2f20.10)') 'SAT ',paropt%itsd,c1,c2,c3,c4,c5,c6,c7,c8,c9,fnrm,paropt%fnrmtolsatur
    if(c1 .and. c2 .and. c3 .and. c4 .and. c5 .and. c6 .and. c7 .and. c8 .and. c9) then
        paropt%isatur=paropt%isatur+1
    else
        paropt%isatur=0
    endif
end subroutine test_saturation
!*****************************************************************************************
subroutine final_sdminimum(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    !local variables
    paropt%iflag=0
    paropt%alpha=-1.d0
    paropt%itsd=0
end subroutine final_sdminimum
!*****************************************************************************************
