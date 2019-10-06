!*****************************************************************************************
subroutine mybfgs(iproc,nr,x,epot,f,nwork,work,paropt)
    use mod_opt, only: typ_paropt, frmt_base
    use yaml_output
    implicit none
    integer, intent(in):: iproc, nr, nwork
    real(8):: x(nr), f(nr), epot, work(nwork)
    integer, allocatable:: ipiv(:)
    type(typ_paropt), intent(inout):: paropt
    !local variables
    real(8):: DDOT, de, fnrm, fmax, beta, fnrmsatur
    real(8):: y_B_y !y_k^T * B_k^-1 * y_k in wiki notation
    real(8):: s_dot_y !dot product of s_k^T and y_k in wiki notation
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    character(14):: filename
    character(56):: comment
    character(51), parameter:: frmt='('//frmt_base//',2es12.4,i3,a)'
    integer:: nrsqtwo, iw1, iw2, iw3, iw4, iw5, iw6, iw7, iw8, iw9, info, i, j, l
    if(nwork/=nr*nr+3*nr+3*nr*nr+3*nr) then
        stop 'ERROR: size of work array is insufficient.'
    endif
    nrsqtwo=nr*nr*2
    iw1=1            !inverse of the quasi-Hessian matrix
    iw2=iw1+nr*nr    !force of previous iteration in wiki notation
    iw3=iw2+nr       !y_k in wiki notation
    iw4=iw3+nr       !s_k in wiki notation
    iw5=iw4+nr       !work array to keep the hessian untouched
    iw6=iw5+nr*nr    !work array DSYEV required by DSYEV
    iw7=iw6+nrsqtwo  !p_k in wiki notation
    iw8=iw7+nr       !position of previous iteration
    iw9=iw8+nr       !eigenvalues of inverse of hessian
    call calnorm(nr,f,fnrm)
    !call calmaxforcecomponent(nr,f,fmax)
    fmax=paropt%fmax
    if(paropt%iflag==0) then
        call init_mybfgs(paropt,epot,fmax)
    else
        paropt%iter=paropt%iter+1
    endif
    !df1=fnrm-paropt%fnrmitm1
    !if(df1<0.d0) then
    if(fnrm<1.05d0*paropt%fnrmitm1) then
        paropt%ifnrminc=0
    else
        paropt%ifnrminc=paropt%ifnrminc+1
    endif
    fnrmsatur=paropt%funits*paropt%prefactor*min(50.d0,max(2.5d0,sqrt(real(nr,8))))
    !fnrmsatur=min(4.d-2,max(2.d-3,8.d-4*sqrt(real(nr,8)))) !I am not sure
    !fnrmsatur=min(30.d0,max(1.5d0,6.d-1*sqrt(real(nr,8)))) !LJ38
    if(paropt%iter>5 .and. fnrm<fnrmsatur) then
        if(paropt%ifnrminc==0) then
            paropt%isatur=min(paropt%isatur+paropt%increment,99)
        endif
    else
        paropt%isatur=0
    endif
    de=epot-paropt%epotold
    !if(iproc==0) write(*,'(a10,i4,es23.15,es11.3,2es12.5,2es12.4,i3)') &
    !    'BFGSMIN   ',paropt%iter,epot,de,fnrm,fmax,paropt%zeta,paropt%alpha,paropt%isatur
    !if(paropt%lprint) write(*,'(a4,i3.3,1x,i5,es23.12,es11.3,2es12.5,2es12.4,i3,a)') &
    !if(paropt%lprint) write(*,'(a4,i3.3,1x,i5,es23.15,es11.3,2es12.5,2es12.4,i3,a)') &
    if(paropt%lprint) call yaml_sequence(advance='no')
    if(paropt%lprint) call yaml_mapping_open('BFGS',flow=.true.)
    !if(paropt%lprint) call yaml_map('iproc',iproc,fmt='(i3.3)')
    if(paropt%lprint) call yaml_map('iter',paropt%iter,fmt='(i5)')
    if(paropt%lprint) call yaml_map('epot',epot,fmt='(es20.12)')
    if(paropt%lprint) call yaml_map('de',de,fmt='(es9.1)')
    if(paropt%lprint) call yaml_map('fmax',fmax,fmt='(es10.3)')
    if(paropt%lprint) call yaml_map('fnrm',fnrm,fmt='(es10.3)')
    if(paropt%lprint) call yaml_map('alpha',paropt%alpha,fmt='(es12.4)')
    if(paropt%lprint) call yaml_map('zeta',paropt%zeta,fmt='(es12.4)')
    if(paropt%lprint) call yaml_map('isatur',paropt%isatur,fmt='(i5)')
    if(paropt%lprint) call yaml_mapping_close()
    !if(paropt%lprint) write(*,frmt) 'MIN:',iproc,paropt%iter,epot,de,fnrm,fmax,paropt%zeta,paropt%alpha,paropt%isatur,' BFGS'
    !if(fmax<paropt%fmaxtol) then
    if(paropt%converged) then
        paropt%iflag=0
        !if(paropt%lprint) write(*,'(a,i4,es23.12,2es12.5)') &
        !if(paropt%lprint) write(*,'(a,i4,es23.15,2es12.5)') &
        !    'BFGS FINISHED: itbfgs,epot,fnrm,fmax ',paropt%iter,epot,fnrm,fmax
        if(paropt%lprint) then
            call yaml_sequence(advance='no')
            call yaml_mapping_open('BFGS FINISHED') !,label='id001')
            call yaml_map('success',.true.)
            call yaml_map('iter',paropt%iter,fmt='(i5)')
            call yaml_map('epot',epot,fmt='(es20.12)')
            call yaml_map('fnrm',fnrm,fmt='(es12.5)')
            call yaml_map('fmax',fmax,fmt='(es12.5)')
            call yaml_mapping_close()
        endif
        return
    endif
    if(paropt%iter==paropt%nit) then
        paropt%iflag=0
        write(*,'(a)') 'BFGS: NO CONVERGENCE '
        return
    endif

    !if(de>0.d0 .and. paropt%zeta>1.d-1) then
    if(de>1.d-4) then
    !if(de>4.d-4) then
        epot=paropt%epotold
        x(1:nr)=work(iw8:iw8-1+nr)
        f(1:nr)=work(iw2:iw2-1+nr)
        paropt%reset_hess=.true.
        !paropt%alpha=max(paropt%alpha*0.5d0/1.1d0,1.d-2)
        paropt%zeta=max(paropt%zeta*2.d-1,1.d-3)
        paropt%isatur=0
    else
        !paropt%zeta=1.d0
        !if(paropt%zeta>1.d-1) paropt%zeta=min(paropt%zeta*1.1d0,1.d0)
        paropt%zeta=min(paropt%zeta*1.1d0,1.d0)
        !paropt%isatur=paropt%isatur+1
        paropt%fnrmitm1=fnrm
    endif
    if(paropt%iter>0) then
        !next line is to calculate s_k in wiki notation
        work(iw4:iw4-1+nr)=x(1:nr)-work(iw8:iw8-1+nr)
        !next line is to calculate y_k in wiki notation
        work(iw3:iw3-1+nr)=work(iw2:iw2-1+nr)-f(1:nr)
        !s_dot_y=DDOT(nr,work(iw3),1,work(iw4),1)
        !previous line is the one from wiki and next line may be a bad choice.
        !the point is that next line is a good for testing whether the last
        !move was along the steepest descent direction
        !this must be changed in future to both take advantage of this point and
        !be consistent to wiki (true) BFGS.
        s_dot_y=DDOT(nr,work(iw2),1,work(iw4),1)
    else
        s_dot_y=1.d0 !just some positive value
    endif
    if(paropt%iter==0 .or. paropt%reset_hess .or. s_dot_y<0.d0 .or. paropt%ifnrminc>3) then
        paropt%reset_hess=.false.
        paropt%ifnrminc=0
        paropt%isatur=0
        !if(paropt%isatur>=10) then
        !    paropt%reset_hess=.false.
        !    !paropt%alpha=5.d-1
        !endif
        work(1:nr*nr)=0.d0
        do i=1,nr
            work(i+(i-1)*nr)=paropt%zeta*paropt%alphax
        enddo
        work(iw7:iw7-1+nr)=paropt%zeta*paropt%alphax*f(1:nr)
    else
        !the following is to calculate B_k^-1 * y_k in wiki notation
        do i=1,nr
            tt2=0.d0
            do j=1,nr
                tt2=tt2+work(i+(j-1)*nr)*work(iw3-1+j)
            enddo
            work(iw6-1+i)=tt2
        enddo
        y_B_y=DDOT(nr,work(iw3),1,work(iw6),1)
        if(paropt%lprint) write(21,*) paropt%iter,s_dot_y,y_B_y
        !s_dot_y=max(s_dot_y,1.d-2)
        !The following loop updates the inverse quasi-Hessian matrix by equation
        !derived by using Sherman-Morrison formula.
        do i=1,nr
            tt1=work(iw4-1+i)
            do j=i,nr
                l=i+(j-1)*nr
                tt2=(s_dot_y+y_B_y)*tt1*work(iw4-1+j)/s_dot_y**2
                tt3=(work(iw6-1+i)*work(iw4-1+j)+work(iw6-1+j)*tt1)/s_dot_y
                work(l)=work(l)+tt2-tt3
                !work(l)=work(l)+(s_dot_y+y_B_y)*tt1*work(iw4-1+j)/s_dot_y**2- &
                !    (work(iw6-1+i)*work(iw4-1+j)+work(iw6-1+j)*tt1)/s_dot_y
                work(j+(i-1)*nr)=work(l)
            enddo
        enddo
        !do i=1,nr
        !    tt2=0.d0
        !    do j=1,nr
        !        tt2=tt2+work(j+(i-1)*nr)*f(j)
        !    enddo
        !    work(iw7-1+i)=tt2
        !enddo
        !write(31,*) paropt%zeta
        work(iw5:iw5-1+nr*nr)=work(1:nr*nr)
        call DSYEV('V','L',nr,work(iw5),nr,work(iw9),work(iw6),nrsqtwo,info)
        if(info/=0) stop 'DSYEV'
        tt1=work(iw9+0)    ; tt2=work(iw9+1)    ; tt3=work(iw9+2)
        tt4=work(iw9+nr-3) ; tt5=work(iw9+nr-2) ; tt6=work(iw9+nr-1)
        if(paropt%lprint) write(41,'(i5,6es15.5)') paropt%iter,tt1,tt2,tt3,tt4,tt5,tt6
        paropt%evalmin=tt1 ; paropt%evalmax=tt6
        work(iw7:iw7-1+nr)=0.d0
        beta=1.d0+paropt%condnum*min((paropt%isatur/6)/15.d0,1.d0)
        beta=1.d0/(beta*paropt%alphax)
        !do j=1,nr
        !    if(work(iw9-1+j)>0.d0) then
        !        tt3=work(iw9-1+j)
        !        exit
        !    enddo
        !enddo
        tt3=paropt%alphax*0.5d0
        do j=1,nr
            tt1=DDOT(nr,work(iw5+nr*(j-1)),1,f,1)
            if(work(iw9-1+j)<tt3) then
                tt4=tt3
            else
                tt4=work(iw9-1+j)
            endif
            tt2=1.d0/sqrt(1.d0/tt4**2+beta**2)
            do i=1,nr
                work(iw7-1+i)=work(iw7-1+i)+tt1*work(iw5-1+i+nr*(j-1))*tt2
            enddo
        enddo
    endif !provided paropt%iter==0 .or. paropt%reset_hess .or. tt1<0.d0
    if(paropt%evalmin<0.d0 .and. abs(paropt%evalmin)>paropt%alphax*1.d-1) then
        paropt%reset_hess=.true.
        epot=paropt%epotold
        x(1:nr)=work(iw8:iw8-1+nr)
        f(1:nr)=work(iw2:iw2-1+nr)
        work(iw7:iw7-1+nr)=0.d0
        paropt%isatur=0
        paropt%evalmin=paropt%alphax ; paropt%evalmax=1.d1*paropt%alphax !just to avoid problem
    endif
    paropt%epotold=epot
    work(iw2:iw2-1+nr)=f(1:nr)
    work(iw8:iw8-1+nr)=x(1:nr)
    !if(fnrm<1.d0) then
    !    paropt%isatur=paropt%isatur+1
    !else
    !    paropt%isatur=0
    !endif
    paropt%alpha=min(paropt%alphamax,paropt%alpha*1.1d0)
    x(1:nr)=x(1:nr)+paropt%alpha*work(iw7:iw7-1+nr)
end subroutine mybfgs
!*****************************************************************************************
subroutine init_mybfgs(paropt,epot,fmax)
    use mod_opt, only: typ_paropt
    use yaml_output
    implicit none
    type(typ_paropt), intent(inout):: paropt
    real(8), intent(in):: epot, fmax
    !local variables
    paropt%iflag=1
    paropt%iter=0
    paropt%epotold=epot
    paropt%alpha=8.d-1
    paropt%reset_hess=.false.
    paropt%alphamax=0.95d0
    paropt%zeta=min(1.d-1/(paropt%alphax*fmax),1.d0)
    if(paropt%zeta<1.d-3) stop 'ERROR: input configuration does not seem fine.'
    paropt%isatur=0
    paropt%evalmin=paropt%alphax ; paropt%evalmax=1.d1*paropt%alphax !just to avoid problem
    paropt%fnrmitm1=1.d50
    paropt%ifnrminc=0
    if(paropt%funits<0.d0) paropt%funits=1.d0
    if(paropt%nit<0) paropt%nit=1000
    if(trim(paropt%precaution)=='high') then
        paropt%prefactor=2.d-3
        paropt%increment=1
    elseif(trim(paropt%precaution)=='normal') then
        paropt%prefactor=2.4d-3
        paropt%increment=2
    elseif(trim(paropt%precaution)=='low') then
        paropt%prefactor=2.8d-3
        paropt%increment=3
    else
        stop 'ERROR: precaution from typ_paropt is not set'
    endif
    call yaml_sequence_open('BFGS optimization iterations')
end subroutine init_mybfgs
!*****************************************************************************************
