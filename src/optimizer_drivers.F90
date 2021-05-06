!*****************************************************************************************
subroutine minimize(parini,iproc,atoms,paropt)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info, update_ratp, update_rat
    use mod_opt, only: typ_paropt
    use mod_bin, only: write_bin_conf
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt
    !local variables
    type(typ_file_info):: file_info
    integer:: nr, istat, m, nwork, iat, info, i
    integer, allocatable:: ipiv(:)
    real(8), allocatable:: diag(:), work(:), xr(:), fr(:)
    real(8):: drift(3)
    real(8):: count_sqnm
    logical:: fail
    character(10):: fnmd
    character(50):: comment
    paropt%converged=.false.
    nr=atoms%ndof
    !if(paropt%lprint) write(*,'(a,a,1x,i3)') 'begin of minimization using ',trim(paropt%approach),iproc
    if(paropt%approach=='unknown') then
        if(iproc==0) write(*,*) 'The minimize routine returns becuase method is not specified.'
        return
    endif
    allocate(xr(nr),fr(nr))
    file_info%filename_positions=trim(paropt%filename)
    file_info%file_position='new'
    file_info%print_force=paropt%print_force
    !if(paropt%lprint) call report_param(paropt)
    paropt%iflag=0
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='SD') then
        nwork=2*nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating diag.'
        do
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            call sdminimum(parini,iproc,nr,xr,fr,atoms%epot,paropt,nwork,work)
    !if(paropt%iflag<=0) then
    !write(73,'(i5)') nr/3
    !do i=1,nr,3
    !    write(73,'(3es23.13)') xr(i),xr(i+1),xr(i+2)
    !enddo
    !    write(*,'(3es23.13)') xr(1),xr(2),xr(3)
    !    stop
    !endif
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0) exit
        enddo
        call yaml_sequence_close()
    !write(72,'(i5)') atoms%nat
    !do iat=1,atoms%nat
    !    write(72,'(3es23.13)') atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
    !    write(*,'(3es23.13)') atoms%rat(1,1),atoms%rat(2,1),atoms%rat(3,1)
    !stop
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='SDCG') then
        nwork=2*nr+nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        if(paropt%nfail<0) paropt%nfail=5
        paropt%ifail=0
        do
        do
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            call sdminimum(parini,iproc,nr,xr,fr,atoms%epot,paropt,nwork,work)
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0 .or. paropt%sdsaturated) exit
        enddo
        do 
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            !write(24,*) paropt%fmax
            call cgminimum(iproc,3*atoms%nat,nr,xr,fr,atoms%epot,paropt,nwork,work)
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0) exit
        enddo
        if(paropt%converged) exit
        paropt%ifail=paropt%ifail+1
        if(paropt%ifail>=paropt%nfail) then
            write(*,'(a)') 'WARNING: too many failures of SDCG'
            exit
        endif
        enddo !end of loop over ifail
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='SDDIIS') then
        nwork=(3*paropt%idsx+3)*nr !2*n+nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        !paropt%sdsaturated=.false.
        !do
        !    call cal_potential_forces(parini,iproc,3*atoms%nat,atoms%rat,atoms%fat,atoms%epot)
        !    if(.not. paropt%sdsaturated) then
        !        call sdminimum(parini,iproc,nr,atoms%rat,atoms%fat,atoms%epot,paropt,nwork,work)
        !    endif
        !    if(paropt%converged) write(*,*) 'converged before starting DIIS'
        !    if(paropt%itsd>paropt%nit .and. .not. paropt%sdsaturated) then
        !        write(*,'(a)') 'SD did not saturate, so cgminimum can not continue.'
        !    endif 
        !    if(paropt%sdsaturated .and. .not. paropt%converged) then
        !        call diisminimum(3*atoms%nat,nr,atoms%rat,atoms%epot,atoms%fat,paropt,nwork,work)
        !    endif
        !    if(paropt%iflag<0 .or. paropt%converged) exit
        !enddo
        if(paropt%nfail<0) paropt%nfail=5
        paropt%ifail=0
        do
        do
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            call sdminimum(parini,iproc,nr,xr,fr,atoms%epot,paropt,nwork,work)
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0 .or. paropt%sdsaturated) exit
        enddo
        do 
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            !write(24,*) paropt%fmax
            !call cgminimum(iproc,3*atoms%nat,nr,xr,fr,atoms%epot,paropt,nwork,work)
            call diisminimum(3*atoms%nat,nr,xr,atoms%epot,fr,paropt,nwork,work)
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0) exit
        enddo
        if(paropt%converged) exit
        paropt%ifail=paropt%ifail+1
        if(paropt%ifail>=paropt%nfail) then
            write(*,'(a)') 'WARNING: too many failures of SDCG'
            exit
        endif
        enddo !end of loop over ifail
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='NLBFGS') then
        m=9
        paropt%DIAGCO=.false.
        allocate(diag(nr),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating diag.'
        allocate(work(nr*(2*m+1)+2*m),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        !do
        !    call cal_potential_forces(parini,iproc,3*atoms%nat,atoms%rat,atoms%fat,atoms%epot)
        !    call nlbfgs(nr,m,atoms%rat,atoms%epot,atoms%fat,diag,work,paropt)
        !    !SUBROUTINE NLBFGS(N,M,X,F,G,DIAG,W,paropt)
        !    !SUBROUTINE NLBFGS(N,M,X,F,G,DIAG,IPRINT,EPS,XTOL,W,IFLAG,paropt)
        !    if(paropt%iflag<=0) exit
        !    !if(paropt%iflag<0 .or. paropt%converged) exit
        !enddo
        do
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            call nlbfgs(nr,m,xr,atoms%epot,fr,diag,work,paropt)
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0) exit
            !if(paropt%iflag<0 .or. paropt%converged) exit
        enddo
        deallocate(diag,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating diag.'
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='BFGS') then
        nwork=nr*nr+3*nr+3*nr*nr+3*nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        !open(unit=1,file='t_posinp.xyz',status='old')
        !do i=1,nr/3
        !    read(1,*) x(i*3-2),x(i*3-1),x(i*3-0)
        !enddo
        !close(1)
        !    comment='produced by nab_pdbwrite'
        !    fnmd='posout.pdb'
        !    call nab_pdbwrite(x,trim(fnmd),len(trim(fnmd)),comment)
        !    !call nab_pdbwrite(x,'posout.xyz',10,'structure_with_wrong_force')
        !    call cal_potential_forces(parini,iproc,3*atoms%nat,x,f,epot)
        !    write(*,*) 'DEB ',epot
        !    x(1:nr)=x(1:nr)+1.d-5*f(1:nr)
        !    call cal_potential_forces(parini,iproc,3*atoms%nat,x,f,epot)
        !    write(*,*) 'DEB ',epot
        !    stop
        do
    !write(72,'(i5)') atoms%nat
    !do iat=1,atoms%nat
    !    write(72,'(3es23.13)') atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
    !write(73,'(i5)') atoms%nat
    !do iat=1,atoms%nat
    !    write(73,'(3es23.13)') atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)
    !enddo

            !if(iproc==0) then
            !drift(1:3)=0.d0
            !do iat=1,atoms%nat
            !    drift(1)=drift(1)+atoms%fat(1,iat)
            !    drift(2)=drift(2)+atoms%fat(2,iat)
            !    drift(3)=drift(3)+atoms%fat(3,iat)
            !enddo
            !write(*,'(a,3es20.10)') 'drift ',drift(1),drift(2),drift(3)
            !endif
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            call mybfgs(iproc,nr,xr,atoms%epot,fr,nwork,work,paropt)
            !if(iproc==0) write(*,*) 'TESET ',paropt%iflag,paropt%converged
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0) exit
            !if(paropt%iflag<0 .or. paropt%converged) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
        call yaml_sequence_close()
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='DFP') then
        stop 'ERROR: DFP not ready to use yet.'
        !nwork=nr*nr+3*nr+2*nr*nr+2*nr
        !allocate(work(nwork),stat=istat)
        !if(istat/=0) stop 'ERROR: failure allocating work.'
        !do
        !    call cal_potential_forces(parini,iproc,3*atoms%nat,x,f,epot)
        !    call mydfp(nr,x,epot,f,nwork,work,paropt)
        !    if(paropt%iflag<=0) exit
        !    !if(paropt%iflag<0 .or. paropt%converged) exit
        !enddo
        !deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='FIRE') then
        allocate(work(3*nr),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        do
            if(paropt%trajectory) then
                !call acf_write(file_info,atoms=atoms,strkey='geopt')
                call write_bin_conf(file_info,atoms,strkey='geopt')
            endif
            file_info%file_position='append'
            call cal_potential_forces(parini,atoms)
            call update_ratp(atoms)
            call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
            call test_convergence(nr,fr,paropt)
            call fire(parini,iproc,nr,xr,atoms%epot,fr,work,paropt)
            call xr_to_x(nr,xr,3*atoms%nat,atoms%bemoved,atoms%ratp)
            call update_rat(atoms)
            if(paropt%iflag<=0) exit
            !if(paropt%iflag<0 .or. paropt%converged) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
        call yaml_sequence_close()
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='SQNM') then
        count_sqnm=0.d0
        fail=.true.
        call sqnm(parini,atoms,paropt,count_sqnm,fail)
        call update_ratp(atoms)
        call x_to_xr(3*atoms%nat,atoms%ratp,atoms%fat,atoms%bemoved,nr,xr,fr)
        call test_convergence(nr,fr,paropt)
        call yaml_sequence_close()
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='GMDFIRE') then
        !allocate(work(5*nr+2*3*atoms%nat*3*atoms%nat),stat=istat)
        !if(istat/=0) stop 'ERROR: failure allocating work.'
        !do
        !    call cal_potential_forces(parini,iproc,3*atoms%nat,atoms%rat,atoms%fat,atoms%epot)
        !    stop 'ERROR: the following line is comment so GMDFIRE cannot be used'
        !    !call calpseudohess(3*atoms%nat,atoms%rat,work(5*nr+1))
        !    do i=1,3*atoms%nat
        !        work(5*nr+(i-1)*3*atoms%nat+i)=work(5*nr+(i-1)*3*atoms%nat+i)+1.d2
        !    enddo
        !    allocate(ipiv(nr))
        !    call DSYTRF('L',3*atoms%nat,work(5*nr+1),3*atoms%nat,ipiv,work(5*nr+3*atoms%nat*3*atoms%nat+1),3*atoms%nat*3*atoms%nat,info)
        !    if(info/=0) then;write(*,*) 'ERROR: DSYTRF failed: info',info;stop;endif
        !    work(3*nr+1:4*nr)=atoms%fat(1:3*atoms%nat)
        !    call DSYTRS('L',3*atoms%nat,1,work(5*nr+1),3*atoms%nat,ipiv,work(3*nr+1),3*atoms%nat,info)
        !    if(info/=0) then;write(*,*) 'ERROR: DSYTRS failed: info',info;stop;endif
        !    deallocate(ipiv)
        !    call gmdfire(nr,atoms%rat,atoms%epot,atoms%fat,work,paropt)
        !    if(paropt%iflag<=0) exit
        !    !if(paropt%iflag<0 .or. paropt%converged) exit
        !enddo
        !deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    deallocate(xr)
    !if(paropt%lprint) write(*,'(a,a,1x,i3)') 'end of minimization using ',trim(paropt%approach),iproc
end subroutine minimize
!*****************************************************************************************
subroutine test_convergence(n,f,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: f(n)
    type(typ_paropt), intent(inout):: paropt
    !local variables
    !real(8):: fmax
    call calmaxforcecomponent(n,f,paropt%fmax)
    if(paropt%fmax<paropt%fmaxtol) then
        paropt%converged=.true.
    else
        paropt%converged=.false.
    endif
end subroutine test_convergence
!*****************************************************************************************
subroutine x_to_xr(n,x,f,bemoved,nr,xr,fr)
    !use mod_atoms, only: 
    implicit none
    integer, intent(in):: n, nr
    real(8), intent(in):: x(n), f(n)
    logical, intent(in):: bemoved(3,n/3)
    real(8), intent(inout):: xr(nr), fr(nr)
    !local variables
    integer:: i, j, ixyz, iat
    j=0
    do i=1,n
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(bemoved(ixyz,iat)) then
            j=j+1
            if(j>nr) stop 'ERROR: j>nr in subroutine x_to_xr'
            xr(j)=x(i)
            fr(j)=f(i)
        endif
    enddo
end subroutine x_to_xr
!*****************************************************************************************
subroutine xr_to_x(nr,xr,n,bemoved,x)
    !use mod_atoms, only: 
    implicit none
    integer, intent(in):: n, nr
    real(8), intent(in):: xr(nr)
    logical, intent(in):: bemoved(3,n/3)
    real(8), intent(inout):: x(n)
    !local variables
    integer:: i, j, ixyz, iat
    j=0
    do i=1,n
        iat=(i-1)/3+1
        ixyz=mod(i-1,3)+1
        if(bemoved(ixyz,iat)) then
            j=j+1
            if(j>nr) stop 'ERROR: j>nr in subroutine xr_to_x'
            x(i)=xr(j)
        endif
    enddo
end subroutine xr_to_x
!*****************************************************************************************
subroutine report_param(paropt)
    use mod_opt, only: typ_paropt
    use yaml_output
    implicit none
    type(typ_paropt), intent(inout):: paropt
    !local variables
    call yaml_mapping_open('optimization parameters') !,label='id001')
    call yaml_map('method',trim(paropt%approach))
    call yaml_map('precaution',trim(paropt%precaution))
    call yaml_map('alphax',paropt%alphax,fmt='(f10.6)')
    call yaml_map('fmaxtol',paropt%fmaxtol,fmt='(f10.6)')
    call yaml_map('dxmax',paropt%dxmax,fmt='(f10.6)')
    call yaml_map('condnum',paropt%condnum,fmt='(f10.6)')
    call yaml_map('fnrmtolsatur',paropt%fnrmtolsatur,fmt='(f10.6)')
    call yaml_map('dt_start',paropt%dt_start,fmt='(f10.6)')
    call yaml_map('nsatur',paropt%nsatur,fmt='(i3)')
    call yaml_map('nit',paropt%nit,fmt='(i5)')
    call yaml_map('lprint',paropt%lprint,fmt='(l)')
    !write(*,'(a,2x,a)'  ) 'optimization parameters:       method ',trim(paropt%approach)
    !write(*,'(a,2x,a)'  ) 'optimization parameters:   precaution ',trim(paropt%precaution)
    !write(*,'(a,es10.2)') 'optimization parameters:       alphax ',paropt%alphax
    !write(*,'(a,es10.2)') 'optimization parameters:      fmaxtol ',paropt%fmaxtol
    !write(*,'(a,es10.2)') 'optimization parameters:        dxmax ',paropt%dxmax
    !write(*,'(a,es10.2)') 'optimization parameters:      condnum ',paropt%condnum
    !write(*,'(a,es10.2)') 'optimization parameters: fnrmtolsatur ',paropt%fnrmtolsatur
    !write(*,'(a,es10.2)') 'optimization parameters:     dt_start ',paropt%dt_start
    !write(*,'(a,1x,i3)' ) 'optimization parameters:       nsatur ',paropt%nsatur
    !write(*,'(a,1x,i5)' ) 'optimization parameters:          nit ',paropt%nit
    !write(*,'(a,1x,l)'  ) 'optimization parameters:       lprint ',paropt%lprint
    paropt%param_reported=.true.
    call yaml_mapping_close()
end subroutine report_param
!*****************************************************************************************
subroutine initminimize(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    !local variables
    integer::istat
    character(2)::tapp1,tapp2
    character(4)::tapp3
    tapp1(1:2)=paropt%approach(1:2)
    if(len(trim(paropt%approach))==4) tapp2(1:2)=paropt%approach(3:4)
    if(len(trim(paropt%approach))==6) tapp3(1:4)=paropt%approach(3:6)
    !write(*,*) tapp1
    !write(*,*) tapp2
    !write(*,*) tapp3
    !write(*,*) len(trim(paropt%approach))
    !paropt%maxforcecall=10000
    !if(paropt%fmaxtol<0.d0) stop 'ERROR: fmaxtol<0, maybe it is not set.'
    !if(tapp1=='SD') then
    !    if(paropt%alphax<0.d0) stop 'ERROR: alphax<0, maybe it is not set.'
    !    paropt%alphamin=1.d-1*paropt%alphax
    !    paropt%alphamax=2.0d0*paropt%alphax
    !    paropt%fnrmtolsatur=0.5d0 !paropt%fmaxtol**0.1d0
    !    paropt%anoise=epsilon(paropt%anoise)
    !    paropt%nitsd=10000
    !    paropt%nsatur=5
    !    if(paropt%feedback<0) paropt%feedback=1
    !endif
    !if(tapp1=='CG' .or. tapp2=='CG') then
    !    paropt%nitcg=500
    !    paropt%nfail=10
    !endif
    if(tapp3=='DIIS') then
        paropt%idsx=3
        allocate(paropt%a(paropt%idsx+1,paropt%idsx+1,3),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating paropt%a.'
        allocate(paropt%b(paropt%idsx+1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating paropt%b.'
        allocate(paropt%ipiv(paropt%idsx+1),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating paropt%ipiv.'
    endif
end subroutine initminimize
!*****************************************************************************************
subroutine finalminimize(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    !local variables
    integer::istat
    character(2)::tapp1,tapp2
    character(4)::tapp3
    tapp1(1:2)=paropt%approach(1:2)
    if(len(trim(paropt%approach))==4) tapp2(1:2)=paropt%approach(3:4)
    if(len(trim(paropt%approach))==6) tapp3(1:4)=paropt%approach(3:6)
    if(tapp3=='DIIS') then
        paropt%idsx=10
        deallocate(paropt%a,stat=istat)
        if(istat/=0) stop 'ERROR: failure deallocating paropt%a.'
        deallocate(paropt%b,stat=istat)
        if(istat/=0) stop 'ERROR: failure deallocating paropt%b.'
        deallocate(paropt%ipiv,stat=istat)
        if(istat/=0) stop 'ERROR: failure deallocating paropt%ipiv.'
    endif
end subroutine finalminimize
!*****************************************************************************************
