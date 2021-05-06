!*****************************************************************************************
subroutine optimizer_saddle(parini,iproc,atoms_s,n,nr,x,f,epot,paropt,uvn)
    use mod_parini, only: typ_parini
    use mod_opt, only: typ_paropt
    use mod_potential, only: fcalls
    use mod_atoms, only: typ_atoms
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, n, nr
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: x(n), f(n), uvn(n), epot
    type(typ_paropt), intent(inout):: paropt
    !local variables
    integer:: icall, istat, m, it, nwork, i, info
    real(8):: t1, curv, curv0, drift(3)
    character(10):: fnmd
    character(50):: comment
    real(8), allocatable:: work(:), xr(:), fr(:), fold(:), feff(:)
    integer, allocatable:: ipiv(:)
    paropt%converged=.false.
    !if(paropt%lprint) write(*,'(a,a,1x,i3)') 'begin of minimization using ',trim(paropt%approach),iproc
    if(paropt%approach=='unknown') then
        if(iproc==0) write(*,*) 'The optimizer_saddle routine returns becuase method is not specified.'
        return
    endif
    allocate(xr(nr),fr(nr),fold(n),feff(n))
    !if(nr/=n) stop 'ERROR: nr/=n'
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='SD') then
        nwork=2*nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        do
            call cal_potential_forces_modified(parini,iproc,atoms_s,n,x,f,epot,nr,uvn,feff,curv0,curv,fold,paropt)
            call x_to_xr(n,x,feff,atoms_s%bemoved,nr,xr,fr)
            call sdminimum(parini,iproc,nr,xr,fr,epot,paropt,nwork,work)
            call xr_to_x(nr,xr,n,atoms_s%bemoved,x)
            if(paropt%iflag<=0 .or. paropt%converged) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='SDCG') then
        nwork=2*nr+nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        if(paropt%nfail<0) paropt%nfail=1
        paropt%ifail=0
        do
        do
            paropt%approach_current='SD'
            !write(*,*) 'REZA-1 ',n
            call cal_potential_forces_modified(parini,iproc,atoms_s,n,x,f,epot,nr,uvn,feff,curv0,curv,fold,paropt)
            !write(*,'(a,i2,i6,i7)') 'TEST-SD ',paropt%ifail,it,int(fcalls)
            call x_to_xr(n,x,feff,atoms_s%bemoved,nr,xr,fr)
            call test_convergence_saddle(nr,fr,curv,paropt)
            call sdminimum(parini,iproc,nr,xr,fr,epot,paropt,nwork,work)
            call xr_to_x(nr,xr,n,atoms_s%bemoved,x)
            if(paropt%iflag<=0 .or. paropt%sdsaturated) exit
        enddo
        do it=1,10000000
            if(it==1) then
                work(1:nwork)=0.d0
            else
            paropt%approach_current='CG'
            call cal_potential_forces_modified(parini,iproc,atoms_s,n,x,f,epot,nr,uvn,feff,curv0,curv,fold,paropt)
            !write(*,'(a,i2,i6,i7)') 'TEST-CG ',paropt%ifail,it,int(fcalls)
            endif
            call x_to_xr(n,x,feff,atoms_s%bemoved,nr,xr,fr)
            call test_convergence_saddle(nr,fr,curv,paropt)
            write(24,*) paropt%fmax
            call cgminimum(iproc,n,nr,xr,fr,epot,paropt,nwork,work)
            call xr_to_x(nr,xr,n,atoms_s%bemoved,x)
            if(paropt%iflag<=0) exit
        enddo
        call yaml_sequence_close()
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
    deallocate(xr,fr,fold,feff)
    !if(paropt%lprint) write(*,'(a,a,1x,i3)') 'end of minimization using ',trim(paropt%approach),iproc
end subroutine optimizer_saddle
!*****************************************************************************************
subroutine cal_potential_forces_modified(parini,iproc,atoms_s,n,x,f,epot,nr,uvn,feff,curv0,curv,fold,paropt)
    use mod_parini, only: typ_parini
    use mod_saddle, only: beta
    use mod_opt, only: typ_paropt
    use mod_atoms, only: typ_atoms, atom_ddot, atom_copy_old, atom_calnorm, atom_deallocate_old
    use mod_atoms, only: set_rat
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc, n, nr
    type(typ_atoms), intent(inout):: atoms_s
    real(8), intent(inout):: x(3,atoms_s%nat)
    real(8), intent(inout):: f(3,atoms_s%nat), epot, uvn(3,atoms_s%nat), feff(3,atoms_s%nat), curv0, curv, fold(n)
    type(typ_paropt), intent(inout):: paropt
    !local variables
    type(typ_atoms):: atoms
    real(8):: t1, DDOT, fnrm, fnrm_small
    call atom_copy_old(atoms_s,atoms,'atoms_s->atoms')
    !atoms%rat(1:3,1:atoms%nat)=x(1:3,1:atoms%nat)
    call set_rat(atoms,x,setall=.true.)
    call cal_potential_forces(parini,atoms)
    f(1:3,1:atoms%nat)=atoms%fat(1:3,1:atoms%nat)
    epot=atoms%epot
    !call calnorm(n,f,fnrm)
    call atom_calnorm(n/3,atoms_s%bemoved,f,fnrm)
    fnrm_small=6.d-3*sqrt(real(nr,8))
    if(.not. (fnrm<fnrm_small .and. paropt%fnrmitm1<fnrm_small .and. curv<0.d0)) then
        call lowestcurvature(parini,iproc,atoms_s,n/3,nr,x,uvn,f,2.d0,100,curv0,curv,1)
    endif

    t1=atom_ddot(n/3,f,uvn,atoms_s%bemoved)
    if(trim(paropt%approach_current)=='SD') then
        beta=16
    elseif(trim(paropt%approach_current)=='CG') then
        beta=4
    endif
    !write(*,*) 'REZA-3 ',n
    feff(1:3,1:atoms%nat)=f(1:3,1:atoms%nat)-beta*t1*uvn(1:3,1:atoms%nat)
    if(curv<0.d0) then
        paropt%optional_control_on_saturation=.true.
    else
        paropt%optional_control_on_saturation=.false.
    endif
    call atom_deallocate_old(atoms)
end subroutine cal_potential_forces_modified
!*****************************************************************************************
subroutine test_convergence_saddle(n,f,curv,paropt)
    use mod_opt, only: typ_paropt
    implicit none
    integer, intent(in):: n
    real(8), intent(in):: f(n), curv
    type(typ_paropt), intent(inout):: paropt
    !local variables
    !real(8):: fmax
    call calmaxforcecomponent(n,f,paropt%fmax)
    if(paropt%fmax<paropt%fmaxtol .and. curv<0.d0) then
        paropt%converged=.true.
    else
        paropt%converged=.false.
    endif
end subroutine test_convergence_saddle
!*****************************************************************************************
