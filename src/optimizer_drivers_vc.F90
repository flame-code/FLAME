!*****************************************************************************************
subroutine vc_minimize(parini,iproc,atoms,paropt)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_file_info
    use mod_atoms, only: update_ratp, update_rat, set_rat, get_rat
    use mod_acf, only: acf_write
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iproc !, n, nr
    !real(8), intent(inout):: x(n), f(n), epot
    type(typ_atoms), intent(inout):: atoms
    type(typ_paropt), intent(inout):: paropt
    !local variables
    type(typ_file_info):: file_info
    integer:: n, nr, istat, m, nwork, iat, info, i
    integer, allocatable:: ipiv(:)
    real(8), allocatable:: diag(:), work(:), xr(:), fr(:)
    real(8), allocatable:: rat_int(:,:), fat_int(:,:)
    real(8):: drift(3), fnrm
    character(10):: fnmd
    character(50):: comment
    paropt%converged=.false.
    n=3*atoms%nat
    nr=atoms%ndof+9
    if(paropt%lprint) write(*,'(a,a,1x,i3)') 'begin of minimization using ',trim(paropt%approach),iproc
    if(paropt%approach=='unknown') then
        if(iproc==0) write(*,*) 'The minimize routine returns becuase method is not specified.'
        return
    endif
    allocate(xr(nr),fr(nr))
    file_info%filename_positions='traj.acf'
    file_info%file_position='new'
    atoms%epot=0.d0
    call acf_write(file_info,atoms=atoms,strkey='traj')
    allocate(rat_int(3,atoms%nat),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating rat_int.'
    allocate(fat_int(3,atoms%nat),stat=istat)
    if(istat/=0) stop 'ERROR: failure allocating fat_int.'
    !atoms%rat(1:3,1:atoms%nat)=rat_int(1:3,1:atoms%nat)
    !if(paropt%lprint) call report_param(paropt)
    paropt%iflag=0
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='SD') then
        nwork=2*nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating diag.'
        do
            !call cal_potential_forces_vc(iproc,atoms%nat,atoms%rat,atoms%cellvec,atoms%pressure,fat_int,atoms%celldv,atoms%stress,atoms%epot,atoms%enth)
            call cal_potential_forces(parini,atoms)
            !do iat=1,atoms%nat
            !    write(*,'(a,i5,3es14.5)') 'FOR ',paropt%itsd,atoms%fat(1,iat),atoms%fat(2,iat),atoms%fat(3,iat)
            !enddo
            call atom_calnorm(atoms%nat,atoms%bemoved,atoms%fat,fnrm)
            write(33,'(i6,10es14.5)') paropt%itsd,fnrm, &
                atoms%celldv(1,1),atoms%celldv(1,2),atoms%celldv(1,3), &
                atoms%celldv(2,1),atoms%celldv(2,2),atoms%celldv(2,3), &
                atoms%celldv(3,1),atoms%celldv(3,2),atoms%celldv(3,3)
            write(61,'(9es14.5)') &
                atoms%stress(1,1),atoms%stress(1,2),atoms%stress(1,3), &
                atoms%stress(2,1),atoms%stress(2,2),atoms%stress(2,3), &
                atoms%stress(3,1),atoms%stress(3,2),atoms%stress(3,3)
            !call fxyz_red2cart(atoms%nat,fat_int,atoms%cellvec,atoms%fat)
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='traj')
            call fxyz_cart2int_alborz(atoms%nat,atoms%fat,atoms%cellvec,fat_int)
            atoms%fat=fat_int
            call update_ratp(atoms)
            call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,rat_int)
            call set_rat(atoms,rat_int,setall=.true.)
            call vc_x_to_xr(atoms,nr,xr,fr)
            call vc_test_convergence(nr,fr,paropt)
            call sdminimum(parini,iproc,nr,xr,fr,atoms%epot,paropt,nwork,work)
            call vc_xr_to_x(nr,xr,atoms)
            call get_rat(atoms,rat_int)
            call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,rat_int,atoms%ratp)
            call update_rat(atoms,upall=.true.)
            if(paropt%iflag<=0) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    if(trim(paropt%approach)=='BFGS') then
        nwork=nr*nr+3*nr+3*nr*nr+3*nr
        allocate(work(nwork),stat=istat)
        if(istat/=0) stop 'ERROR: failure allocating work.'
        do
            !call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%rat,rat_int)
            !do iat=1,atoms%nat
            !    write(*,'(3es14.5)') atoms%rat(1,iat),atoms%rat(2,iat),atoms%rat(3,iat)*1.1d0
            !enddo
            !stop
            !call cal_potential_forces_vc(iproc,atoms%nat,rat_int,atoms%cellvec,atoms%pressure,atoms%fat,atoms%celldv,atoms%stress,atoms%epot,atoms%enth)
            call cal_potential_forces(parini,atoms)
            call atom_calnorm(atoms%nat,atoms%bemoved,atoms%fat,fnrm)
            write(33,'(i6,10es14.5)') paropt%itsd,fnrm, &
                atoms%celldv(1,1),atoms%celldv(1,2),atoms%celldv(1,3), &
                atoms%celldv(2,1),atoms%celldv(2,2),atoms%celldv(2,3), &
                atoms%celldv(3,1),atoms%celldv(3,2),atoms%celldv(3,3)
            write(61,'(9es14.5)') &
                atoms%stress(1,1),atoms%stress(1,2),atoms%stress(1,3), &
                atoms%stress(2,1),atoms%stress(2,2),atoms%stress(2,3), &
                atoms%stress(3,1),atoms%stress(3,2),atoms%stress(3,3)
            !call fxyz_cart2int_alborz(atoms%nat,atoms%rat,atoms%cellvec,rat_int)
            !call fxyz_red2cart(atoms%nat,fat_int,atoms%cellvec,atoms%fat)
            file_info%file_position='append'
            call acf_write(file_info,atoms=atoms,strkey='traj')
            call fxyz_cart2int_alborz(atoms%nat,atoms%fat,atoms%cellvec,fat_int)
            atoms%fat=fat_int
            call update_ratp(atoms)
            call rxyz_cart2int_alborz(atoms%nat,atoms%cellvec,atoms%ratp,rat_int)
            call set_rat(atoms,rat_int,setall=.true.)
            !call fxyz_cart2int_alborz(atoms%nat,fat_int,atoms%cellvec,atoms%fat)
            !atoms%fat=fat_int
            call vc_x_to_xr(atoms,nr,xr,fr)
            call vc_test_convergence(nr,fr,paropt)
            call mybfgs(iproc,nr,xr,atoms%epot,fr,nwork,work,paropt)
            call vc_xr_to_x(nr,xr,atoms)
            call get_rat(atoms,rat_int)
            call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,rat_int,atoms%ratp)
            call update_rat(atoms,upall=.true.)
            !call xr_to_x(nr,xr,n,atoms%bemoved,atoms%rat)
            if(paropt%iflag<=0) exit
        enddo
        deallocate(work,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating work.'
    endif
    !-------------------------------------------------------------------------------------
    call get_rat(atoms,rat_int)
    call rxyz_int2cart_alborz(atoms%nat,atoms%cellvec,rat_int,atoms%ratp)
    call update_rat(atoms,upall=.true.)
    deallocate(rat_int,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating rat_int.'
    deallocate(fat_int,stat=istat);if(istat/=0) stop 'ERROR: failure deallocating fat_int.'
    deallocate(xr)
    if(paropt%lprint) write(*,'(a,a,1x,i3)') 'end of minimization using ',trim(paropt%approach),iproc
end subroutine vc_minimize
!*****************************************************************************************
subroutine vc_test_convergence(n,f,paropt)
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
end subroutine vc_test_convergence
!*****************************************************************************************
subroutine vc_x_to_xr(atoms,nr,xr,fr)
    use mod_atoms, only: typ_atoms, get_rat
    implicit none
    type(typ_atoms), intent(in):: atoms
    integer, intent(in):: nr
    !real(8), intent(in):: x(n), f(n)
    !logical, intent(in):: bemoved(3,n/3)
    real(8), intent(inout):: xr(nr), fr(nr)
    !local variables
    integer:: i, j, ixyz, iat
    real(8), allocatable:: rat(:,:)
    allocate(rat(3,atoms%nat))
    call get_rat(atoms,rat)
    j=0
    do iat=1,atoms%nat
    do ixyz=1,3
        !iat=(i-1)/3+1
        !ixyz=mod(i-1,3)+1
            !write(*,*) 'AAAAAAAAAA ',allocated(atoms%bemoved)
        if(atoms%bemoved(ixyz,iat)) then
            j=j+1
            if(j>nr) stop 'ERROR: j>nr in subroutine x_to_xr'
            xr(j)=rat(ixyz,iat)
            fr(j)=atoms%fat(ixyz,iat)
        endif
    enddo
    enddo
    j=j+1 ; xr(j)=atoms%cellvec(1,1) ; fr(j)=atoms%celldv(1,1)
    j=j+1 ; xr(j)=atoms%cellvec(2,1) ; fr(j)=atoms%celldv(2,1)
    j=j+1 ; xr(j)=atoms%cellvec(3,1) ; fr(j)=atoms%celldv(3,1)
    j=j+1 ; xr(j)=atoms%cellvec(1,2) ; fr(j)=atoms%celldv(1,2)
    j=j+1 ; xr(j)=atoms%cellvec(2,2) ; fr(j)=atoms%celldv(2,2)
    j=j+1 ; xr(j)=atoms%cellvec(3,2) ; fr(j)=atoms%celldv(3,2)
    j=j+1 ; xr(j)=atoms%cellvec(1,3) ; fr(j)=atoms%celldv(1,3)
    j=j+1 ; xr(j)=atoms%cellvec(2,3) ; fr(j)=atoms%celldv(2,3)
    j=j+1 ; xr(j)=atoms%cellvec(3,3) ; fr(j)=atoms%celldv(3,3)
    deallocate(rat)
end subroutine vc_x_to_xr
!*****************************************************************************************
subroutine vc_xr_to_x(nr,xr,atoms)
    use mod_atoms, only: typ_atoms, update_rat
    implicit none
    integer, intent(in):: nr
    real(8), intent(in):: xr(nr)
    type(typ_atoms), intent(inout):: atoms
    !local variables
    integer:: i, j, ixyz, iat
    j=0
    do iat=1,atoms%nat
    do ixyz=1,3
        !iat=(i-1)/3+1
        !ixyz=mod(i-1,3)+1
        if(atoms%bemoved(ixyz,iat)) then
            j=j+1
            if(j>nr) stop 'ERROR: j>nr in subroutine vc_xr_to_x'
            atoms%ratp(ixyz,iat)=xr(j)
        endif
    enddo
    enddo
    call update_rat(atoms)
    j=j+1 ; atoms%cellvec(1,1)=xr(j)
    j=j+1 ; atoms%cellvec(2,1)=xr(j)
    j=j+1 ; atoms%cellvec(3,1)=xr(j)
    j=j+1 ; atoms%cellvec(1,2)=xr(j)
    j=j+1 ; atoms%cellvec(2,2)=xr(j)
    j=j+1 ; atoms%cellvec(3,2)=xr(j)
    j=j+1 ; atoms%cellvec(1,3)=xr(j)
    j=j+1 ; atoms%cellvec(2,3)=xr(j)
    j=j+1 ; atoms%cellvec(3,3)=xr(j)
end subroutine vc_xr_to_x
!*****************************************************************************************
subroutine vc_report_param(paropt)
    use mod_opt, only: typ_paropt
    implicit none
    type(typ_paropt), intent(inout):: paropt
    !local variables
    write(*,'(a,2x,a)'  ) 'optimization parameters:       method ',trim(paropt%approach)
    write(*,'(a,2x,a)'  ) 'optimization parameters:   precaution ',trim(paropt%precaution)
    write(*,'(a,es10.2)') 'optimization parameters:       alphax ',paropt%alphax
    write(*,'(a,es10.2)') 'optimization parameters:      fmaxtol ',paropt%fmaxtol
    write(*,'(a,es10.2)') 'optimization parameters:        dxmax ',paropt%dxmax
    write(*,'(a,es10.2)') 'optimization parameters:      condnum ',paropt%condnum
    write(*,'(a,es10.2)') 'optimization parameters: fnrmtolsatur ',paropt%fnrmtolsatur
    write(*,'(a,es10.2)') 'optimization parameters:     dt_start ',paropt%dt_start
    write(*,'(a,1x,i3)' ) 'optimization parameters:       nsatur ',paropt%nsatur
    write(*,'(a,1x,i5)' ) 'optimization parameters:          nit ',paropt%nit
    write(*,'(a,1x,l)'  ) 'optimization parameters:       lprint ',paropt%lprint
    paropt%param_reported=.true.
end subroutine vc_report_param
!*****************************************************************************************
