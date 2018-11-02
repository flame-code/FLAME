!*****************************************************************************************
subroutine ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_ekf, only: typ_ekf
    use mod_parlm, only: typ_parlm
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr):: ann_arr
    type(typ_ekf):: ekf
    type(typ_atoms_arr):: atoms_train
    type(typ_atoms_arr):: atoms_valid
    type(typ_symfunc_arr):: symfunc_train
    type(typ_symfunc_arr):: symfunc_valid
    type(typ_parlm):: parlm
    real(8):: tt, epot
    integer:: m, iann
    allocate(ekf%g(ekf%n)) !,v1(ekf%n),ekf%epotd(ekf%num(1)))
    if(parini%fit_hoppint) then
        call fit_hgen(parini,atoms_valid,ann_arr,ekf)
    endif
    parlm%xtol=1.d-8
    parlm%ftol=1.d-8
    parlm%gtol=1.d-8
    m=atoms_train%nconf
    parlm%n=ekf%n
    call init_lmder_modified(parlm,m,m)
    parlm%x(1:parlm%n)=ekf%x(1:parlm%n)
    do
        call lmder_modified(parlm,m,m)
        write(*,*) 'nfev,njev: ',parlm%nfev,parlm%njev
        if(parlm%finish) exit
        !do i=1,parlm%n
        !    write(500+istep,'(i4,2es19.10)') i,parlm%x(i),parlm%wa2(i)
        !enddo
        !iann=1
        if(parlm%icontinue==700) then
            do iann=1,ann_arr%n
                call convert_x_ann(ekf%num(iann),parlm%wa2(ekf%loc(iann)),ann_arr%ann(iann))
            enddo
            !ekf%x(1:parlm%n)=parlm%wa2(1:parlm%n)
            call fcn_epot(m,parlm%n,parlm%wa2,parlm%wa4,parlm%fjac,m,parlm%iflag, &
                parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
        else
            do iann=1,ann_arr%n
                call convert_x_ann(ekf%num(iann),parlm%x(ekf%loc(iann)),ann_arr%ann(iann))
            enddo
            !ekf%x(1:parlm%n)=parlm%x(1:parlm%n)
            call fcn_epot(m,parlm%n,parlm%x,parlm%fvec,parlm%fjac,m,parlm%iflag, &
                parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
        endif
        !istep=istep+1
    enddo
    ekf%x(1:parlm%n)=parlm%x(1:parlm%n)
    write(*,*) 'info= ',parlm%info
    call final_lmder_modified(parlm)
end subroutine ann_lm
!*****************************************************************************************
subroutine fcn_epot(m,n,x,fvec,fjac,ldfjac,iflag,parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_ekf, only: typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    use mod_ekf, only: ann_evaluate
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_ekf), intent(inout):: ekf
    integer:: m, n, ldfjac, iflag
    real(8):: x(n), fvec(m), fjac(ldfjac,n)
    !local variables
    type(typ_atoms):: atoms
    integer:: iconf, iter, ia
    integer, save:: icall=0
    integer, save:: icall0=0
    icall=icall+1
    !ekf%x(1:n)=x(1:n)
    !do ia=1,ann_arr%n
    !    call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
    !enddo
    write(*,'(a,i,a,i,a)') '**************** icall= ',icall,'  iflag= ',iflag,'  ************'
    if(iflag==1) then
        ann_arr%event='evalu'
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            fvec(iconf)=(atoms%epot-symfunc_train%symfunc(iconf)%epot)**2
        enddo
    elseif(iflag==2) then
        ann_arr%event='train'
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            fjac(iconf,1:ekf%n)=ekf%g(1:ekf%n)*(atoms%epot-symfunc_train%symfunc(iconf)%epot)*2.d0
        enddo
    elseif(iflag==0) then
        iter=icall0
        call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
        call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
        icall0=icall0+1
    endif
end subroutine fcn_epot
!*****************************************************************************************
