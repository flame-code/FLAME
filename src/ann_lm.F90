!*****************************************************************************************
subroutine ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
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
    real(8):: tt, epot

    integer:: m, n
    real(8), allocatable:: fvec(:), fjac(:,:), diag(:), qtf(:), wa1(:), wa2(:), wa3(:), wa4(:), x(:)
    integer, allocatable:: ipvt(:)
    real(8), parameter:: ftol=1.d-10
    real(8), parameter:: xtol=1.d-10
    real(8), parameter:: gtol=1.d-10
    integer, parameter:: maxfev=500
    real(8), parameter:: factor=100.d0
    integer, parameter:: nprint=1
    integer:: info, nfev, njev
    allocate(ekf%g(ekf%n))
    m=atoms_train%nconf
    n=ekf%n
    allocate(fvec(m),fjac(m,n),diag(n),qtf(n),wa1(n),wa2(n),wa3(n),wa4(m),ipvt(n),x(n))
    x(1:n)=ekf%x(1:n)
    fvec(1:m)=0.d0 ; fjac(1:m,1:n)=0.d0 ; diag(1:n)=0.d0 ; qtf(1:n)=0.d0
    wa1(1:n)=0.d0 ; wa2(1:n)=0.d0 ; wa3(1:n)=0.d0 ; wa4(1:m)=0.d0 ; ipvt(1:n)=0
    call lmder_modified(fcn_least_squares,m,n,x,fvec,fjac,m,ftol,xtol,gtol,maxfev,diag,1,factor, &
        nprint,info,nfev,njev,ipvt,qtf,wa1,wa2,wa3,wa4, &
        parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    write(*,'(a,3i6)') '#nfev,njev,info ',nfev,njev,info
    deallocate(fvec,fjac,diag,qtf,wa1,wa2,wa3,wa4,ipvt,x)
    deallocate(ekf%g)
end subroutine ann_lm
!*****************************************************************************************
subroutine fcn_least_squares(m,n,x,fvec,fjac,ldfjac,iflag,parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,ekf)
    use mod_interface
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_ekf
    use mod_atoms, only: typ_atoms, typ_atoms_arr
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
    ekf%x(1:n)=x(1:n)
    do ia=1,ann_arr%n
        call convert_x_ann(ekf%num(ia),ekf%x(ekf%loc(ia)),ann_arr%ann(ia))
    enddo
    write(*,'(a,i,a,i,a)') '**************** icall= ',icall,'  iflag= ',iflag,'  ************'
    if(iflag==1) then
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,ekf)
            fvec(iconf)=(atoms%epot-symfunc_train%symfunc(iconf)%epot)**2
        enddo
    elseif(iflag==2) then
        do iconf=1,atoms_train%nconf
            call cal_ann_main(parini,atoms_train%atoms(iconf),symfunc_train%symfunc(iconf),ann_arr,ekf)
            fjac(iconf,1:ekf%n)=ekf%g(1:ekf%n)*(atoms%epot-symfunc_train%symfunc(iconf)%epot)*2.d0
        enddo
    elseif(iflag==0) then
        iter=icall0
        call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,11)
        call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,12)
        icall0=icall0+1
    endif
end subroutine fcn_least_squares
!*****************************************************************************************
