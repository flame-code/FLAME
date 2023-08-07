!*****************************************************************************************
module mod_opt_ann
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, convert_x_ann_arr
    use mod_flm_futile
    implicit none
    private
    public:: ekf_rivals
    public:: init_x_opt_ann
    public:: set_opt_ann_x !delete in future: both are used only for TB
    public:: ann_lm
    public:: typ_cost_object
    type, public:: typ_opt_ann
        private
        integer, public:: n=-1
        integer:: ndp_train=-1
        integer:: iunit=6
        logical:: x_initialized=.false.
        real(8), allocatable:: x(:)
        !IMPORTANT: g(:) must be changed to private once cal_ann_tb does not need it anymore.
        real(8), public, allocatable:: g(:) !gradient of neural artificial neural network output
        contains
        procedure, public, pass(self):: init_opt_ann
        procedure, public, pass(self):: fini_opt_ann
        procedure, public, pass(self):: init_x_opt_ann
        procedure, public, pass(self):: get_opt_ann_x
    end type typ_opt_ann
    type, abstract:: typ_cost_object
        contains
        procedure(routine_value), public, pass(self), deferred:: func_value
        procedure(routine_write), public, pass(self), deferred:: func_write
        procedure(routine_evaluate), public, pass(self), deferred:: func_evaluate
    end type
    abstract interface
    subroutine routine_value(self,parini,idp,opt_ann,fcn_ann,fcn_ref,g)
        import:: typ_cost_object, typ_opt_ann, typ_parini
        implicit none
        class(typ_cost_object), intent(inout):: self
        type(typ_parini), intent(in):: parini
        integer, intent(in):: idp
        type(typ_opt_ann), intent(inout):: opt_ann
        real(8), intent(out):: fcn_ann
        real(8), intent(out):: fcn_ref
        real(8), intent(out):: g(opt_ann%n)
    end subroutine routine_value
    subroutine routine_write(self,parini,iter)
        import:: typ_cost_object, typ_parini
        implicit none
        class(typ_cost_object), intent(inout):: self
        type(typ_parini), intent(in):: parini
        integer, intent(in):: iter
    end subroutine routine_write
    subroutine routine_evaluate(self,parini,iter)
        import:: typ_cost_object, typ_parini
        implicit none
        class(typ_cost_object), intent(inout):: self
        type(typ_parini), intent(in):: parini
        integer, intent(in):: iter
    end subroutine routine_evaluate
    end interface
contains
!*****************************************************************************************
subroutine init_opt_ann(self,ndp_train,n,iunit)
    implicit none
    class(typ_opt_ann), intent(inout):: self
    integer, intent(in):: ndp_train, n, iunit
    !local variables
    self%n=n
    self%ndp_train=ndp_train
    self%g=f_malloc0([1.to.self%n],id='opt_ann%g')
    self%x=f_malloc0([1.to.self%n],id='opt_ann%x')
    self%iunit=iunit
end subroutine init_opt_ann
!*****************************************************************************************
subroutine fini_opt_ann(self)
    implicit none
    class(typ_opt_ann), intent(inout):: self
    !local variables
    call f_free(self%g)
    call f_free(self%x)
end subroutine fini_opt_ann
!*****************************************************************************************
subroutine init_x_opt_ann(self,n,x)
    implicit none
    class(typ_opt_ann), intent(inout):: self
    integer, intent(in):: n
    real(8), intent(in):: x(n)
    !local variables
    if(n/=self%n) then
        write(*,'(a)') 'ERROR: n/=self%n in init_x_opt_ann'
        write(*,'(a,2i8)') 'n,self%n= ',n,self%n
        stop
    endif
    self%x=x
    self%x_initialized=.true.
end subroutine init_x_opt_ann
!*****************************************************************************************
subroutine set_opt_ann_x(ann_arr,x_arr,opt_ann)
    implicit none
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    real(8), intent(in):: x_arr(ann_arr%nweight_max,ann_arr%nann)
    !local variables
    integer:: i, j
    do i=1,ann_arr%nann
        do j=1,ann_arr%num(i)
            opt_ann%x(ann_arr%loc(i)+j-1)=x_arr(j,i)
        enddo
    enddo
end subroutine set_opt_ann_x
!*****************************************************************************************
subroutine get_opt_ann_x(self,n,x)
    implicit none
    class(typ_opt_ann), intent(in):: self
    integer, intent(in):: n
    real(8), intent(out):: x(n)
    !local variables
    integer:: i, j
    if(n/=self%n) then
        write(*,'(a)') 'ERROR: n/=self%n in get_opt_ann_x'
        write(*,'(a,2i8)') 'n,self%n= ',n,self%n
        stop
    endif
    x=self%x
end subroutine get_opt_ann_x
!*****************************************************************************************
subroutine ekf_rivals(cost_object,parini,opt_ann)
    use mod_parini, only: typ_parini
    use mod_processors, only: iproc
    implicit none
    class(typ_cost_object), intent(inout) :: cost_object
    type(typ_parini), intent(in):: parini
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    real(8), allocatable:: f(:) !Kalman gain matrix
    real(8), allocatable:: p(:,:) !covariance matrix
    real(8), allocatable:: v1(:) !work array
    integer:: i, j, iter, idp, ios, ia
    real(8):: DDOT, tt, den, fcn_ann, fcn_ref
    real(8):: r, rinv, r0, rf, alpha
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    integer:: is, ie, im, mproc, jproc, ierr
    integer, allocatable:: ndispls(:) ,ncounts(:)
#if defined(MPI)
    include 'mpif.h'
#endif
    if(.not. opt_ann%x_initialized) then
        write(*,'(a)') 'ERROR: opt_ann%x is not initialized.'
        stop
    endif
    if(parini%mpi_env%nproc>1) then
        allocate(ndispls(parini%mpi_env%nproc),ncounts(parini%mpi_env%nproc))
        do jproc=0,parini%mpi_env%nproc-1
            im=opt_ann%n/parini%mpi_env%nproc
            is=jproc*im+1
            mproc=mod(opt_ann%n,parini%mpi_env%nproc)
            is=is+max(0,jproc-parini%mpi_env%nproc+mproc)
            if(jproc>parini%mpi_env%nproc-mproc-1) im=im+1
            ie=is+im-1
            ncounts(jproc+1)=ie-is+1
            ndispls(jproc+1)=is-1
        enddo
        im=opt_ann%n/parini%mpi_env%nproc
        is=parini%mpi_env%iproc*im+1
        mproc=mod(opt_ann%n,parini%mpi_env%nproc)
        is=is+max(0,parini%mpi_env%iproc-parini%mpi_env%nproc+mproc)
        if(parini%mpi_env%iproc>parini%mpi_env%nproc-mproc-1) im=im+1
        ie=is+im-1
    else
        is=1
        ie=opt_ann%n
    endif
    allocate(f(opt_ann%n),source=0.d0)
    allocate(v1(opt_ann%n),source=0.d0)
    allocate(p(opt_ann%n,is:ie))
    p(1:opt_ann%n,is:ie)=0.d0
    do i=is,ie
        p(i,i)=1.d-2
    enddo
    if(trim(parini%approach_ann)=='eem1' .or. trim(parini%approach_ann)=='cent1') then
        r0=10.d0
        alpha=100.d-2
        rf=1.d-6
    elseif(trim(parini%approach_ann)=='cent2') then
        r0=10.d0
        alpha=100.d-2
        rf=1.d-6
    elseif(trim(parini%approach_ann)=='centt') then
        r0=10.d0
        alpha=100.d-2
        rf=1.d-6
    elseif(trim(parini%approach_ann)=='cent3') then
        r0=10.d0
        alpha=100.d-2
        rf=1.d-6
    elseif(trim(parini%approach_ann)=='tb') then
        r0=100.d0
        alpha=100.d-2
        rf=1.d-5
        !alpha=20.d-2
        !rf=1.d-10
    else
        !r0=1.d0
        !alpha=5.d-1
        !rf=1.d-8
        r0=100.d0
        alpha=100.d-2
        rf=1.d-6
    endif
    do iter=0,parini%nstep_opt_ann
        call yaml_sequence(advance='no',unit=opt_ann%iunit)
        call cpu_time(time_s)
        if(iproc==0) then
            call cost_object%func_write(parini,iter)
        endif
        if(mod(iter,1)==0) then
            call cost_object%func_evaluate(parini,iter)
        endif
        if(iter==parini%nstep_opt_ann) exit
        call cpu_time(time1)
        dtime1=time1-time_s !time to get ANN energies for train and valid
        dtime2=0.d0 !time to convert ANN 1D array to ANN typ_ann
        dtime3=0.d0 !time to calculate ANN and its derivatives w.r.t. weights
        dtime4=0.d0 !time to convert derivative of ANN in typ_ann to 1D array
        dtime5=0.d0 !time to matrix-vector multiplication in Kalman filter
        dtime6=0.d0 !time of the rest of Kalman filter algorithm
       ! if (.not. parini%restart_param) then
       !     r=(r0-rf)*exp(-alpha*(iter))+rf
       ! else
            r=(r0-rf)*exp(-alpha*(iter+parini%restart_iter))+rf
       ! endif
        rinv=1.d0/r
        write(31,'(i6,es14.5)') iter,r
        do idp=1,opt_ann%ndp_train
            !The following is allocated with ann_arr%num(1), this means number of
            !nodes in the input layer is the same for all atom types.
            !Therefore, it must be fixed later.
            !g_per_atom=f_malloc([1.to.ann_arr%num(1),1.to.atoms%nat],id='g_per_atom') !HERE
            call cost_object%func_value(parini,idp,opt_ann,fcn_ann,fcn_ref,opt_ann%g)
            !if(trim(parini%approach_ann)=='tb') then
                do j=1,opt_ann%n
                    opt_ann%g(j)=opt_ann%g(j)+parini%weight_hardness*opt_ann%x(j)
                enddo
            !endif
            call cpu_time(time1)
            call DGEMV('T',opt_ann%n,ie-is+1,1.d0,p(1,is),opt_ann%n,opt_ann%g,1,0.d0,v1(is),1)
            if(parini%mpi_env%nproc>1) then
            call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,v1,ncounts,ndispls,MPI_DOUBLE_PRECISION,parini%mpi_env%mpi_comm,ierr)
            endif
            !if(parini%mpi_env%iproc==0) then
            !    do i=1,opt_ann%n
            !        write((parini%mpi_env%iproc+1)*1000+iter,'(i6,es19.10)') i,v1(i)
            !    enddo
            !endif
            !call cal_matvec_mpi(opt_ann%n,p,opt_ann%g,v1)
            call cpu_time(time2)
            tt=DDOT(opt_ann%n,opt_ann%g,1,v1,1)
            den=1.d0/(tt+r)
            do j=is,ie
                do i=1,opt_ann%n
                    p(i,j)=p(i,j)-v1(i)*v1(j)*den
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                    !write(21,'(2i5,es20.10)') i,j,p(i,j)
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,idp,r,den
            call DGEMV('T',opt_ann%n,ie-is+1,rinv,p(1,is),opt_ann%n,opt_ann%g,1,0.d0,f(is),1)
            if(parini%mpi_env%nproc>1) then
            call MPI_ALLGATHERV(MPI_IN_PLACE,0,MPI_DATATYPE_NULL,f,ncounts,ndispls,MPI_DOUBLE_PRECISION,parini%mpi_env%mpi_comm,ierr)
            endif
            do i=1,opt_ann%n
                !write(81,*) iter,idp,i,f(i)*(epotall(idp)-opt_ann%epot)
                opt_ann%x(i)=opt_ann%x(i)+f(i)*(fcn_ref-fcn_ann)
            enddo
            call cpu_time(time3)
            dtime5=dtime5+time2-time1
            dtime6=dtime6+time3-time2
        enddo
        call cpu_time(time_e)
        dtime=time_e-time_s
        !dtime=dtime1+dtime2+dtime3+dtime4+dtime5+dtime6
        tt1=dtime1/dtime*100.d0
        tt2=dtime2/dtime*100.d0
        tt3=dtime3/dtime*100.d0
        tt4=dtime4/dtime*100.d0
        tt5=dtime5/dtime*100.d0
        tt6=dtime6/dtime*100.d0
        tt=tt1+tt2+tt3+tt4+tt5+tt6
        if(iter==0) then
            open(unit=41,file='time_real.prc',status='replace',iostat=ios)
            if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
            write(41,'(a1,4x,7a10)') '#',' dtime1 ',' dtime2 ',' dtime3 ',' dtime4 ',' dtime5 ',' dtime6 ',' sum '
            open(unit=42,file='time_frac.prc',status='replace',iostat=ios)
            if(ios/=0) then;write(*,'(a)') 'ERROR: failure openning output file';stop;endif
            write(42,'(a1,4x,7a10)') '#','    tt1 ','    tt2 ','    tt3 ','    tt4 ','    tt5 ','    tt6 ', 'sum '
        endif
        write(41,'(i5,7f10.2)') iter,dtime1,dtime2,dtime3,dtime4,dtime5,dtime6,dtime
        write(42,'(i5,7f10.1)') iter,tt1,tt2,tt3,tt4,tt5,tt6,tt
    enddo
    close(41)
    close(42)
    deallocate(f,p,v1)
    if(parini%mpi_env%nproc>1) then
        deallocate(ndispls,ncounts)
    endif
end subroutine ekf_rivals
!*****************************************************************************************
subroutine ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_symfunc, only: typ_symfunc_arr
    use mod_parlm, only: typ_parlm
    use mod_atoms, only: typ_atoms_arr
    implicit none
    type(typ_parini), intent(in):: parini
    !local variables
    type(typ_ann_arr):: ann_arr
    type(typ_opt_ann):: opt_ann
    type(typ_atoms_arr):: atoms_train
    type(typ_atoms_arr):: atoms_valid
    type(typ_symfunc_arr):: symfunc_train
    type(typ_symfunc_arr):: symfunc_valid
    type(typ_parlm):: parlm
    real(8):: tt, epot
    integer:: m, iann
    if(parini%fit_hoppint) then
        call fit_hgen(parini,ann_arr,opt_ann)
    endif
    parlm%xtol=1.d-8
    parlm%ftol=1.d-8
    parlm%gtol=1.d-8
    m=atoms_train%nconf
    parlm%n=opt_ann%n
    call init_lmder_modified(parlm,m,m)
    parlm%x(1:parlm%n)=opt_ann%x(1:parlm%n)
    do
        call lmder_modified(parlm,m,m)
        write(*,*) 'nfev,njev: ',parlm%nfev,parlm%njev
        if(parlm%finish) exit
        !do i=1,parlm%n
        !    write(500+istep,'(i4,2es19.10)') i,parlm%x(i),parlm%wa2(i)
        !enddo
        !iann=1
        if(parlm%icontinue==700) then
            call convert_opt_x_ann_arr(opt_ann,ann_arr)
            !opt_ann%x(1:parlm%n)=parlm%wa2(1:parlm%n)
            call fcn_epot(m,parlm%n,parlm%wa2,parlm%wa4,parlm%fjac,m,parlm%iflag, &
                parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,opt_ann)
        else
            call convert_opt_x_ann_arr(opt_ann,ann_arr)
            !opt_ann%x(1:parlm%n)=parlm%x(1:parlm%n)
            call fcn_epot(m,parlm%n,parlm%x,parlm%fvec,parlm%fjac,m,parlm%iflag, &
                parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,opt_ann)
        endif
        !istep=istep+1
    enddo
    opt_ann%x(1:parlm%n)=parlm%x(1:parlm%n)
    write(*,*) 'info= ',parlm%info
    call final_lmder_modified(parlm)
end subroutine ann_lm
!*****************************************************************************************
subroutine fcn_epot(m,n,x,fvec,fjac,ldfjac,iflag,parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_symfunc, only: typ_symfunc_arr
    !use mod_opt_ann, only: typ_opt_ann
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old, atom_deallocate_old
    !use mod_opt_ann, only: ann_evaluate
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_atoms_arr), intent(inout):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_opt_ann), intent(inout):: opt_ann
    integer:: m, n, ldfjac, iflag
    real(8):: x(n), fvec(m), fjac(ldfjac,n)
    !local variables
    type(typ_atoms):: atoms
    integer:: iconf, iter, ia
    integer, save:: icall=0
    integer, save:: icall0=0
    icall=icall+1
    !opt_ann%x(1:n)=x(1:n)
    !call convert_opt_x_ann_arr(opt_ann,ann_arr)
    write(*,'(a24,i8,a9,i8,a14)') '**************** icall= ',icall,'  iflag= ',iflag,'  ************'
    if(iflag==1) then
        ann_arr%event='evalu'
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,opt_ann)
            fvec(iconf)=(atoms%epot-atoms_train%atoms(iconf)%epot)**2
            call atom_deallocate_old(atoms)
        enddo
    elseif(iflag==2) then
        ann_arr%event='train'
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,opt_ann)
            fjac(iconf,1:opt_ann%n)=opt_ann%g(1:opt_ann%n)*(atoms%epot-atoms_train%atoms(iconf)%epot)*2.d0
            call atom_deallocate_old(atoms)
        enddo
    elseif(iflag==0) then
        iter=icall0
        call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
        call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
        icall0=icall0+1
    endif
end subroutine fcn_epot
!*****************************************************************************************
end module mod_opt_ann
!*****************************************************************************************
