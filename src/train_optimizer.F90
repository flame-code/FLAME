!*****************************************************************************************
module mod_opt_ann
    implicit none
    !private
    public:: ekf_rivals, ekf_behler, ann_evaluate, init_opt_ann
    type, public:: typ_opt_ann
        integer:: n=-1
        integer:: loc(10)
        integer:: num(10)
        integer:: ndp_train=-1
        integer:: ndp_valid=-1
        real(8), allocatable:: x(:)
        real(8), allocatable:: epotd(:)
        real(8), allocatable:: g(:) !gradient of neural artificial neural network output
        real(8), allocatable:: gc(:,:)
        real(8), allocatable:: gs(:,:)
    end type typ_opt_ann
contains
!*****************************************************************************************
subroutine init_opt_ann(ndp_train,ndp_valid,opt_ann,ann_arr)
    use mod_ann, only: typ_ann_arr, ann_allocate
    use dynamic_memory
    implicit none
    integer, intent(in):: ndp_train
    integer, intent(in):: ndp_valid
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: istat, ng
    opt_ann%ndp_train=ndp_train
    opt_ann%ndp_valid=ndp_valid
    call ann_allocate(opt_ann%n,opt_ann%num,ann_arr)
end subroutine init_opt_ann
!*****************************************************************************************
subroutine ekf_rivals(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    !use mod_train, only: convert_x_ann
    use mod_processors, only: iproc, mpi_comm_abz
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    real(8), allocatable:: f(:) !Kalman gain matrix
    real(8), allocatable:: p(:,:) !covariance matrix
    real(8), allocatable:: v1(:) !work array
    integer:: i, j, iter, iconf, ios, ia
    real(8):: DDOT, tt, den, fcn_ann, fcn_ref
    real(8):: r, rinv, r0, rf, alpha
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    allocate(f(opt_ann%n),p(opt_ann%n,opt_ann%n),opt_ann%g(opt_ann%n),v1(opt_ann%n),opt_ann%epotd(opt_ann%num(1)))
    p(1:opt_ann%n,1:opt_ann%n)=0.d0
    do i=1,opt_ann%n
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
        r0=10.d0
        alpha=100.d-2
        rf=1.d-10
    endif
    if(parini%fit_hoppint) then
        call fit_hgen(parini,atoms_valid,ann_arr,opt_ann)
    endif
    do iter=0,parini%nstep_opt_ann
        call yaml_sequence(advance='no',unit=ann_arr%iunit)
        call cpu_time(time_s)
        !call randomize_data_order(atoms_train)
        do ia=1,ann_arr%n
            call convert_x_ann(opt_ann%num(ia),opt_ann%x(opt_ann%loc(ia)),ann_arr%ann(ia))
        enddo
        if(iproc==0) then
            if( ann_arr%exists_yaml_file) then
                call write_ann_all_yaml(parini,ann_arr,iter)
            else
                call write_ann_all(parini,ann_arr,iter)
            endif
        endif
        if(mod(iter,1)==0) then
            call analyze_epoch_init(parini,atoms_train,ann_arr)
            call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
            call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
            call analyze_epoch_print(parini,iter,atoms_train,ann_arr)
        endif
        if(iter==parini%nstep_opt_ann) exit
        call cpu_time(time1)
        dtime1=time1-time_s !time to get ANN energies for atoms_train and atoms_valid
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
        do iconf=1,atoms_train%nconf
            ann_arr%event='train'
            call get_fcn_ann(parini,iconf,'train',ann_arr,opt_ann,fcn_ann,fcn_ref)
            if(trim(parini%approach_ann)=='tb') then
                do j=1,opt_ann%n
                    opt_ann%g(j)=opt_ann%g(j)+parini%weight_hardness*opt_ann%x(j)
                enddo
            endif
            call cpu_time(time1)
            call DGEMV('T',opt_ann%n,opt_ann%n,1.d0,p,opt_ann%n,opt_ann%g,1,0.d0,v1,1)
            !call cal_matvec_mpi(opt_ann%n,p,opt_ann%g,v1)
            call cpu_time(time2)
            tt=DDOT(opt_ann%n,opt_ann%g,1,v1,1)
            den=1.d0/(tt+r)
            do j=1,opt_ann%n
                do i=1,opt_ann%n
                    p(i,j)=p(i,j)-v1(i)*v1(j)*den
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                    !write(21,'(2i5,es20.10)') i,j,p(i,j)
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,iconf,r,den
            call DGEMV('N',opt_ann%n,opt_ann%n,rinv,p,opt_ann%n,opt_ann%g,1,0.d0,f,1)
            do i=1,opt_ann%n
                !write(81,*) iter,iconf,i,f(i)*(epotall(iconf)-opt_ann%epot)
                !opt_ann%x(i)=opt_ann%x(i)+f(i)*(symfunc_train%symfunc(iconf)%epot-fcn_ann)
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
    deallocate(f,p,opt_ann%g,v1,opt_ann%epotd)
end subroutine ekf_rivals
!*****************************************************************************************
subroutine set_ref_energy(parini,atoms_train,atoms_ref,ind)
    use mod_parini, only: typ_parini
    !use mod_ann, only: typ_symfunc_arr
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    implicit none
    type(typ_parini), intent(in):: parini
    !type(typ_symfunc_arr), intent(in):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(in):: atoms_train
    !type(typ_atoms_arr), intent(in):: atoms_valid
    type(typ_atoms_arr), intent(inout):: atoms_ref
    integer, intent(out):: ind(200)
    !local variables
    !real(8), allocatable:: 
    integer:: iconf, mat
    atoms_ref%nconf=200
    allocate(atoms_ref%atoms(atoms_ref%nconf))
    do mat=1,200
        atoms_ref%atoms(mat)%epot=huge(1.d0)
    enddo
    write(*,*) 'atoms_train%nconf ',atoms_train%nconf
    ind(1:200)=0
    do iconf=1,atoms_train%nconf
        write(*,*) 'nat ',atoms_train%atoms(iconf)%nat
        mat=atoms_train%atoms(iconf)%nat
        if(atoms_train%atoms(iconf)%epot<atoms_ref%atoms(mat)%epot) then
            ind(mat)=iconf
            atoms_ref%atoms(mat)%epot=atoms_train%atoms(iconf)%epot
        endif
    enddo
    do mat=1,200
        if(ind(mat)==0) cycle
        iconf=ind(mat)
        call atom_copy_old(atoms_train%atoms(iconf),atoms_ref%atoms(mat),'copy to atoms_ref')
    enddo
end subroutine set_ref_energy
!*****************************************************************************************
subroutine analyze_epoch_init(parini,atoms_train,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms_arr), intent(in):: atoms_train
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    if(.not. (trim(ann_arr%approach)=='eem1' .or. trim(parini%approach_ann)=='cent1' .or. trim(ann_arr%approach)=='cent2')) return
    ann_arr%natsum(1:10)=0
    ann_arr%qmin(1:10)=huge(1.d0)
    ann_arr%qmax(1:10)=-huge(1.d0)
    ann_arr%qsum(1:10)=0.d0
    ann_arr%chi_min(1:10)=huge(1.d0)
    ann_arr%chi_max(1:10)=-huge(1.d0)
    ann_arr%chi_sum(1:10)=0.d0
    ann_arr%chi_delta(1:10)=0.d0
end subroutine analyze_epoch_init
!*****************************************************************************************
subroutine analyze_epoch_print(parini,iter,atoms_train,ann_arr)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann_arr
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_atoms_arr), intent(in):: atoms_train
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: i, ios
    real(8):: ttavg, ttmin, ttmax, ssavg, ssmin, ssmax
    character(50):: fn_charge, fn_chi
    character(20):: str_key
    if(.not. (trim(ann_arr%approach)=='eem1' .or. trim(parini%approach_ann)=='cent1' .or. trim(ann_arr%approach)=='cent2')) return
    do i=1,parini%ntypat
        !fn_charge='charge.'//trim(parini%stypat(i))
        !fn_chi='chi.'//trim(parini%stypat(i))
        !if(iter==0) then
        !    open(unit=61,file=trim(fn_charge),status='replace',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_charge)
        !        stop
        !    endif
        !    open(unit=71,file=trim(fn_chi),status='replace',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_chi)
        !        stop
        !    endif
        !else
        !    open(unit=61,file=trim(fn_charge),status='old',position='append',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_charge)
        !        stop
        !    endif
        !    open(unit=71,file=trim(fn_chi),status='old',position='append',iostat=ios)
        !    if(ios/=0) then
        !        write(*,'(2a)') 'ERROR: failure openning ',trim(fn_chi)
        !        stop
        !    endif
        !endif
        ttavg=ann_arr%qsum(i)/real(ann_arr%natsum(i),8)
        ttmin=ann_arr%qmin(i)
        ttmax=ann_arr%qmax(i)
        !write(61,'(i6,4f8.3)') iter,ttavg,ttmin,ttmax,ttmax-ttmin
        
        write(str_key,'(2a)') 'charge_',trim(parini%stypat(i))
        call yaml_mapping_open(trim(str_key),flow=.true.,unit=ann_arr%iunit)
        call yaml_map('iter',iter,unit=ann_arr%iunit)
        call yaml_map('qavg',ttavg,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('qmin',ttmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('qmax',ttmax,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('qvar',ttmax-ttmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_mapping_close(unit=ann_arr%iunit)

        !write(61,'(i6,4es14.5)') iter,ttavg,ttmin,ttmax,ttmax-ttmin
        ssavg=ann_arr%chi_sum(i)/real(ann_arr%natsum(i),8)
        ssmin=ann_arr%chi_min(i)
        ssmax=ann_arr%chi_max(i)
        !write(71,'(i6,5f8.3)') iter,ssavg,ssmin,ssmax,ssmax-ssmin,ann_arr%chi_delta(i)
        write(str_key,'(2a)') 'chi_',trim(parini%stypat(i))
        call yaml_mapping_open(trim(str_key),flow=.true.,unit=ann_arr%iunit)
        call yaml_map('iter',iter,unit=ann_arr%iunit)
        call yaml_map('chiavg',ssavg,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('chimin',ssmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('chimax',ssmax,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_map('chivar',ssmax-ssmin,fmt='(f8.3)',unit=ann_arr%iunit)
        call yaml_mapping_close(unit=ann_arr%iunit)
        !write(71,'(i6,4es14.5)') iter,ssavg,ssmin,ssmax,ssmax-ssmin
        !if (trim(parini%stypat(i))=='O' .and. ssmax-ssmin> 0.01) stop
        !close(61)
        !close(71)
    enddo
end subroutine analyze_epoch_print
!*****************************************************************************************
subroutine ekf_behler(parini,ann_arr,symfunc_train,symfunc_valid,atoms_train,atoms_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_processors, only: iproc, mpi_comm_abz
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_train, symfunc_valid
    type(typ_atoms_arr), intent(inout):: atoms_train
    type(typ_atoms_arr), intent(inout):: atoms_valid
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    real(8), allocatable:: f(:) !Kalman gain matrix
    real(8), allocatable:: p(:,:) !covariance matrix
    real(8), allocatable:: v1(:) !work array
    integer:: i, j, iter, iconf, ios, ia
    real(8):: DDOT, tt, den, alambda, alambdainv, alambda0
    real(8):: fcn_ann, fcn_ref
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    allocate(f(opt_ann%n),p(opt_ann%n,opt_ann%n),opt_ann%g(opt_ann%n),v1(opt_ann%n),opt_ann%epotd(opt_ann%num(1)))
    p(1:opt_ann%n,1:opt_ann%n)=0.d0
    do i=1,opt_ann%n
        p(i,i)=1.d-2
    enddo
    !alambda0=0.9994d0
    !alambda=0.99d0
    alambda0=0.9998d0
    alambda=0.997d0
    do iter=0,parini%nstep_opt_ann
        call cpu_time(time_s)
        !call randomize_data_order(atoms_train)
        do ia=1,ann_arr%n
            call convert_x_ann(opt_ann%num(ia),opt_ann%x(opt_ann%loc(ia)),ann_arr%ann(ia))
        enddo
        if(iproc==0) then
            if( ann_arr%exists_yaml_file) then
                call write_ann_all_yaml(parini,ann_arr,iter)
            else
                call write_ann_all(parini,ann_arr,iter) 
            endif
        endif
        if(mod(iter,1)==0) then
            call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
            call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
        endif
        if(iter==parini%nstep_opt_ann) exit
        call cpu_time(time1)
        dtime1=time1-time_s !time to get ANN energies for atoms_train and atoms_valid
        dtime2=0.d0 !time to convert ANN 1D array to ANN typ_ann
        dtime3=0.d0 !time to calculate ANN and its derivatives w.r.t. weights
        dtime4=0.d0 !time to convert derivative of ANN in typ_ann to 1D array
        dtime5=0.d0 !time to matrix-vector multiplication in Kalman filter
        dtime6=0.d0 !time of the rest of Kalman filter algorithm
        do iconf=1,atoms_train%nconf
            alambda=min(alambda0*alambda+1.d0-alambda0,0.999999d0)
            alambdainv=1.d0/alambda
            !write(31,'(2i6,f14.10)') iter,iconf,alambda
            ann_arr%event='train'
            call get_fcn_ann(parini,iconf,'train',ann_arr,opt_ann,fcn_ann,fcn_ref)
            call cpu_time(time1)
            !opt_ann%g(1:opt_ann%n)=opt_ann%g(1:opt_ann%n)*(opt_ann%epot-epotall(iconf))
            call DGEMV('T',opt_ann%n,opt_ann%n,1.d0,p,opt_ann%n,opt_ann%g,1,0.d0,v1,1)
            !call cal_matvec_mpi(opt_ann%n,p,opt_ann%g,v1)
            call cpu_time(time2)
            tt=DDOT(opt_ann%n,opt_ann%g,1,v1,1)
            den=alambdainv/(1.d0+tt*alambdainv)
            !write(*,'(i7,i5,f14.10,es11.2)') iter,iconf,alambda,den
            !call DGEMV('T',opt_ann%n,opt_ann%n,den,p,opt_ann%n,opt_ann%g,1,0.d0,f,1)
            f(1:opt_ann%n)=v1(1:opt_ann%n)*den
            do j=1,opt_ann%n
                do i=1,opt_ann%n
                    p(i,j)=(p(i,j)-f(i)*v1(j))*alambdainv
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,iconf,alambda,den
            do i=1,opt_ann%n
                !write(81,*) iter,iconf,i,f(i)*(epotall(iconf)-opt_ann%epot)
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
    deallocate(f,p,opt_ann%g,v1,opt_ann%epotd)
end subroutine ekf_behler
!*****************************************************************************************
subroutine ann_evaluate(parini,iter,ann_arr,symfunc_arr,atoms_arr,data_set,partb)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, typ_symfunc_arr, typ_symfunc
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
    use mod_processors, only: iproc
    use mod_tightbinding, only: typ_partb
    use futile
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_symfunc_arr), intent(inout):: symfunc_arr
    type(typ_atoms_arr), intent(inout):: atoms_arr
    character(*), intent(in):: data_set
    type(typ_partb), optional, intent(inout):: partb
    !local variables
    type(typ_atoms):: atoms
    type(typ_symfunc):: symfunc
    type(typ_opt_ann):: opt_ann
    real(8):: rmse, errmax, tt, pi
    real(8):: frmse, ttx, tty, ttz, ppx, ppy, ppz, tt1, tt2, tt3, ttn, tta, ttmax
    integer:: iconf, ierrmax, iat, nat_tot, nconf_force
    real(8):: time1=0.d0
    real(8):: time2=0.d0
    real(8), save:: time_p=0.d0
    real(8):: dtime1, dtime2
    integer:: ilarge1, ilarge2, ilarge3, iunit, ios
    !character(28):: frmt1='(i6,5f10.3,i7,i5,3i6,a40,i6)'
    !character(28):: frmt2='(i6,5e10.1,i7,i5,3i6,a40,i6)'
    !character(28):: frmt
    character(28):: fmt_main
    character(28):: fmt1='(f10.3)'
    character(28):: fmt2='(e10.1)'
    character(15):: filename
    logical:: file_exists
    character(len=8):: str_key
    call cpu_time(time1)
    pi=4.d0*atan(1.d0)
    rmse=0.d0
    frmse=0.d0
    ttn=0.d0
    tta=0.d0
    nat_tot=0
    nconf_force=0
    errmax=0.d0
    ierrmax=0
    ilarge1=0
    ilarge2=0
    ilarge3=0
    ann_arr%event='evalu'
    if(parini%save_symfunc_behnam) then
        ann_arr%compute_symfunc=.false.
    else
        ann_arr%compute_symfunc=.true.
    endif
    if(parini%print_energy) then
        write(filename,'(a12,i3.3)') 'detailed_err',iter
        iunit=f_get_free_unit(10**5)
        if(trim(data_set)=="train") then
            open(unit=iunit,file=trim(filename),status='unknown',iostat=ios)
        elseif(trim(data_set)=="valid") then
            open(unit=iunit,file=trim(filename),status='old',access='append',iostat=ios)
        endif
        if(ios/=0) then
            write(*,'(a,a)') 'ERROR: failure openning file: ',trim(filename)
            stop
        endif
        !write(iunit,'(a2,a44,4a23)') "#", " ","E_dft","E_ann","E_dft-E_ann/atom (Ha)","E_dft-E_ann (eV)" 
    endif
    configuration: do iconf=1,atoms_arr%nconf
        if(.not. atoms_arr%conf_inc(iconf)) cycle
        call atom_copy_old(atoms_arr%atoms(iconf),atoms,'atoms_arr%atoms(iconf)->atoms')
        if(parini%save_symfunc_behnam) then
            call cal_ann_main(parini,atoms,symfunc_arr%symfunc(iconf),ann_arr,opt_ann)
        else
            call cal_ann_main(parini,atoms,symfunc,ann_arr,opt_ann)
        endif
        !if(iter==parini%nstep_opt_ann) then
        !    write(40+ifile,'(2i6,2es24.15,es14.5)') iconf,atoms_arr%atoms(iconf)%nat, &
        !        atoms_arr%atoms(iconf)%epot/atoms_arr%atoms(iconf)%nat,atoms%epot/atoms_arr%atoms(iconf)%nat, &
        !        (atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        !endif
        tt=abs(atoms%epot-atoms_arr%atoms(iconf)%epot)/atoms_arr%atoms(iconf)%nat
        !HERE
        if(parini%print_energy) then
            write(iunit,'(i7,es14.5,a40,i6,a)') iconf,tt,trim(atoms_arr%fn(iconf)),atoms_arr%lconf(iconf),trim(data_set)
        endif
        if(tt>1.d-2) ilarge1=ilarge1+1
        if(tt>1.d-3) ilarge2=ilarge2+1
        if(tt>1.d-4) ilarge3=ilarge3+1
        if(tt>errmax) then
            errmax=tt
            ierrmax=iconf
        endif
        rmse=rmse+tt**2
        !if(tt>1.d-2) then
        !    atoms_arr%inclusion(iconf)=0
        !else
        !    atoms_arr%inclusion(iconf)=1
        !endif
        !write(22,'(a,i5.5)') 'configuration ',iconf
        !if(atoms%nat<=parini%nat_force) then
        nat_tot=nat_tot+atoms%nat
        nconf_force=nconf_force+1
        do iat=1,atoms%nat
            ttx=atoms_arr%atoms(iconf)%fat(1,iat)
            tty=atoms_arr%atoms(iconf)%fat(2,iat)
            ttz=atoms_arr%atoms(iconf)%fat(3,iat)
            ppx=atoms%fat(1,iat)
            ppy=atoms%fat(2,iat)
            ppz=atoms%fat(3,iat)
            !write(41,'(2i6,6f7.3)') iconf,iat,ttx,tty,ttz,ppx,ppy,ppz
            tt1=sqrt(ttx**2+tty**2+ttz**2)
            tt2=sqrt(ppx**2+ppy**2+ppz**2)
            tt3=(ppx-ttx)**2+(ppy-tty)**2+(ppz-ttz)**2
            frmse=frmse+tt3
            tt3=sqrt(tt3)
        enddo
        ttn=ttn+ann_arr%fchi_norm
        tta=tta+ann_arr%fchi_angle
        !write(44,'(2i7,4es14.5)') iter,iconf,ann_arr%fchi_norm,ann_arr%fchi_angle,ttn/nconf_force,tta/nconf_force
        !endif
    enddo configuration
    rmse=sqrt(rmse/real(atoms_arr%nconf_inc,8))
    if(nconf_force==0) nconf_force=1
    ttn=ttn/real(nconf_force,8)
    tta=tta/real(nconf_force,8)
    if(nat_tot==0) nat_tot=1
    frmse=sqrt(frmse/real(3*nat_tot,8))
    if(iproc==0) then
        rmse=rmse*1.d3
        errmax=errmax*1.d3
        if(rmse>99999.d0) then
            !frmt=frmt2
            fmt_main=fmt2
        else
            !frmt=frmt1
            fmt_main=fmt1
        endif

        !write(ifile,frmt) iter,rmse,ttn,tta,frmse,errmax,ierrmax,atoms_arr%atoms(ierrmax)%nat, &
        !    ilarge1,ilarge2,ilarge3,trim(atoms_arr%fn(ierrmax)),atoms_arr%lconf(ierrmax)
        write(str_key,'(a)') trim(data_set)
        call yaml_mapping_open(trim(str_key),flow=.true.,unit=ann_arr%iunit)
        call yaml_map('iter',iter,unit=ann_arr%iunit)
        call yaml_map('rmse',rmse,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('ttn',ttn,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('tta',tta,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('frmse',frmse,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('errmax',errmax,fmt=trim(fmt_main),unit=ann_arr%iunit)
        call yaml_map('ierrmax',ierrmax,unit=ann_arr%iunit)
        call yaml_map('nat',atoms_arr%atoms(ierrmax)%nat,unit=ann_arr%iunit)
        call yaml_map('ilarge1',ilarge1,unit=ann_arr%iunit)
        call yaml_map('ilarge2',ilarge2,unit=ann_arr%iunit)
        call yaml_map('ilarge3',ilarge3,unit=ann_arr%iunit)
        call yaml_map('fn_ierrmax',trim(atoms_arr%fn(ierrmax)),unit=ann_arr%iunit)
        call yaml_map('lconf_ierrmax',atoms_arr%lconf(ierrmax),unit=ann_arr%iunit)
        call yaml_mapping_close(unit=ann_arr%iunit)
    endif
    call cpu_time(time2)
    dtime1=time1-time_p
    dtime2=time2-time1
    !write(*,'(a,2f20.2)') 'TIME ',dtime1,dtime2
    time_p=time2
    if(parini%print_energy) then
        close(iunit)
    endif
    ann_arr%compute_symfunc=.false.
end subroutine ann_evaluate
!*****************************************************************************************
end module mod_opt_ann
!*****************************************************************************************
