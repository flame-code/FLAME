!*****************************************************************************************
module mod_opt_ann
    implicit none
    !private
    public:: ekf_rivals, ekf_behler, init_opt_ann, fini_opt_ann
    public:: set_opt_ann_grad, set_annweights
    public:: get_opt_ann_x, set_opt_ann_x !delete in future: both are used only for TB
    public:: ann_lm, convert_opt_x_ann_arr
    type, public:: typ_opt_ann
        private
        integer, public:: n=-1
        integer:: ndp_train=-1
        real(8), allocatable:: x(:)
        real(8), allocatable:: g(:) !gradient of neural artificial neural network output
    end type typ_opt_ann
contains
!*****************************************************************************************
subroutine init_opt_ann(ndp_train,opt_ann,ann_arr)
    use mod_ann, only: typ_ann_arr
    use dynamic_memory
    implicit none
    integer, intent(in):: ndp_train
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    opt_ann%n=sum(ann_arr%num(1:ann_arr%nann))
    opt_ann%ndp_train=ndp_train
    opt_ann%g=f_malloc0([1.to.opt_ann%n],id='opt_ann%g')
    opt_ann%x=f_malloc0([1.to.opt_ann%n],id='opt_ann%x')
end subroutine init_opt_ann
!*****************************************************************************************
subroutine fini_opt_ann(opt_ann)
    use mod_ann, only: typ_ann_arr
    use dynamic_memory
    implicit none
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    call f_free(opt_ann%g)
    call f_free(opt_ann%x)
end subroutine fini_opt_ann
!*****************************************************************************************
subroutine set_opt_ann_grad(ann_arr,grad,opt_ann)
    use mod_ann, only: typ_ann_arr
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    real(8), intent(in):: grad(ann_arr%nweight_max,ann_arr%nann)
    !local variables
    integer:: i, j
    do i=1,ann_arr%nann
        do j=1,ann_arr%num(i)
            opt_ann%g(ann_arr%loc(i)+j-1)=grad(j,i)
        enddo
    enddo
end subroutine set_opt_ann_grad
!*****************************************************************************************
subroutine set_opt_ann_x(ann_arr,x_arr,opt_ann)
    use mod_ann, only: typ_ann_arr
    use dynamic_memory
    use yaml_output
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
subroutine get_opt_ann_x(ann_arr,opt_ann,x_arr)
    use mod_ann, only: typ_ann_arr
    use dynamic_memory
    use yaml_output
    implicit none
    type(typ_ann_arr), intent(in):: ann_arr
    type(typ_opt_ann), intent(in):: opt_ann
    real(8), intent(out):: x_arr(ann_arr%nweight_max,ann_arr%nann)
    !local variables
    integer:: i, j
    do i=1,ann_arr%nann
        do j=1,ann_arr%num(i)
            x_arr(j,i)=opt_ann%x(ann_arr%loc(i)+j-1)
        enddo
    enddo
end subroutine get_opt_ann_x
!*****************************************************************************************
subroutine convert_opt_x_ann_arr(opt_ann,ann_arr)
    use mod_ann, only: typ_ann_arr, convert_x_ann_arr
    implicit none
    type(typ_opt_ann), intent(in):: opt_ann
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    call convert_x_ann_arr(opt_ann%n,opt_ann%x,ann_arr)
end subroutine convert_opt_x_ann_arr
!*****************************************************************************************
subroutine ekf_rivals(parini,ann_arr,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_processors, only: iproc
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
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
    allocate(f(opt_ann%n),p(opt_ann%n,opt_ann%n),v1(opt_ann%n))
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
    if(parini%fit_hoppint) then
        call fit_hgen(parini,ann_arr,opt_ann)
    endif
    do iter=0,parini%nstep_opt_ann
        call yaml_sequence(advance='no',unit=ann_arr%iunit)
        call cpu_time(time_s)
        call convert_opt_x_ann_arr(opt_ann,ann_arr)
        if(iproc==0) then
            if( ann_arr%exists_yaml_file) then
                call write_ann_all_yaml(parini,ann_arr,iter)
            else
                call write_ann_all(parini,ann_arr,iter)
            endif
        endif
        if(mod(iter,1)==0) then
            call analyze_epoch_init(parini,ann_arr)
            call ann_evaluate_all(parini,iter,ann_arr)
            call analyze_epoch_print(parini,iter,ann_arr)
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
        ann_arr%dpm_rmse=0
        do idp=1,opt_ann%ndp_train
            ann_arr%event='train'
            !The following is allocated with ann_arr%num(1), this means number of
            !nodes in the input layer is the same for all atom types.
            !Therefore, it must be fixed later.
            !g_per_atom=f_malloc([1.to.ann_arr%num(1),1.to.atoms%nat],id='g_per_atom') !HERE
            call convert_opt_x_ann_arr(opt_ann,ann_arr)
            call get_fcn_ann(parini,idp,'train',ann_arr,opt_ann,fcn_ann,fcn_ref)
            ann_arr%dpm_rmse=ann_arr%dpm_rmse+ann_arr%dpm_err
            if(parini%iverbose>=2) then 
                write(1370,'(2i4,2es18.8)') iter,idp,ann_arr%dpm_err, ann_arr%dpm_rmse
            end if
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
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,idp,r,den
            call DGEMV('N',opt_ann%n,opt_ann%n,rinv,p,opt_ann%n,opt_ann%g,1,0.d0,f,1)
            do i=1,opt_ann%n
                !write(81,*) iter,idp,i,f(i)*(epotall(idp)-opt_ann%epot)
                opt_ann%x(i)=opt_ann%x(i)+f(i)*(fcn_ref-fcn_ann)
            enddo
            call cpu_time(time3)
            dtime5=dtime5+time2-time1
            dtime6=dtime6+time3-time2
        enddo
        ann_arr%dpm_rmse=sqrt(ann_arr%dpm_rmse/(3.d0*opt_ann%ndp_train))
        write(88,'(i4,es18.8)') iter,ann_arr%dpm_rmse
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
end subroutine ekf_rivals
!*****************************************************************************************
subroutine analyze_epoch_init(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    if(.not. (trim(ann_arr%approach)=='eem1' .or. trim(parini%approach_ann)=='cent1' &
        .or. trim(ann_arr%approach)=='centt' .or. trim(ann_arr%approach)=='cent3' &
        .or. trim(ann_arr%approach)=='cent2')) return
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
subroutine analyze_epoch_print(parini,iter,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use yaml_output
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: i, ios
    real(8):: ttavg, ttmin, ttmax, ssavg, ssmin, ssmax
    character(50):: fn_charge, fn_chi
    character(20):: str_key
    if(.not. (trim(ann_arr%approach)=='eem1' .or. trim(parini%approach_ann)=='cent1' &
        .or. trim(ann_arr%approach)=='centt' .or. trim(ann_arr%approach)=='cent3' &
        .or. trim(ann_arr%approach)=='cent2' )) return
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
subroutine ekf_behler(parini,ann_arr,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_processors, only: iproc
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    !local variables
    real(8), allocatable:: f(:) !Kalman gain matrix
    real(8), allocatable:: p(:,:) !covariance matrix
    real(8), allocatable:: v1(:) !work array
    integer:: i, j, iter, idp, ios, ia
    real(8):: DDOT, tt, den, alambda, alambdainv, alambda0
    real(8):: fcn_ann, fcn_ref
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    allocate(f(opt_ann%n),p(opt_ann%n,opt_ann%n),v1(opt_ann%n))
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
        call convert_opt_x_ann_arr(opt_ann,ann_arr)
        if(iproc==0) then
            if( ann_arr%exists_yaml_file) then
                call write_ann_all_yaml(parini,ann_arr,iter)
            else
                call write_ann_all(parini,ann_arr,iter) 
            endif
        endif
        if(mod(iter,1)==0) then
            call ann_evaluate_all(parini,iter,ann_arr)
        endif
        if(iter==parini%nstep_opt_ann) exit
        call cpu_time(time1)
        dtime1=time1-time_s !time to get ANN energies for train and valid
        dtime2=0.d0 !time to convert ANN 1D array to ANN typ_ann
        dtime3=0.d0 !time to calculate ANN and its derivatives w.r.t. weights
        dtime4=0.d0 !time to convert derivative of ANN in typ_ann to 1D array
        dtime5=0.d0 !time to matrix-vector multiplication in Kalman filter
        dtime6=0.d0 !time of the rest of Kalman filter algorithm
        do idp=1,opt_ann%ndp_train
            alambda=min(alambda0*alambda+1.d0-alambda0,0.999999d0)
            alambdainv=1.d0/alambda
            !write(31,'(2i6,f14.10)') iter,idp,alambda
            ann_arr%event='train'
            call get_fcn_ann(parini,idp,'train',ann_arr,opt_ann,fcn_ann,fcn_ref)
            call cpu_time(time1)
            !opt_ann%g(1:opt_ann%n)=opt_ann%g(1:opt_ann%n)*(opt_ann%epot-epotall(idp))
            call DGEMV('T',opt_ann%n,opt_ann%n,1.d0,p,opt_ann%n,opt_ann%g,1,0.d0,v1,1)
            !call cal_matvec_mpi(opt_ann%n,p,opt_ann%g,v1)
            call cpu_time(time2)
            tt=DDOT(opt_ann%n,opt_ann%g,1,v1,1)
            den=alambdainv/(1.d0+tt*alambdainv)
            !write(*,'(i7,i5,f14.10,es11.2)') iter,idp,alambda,den
            !call DGEMV('T',opt_ann%n,opt_ann%n,den,p,opt_ann%n,opt_ann%g,1,0.d0,f,1)
            f(1:opt_ann%n)=v1(1:opt_ann%n)*den
            do j=1,opt_ann%n
                do i=1,opt_ann%n
                    p(i,j)=(p(i,j)-f(i)*v1(j))*alambdainv
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,idp,alambda,den
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
end subroutine ekf_behler
!*****************************************************************************************
subroutine ann_lm(parini,ann_arr,atoms_train,atoms_valid,symfunc_train,symfunc_valid,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
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
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    !use mod_opt_ann, only: typ_opt_ann
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_copy_old
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
        enddo
    elseif(iflag==2) then
        ann_arr%event='train'
        do iconf=1,atoms_train%nconf
            call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
            call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,opt_ann)
            fjac(iconf,1:opt_ann%n)=opt_ann%g(1:opt_ann%n)*(atoms%epot-atoms_train%atoms(iconf)%epot)*2.d0
        enddo
    elseif(iflag==0) then
        iter=icall0
        call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
        call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
        icall0=icall0+1
    endif
end subroutine fcn_epot
!*****************************************************************************************
subroutine set_annweights(parini,opt_ann,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, convert_ann_x
    use mod_processors, only: iproc
    use mod_utils
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_ann_arr), intent(in):: ann_arr
    !local variables
    integer:: ierr, i, ia
    real(8):: tt
    character(50):: approach
#if defined(MPI)
    include 'mpif.h'
#endif
    if (parini%restart_param) then
        do ia=1,ann_arr%nann
            call convert_ann_x(ann_arr%num(ia),opt_ann%x(ann_arr%loc(ia)),ann_arr%ann(ia))
        enddo
        return
    endif
    !approach='uniform'
    approach='pure_electrostatic'
    !approach='type_dependent'
    if(iproc==0) then
        if(trim(approach)=='uniform') then
            if(trim(parini%rng_type)=='only_for_tests') then
                call random_number_generator_simple(opt_ann%n,opt_ann%x)
            else
                call random_number(opt_ann%x)
            endif
            opt_ann%x(1:opt_ann%n)=(opt_ann%x(1:opt_ann%n)-0.5d0)*2.d0*parini%ampl_rand
        elseif(trim(approach)=='pure_electrostatic') then
            do i=1,opt_ann%n
                if(i<=opt_ann%n/2) then
                if(trim(parini%rng_type)=='only_for_tests') then
                    call random_number_generator_simple(opt_ann%x(i))
                else
                    call random_number(opt_ann%x(i))
                endif
                    tt=-2.d0*parini%ampl_rand
                    opt_ann%x(i)=(opt_ann%x(i)-0.5d0)*tt
                else
                    opt_ann%x(i)=-opt_ann%x(i-opt_ann%n/2)
                endif
            enddo
            !opt_ann%x(opt_ann%n/2)=0.d0
            !opt_ann%x(opt_ann%n)=0.d0
        elseif(trim(approach)=='type_dependent') then
            do i=1,opt_ann%n
                if(i<=opt_ann%n/4) then
                    if(trim(parini%rng_type)=='only_for_tests') then
                        call random_number_generator_simple(opt_ann%x(i))
                    else
                        call random_number(opt_ann%x(i))
                    endif
                    tt=2.d0*parini%ampl_rand
                    opt_ann%x(i)=(opt_ann%x(i)-0.5d0)*tt
                elseif(i<=opt_ann%n/2) then
                    opt_ann%x(i)=-opt_ann%x(i-opt_ann%n/4)
                elseif(i<=3*opt_ann%n/4) then
                    if(trim(parini%rng_type)=='only_for_tests') then
                        call random_number_generator_simple(opt_ann%x(i))
                    else
                        call random_number(opt_ann%x(i))
                    endif
                    tt=2.d0*parini%ampl_rand
                    opt_ann%x(i)=-(opt_ann%x(i)-0.5d0)*tt
                else
                    opt_ann%x(i)=-opt_ann%x(i-opt_ann%n/4)
                endif
            enddo
        else
            stop 'ERROR: unknown approach in set_annweights'
        endif

    endif
#if defined(MPI)
    call MPI_BCAST(opt_ann%x,opt_ann%n,MPI_DOUBLE_PRECISION,0,parini%mpi_env%mpi_comm,ierr)
#endif
end subroutine set_annweights
!*****************************************************************************************
end module mod_opt_ann
!*****************************************************************************************
