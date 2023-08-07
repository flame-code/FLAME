module mod_callback_fit_dvec
    use mod_atoms, only: typ_atoms_arr
    use mod_parini, only: typ_parini
    !use mod_ann, only: typ_ann_arr
    !use mod_symfunc, only: typ_symfunc_arr
    !use mod_opt_ann, only: typ_opt_ann
    implicit none
    type(typ_parini), pointer:: parini_t
    !type(typ_atoms_arr), pointer:: atoms_train_t
    !type(typ_atoms_arr), pointer:: atoms_valid_t
    type(typ_atoms_arr), pointer:: atoms_smplx_t
    !type(typ_ann_arr), pointer:: ann_arr_t
    !type(typ_opt_ann), pointer:: opt_ann_t
    !type(typ_symfunc_arr), pointer:: symfunc_train_t
    !type(typ_symfunc_arr), pointer:: symfunc_valid_t
    real(8), pointer:: dvec_all_t(:,:,:)
end module mod_callback_fit_dvec
!*****************************************************************************************
subroutine ekf_rivals_fitdvec(parini,ann_arr_main,opt_ann_main,atoms_train,atoms_valid,symfunc_train,symfunc_valid) !,opt_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, convert_x_ann_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_ann_io_yaml, only: write_ann_yaml
    use mod_opt_ann, only: typ_opt_ann
    use mod_ann_io_yaml, only: get_symfunc_parameters_yaml
    use mod_processors, only: iproc
    use yaml_output
    use futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr_main
    type(typ_opt_ann), intent(inout):: opt_ann_main
    type(typ_atoms_arr), intent(in):: atoms_train, atoms_valid
    type(typ_symfunc_arr), intent(in):: symfunc_train, symfunc_valid
    !local variables
    real(8), allocatable:: f(:) !Kalman gain matrix
    real(8), allocatable:: p(:,:) !covariance matrix
    real(8), allocatable:: v1(:) !work array
    real(8), allocatable:: x(:), g(:)
    integer:: n, i, j, iter, idp, ios, iann, ndp_train, ndp_valid
    integer:: iconf, iat, ig, ita, itt
    real(8):: DDOT, tt, den, fcn_ann, fcn_ref
    real(8):: r, rinv, r0, rf, alpha
    real(8):: time_s, time_e, time1, time2, time3 !, time4
    real(8):: dtime, dtime1, dtime2, dtime3, dtime4, dtime5, dtime6
    real(8):: tt1, tt2, tt3, tt4, tt5, tt6
    type(typ_ann_arr):: ann_arr
    type(typ_opt_ann):: opt_ann
    character(5):: stypat
    real(8):: rcut, chi_rmse
    character(50):: fname
    real(8):: ssavg, ssmin, ssmax
    real(8), allocatable:: features_train(:,:)
    real(8), allocatable:: features_valid(:,:)
    real(8), allocatable:: chi_ref_all_train(:,:)
    real(8), allocatable:: chi_ref_all_valid(:,:)
    real(8), allocatable:: chi_ref_train(:)
    real(8), allocatable:: chi_ref_valid(:)
    character(20):: str1, str2, str3, str4
    character(21):: fn
    character(50):: filename

    call ann_arr%set_number_of_ann(1)
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in ann_train'
    allocate(ann_arr%ann(ann_arr%nann))
    ann_arr%approach=trim(parini%approach_ann)
    do iann=1,ann_arr%nann
        stypat=parini%stypat(iann)
        fname = trim(stypat)//'.ann.input.yaml'
        call get_symfunc_parameters_yaml(parini,iproc,fname,ann_arr%ann(iann),rcut)
        ann_arr%rcut=rcut
        call dict_free(ann_arr%ann(iann)%dict)
        nullify(ann_arr%ann(iann)%dict)
    enddo

    call ann_arr%init_ann_arr()
    n=sum(ann_arr%num(1:ann_arr%nann))
    write(*,*) n,ann_arr%natmax
    allocate(g(n))
    allocate(x(n))
    allocate(f(n),p(n,n),v1(n)) !,opt_ann%epotd(ann_arr%num(1)))
    allocate(chi_ref_all_train(ann_arr%natmax,atoms_train%nconf))
    allocate(chi_ref_all_valid(ann_arr%natmax,atoms_valid%nconf))
    !open(unit=12,file='chi_ref_train.dat',status='old',iostat=ios)
    do iconf=1,atoms_train%nconf
    do iat=1,atoms_train%atoms(iconf)%nat
        !read(12,*) str1,str2,str3,chi_ref_all_train(iat,iconf)
        chi_ref_all_train(iat,iconf)=ann_arr_main%chi_ref_train(iconf)%chis(iat)
        !write(*,*) str1,str2,str3,chi_ref_all_train(iat,iconf)
    enddo
    enddo
    !close(12)
    !open(unit=12,file='chi_ref_valid.dat',status='old',iostat=ios)
    do iconf=1,atoms_valid%nconf
    do iat=1,atoms_valid%atoms(iconf)%nat
        !read(12,*) str1,str2,str3,chi_ref_all_valid(iat,iconf)
        chi_ref_all_valid(iat,iconf)=ann_arr_main%chi_ref_valid(iconf)%chis(iat)
        !write(*,*) str1,str2,str3,chi_ref_all_valid(iat,iconf)
    enddo
    enddo
    !close(12)

    do ita=1,parini%ntypat
    ndp_train=0
    do iconf=1,atoms_train%nconf
        do iat=1,atoms_train%atoms(iconf)%nat
            i=atoms_train%atoms(iconf)%itypat(iat)
            if(i==ita) ndp_train=ndp_train+1
        enddo
    enddo
    ndp_valid=0
    do iconf=1,atoms_valid%nconf
        do iat=1,atoms_valid%atoms(iconf)%nat
            i=atoms_valid%atoms(iconf)%itypat(iat)
            if(i==ita) ndp_valid=ndp_valid+1
        enddo
    enddo
    call opt_ann_main%get_opt_ann_x(opt_ann_main%n,x)
    !write(*,*) x(:)
    if(ann_arr%ann(1)%nn(0)/=symfunc_train%symfunc(1)%ng .or. ann_arr%ann(1)%nn(0)/=symfunc_valid%symfunc(1)%ng) then
        stop 'ERROR: ann_arr%ann(1)%nn(0) and symfunc_train%symfunc(1)%ng differ!'
    endif
    allocate(features_train(symfunc_train%symfunc(1)%ng,ndp_train))
    allocate(chi_ref_train(ndp_train))
    idp=0
    do iconf=1,atoms_train%nconf
        do iat=1,atoms_train%atoms(iconf)%nat
            i=atoms_train%atoms(iconf)%itypat(iat)
            if(i==ita) then
            idp=idp+1
            chi_ref_train(idp)=chi_ref_all_train(iat,iconf)
            do ig=1,symfunc_train%symfunc(iconf)%ng
                features_train(ig,idp)=symfunc_train%symfunc(iconf)%y(ig,iat)
                !write(21,'(es14.5,3i5)') symfunc_train%symfunc(iconf)%y(ig,iat),iconf,iat,ig
            enddo
            !write(22,'(es14.5,2i5)') chi_ref_all_train(iat,iconf),iconf,iat
            endif
        enddo
    enddo
    allocate(features_valid(symfunc_valid%symfunc(1)%ng,ndp_valid))
    allocate(chi_ref_valid(ndp_valid))
    idp=0
    do iconf=1,atoms_valid%nconf
        do iat=1,atoms_valid%atoms(iconf)%nat
            i=atoms_valid%atoms(iconf)%itypat(iat)
            if(i==ita) then
            idp=idp+1
            chi_ref_valid(idp)=chi_ref_all_valid(iat,iconf)
            do ig=1,symfunc_valid%symfunc(iconf)%ng
                features_valid(ig,idp)=symfunc_valid%symfunc(iconf)%y(ig,iat)
                !write(21,'(es14.5,3i5)') symfunc_train%symfunc(iconf)%y(ig,iat),iconf,iat,ig
            enddo
            !write(22,'(es14.5,2i5)') chi_ref_all_train(iat,iconf),iconf,iat
            endif
        enddo
    enddo
    !do idp=1,ndp_train
    !    do ig=1,ann_arr%ann(1)%nn(0)
    !        write(ita*100+1,'(es14.5)') features_train(ig,idp)
    !    enddo
    !    write(ita*100+3,'(es14.5)') chi_ref_train(idp)
    !enddo
    !do idp=1,ndp_valid
    !    do ig=1,ann_arr%ann(1)%nn(0)
    !        write(ita*100+2,'(es14.5)') features_valid(ig,idp)
    !    enddo
    !    write(ita*100+4,'(es14.5)') chi_ref_valid(idp)
    !enddo
    !stop 'PPPPPPPPPPPPPPPP'
    tt1=0.d0
    do idp=1,ndp_train
        tt1=tt1+chi_ref_train(idp)
    enddo
    tt1=tt1/real(ndp_train,kind=8)
    tt1=real(int(tt1*1.d2),kind=8)*1.d-2
    ann_arr%ann(1)%chi0=tt1
    ann_arr_main%ann(ita)%chi0=tt1



    p(1:n,1:n)=0.d0
    do i=1,n
        p(i,i)=1.d-2
    enddo
    if(trim(parini%approach_ann)=='cent2') then
        r0=10.d0
        alpha=100.d-2
        rf=1.d-6
    else
        write(*,*) 'ERROR: ekf_rivals_fitdvec to be used only with CENT2'
        stop
    endif
    do iter=0,parini%nstep_opt_ann
        !call yaml_sequence(advance='no',unit=ann_arr%iunit)
        !call convert_opt_x_ann_arr(opt_ann,ann_arr)
        call convert_x_ann_arr(n,x,ann_arr)

        ann_arr_main%ann(ita)%a(:,:,:)=ann_arr%ann(1)%a(:,:,:)
        ann_arr_main%ann(ita)%b(:,:)=ann_arr%ann(1)%b(:,:)
        write(fn,'(a16,i5.5)') '.ann.param.yaml.',iter
        filename=trim(parini%stypat(ita))//trim(fn)
        !write(*,'(a)') trim(filename)
        call yaml_comment(trim(filename))
        call write_ann_yaml(parini,filename,ann_arr_main%ann(ita),ann_arr_main%rcut)

        if(mod(iter,1)==0) then
            !call analyze_epoch_init(parini,ann_arr)
            call ann_evaluate_fitchi(parini,ita,iter,ann_arr,n,g, &
                ndp_train,features_train,chi_ref_train, &
                ndp_valid,features_valid,chi_ref_valid)
            !call analyze_epoch_print(parini,iter,ann_arr)
        endif
        if(iter==parini%nstep_opt_ann) exit
        r=(r0-rf)*exp(-alpha*(iter+parini%restart_iter))+rf
        rinv=1.d0/r
        write(31,'(i6,es14.5)') iter,r
        chi_rmse=0
        do idp=1,ndp_train
            ann_arr%event='train'
            !The following is allocated with ann_arr%num(1), this means number of
            !nodes in the input layer is the same for all atom types.
            !Therefore, it must be fixed later.
            !g_per_atom=f_malloc([1.to.ann_arr%num(1),1.to.atoms%nat],id='g_per_atom') !HERE
            !call convert_opt_x_ann_arr(opt_ann,ann_arr)
            call convert_x_ann_arr(n,x,ann_arr)
            call cal_chi_from_features(parini,ann_arr,features_train(1,idp),'train',n,g,fcn_ann)
            fcn_ref=chi_ref_train(idp)
            chi_rmse=chi_rmse+(fcn_ann-fcn_ref)**2
            !if(parini%iverbose>=2) then
            !    write(1370,'(2i4,2es18.8)') iter,idp,ann_arr%dpm_err, ann_arr%dpm_rmse
            !end if
            call DGEMV('T',n,n,1.d0,p,n,g,1,0.d0,v1,1)
            !call cal_matvec_mpi(n,p,opt_ann%g,v1)
            tt=DDOT(n,g,1,v1,1)
            den=1.d0/(tt+r)
            do j=1,n
                do i=1,n
                    p(i,j)=p(i,j)-v1(i)*v1(j)*den
                    !if(i==j) p(i,j)=p(i,j)+1.d-1
                    !write(21,'(2i5,es20.10)') i,j,p(i,j)
                enddo
            enddo
                    !write(21,'(a)') '----------------------------------------'
            !write(*,'(a,i7,i6,2es15.5)') 'forgetting ',iter,idp,r,den
            call DGEMV('N',n,n,rinv,p,n,g,1,0.d0,f,1)
            do i=1,n
                !write(81,*) iter,idp,i,f(i)*(epotall(idp)-opt_ann%epot)
                x(i)=x(i)+f(i)*(fcn_ref-fcn_ann)
            enddo
        enddo
        chi_rmse=sqrt(chi_rmse/(3.d0*ndp_train))
        write(45,'(i4,es18.8)') iter,chi_rmse
    enddo
    deallocate(features_train)
    deallocate(chi_ref_train)
    deallocate(features_valid)
    deallocate(chi_ref_valid)
    enddo !end of loop over ita
    close(41)
    close(42)
    deallocate(f,p,v1) !,opt_ann%epotd)
end subroutine ekf_rivals_fitdvec
!*****************************************************************************************
subroutine cal_dvec_from_features(parini,ann_arr,features,str_dataset,n,g,fcn_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: features(ann_arr%ann(1)%nn(0))
    character(*), intent(in):: str_dataset
    integer, intent(in):: n
    real(8), intent(out):: fcn_ann, g(n)
    !local variables
    integer:: ng
    real(8):: out_ann, tt1
    ng=ann_arr%ann(1)%nn(0)
    ann_arr%ann(1)%y(1:ng,0)=features(1:ng)
    if(trim(ann_arr%event)=='potential' .or. trim(ann_arr%event)=='evalu') then
        call cal_architecture(ann_arr%ann(1),out_ann)
        tt1=tanh(ann_arr%ann(1)%prefactor_chi*out_ann)
        fcn_ann=ann_arr%ann(1)%ampl_chi*tt1+ann_arr%ann(1)%chi0
    elseif(trim(ann_arr%event)=='train') then
        call cal_architecture_der(ann_arr%ann(1),out_ann)
        call convert_ann_epotd(ann_arr%ann(1),ann_arr%num(1),g)
        tt1=tanh(ann_arr%ann(1)%prefactor_chi*out_ann)
        fcn_ann=ann_arr%ann(1)%ampl_chi*tt1+ann_arr%ann(1)%chi0
        g(1:n)=g(1:n)*ann_arr%ann(1)%ampl_chi*ann_arr%ann(1)%prefactor_chi*(1.d0-tt1**2)
    else
        stop 'ERROR: undefined content for ann_arr%event'
    endif
end subroutine cal_dvec_from_features
!*****************************************************************************************
subroutine ann_evaluate_fitdvec(parini,ita,iter,ann_arr,n,g,ndp_train,ft_train,chi_ref_train,ndp_valid,ft_valid,chi_ref_valid)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: ita, iter
    type(typ_ann_arr), intent(inout):: ann_arr
    integer, intent(in):: n
    real(8), intent(inout):: g(n)
    integer, intent(in):: ndp_train, ndp_valid
    real(8), intent(in):: ft_train(ann_arr%ann(1)%nn(0),ndp_train), chi_ref_train(ndp_train)
    real(8), intent(in):: ft_valid(ann_arr%ann(1)%nn(0),ndp_valid), chi_ref_valid(ndp_valid)
    !local variables
    real(8):: rmse_train, rmse_valid, chi, tt1, tt2, tt3, ss1, ss2
    integer:: idp, iconf, jdp, jconf
    !allocate(chiall(ndp_train))
    ann_arr%event='evalu'
    !iconf=0
    rmse_train=0.d0
    do idp=1,ndp_train
        call cal_chi_from_features(parini,ann_arr,ft_train(1,idp),'train',n,g,chi)
        !chiall(idp)=chi
        !if(mod(idp-1,40)==0) iconf=iconf+1
        !write(1000+iter,*) iconf,abs(chi-chi_ref_train(idp))
        rmse_train=rmse_train+(chi-chi_ref_train(idp))**2
    enddo
    !close(1000+iter)
    rmse_train=sqrt(rmse_train/real(ndp_train,kind=8))
    !iconf=0
    rmse_valid=0.d0
    do idp=1,ndp_valid
        call cal_chi_from_features(parini,ann_arr,ft_valid(1,idp),'valid',n,g,chi)
        !chiall(idp)=chi
        !if(mod(idp-1,40)==0) iconf=iconf+1
        !write(2000+iter,*) iconf,abs(chi-chi_ref_valid(idp))
        rmse_valid=rmse_valid+(chi-chi_ref_valid(idp))**2
    enddo
    !close(2000+iter)
    rmse_valid=sqrt(rmse_valid/real(ndp_valid,kind=8))
    write(80+ita,'(a,i5,2es14.5)') 'RMSE ',iter,rmse_train,rmse_valid
    !-------------------------------------------------------
    !if(iter==0) then
    !tt1=0.d0
    !tt2=0.d0
    !tt3=0.d0
    !iconf=0
    !do idp=1,ndp_train
    !    if(mod(idp-1,40)==0) iconf=iconf+1
    !    jconf=0
    !    do jdp=1,idp-1
    !        if(mod(jdp-1,40)==0) jconf=jconf+1
    !        ss1=sqrt(sum((ft_train(1:ann_arr%ann(1)%nn(0),idp)-ft_train(1:ann_arr%ann(1)%nn(0),jdp))**2))
    !        ss2=abs(chi_ref_train(idp)-chi_ref_train(jdp))
    !        tt1=tt1+ss1
    !        tt2=tt2+ss2
    !        tt3=tt3+1.d0
    !        if(ss1<0.05 .and. ss2>0.005) write(65,'(i4,2i6,2f15.5)') iter,iconf,jconf,ss1,ss2
    !    enddo
    !enddo
    !tt1=tt1/tt3
    !tt2=tt2/tt3
    !write(*,'(a,2f15.5)') 'TT1,TT2 ',tt1,tt2
    !!deallocate(chiall)
    !endif
end subroutine ann_evaluate_fitdvec
!*****************************************************************************************
subroutine fit_dvec(parini)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, typ_atoms_arr, atom_deallocate_old !, atom_copy
    use mod_atoms, only: update_ratp
    !use mod_potential, only: fcalls, perfstatus, potcode
    !use mod_acf, only: acf_read_new, acf_write
    use mod_callback_fit_dvec, only: atoms_smplx_t, parini_t
    use mod_callback_fit_dvec, only: dvec_all_t
    use mod_yaml_conf, only: read_yaml_conf
    use mod_ann_io_yaml, only: read_input_ann_yaml, read_data_yaml
    use mod_const, only: ha2ev
    !use mod_processors, only: iproc
    !use mod_const, only: ev2ha, ang2bohr
    !use yaml_output
    implicit none
    type(typ_parini), target, intent(in):: parini
    !local variables
    type(typ_atoms_arr), target:: atoms_arr
    !type(typ_atoms):: atoms
    integer:: iconf, iat, iter, i, id, ii, jj, ios, itypat, info, jat, kat, istep, ityp
    integer:: ndim, nfiles_max, nfiles
    real(8):: tt, tt1, tt2, tt3, tt4, tt5, tt6, ss, step, ftol
    real(8):: dx, dy, dz, ccc(3), vol, eval(3), qt(2,2), gwall(3,2,2), pref(3,2)
    character(len=2):: str
    real(8), allocatable, target:: dvec_all(:,:,:)
    real(8), allocatable:: vertices(:,:), fval(:)
    real(8), allocatable:: refvec(:,:)
    real(8), allocatable:: qat(:,:)
    character(256), allocatable:: fn_list(:)
    external cal_rmse_dvec
    !call read_yaml_conf(parini,'posinp.yaml',10000,atoms_arr)
    call read_data_yaml(parini,'list_posinp.yaml',atoms_arr)
    do iconf=1,atoms_arr%nconf
    call update_ratp(atoms_arr%atoms(iconf))
    !    call atom_copy(atoms_arr%atoms(iconf),atoms,'atoms_arr%atoms->atoms_arr')
    do iat=1,atoms_arr%atoms(iconf)%nat
        do itypat=1,parini%ntypat
            if(trim(atoms_arr%atoms(iconf)%sat(iat))==trim(parini%stypat(itypat))) then
                atoms_arr%atoms(iconf)%itypat(iat)=parini%ltypat(itypat)
                exit
            endif
        enddo
    enddo
    enddo
    !write(*,*) 'NCONF ',atoms_arr%nconf
    !do iconf=1,atoms_arr%nconf
    !do iat=1,atoms_arr%atoms(iconf)%nat
    !write(51,'(2i6,a5,3f15.6)') iconf,atoms_arr%atoms(iconf)%itypat(iat),atoms_arr%atoms(iconf)%sat(iat), &
    !    atoms_arr%atoms(iconf)%rat(1,iat), &
    !    atoms_arr%atoms(iconf)%rat(2,iat), &
    !    atoms_arr%atoms(iconf)%rat(3,iat)
    !enddo
    !enddo
    allocate(dvec_all(3,atoms_arr%atoms(1)%nat,atoms_arr%nconf))
    nfiles_max=10**5
    allocate(fn_list(nfiles_max))
    call read_list_files_yaml('list_dvec.yaml',nfiles_max,fn_list,nfiles)
    if(nfiles/=atoms_arr%nconf) then
        write(*,*) 'ERROR: nfiles/=atoms_arr%nconf in fit_dvec'
        stop
    endif
    !-----------------------------------------------------------------
    do iconf=1,atoms_arr%nconf
        open(unit=21,file=trim(fn_list(iconf)),status='old',iostat=ios)
        do iat=1,atoms_arr%atoms(iconf)%nat
            read(21,*) ii,str,tt1,tt2,tt3,tt4,tt5,tt6
            dvec_all(1,iat,iconf)=tt4
            dvec_all(2,iat,iconf)=tt5
            dvec_all(3,iat,iconf)=tt6
            tt1=sqrt(sum(dvec_all(:,iat,iconf)**2))
            write(*,'(a,i6,f10.2)') 'DDD ',iat,tt1
            !dvec(:,iat)=dvec(:,iat)/tt1
        enddo
        close(21)
    enddo
    deallocate(fn_list)
    !-----------------------------------------------------------------
    atoms_smplx_t=>atoms_arr
    parini_t=>parini
    dvec_all_t=>dvec_all
    ndim=18
    ftol=1.d-10 !parini%ftol_ann
    step=0.d0
    allocate(vertices(ndim,ndim+1),fval(ndim+1))
    !vertices(1,1)=3.01d0/ha2ev*2.0d0
    !vertices(2,1)=3.01d0/ha2ev*0.5d0
    !vertices(3,1)=6.22d0/ha2ev*2.0d0
    !vertices(4,1)=6.22d0/ha2ev*0.5d0

    vertices( 1,1)= 3.d0
    vertices( 2,1)= 3.d0
    vertices( 3,1)= 5.d0
    vertices( 4,1)= 3.d0
    vertices( 5,1)= 3.d0
    vertices( 6,1)= 5.d0

    vertices( 7,1)= 3.d0
    vertices( 8,1)= 3.d0
    vertices( 9,1)= 5.d0
    vertices(10,1)= 3.d0
    vertices(11,1)= 3.d0
    vertices(12,1)= 5.d0

    vertices(13,1)= 3.d0
    vertices(14,1)= 3.d0
    vertices(15,1)= 5.d0
    vertices(16,1)= 3.d0
    vertices(17,1)= 3.d0
    vertices(18,1)= 5.d0
    !vertices( 7,1)= 1.d-1
    !vertices( 8,1)= 1.d-1
    !vertices( 9,1)= 1.d-1
    !vertices(10,1)= 1.d-1

    !vertices(1,1)= 4.d0
    !vertices(2,1)= 1.d0
    !vertices(3,1)=-2.d0
    !vertices(4,1)= 1.d0
    !vertices(5,1)= 5.d0
    !vertices(6,1)= 4.d0
    !vertices(7,1)= 5.d0
    !vertices(8,1)= 4.d0

    do id=2,ndim+1
    do i=1,ndim
        call random_number(tt)
        tt=2.d0*(tt-0.5d0)*0.03d0
        vertices(i,id)=vertices(i,id-1)+tt
    enddo
    enddo

    !do istep=1,1
    !call get_ref_vector(parini,atoms,refvec,dvec,qat,pref,qt)
    call simplex(vertices,fval,step,ndim,ftol,cal_rmse_dvec,iter)
    do iconf=1,atoms_arr%nconf
    associate(atoms=>atoms_arr%atoms(iconf))
    allocate(refvec(3,atoms%nat))
    allocate(qat(3,atoms%nat))
!    if(parini%mpi_env%iproc==0) then
!    !write(*,'(a,i6,3es14.5)') 'COST ',istep,sqrt(sum((refvec-dvec)**2)),maxval(abs(refvec-dvec)),sqrt(sum((dvec)**2))
!    write(*,'(a,i6,3es14.5)') 'COST ',1,sqrt(sum((refvec-dvec)**2)),maxval(abs(refvec-dvec)),sqrt(sum((dvec)**2))
!    endif
    !qtarget_der=0.d0
    !do ityp=1,2
    !    qtarget(1,ityp)=qtarget(1,ityp)-1.d-2*qtarget_der(1,ityp)
    !    qtarget(2,ityp)=qtarget(2,ityp)-1.d-2*qtarget_der(2,ityp)
    !    qtarget(3,ityp)=qtarget(3,ityp)-1.d-2*qtarget_der(3,ityp)
    !enddo
    !write(22,'(i6,6f6.2)') istep,qtarget(1,1),qtarget(2,1),qtarget(3,1),qtarget(1,2),qtarget(2,2),qtarget(3,2)
    !enddo !end of loop over istep

    pref(1,1)=vertices(1,1)
    pref(1,2)=vertices(2,1)
    gwall(1,1,1)=vertices(3,1)
    gwall(1,2,1)=vertices(4,1)
    gwall(1,1,2)=vertices(5,1)
    gwall(1,2,2)=vertices(6,1)
    pref(2,1)=vertices(7,1)
    pref(2,2)=vertices(8,1)
    gwall(2,1,1)=vertices(9,1)
    gwall(2,2,1)=vertices(10,1)
    gwall(2,1,2)=vertices(11,1)
    gwall(2,2,2)=vertices(12,1)
    pref(3,1)=vertices(13,1)
    pref(3,2)=vertices(14,1)
    gwall(3,1,1)=vertices(15,1)
    gwall(3,2,1)=vertices(16,1)
    gwall(3,1,2)=vertices(17,1)
    gwall(3,2,2)=vertices(18,1)

    !qt(1,1)=vertices(7,1)
    !qt(2,1)=vertices(8,1)
    !qt(1,2)=vertices(9,1)
    !qt(2,2)=vertices(10,1)
    call get_ref_vector(parini,atoms,refvec,dvec_all(1,1,iconf),qat,pref,qt,gwall)
    if(parini%mpi_env%iproc==0) then
    do iat=1,atoms%nat
        write(*,'(a,i5,4f10.1)') 'QAT ',iat,qat(1,iat),qat(2,iat),qat(3,iat),pref(1,atoms%itypat(iat))
        write(77,'(i5,2x,a2,2x,6es14.5)') iat,trim(atoms%sat(iat)), &
            refvec(1,iat),refvec(2,iat),refvec(3,iat),dvec_all(1,iat,iconf),dvec_all(2,iat,iconf),dvec_all(3,iat,iconf)
    enddo
    endif
    deallocate(refvec)
    deallocate(qat)
    end associate
    enddo !end of loop over iconf
    nullify(atoms_smplx_t)
    nullify(parini_t)
    nullify(dvec_all_t)
    do iconf=1,atoms_arr%nconf
        call atom_deallocate_old(atoms_arr%atoms(iconf))
    enddo
    deallocate(atoms_arr%atoms)
    deallocate(dvec_all)
end subroutine fit_dvec
!*****************************************************************************************
subroutine get_ref_vector_new(parini,atoms,refvec,refvec_qder,qat,pref,icycle,qpm,qt,gwall)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms, update_ratp
    use mod_ann, only: typ_ann_arr
    use mod_const, only: bohr2ang, ha2ev
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    !type(typ_ann_arr), intent(in):: ann_arr
    real(8), intent(out):: refvec(3,atoms%nat)
    real(8), intent(out):: refvec_qder(3,atoms%nat,3,atoms%nat)
    real(8), intent(in):: qat(3,atoms%nat)
    real(8), intent(in):: pref(3,2)
    integer, intent(in):: icycle
    real(8), intent(inout):: qpm(3,3,atoms%nat)
    real(8), intent(in):: qt(2,2), gwall(3,2,2)
    !local variables
    integer:: iat, jat, kat, nat, info, ii
    real(8):: dx, dy, dz, r, rsq, pi, beta_iat, beta_jat, gama
    real(8):: ri, rj, fc, fcj, tt1, tt2, tt3, gg, shift, ee1, ee2
    real(8):: gw(2), chi(2), hn(2), rcut, gg_kat, gg_iat, tta(3)
    real(8):: txx, tyx, tzx, txy, tyy, tzy, txz, tyz, tzz
    real(8):: cutoff_function
    real(8), allocatable:: a(:,:)
    real(8), allocatable:: rhs(:,:)
    real(8), allocatable:: rat(:,:)
    real(8), allocatable:: qat_t(:,:)
    integer, allocatable:: itypat(:)
    integer, allocatable:: nn(:)
    integer, allocatable:: ind(:)
    pi=4.d0*atan(1.d0)
    gw(1)=1.28d0/bohr2ang*1.d0
    gw(2)=1.05d0/bohr2ang*1.d0
    chi(1)=3.01d0/ha2ev*10.d0
    chi(2)=6.22d0/ha2ev*10.d0
    hn(1)=0.4d0 !2.39d0/ha2ev
    hn(2)=0.2d0 !4.14d0/ha2ev
    rcut=15.d0
    !qt(1)= 1.0d-1
    !qt(2)=-4.d-1
    !if(icycle==1) rcut=6.d0
    !if(icycle==2) rcut=9.d0
    !if(icycle==3) rcut=11.d0
    allocate(nn(atoms%nat),source=1)
    do iat=1,atoms%nat
        do jat=iat+1,atoms%nat
            dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
            dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
            dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
            rsq=dx*dx+dy*dy+dz*dz
            if(rsq<rcut**2) then
                nn(iat)=nn(iat)+1
                nn(jat)=nn(jat)+1
            endif
        enddo
    enddo
    !write(*,*) 'NN= ',nn
    !-----------------------------------------------------------------
    refvec_qder=0.d0
    do kat=1,atoms%nat
    nat=nn(kat)
    allocate(a(nat+1,nat+1))
    allocate(rhs(3,1:nat+1))
    allocate(rat(3,nat))
    allocate(qat_t(3,nat))
    allocate(itypat(nat))
    allocate(ind(nat))
    rat(1:3,1)=atoms%ratp(1:3,kat)
    itypat(1)=atoms%itypat(kat)
    ind(1)=kat
    jat=1
    do iat=1,atoms%nat
        if(iat==kat) cycle
        dx=atoms%ratp(1,kat)-atoms%ratp(1,iat)
        dy=atoms%ratp(2,kat)-atoms%ratp(2,iat)
        dz=atoms%ratp(3,kat)-atoms%ratp(3,iat)
        rsq=dx*dx+dy*dy+dz*dz
        if(rsq<rcut**2) then
            r=sqrt(rsq)
            jat=jat+1
            rat(1:3,jat)=atoms%ratp(1:3,iat)
            itypat(jat)=atoms%itypat(iat)
            ind(jat)=iat
            qat_t(1,jat)=qat(1,iat)
            qat_t(2,jat)=qat(2,iat)
            qat_t(3,jat)=qat(3,iat)
            !write(*,'(a,4f10.3,i4)') 'TTT ',qat_t(jat),chi_ends(1,itypat(jat)),chi_ends(2,itypat(jat)),r,itypat(jat)
            !write(*,'(a,3i5,f10.3,i8,3f10.3)') 'KKK ',kat,iat,jat,sqrt(r),itypat(jat),rat(1:3,jat)
        endif
    enddo
    qat_t(1,1)=qat(1,kat)
    qat_t(2,1)=qat(2,kat)
    qat_t(3,1)=qat(3,kat)
    if(parini%mpi_env%iproc==0) then
    if(jat/=nat) write(*,*) 'ERROR: ',kat,nat,jat
    endif
    do iat=1,nat
        rhs(1,iat)=qat_t(1,iat)
        rhs(2,iat)=qat_t(2,iat)
        rhs(3,iat)=qat_t(3,iat)
        !write(33,'(i5,f10.3)') iat,rhs(iat)
    enddo
    refvec(1:3,kat)=0.d0
    qpm(1:3,1:3,kat)=0.d0
    gg_kat=gw(itypat(1))
    do iat=2,nat
        dx=rat(1,iat)-rat(1,1)
        dy=rat(2,iat)-rat(2,1)
        dz=rat(3,iat)-rat(3,1)
        r=sqrt(dx*dx+dy*dy+dz*dz)
        gg=gw(itypat(iat))
        gg_iat=gw(itypat(iat))
        fc=cutoff_function(r,rcut)
        !if(trim(parini%stypat(itypat(iat)))=='Li') then
        !tt=-1.d0*fc*(2.d0/(gg*exp(r**2/gg**2)*sqrt(pi)*r)-erf(r/gg)/r**2)
        !else
        !tt=-2.d0*fc*(2.d0/(gg*exp(r**2/gg**2)*sqrt(pi)*r)-erf(r/gg)/r**2)
        !endif
        !tt=-fc*(2.d0/(gg*exp(r**2/gg**2)*sqrt(pi)*r)-erf(r/gg)/r**2)
        !tt=-fc*(-exp(-r/gg)/r**3)
        !tt1=-fc*(-1.d0/r**3)
        !tt2=-fc*(-1.d0/r**3.5d0)
        !tt3=-fc*(-1.d0/r**4)
        !if(icycle==1) tt1=-fc*(-1.d0/r**3)*1.0d0
        !if(icycle==1) tt1=-fc*(-1.d0)*1.0d0

        !if(icycle==1) tt1=-fc*(-1.d0/r**5)*1.0d0

        !if(icycle==1) tt1=-fc*(-exp(-(r/(0.80d0*(gg_kat+gg_iat)))**2))*1.0d0
        !if(icycle==2) tt1=-fc*(-exp(-(r/(0.85d0*(gg_kat+gg_iat)))**2))*1.0d0
        !if(icycle==3) tt1=-fc*(-exp(-(r/(0.90d0*(gg_kat+gg_iat)))**2))*1.0d0
        tt1=-fc*(-exp(-(r/gwall(1,itypat(iat),itypat(1)))**2))*1.0d0
        !tt1=-fc*(-exp(-(r/gwall(itypat(iat),itypat(1)))**1))*1.0d0
        tt2=-fc*(-exp(-(r/gwall(2,itypat(iat),itypat(1)))**2))*1.0d0
        tt3=-fc*(-exp(-(r/gwall(3,itypat(iat),itypat(1)))**2))*1.0d0

        !tt1=-fc*(-exp(-(r/(0.7d0*(gg_kat+gg_iat)))**2))*1.0d0

        !if(icycle==2) tt1=-fc*(-exp(-(r/7.d0)**2))*1.0d0
        !if(icycle==3) tt1=-fc*(-exp(-(r/9.d0)**2))*1.0d0
        !write(33,'(2i6,2f10.3)') kat,iat,fc,exp(-(r/(2.d0*gg))**2)
        !tt1=fc*1.d0 !/r**2
        shift=r !0.d0
        ee1=1.d0
        ee2=1.d-2
        txx=(3.d0*dx*dx-shift**2)*tt1*ee1
        tyx=(3.d0*dy*dx         )*tt1*ee1
        tzx=(3.d0*dz*dx         )*tt1*ee1
        txy=(3.d0*dx*dy         )*tt1*ee1
        tyy=(3.d0*dy*dy-shift**2)*tt1*ee1
        tzy=(3.d0*dz*dy         )*tt1*ee1
        txz=(3.d0*dx*dz         )*tt1*ee1
        tyz=(3.d0*dy*dz         )*tt1*ee1
        tzz=(3.d0*dz*dz-shift**2)*tt1*ee1
        !write(21,'(2(2i6,f10.3))') kat,itypat(1),qt(itypat(1)),iat,itypat(iat),qt(itypat(iat))
        qpm(1,1,kat)=qpm(1,1,kat)+qt(itypat(iat),itypat(1))*txx
        qpm(2,1,kat)=qpm(2,1,kat)+qt(itypat(iat),itypat(1))*tyx
        qpm(3,1,kat)=qpm(3,1,kat)+qt(itypat(iat),itypat(1))*tzx
        qpm(1,2,kat)=qpm(1,2,kat)+qt(itypat(iat),itypat(1))*txy
        qpm(2,2,kat)=qpm(2,2,kat)+qt(itypat(iat),itypat(1))*tyy
        qpm(3,2,kat)=qpm(3,2,kat)+qt(itypat(iat),itypat(1))*tzy
        qpm(1,3,kat)=qpm(1,3,kat)+qt(itypat(iat),itypat(1))*txz
        qpm(2,3,kat)=qpm(2,3,kat)+qt(itypat(iat),itypat(1))*tyz
        qpm(3,3,kat)=qpm(3,3,kat)+qt(itypat(iat),itypat(1))*tzz
        !refvec(1,kat)=refvec(1,kat)+pref(kat)*rhs(1,iat)*(txx+txy+txz)*tt1*ee2
        !refvec(2,kat)=refvec(2,kat)+pref(kat)*rhs(1,iat)*(tyx+tyy+tyz)*tt1*ee2
        !refvec(3,kat)=refvec(3,kat)+pref(kat)*rhs(1,iat)*(tzx+tzy+tzz)*tt1*ee2
        !refvec_qder(1,kat,1,ind(iat))=pref(kat)*(txx+txy+txz)*tt1*ee2
        !refvec_qder(2,kat,1,ind(iat))=pref(kat)*(tyx+tyy+tyz)*tt1*ee2
        !refvec_qder(3,kat,1,ind(iat))=pref(kat)*(tzx+tzy+tzz)*tt1*ee2
        !tt2=0.d0
        !tt3=0.d0
        refvec(1,kat)=refvec(1,kat)+pref(1,itypat(1))*rhs(1,iat)*tt1*dx
        refvec(2,kat)=refvec(2,kat)+pref(1,itypat(1))*rhs(1,iat)*tt1*dy
        refvec(3,kat)=refvec(3,kat)+pref(1,itypat(1))*rhs(1,iat)*tt1*dz
        refvec(1,kat)=refvec(1,kat)+pref(2,itypat(1))*rhs(2,iat)*tt2*dx
        refvec(2,kat)=refvec(2,kat)+pref(2,itypat(1))*rhs(2,iat)*tt2*dy
        refvec(3,kat)=refvec(3,kat)+pref(2,itypat(1))*rhs(2,iat)*tt2*dz
        refvec(1,kat)=refvec(1,kat)+pref(3,itypat(1))*rhs(3,iat)*tt3*dx
        refvec(2,kat)=refvec(2,kat)+pref(3,itypat(1))*rhs(3,iat)*tt3*dy
        refvec(3,kat)=refvec(3,kat)+pref(3,itypat(1))*rhs(3,iat)*tt3*dz
        refvec_qder(1,kat,1,ind(iat))=pref(1,itypat(1))         *tt1*dx
        refvec_qder(2,kat,1,ind(iat))=pref(1,itypat(1))         *tt1*dy
        refvec_qder(3,kat,1,ind(iat))=pref(1,itypat(1))         *tt1*dz
        refvec_qder(1,kat,2,ind(iat))=pref(2,itypat(1))         *tt2*dx
        refvec_qder(2,kat,2,ind(iat))=pref(2,itypat(1))         *tt2*dy
        refvec_qder(3,kat,2,ind(iat))=pref(2,itypat(1))         *tt2*dz
        refvec_qder(1,kat,3,ind(iat))=pref(3,itypat(1))         *tt3*dx
        refvec_qder(2,kat,3,ind(iat))=pref(3,itypat(1))         *tt3*dy
        refvec_qder(3,kat,3,ind(iat))=pref(3,itypat(1))         *tt3*dz
    enddo
    if(icycle>1) then
    do ii=1,icycle-1
    tta=refvec(1:3,kat)
    refvec(1,kat)=dot_product(qpm(1:3,1,kat),tta)
    refvec(2,kat)=dot_product(qpm(1:3,2,kat),tta)
    refvec(3,kat)=dot_product(qpm(1:3,3,kat),tta)
    do iat=2,nat
    tta=refvec_qder(1:3,kat,1,ind(iat))
    refvec_qder(1,kat,1,ind(iat))=dot_product(qpm(1:3,1,kat),tta)
    refvec_qder(2,kat,1,ind(iat))=dot_product(qpm(1:3,2,kat),tta)
    refvec_qder(3,kat,1,ind(iat))=dot_product(qpm(1:3,3,kat),tta)
    enddo
    enddo
    endif
    !refvec(1,kat)=refvec(1,kat)*rhs(1)
    !refvec(2,kat)=refvec(2,kat)*rhs(1)
    !refvec(3,kat)=refvec(3,kat)*rhs(1)
    if(parini%mpi_env%iproc==0) then
        !write(*,'(a,4f10.3)') 'ENDs',chi_ends(1,1),chi_ends(2,1),chi_ends(1,2),chi_ends(2,2)
        !write(*,'(a,4i4)') 'ENDs',itypat(1),itypat(2),itypat(nat-1),itypat(nat)
    !write(*,'(a,i5,3es14.5)') 'EW= ',kat,refvec(1,kat),refvec(2,kat),refvec(3,kat)
    !do iat=1,nat
    !    write(*,'(a,2i5,f10.3,2x,a2)') 'QW= ',kat,iat,rhs(iat),trim(parini%stypat(itypat(iat)))
    !enddo
    endif
    deallocate(a,rhs,rat,itypat,ind,qat_t)
    enddo
end subroutine get_ref_vector_new
!*****************************************************************************************
subroutine get_ref_vector(parini,atoms,refvec,dvec,qat,pref,qt,gwall)
!,vertex)
    use mod_parini, only: typ_parini
    use mod_atoms, only: typ_atoms !, update_ratp
    !use mod_ann, only: typ_ann_arr
    !use mod_const, only: bohr2ang, ha2ev
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_atoms), intent(in):: atoms
    !type(typ_ann_arr), intent(in):: ann_arr
    real(8), intent(out):: refvec(3,atoms%nat)
    real(8), intent(in):: dvec(3,atoms%nat)
    real(8), intent(inout):: qat(3,atoms%nat)
    real(8), intent(in):: pref(3,2)
    real(8), intent(in):: qt(2,2), gwall(3,2,2)
    !local variables
    integer:: icycle, iat, jat, kat, ii, jj, info
    integer:: nwork
    real(8):: tt, ss, qtarget(3,2), regcoeff(3)
    real(8), allocatable:: qpm(:,:,:)
    real(8), allocatable:: refvec_t(:,:)
    real(8), allocatable:: refvec_qder(:,:,:,:)
    real(8), allocatable:: amat(:,:)
    real(8), allocatable:: amat_t(:,:)
    real(8), allocatable:: rhs(:)
    real(8), allocatable:: work(:)
    integer, allocatable:: ipiv(:)
    real(8), allocatable:: dvec_t(:,:)
    allocate(dvec_t(3,atoms%nat))
    allocate(qpm(3,3,atoms%nat))
    allocate(refvec_t(3,atoms%nat))
    allocate(refvec_qder(3,atoms%nat,3,atoms%nat))
    allocate(amat(3*atoms%nat,3*atoms%nat),source=0.d0)
    allocate(amat_t(3*atoms%nat,3*atoms%nat),source=0.d0)
    allocate(rhs(3*atoms%nat))
    allocate(ipiv(3*atoms%nat))
    nwork=10*(3*atoms%nat)**2
    allocate(work(nwork))
    regcoeff(1)=1.d-7
    regcoeff(2)=1.d-7
    regcoeff(3)=1.d-7
    qtarget(1,1)= 0.0d0 ! 0.1d0
    qtarget(2,1)= 0.0d0 ! 0.2d0
    qtarget(3,1)= 0.0d0 ! 0.3d0
    qtarget(1,2)=-0.0d0 !-0.1d0
    qtarget(2,2)=-0.0d0 !-0.2d0
    qtarget(3,2)=-0.0d0 !-0.3d0
    dvec_t=0.d0
    refvec=0.d0
    refvec_t=0.d0
    do icycle=1,1 !3
    dvec_t=dvec-refvec_t
    !write(*,'(a,f10.2)') 'GGG ',sqrt(sum(dvec_t**2))
    !do iat=1,atoms%nat
    !    if(sqrt(sum(dvec_t(:,iat)**2))>0.05d0) then
    !    write(*,'(a,i5,3f10.2)') 'EEE ',iat,dvec_t(1,iat),dvec_t(2,iat),dvec_t(3,iat)
    !    endif
    !enddo
    call get_ref_vector_new(parini,atoms,refvec,refvec_qder,qat,pref,icycle,qpm,qt,gwall)
    do iat=1,atoms%nat
    do ii=1,3
        do jat=1,atoms%nat
        do jj=1,3
            tt=0.d0
            do kat=1,atoms%nat
                ss=refvec_qder(1,kat,ii,iat)*refvec_qder(1,kat,jj,jat)+ &
                   refvec_qder(2,kat,ii,iat)*refvec_qder(2,kat,jj,jat)+ &
                   refvec_qder(3,kat,ii,iat)*refvec_qder(3,kat,jj,jat)
                tt=tt+ss !*(dvec(1,kat)**2+dvec(2,kat)**2+dvec(3,kat)**2)
            enddo
            amat(3*(iat-1)+ii,3*(jat-1)+jj)=2.d0*tt
            !dx=atoms%ratp(1,jat)-atoms%ratp(1,iat)
            !dy=atoms%ratp(2,jat)-atoms%ratp(2,iat)
            !dz=atoms%ratp(3,jat)-atoms%ratp(3,iat)
            !!if(iat==jat .and. ii==jj) then
            !write(*,'(4i6,2es14.5)') iat,ii,jat,jj,amat(3*(iat-1)+ii,3*(jat-1)+jj),sqrt(dx*dx+dy*dy+dz*dz)
            !!endif
        enddo
        enddo
    enddo
    enddo
    
    amat_t=amat
    call DSYEV('N','U',3*atoms%nat,amat_t,3*atoms%nat,rhs,work,nwork,info)
    write(*,'(a,i6,es14.5)') 'EVAL ',1,rhs(1)
    write(*,'(a,i6,es14.5)') 'EVAL ',2,rhs(2)
    write(*,'(a,i6,es14.5)') 'EVAL ',3,rhs(3)
    write(*,'(a,i6,es14.5)') 'EVAL ',3*atoms%nat-2,rhs(3*atoms%nat-2)
    write(*,'(a,i6,es14.5)') 'EVAL ',3*atoms%nat-1,rhs(3*atoms%nat-1)
    write(*,'(a,i6,es14.5)') 'EVAL ',3*atoms%nat-0,rhs(3*atoms%nat-0)
    !do ii=1,3*atoms%nat
    !!do ii=2*atoms%nat+1,3*atoms%nat
    !    write(*,'(a,i6,es14.5)') 'EVAL ',ii,rhs(ii)
    !enddo
    !stop 'WWWWWWWWWWWWWW'
    regcoeff(1)=rhs(3*atoms%nat)*1.d-7
    regcoeff(2)=rhs(3*atoms%nat)*1.d-7
    regcoeff(3)=rhs(3*atoms%nat)*1.d-7


    !do iat=1,atoms%nat
    !    do jat=iat+1,atoms%nat
    !        amat(jat,iat)=amat(iat,jat)
    !    enddo
    !enddo
    do iat=1,atoms%nat
    do ii=1,3
        tt=0.d0
        do kat=1,atoms%nat
            tt=tt+refvec_qder(1,kat,ii,iat)*dvec_t(1,kat)+ &
                  refvec_qder(2,kat,ii,iat)*dvec_t(2,kat)+ &
                  refvec_qder(3,kat,ii,iat)*dvec_t(3,kat)
        enddo
        !if(trim(atoms%sat(iat))=='Li') then
        !    qtarget(1)= 0.1d0
        !    qtarget(2)= 0.2d0
        !    qtarget(3)= 0.3d0
        !elseif(trim(atoms%sat(iat))=='S') then
        !    qtarget(1)=-0.1d0
        !    qtarget(2)=-0.2d0
        !    qtarget(3)=-0.3d0
        !else
        !    stop 'ERROR: WHAT TYPE?'
        !endif
        rhs(3*(iat-1)+ii)=2.d0*tt+regcoeff(ii)*2.d0*qtarget(ii,atoms%itypat(iat))
        amat(3*(iat-1)+ii,3*(iat-1)+ii)=amat(3*(iat-1)+ii,3*(iat-1)+ii)+regcoeff(ii)*2.d0
    enddo
    enddo


    call DGETRF(3*atoms%nat,3*atoms%nat,amat,3*atoms%nat,ipiv,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRF info=',info
        stop
    endif
    call DGETRS('N',3*atoms%nat,1,amat,3*atoms%nat,ipiv,rhs,3*atoms%nat,info)
    if(info/=0) then
        write(*,'(a19,i8)') 'ERROR: DGETRS info=',info
        stop
    endif
    !qat=rhs
    do iat=1,atoms%nat
        qat(1,iat)=rhs(3*(iat-1)+1)
        qat(2,iat)=rhs(3*(iat-1)+2)
        qat(3,iat)=rhs(3*(iat-1)+3)
    enddo
    !if(icycle==1) write(*,*) 'AAA',qat(1,5),qat(1,105)
    !if(parini%mpi_env%iproc==0) then
    !do iat=1,atoms%nat
    !    write(*,'(a,3i5,6f10.1)') 'QQQ ',1,icycle,iat,qat(1,iat),qat(2,iat),qat(3,iat), &
    !        pref(1,atoms%itypat(iat)),pref(2,atoms%itypat(iat)),pref(3,atoms%itypat(iat))
    !enddo
    !endif
    call get_ref_vector_new(parini,atoms,refvec,refvec_qder,qat,pref,icycle,qpm,qt,gwall)
    !write(*,'(a,f10.2)') 'FFF ',sqrt(sum(refvec**2))
    refvec_t=refvec_t+refvec
    enddo !end of loop over icycle
    dvec_t=dvec-refvec_t
    !write(*,'(a,f10.2)') 'GGG ',sqrt(sum(dvec_t**2))
    refvec=refvec_t
end subroutine get_ref_vector
!*****************************************************************************************
subroutine cal_rmse_dvec(ndim,vertex,rmse)
    use mod_callback_fit_dvec, only: atoms_smplx=>atoms_smplx_t, parini=>parini_t
    use mod_callback_fit_dvec, only: dvec_all=>dvec_all_t
    use mod_atoms, only: typ_atoms, atom_copy_old, atom_deallocate_old
    implicit none
    integer, intent(in) :: ndim
    real(8), intent(in) :: vertex(ndim)
    real(8), intent(out) :: rmse
    !local variables
    type(typ_atoms):: atoms
    real(8):: errmax, tot_nat, tt
    integer:: iat, iconf
    real(8), allocatable:: refvec(:,:)
    real(8), allocatable:: qat(:,:)
    real(8):: qt(2,2), gwall(3,2,2), pref(3,2)
    !write(*,'(a,4f7.2)') 'vertex ',vertex(1),vertex(2),vertex(3),vertex(4)
    !write(*,'(a,10f7.2)') 'vertex ',vertex(1),vertex(2),vertex(3),vertex(4),vertex(5),vertex(6),vertex(7),vertex(8),vertex(9),vertex(10)
    !write(*,'(a,6f7.2)') 'vertex ',vertex(1),vertex(2),vertex(3),vertex(4),vertex(5),vertex(6)
    write(*,'(a,3(4x,6f7.2))') 'vertex ',vertex(1),vertex(2),vertex(3),vertex(4),vertex(5),vertex(6), &
                                    vertex(7),vertex(8),vertex(9),vertex(10),vertex(11),vertex(12), &
                                    vertex(13),vertex(14),vertex(15),vertex(16),vertex(17),vertex(18)
    rmse=0.d0
    errmax=0.d0
    tot_nat=0.d0
    do iconf=1,atoms_smplx%nconf
        call atom_copy_old(atoms_smplx%atoms(iconf),atoms,'atoms_smplx%atoms(iconf)->atoms')
        allocate(refvec(3,atoms%nat))
        allocate(qat(3,atoms%nat))
        !call get_ref_vector(parini,atoms,refvec,vertex)
        pref(1,1)=vertex(1)
        pref(1,2)=vertex(2)
        gwall(1,1,1)=vertex(3)
        gwall(1,2,1)=vertex(4)
        gwall(1,1,2)=vertex(5)
        gwall(1,2,2)=vertex(6)
        pref(2,1)=vertex(7)
        pref(2,2)=vertex(8)
        gwall(2,1,1)=vertex(9)
        gwall(2,2,1)=vertex(10)
        gwall(2,1,2)=vertex(11)
        gwall(2,2,2)=vertex(12)
        pref(3,1)=vertex(13)
        pref(3,2)=vertex(14)
        gwall(3,1,1)=vertex(15)
        gwall(3,2,1)=vertex(16)
        gwall(3,1,2)=vertex(17)
        gwall(3,2,2)=vertex(18)
        !qt(1,1)=vertex(7)
        !qt(2,1)=vertex(8)
        !qt(1,2)=vertex(9)
        !qt(2,2)=vertex(10)
        call get_ref_vector(parini,atoms,refvec,dvec_all(1,1,iconf),qat,pref,qt,gwall)
        tt=0.d0
        do iat=1,atoms%nat
            tt=tt+(refvec(1,iat)-dvec_all(1,iat,iconf))**2+ &
                      (refvec(2,iat)-dvec_all(2,iat,iconf))**2+ &
                      (refvec(3,iat)-dvec_all(3,iat,iconf))**2
        !write(88,'(i5,2x,a2,2x,6es14.5)') iat,trim(atoms%sat(iat)), &
        !    refvec(1,iat),refvec(2,iat),refvec(3,iat),dvec(1,iat),dvec(2,iat),dvec(3,iat)
        enddo
        write(*,'(a,i6,f10.3)') 'RMSE ',iconf,tt
        rmse=rmse+tt
        errmax=max(errmax,maxval(abs(refvec-dvec_all(:,:,iconf))))
        deallocate(refvec,qat)
        tot_nat=tot_nat+real(atoms%nat,kind=8)
        call atom_deallocate_old(atoms)
    enddo
    !rmse=sqrt(rmse/(3.d0*tot_nat))
    write(*,'(a,2f10.3)') 'rmse,errmax ',rmse,errmax
end subroutine cal_rmse_dvec
!*****************************************************************************************
!function cutoff_function(r, rc) result(fc)
!    implicit none
!    real(8), intent(in):: r, rc
!    !local variables
!    real(8):: fc, pi
!    if(r<rc) then
!        fc=(1.d0-(r/rc)**2)**3
!    else
!        fc=0.d0
!    endif
!    !pi=4.d0*atan(1.d0)
!    !if(r<rc*.5d0) then
!    !    fc=1.d0
!    !elseif(r<rc) then
!    !    fc=cos((1.d0-((1.d0-((r-rc*0.5d0)/(rc*0.5d0))**2)**3))*pi*0.5d0)
!    !else
!    !    fc=0.d0
!    !endif
!    !Compare it with: PHYSICAL REVIEW B 93, 155203 (2016) -> comparison is done:
!    !the second derivative in the cutoff function in PRB paper does not vanish.
!end function cutoff_function
!!*****************************************************************************************
!function cutoff_function_der(r, rc) result(fcd)
!    implicit none
!    real(8), intent(in):: r, rc
!    !local variables
!    real(8):: fcd, pi
!    if(r<rc) then
!        fcd=-6.d0*(r/rc**2)*(1.d0-(r/rc)**2)**2
!    else
!        fcd=0.d0
!    endif
!end function cutoff_function_der
!!*****************************************************************************************
