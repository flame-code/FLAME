!*****************************************************************************************
module mod_train_chi
    use mod_opt_ann, only: typ_opt_ann, typ_cost_object
    use mod_atoms, only: typ_atoms_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_ann, only: typ_ann_arr
    use mod_refdata, only: typ_refdata
    implicit none
    private
    public:: typ_cost_object_chi, ekf_rivals_fitchi, cal_chi_from_features
    type, extends(typ_cost_object):: typ_cost_object_chi
        logical:: initialized=.false.
        integer, pointer:: ita=>null()
        integer:: ndp_train
        integer:: ndp_valid
        type(typ_ann_arr), pointer:: ann_arr=>null()
        type(typ_ann_arr), pointer:: ann_arr_main=>null()
        type(typ_opt_ann), pointer:: opt_ann=>null()
        real(8), allocatable:: features_train(:,:)
        real(8), allocatable:: features_valid(:,:)
        real(8), allocatable:: chi_ref_train(:)
        real(8), allocatable:: chi_ref_valid(:)
        contains
        procedure, public, pass(self):: init_cost_object_chi
        procedure, public, pass(self):: fini_cost_object_chi
        procedure, public, pass(self):: func_value => get_fcn_chi
        procedure, public, pass(self):: func_write => export_weights_chi
        procedure, public, pass(self):: func_evaluate => ann_evaluate_all_chi
    end type typ_cost_object_chi
contains
!*****************************************************************************************
subroutine init_cost_object_chi(self,ann_arr,ann_arr_main,opt_ann,ita,ng)
    implicit none
    class(typ_cost_object_chi), intent(inout):: self
    integer, intent(in), target:: ita
    integer, intent(in), target:: ng
    type(typ_ann_arr), intent(in), target:: ann_arr, ann_arr_main
    type(typ_opt_ann), intent(in), target:: opt_ann
    !local variables
    self%ita=>ita
    self%ann_arr=>ann_arr
    self%ann_arr_main=>ann_arr_main
    self%opt_ann=>opt_ann
    allocate(self%features_train(ng,self%ndp_train))
    allocate(self%features_valid(ng,self%ndp_valid))
    allocate(self%chi_ref_train(self%ndp_train))
    allocate(self%chi_ref_valid(self%ndp_valid))
    self%initialized=.true.
end subroutine init_cost_object_chi
!*****************************************************************************************
subroutine fini_cost_object_chi(self)
    implicit none
    class(typ_cost_object_chi), intent(inout):: self
    !local variables
    nullify(self%ann_arr)
    nullify(self%ann_arr_main)
    deallocate(self%features_train)
    deallocate(self%features_valid)
    deallocate(self%chi_ref_train)
    deallocate(self%chi_ref_valid)
    self%initialized=.false.
end subroutine fini_cost_object_chi
!*****************************************************************************************
subroutine get_fcn_chi(self,parini,idp,opt_ann,fcn_ann,fcn_ref,g)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann
    use mod_atoms, only: typ_atoms, atom_copy_old, update_ratp, atom_deallocate_old
    implicit none
    class(typ_cost_object_chi), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: idp
    type(typ_opt_ann), intent(inout):: opt_ann
    real(8), intent(out):: fcn_ann
    real(8), intent(out):: fcn_ref
    real(8), intent(out):: g(opt_ann%n)
    !local variables
    self%ann_arr%event='train'
    call convert_opt_x_ann_arr(opt_ann,self%ann_arr)
    associate(features_train=>self%features_train)
    call cal_chi_from_features(parini,self%ann_arr,features_train(1,idp),'train',opt_ann%n,g,fcn_ann)
    end associate
    fcn_ref=self%chi_ref_train(idp)
end subroutine get_fcn_chi
!*****************************************************************************************
subroutine export_weights_chi(self,parini,iter)
    use mod_parini, only: typ_parini
    !use mod_ann_io_yaml, only: write_ann_all_yaml
    use mod_ann_io_yaml, only: write_ann_yaml
    use mod_flm_futile
    implicit none
    class(typ_cost_object_chi), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    !local variables
    character(21):: fn
    character(50):: filename
    !call convert_x_ann_arr(n,x,ann_arr)
    if(.not. self%initialized) then
        write(*,'(a)') 'ERROR: typ_cost_object_chi is not initialized!'
        stop
    endif
    associate(ann_arr=>self%ann_arr,ann_arr_main=>self%ann_arr_main,opt_ann=>self%opt_ann)
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
    if(ann_arr_main%exists_yaml_file) then
        ann_arr_main%ann(self%ita)%a(:,:,:)=ann_arr%ann(1)%a(:,:,:)
        ann_arr_main%ann(self%ita)%b(:,:)=ann_arr%ann(1)%b(:,:)
        write(fn,'(a16,i5.5)') '.ann.param.yaml.',iter
        filename=trim(parini%stypat(self%ita))//trim(fn)
        !write(*,'(a)') trim(filename)
        call yaml_comment(trim(filename))
        call write_ann_yaml(parini,filename,ann_arr_main%ann(self%ita),ann_arr_main%rcut)
    else
        stop 'ERROR: export_weights_chi'
    endif
    end associate
end subroutine export_weights_chi
!*****************************************************************************************
subroutine ann_evaluate_all_chi(self,parini,iter)
    use mod_parini, only: typ_parini
    implicit none
    class(typ_cost_object_chi), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    !local variables
    real(8), allocatable:: g_tmp(:)
    associate(ann_arr=>self%ann_arr,opt_ann=>self%opt_ann)
    associate(chi_ref_train=>self%chi_ref_train,chi_ref_valid=>self%chi_ref_valid)
    associate(features_train=>self%features_train,features_valid=>self%features_valid)
    allocate(g_tmp(opt_ann%n))
    call ann_evaluate_fitchi(parini,self%ita,iter,ann_arr,opt_ann%n,g_tmp, &
        self%ndp_train,features_train,chi_ref_train, &
        self%ndp_valid,features_valid,chi_ref_valid)
    deallocate(g_tmp)
    end associate
    end associate
    end associate
end subroutine ann_evaluate_all_chi
!*****************************************************************************************
subroutine ekf_rivals_fitchi(parini,ann_arr_main,opt_ann_main,refdata)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, convert_x_ann_arr
    use mod_atoms, only: typ_atoms_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_ann_io_yaml, only: write_ann_yaml
    use mod_opt_ann, only: typ_opt_ann, ekf_rivals
    use mod_ann_io_yaml, only: get_symfunc_parameters_yaml
    use mod_processors, only: iproc
    use yaml_output
    use futile
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr_main
    type(typ_opt_ann), intent(inout):: opt_ann_main
    type(typ_refdata), intent(inout):: refdata
    !local variables
    integer:: n, i, j, iter, idp, ios, iann
    integer:: iconf, iat, ig, ita, itt
    real(8):: DDOT, tt, den, fcn_ann, fcn_ref
    real(8):: tt1
    type(typ_cost_object_chi):: cost_object_chi
    type(typ_ann_arr):: ann_arr
    type(typ_opt_ann):: opt_ann
    character(5):: stypat
    real(8):: rcut, chi_rmse
    character(50):: fname
    real(8):: ssavg, ssmin, ssmax
    real(8), allocatable:: chi_ref_all_train(:,:)
    real(8), allocatable:: chi_ref_all_valid(:,:)
    character(20):: str1, str2, str3, str4
    character(21):: fn
    character(50):: filename
    associate(atoms_train=>refdata%atoms_train,atoms_valid=>refdata%atoms_valid)
    associate(symfunc_train=>refdata%symfunc_train,symfunc_valid=>refdata%symfunc_valid)
    call refdata%read_refdata(parini,ann_arr_main,parini%ntypat)
    call ann_arr_main%init_ann_arr()
    call opt_ann_main%init_opt_ann(refdata%atoms_train%nconf,ann_arr_main%nwtot,ann_arr_main%iunit)
    call refdata%prepare_refdata(parini,ann_arr_main)
    call set_annweights(parini,opt_ann_main,ann_arr_main)
    ann_arr_main%compute_symfunc=.false.
    write(*,*) opt_ann_main%n
    call ann_arr%set_number_of_ann(1)
    if(ann_arr%nann==0) stop 'ERROR: number of type of atoms zero in ekf_rivals_fitchi'
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
    cost_object_chi%ndp_train=0
    do iconf=1,atoms_train%nconf
        do iat=1,atoms_train%atoms(iconf)%nat
            i=atoms_train%atoms(iconf)%itypat(iat)
            if(i==ita) cost_object_chi%ndp_train=cost_object_chi%ndp_train+1
        enddo
    enddo
    cost_object_chi%ndp_valid=0
    do iconf=1,atoms_valid%nconf
        do iat=1,atoms_valid%atoms(iconf)%nat
            i=atoms_valid%atoms(iconf)%itypat(iat)
            if(i==ita) cost_object_chi%ndp_valid=cost_object_chi%ndp_valid+1
        enddo
    enddo
    call opt_ann%init_opt_ann(cost_object_chi%ndp_train,n,ann_arr_main%iunit)
    call set_annweights(parini,opt_ann,ann_arr)
    if(ann_arr%ann(1)%nn(0)/=symfunc_train%symfunc(1)%ng .or. ann_arr%ann(1)%nn(0)/=symfunc_valid%symfunc(1)%ng) then
        stop 'ERROR: ann_arr%ann(1)%nn(0) and symfunc_train%symfunc(1)%ng differ!'
    endif
    call cost_object_chi%init_cost_object_chi(ann_arr,ann_arr_main,opt_ann,ita, &
        symfunc_train%symfunc(1)%ng)
    idp=0
    do iconf=1,atoms_train%nconf
        do iat=1,atoms_train%atoms(iconf)%nat
            i=atoms_train%atoms(iconf)%itypat(iat)
            if(i==ita) then
            idp=idp+1
            cost_object_chi%chi_ref_train(idp)=chi_ref_all_train(iat,iconf)
            do ig=1,symfunc_train%symfunc(iconf)%ng
                cost_object_chi%features_train(ig,idp)=symfunc_train%symfunc(iconf)%y(ig,iat)
                !write(21,'(es14.5,3i5)') symfunc_train%symfunc(iconf)%y(ig,iat),iconf,iat,ig
            enddo
            !write(22,'(es14.5,2i5)') chi_ref_all_train(iat,iconf),iconf,iat
            endif
        enddo
    enddo
    idp=0
    do iconf=1,atoms_valid%nconf
        do iat=1,atoms_valid%atoms(iconf)%nat
            i=atoms_valid%atoms(iconf)%itypat(iat)
            if(i==ita) then
            idp=idp+1
            cost_object_chi%chi_ref_valid(idp)=chi_ref_all_valid(iat,iconf)
            do ig=1,symfunc_valid%symfunc(iconf)%ng
                cost_object_chi%features_valid(ig,idp)=symfunc_valid%symfunc(iconf)%y(ig,iat)
                !write(21,'(es14.5,3i5)') symfunc_train%symfunc(iconf)%y(ig,iat),iconf,iat,ig
            enddo
            !write(22,'(es14.5,2i5)') chi_ref_all_train(iat,iconf),iconf,iat
            endif
        enddo
    enddo
    tt1=0.d0
    do idp=1,cost_object_chi%ndp_train
        tt1=tt1+cost_object_chi%chi_ref_train(idp)
    enddo
    tt1=tt1/real(cost_object_chi%ndp_train,kind=8)
    tt1=real(int(tt1*1.d2),kind=8)*1.d-2
    ann_arr%ann(1)%chi0=tt1
    ann_arr_main%ann(ita)%chi0=tt1
    call ekf_rivals(cost_object_chi,parini,opt_ann)
    call opt_ann%fini_opt_ann()
    call cost_object_chi%fini_cost_object_chi()
    enddo !end of loop over ita
    do iconf=1,atoms_train%nconf
        deallocate(ann_arr_main%chi_ref_train(iconf)%chis)
    enddo
    deallocate(ann_arr_main%chi_ref_train)
    do iconf=1,atoms_valid%nconf
        deallocate(ann_arr_main%chi_ref_valid(iconf)%chis)
    enddo
    deallocate(ann_arr_main%chi_ref_valid)
    call ann_arr%fini_ann_arr()
    end associate
    end associate
end subroutine ekf_rivals_fitchi
!*****************************************************************************************
subroutine cal_chi_from_features(parini,ann_arr,features,str_dataset,n,g,fcn_ann)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr, convert_ann_epotd
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    real(8), intent(in):: features(1:ann_arr%ann(1)%nn(0))
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
end subroutine cal_chi_from_features
!*****************************************************************************************
subroutine ann_evaluate_fitchi(parini,ita,iter,ann_arr,n,g,ndp_train,ft_train,chi_ref_train,ndp_valid,ft_valid,chi_ref_valid)
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
end subroutine ann_evaluate_fitchi
!*****************************************************************************************
end module mod_train_chi
!*****************************************************************************************
