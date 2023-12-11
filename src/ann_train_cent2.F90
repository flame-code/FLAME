!*****************************************************************************************
module mod_train_cent2
    use mod_opt_ann, only: typ_opt_ann, typ_cost_object
    use mod_atoms, only: typ_atoms_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_ann, only: typ_ann_arr
    use mod_refdata, only: typ_refdata
    implicit none
    private
    public:: ekf_rivals_cent2
    type, extends(typ_cost_object):: typ_cost_object_cent2
        logical:: initialized=.false.
        type(typ_ann_arr), pointer:: ann_arr=>null()
        type(typ_opt_ann), pointer:: opt_ann=>null()
        type(typ_atoms_arr), pointer:: atoms_train=>null()
        type(typ_atoms_arr), pointer:: atoms_valid=>null()
        type(typ_symfunc_arr), pointer:: symfunc_train=>null()
        type(typ_symfunc_arr), pointer:: symfunc_valid=>null()
        contains
        procedure, public, pass(self):: init_cost_object_cent2
        procedure, public, pass(self):: fini_cost_object_cent2
        procedure, public, pass(self):: func_value => get_fcn_ann_cent2
        procedure, public, pass(self):: func_write => export_weights_cent2
        procedure, public, pass(self):: func_evaluate => ann_evaluate_all_cent2
    end type typ_cost_object_cent2
contains
!*****************************************************************************************
subroutine init_cost_object_cent2(self,ann_arr,opt_ann,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    implicit none
    class(typ_cost_object_cent2), intent(inout):: self
    type(typ_ann_arr), intent(in), target:: ann_arr
    type(typ_opt_ann), intent(in), target:: opt_ann
    type(typ_atoms_arr), intent(in), target:: atoms_train
    type(typ_atoms_arr), intent(in), target:: atoms_valid
    type(typ_symfunc_arr), intent(in), target:: symfunc_train
    type(typ_symfunc_arr), intent(in), target:: symfunc_valid
    !local variables
    self%ann_arr=>ann_arr
    self%opt_ann=>opt_ann
    self%atoms_train=>atoms_train
    self%atoms_valid=>atoms_valid
    self%symfunc_train=>symfunc_train
    self%symfunc_valid=>symfunc_valid
    self%initialized=.true.
end subroutine init_cost_object_cent2
!*****************************************************************************************
subroutine fini_cost_object_cent2(self)
    implicit none
    class(typ_cost_object_cent2), intent(inout):: self
    !local variables
    nullify(self%ann_arr)
    nullify(self%atoms_train)
    nullify(self%atoms_valid)
    nullify(self%symfunc_train)
    nullify(self%symfunc_valid)
    self%initialized=.false.
end subroutine fini_cost_object_cent2
!*****************************************************************************************
subroutine ekf_rivals_cent2(parini,ann_arr,opt_ann,refdata)
    use mod_parini, only: typ_parini
    use mod_opt_ann, only: typ_opt_ann, ekf_rivals
    use mod_atoms, only: typ_atoms_arr
    use mod_ann, only: typ_ann_arr
    use mod_symfunc, only: typ_symfunc_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    type(typ_opt_ann), intent(inout):: opt_ann
    type(typ_refdata), intent(inout):: refdata
    !local variables
    type(typ_cost_object_cent2):: cost_object_cent2
    associate(atoms_train=>refdata%atoms_train,atoms_valid=>refdata%atoms_valid)
    associate(symfunc_train=>refdata%symfunc_train,symfunc_valid=>refdata%symfunc_valid)
    call refdata%read_refdata(parini,ann_arr,parini%ntypat)
    call ann_arr%init_ann_arr()
    call opt_ann%init_opt_ann(refdata%atoms_train%nconf,ann_arr%nwtot,ann_arr%iunit)
    call refdata%prepare_refdata(parini,ann_arr)
    call set_annweights(parini,opt_ann,ann_arr)
    ann_arr%compute_symfunc=.false.
    call cost_object_cent2%init_cost_object_cent2(ann_arr,opt_ann,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    call ekf_rivals(cost_object_cent2,parini,opt_ann)
    call cost_object_cent2%fini_cost_object_cent2()
    end associate
    end associate
end subroutine ekf_rivals_cent2
!*****************************************************************************************
subroutine get_fcn_ann_cent2(self,parini,idp,opt_ann,fcn_ann,fcn_ref,g)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann
    use mod_atoms, only: typ_atoms, atom_copy_old, update_ratp, atom_deallocate_old
    implicit none
    class(typ_cost_object_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: idp
    type(typ_opt_ann), intent(inout):: opt_ann
    real(8), intent(out):: fcn_ann
    real(8), intent(out):: fcn_ref
    real(8), intent(out):: g(opt_ann%n)
    !local variables
    type(typ_atoms):: atoms
    real(8), allocatable:: ann_grad(:,:)
    integer:: iat, i, j, iconf, ixyz
    if(.not. self%initialized) then
        write(*,'(a)') 'ERROR: typ_cost_object_cent2 is not initialized!'
        stop
    endif
    associate(ann_arr=>self%ann_arr)
    associate(atoms_train=>self%atoms_train,symfunc_train=>self%symfunc_train)
    ann_arr%event='train'
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
    !-----------------------------------------------------------------
    iconf=idp
    !-----------------------------------------------------------------
    call atom_copy_old(atoms_train%atoms(iconf),atoms,'atoms_train%atoms(iconf)->atoms')
    !IMPORTANT: it will be fixed: opt_ann will not be passed to cal_ann_main, for the
    !           moment it is passed because of cal_ann_tb 
    call cal_ann_main(parini,atoms,symfunc_train%symfunc(iconf),ann_arr,opt_ann)
    !-----------------------------------------------------------------
    allocate(ann_grad(ann_arr%nweight_max,ann_arr%nann),source=0.d0)
        do iat=1,atoms%nat
            i=atoms%itypat(iat)
            do j=1,ann_arr%nweight_max
                !ann_grad(j,i)=ann_grad(j,i)+atoms%qat(iat)*ann_arr%g_per_atom(j,iat)
                ann_grad(j,i)=ann_grad(j,i)+ann_arr%g_per_atom(j,iat) !WHY IS THIS NOT MULTIPLIED BY CHARGE??!!
            enddo
        enddo
        call set_opt_ann_grad(ann_arr,ann_grad,opt_ann%n,g)
    deallocate(ann_grad)
    !-----------------------------------------------------------------
    fcn_ann=atoms%epot
    fcn_ref=atoms_train%atoms(iconf)%epot
    call atom_deallocate_old(atoms)
    end associate
    end associate
end subroutine get_fcn_ann_cent2
!*****************************************************************************************
subroutine export_weights_cent2(self,parini,iter)
    use mod_parini, only: typ_parini
    use mod_ann_io_yaml, only: write_ann_all_yaml
    implicit none
    class(typ_cost_object_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    !local variables
    if(.not. self%initialized) then
        write(*,'(a)') 'ERROR: typ_cost_object_cent2 is not initialized!'
        stop
    endif
    associate(ann_arr=>self%ann_arr,opt_ann=>self%opt_ann)
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
    if(ann_arr%exists_yaml_file) then
        call write_ann_all_yaml(parini,ann_arr,iter)
    else
        call write_ann_all(parini,ann_arr,iter)
    endif
    end associate
end subroutine export_weights_cent2
!*****************************************************************************************
subroutine ann_evaluate_all_cent2(self,parini,iter)
    use mod_parini, only: typ_parini
    use mod_symfunc, only: typ_symfunc
    use mod_symfunc, only: typ_symfunc_arr
    implicit none
    class(typ_cost_object_cent2), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    !local variables
    if(.not. self%initialized) then
        write(*,'(a)') 'ERROR: typ_cost_object_cent2 is not initialized!'
        stop
    endif
    associate(ann_arr=>self%ann_arr,opt_ann=>self%opt_ann)
    associate(atoms_train=>self%atoms_train,atoms_valid=>self%atoms_valid)
    associate(symfunc_train=>self%symfunc_train,symfunc_valid=>self%symfunc_valid)
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
    call analyze_epoch_init_cent2(parini,ann_arr)
    call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
    call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
    call analyze_epoch_print_cent2(parini,iter,ann_arr)
    end associate
    end associate
    end associate
end subroutine ann_evaluate_all_cent2
!*****************************************************************************************
subroutine analyze_epoch_init_cent2(parini,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    implicit none
    type(typ_parini), intent(in):: parini
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    ann_arr%natsum(1:10)=0
    ann_arr%qmin(1:10)=huge(1.d0)
    ann_arr%qmax(1:10)=-huge(1.d0)
    ann_arr%qsum(1:10)=0.d0
    ann_arr%chi_min(1:10)=huge(1.d0)
    ann_arr%chi_max(1:10)=-huge(1.d0)
    ann_arr%chi_sum(1:10)=0.d0
    ann_arr%chi_delta(1:10)=0.d0
end subroutine analyze_epoch_init_cent2
!*****************************************************************************************
subroutine analyze_epoch_print_cent2(parini,iter,ann_arr)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_flm_futile
    implicit none
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    type(typ_ann_arr), intent(inout):: ann_arr
    !local variables
    integer:: i, ios
    real(8):: ttavg, ttmin, ttmax, ssavg, ssmin, ssmax
    character(50):: fn_charge, fn_chi
    character(20):: str_key
    do i=1,parini%ntypat
        ttavg=ann_arr%qsum(i)/real(ann_arr%natsum(i),8)
        ttmin=ann_arr%qmin(i)
        ttmax=ann_arr%qmax(i)

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
end subroutine analyze_epoch_print_cent2
!*****************************************************************************************
end module mod_train_cent2
!*****************************************************************************************
