!*****************************************************************************************
module mod_train_atombased
    use mod_opt_ann, only: typ_opt_ann, typ_cost_object
    use mod_atoms, only: typ_atoms_arr
    use mod_symfunc, only: typ_symfunc_arr
    use mod_ann, only: typ_ann_arr
    use mod_refdata, only: typ_refdata
    implicit none
    private
    public:: ekf_rivals_atombased
    type, extends(typ_cost_object):: typ_cost_object_atombased
        logical:: initialized=.false.
        type(typ_ann_arr), pointer:: ann_arr=>null()
        type(typ_opt_ann), pointer:: opt_ann=>null()
        type(typ_atoms_arr), pointer:: atoms_train=>null()
        type(typ_atoms_arr), pointer:: atoms_valid=>null()
        type(typ_symfunc_arr), pointer:: symfunc_train=>null()
        type(typ_symfunc_arr), pointer:: symfunc_valid=>null()
        contains
        procedure, public, pass(self):: init_cost_object_atombased
        procedure, public, pass(self):: fini_cost_object_atombased
        procedure, public, pass(self):: func_value => get_fcn_ann_atombased
        procedure, public, pass(self):: func_write => export_weights_atombased
        procedure, public, pass(self):: func_evaluate => ann_evaluate_all_atombased
    end type typ_cost_object_atombased
contains
!*****************************************************************************************
subroutine init_cost_object_atombased(self,ann_arr,opt_ann,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    implicit none
    class(typ_cost_object_atombased), intent(inout):: self
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
end subroutine init_cost_object_atombased
!*****************************************************************************************
subroutine fini_cost_object_atombased(self)
    implicit none
    class(typ_cost_object_atombased), intent(inout):: self
    !local variables
    nullify(self%ann_arr)
    nullify(self%atoms_train)
    nullify(self%atoms_valid)
    nullify(self%symfunc_train)
    nullify(self%symfunc_valid)
    self%initialized=.false.
end subroutine fini_cost_object_atombased
!*****************************************************************************************
subroutine ekf_rivals_atombased(parini,ann_arr,opt_ann,refdata)
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
    type(typ_cost_object_atombased):: cost_object_atombased
    associate(atoms_train=>refdata%atoms_train,atoms_valid=>refdata%atoms_valid)
    associate(symfunc_train=>refdata%symfunc_train,symfunc_valid=>refdata%symfunc_valid)
    call refdata%read_refdata(parini,ann_arr,parini%ntypat)
    call ann_arr%init_ann_arr()
    call opt_ann%init_opt_ann(refdata%atoms_train%nconf,ann_arr%nwtot,ann_arr%iunit)
    call refdata%prepare_refdata(parini,ann_arr)
    call set_annweights(parini,opt_ann,ann_arr)
    ann_arr%compute_symfunc=.false.
    call cost_object_atombased%init_cost_object_atombased(ann_arr,opt_ann,atoms_train,atoms_valid,symfunc_train,symfunc_valid)
    call ekf_rivals(cost_object_atombased,parini,opt_ann)
    call cost_object_atombased%fini_cost_object_atombased()
    end associate
    end associate
end subroutine ekf_rivals_atombased
!*****************************************************************************************
subroutine get_fcn_ann_atombased(self,parini,idp,opt_ann,fcn_ann,fcn_ref,g)
    use mod_parini, only: typ_parini
    use mod_ann, only: typ_ann_arr
    use mod_opt_ann, only: typ_opt_ann
    use mod_atoms, only: typ_atoms, atom_copy_old, update_ratp, atom_deallocate_old
    implicit none
    class(typ_cost_object_atombased), intent(inout):: self
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
        write(*,'(a)') 'ERROR: typ_cost_object_atombased is not initialized!'
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
                ann_grad(j,i)=ann_grad(j,i)+ann_arr%g_per_atom(j,iat)
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
end subroutine get_fcn_ann_atombased
!*****************************************************************************************
subroutine export_weights_atombased(self,parini,iter)
    use mod_parini, only: typ_parini
    use mod_ann_io_yaml, only: write_ann_all_yaml
    implicit none
    class(typ_cost_object_atombased), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    !local variables
    if(.not. self%initialized) then
        write(*,'(a)') 'ERROR: typ_cost_object_atombased is not initialized!'
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
end subroutine export_weights_atombased
!*****************************************************************************************
subroutine ann_evaluate_all_atombased(self,parini,iter)
    use mod_parini, only: typ_parini
    use mod_symfunc, only: typ_symfunc
    implicit none
    class(typ_cost_object_atombased), intent(inout):: self
    type(typ_parini), intent(in):: parini
    integer, intent(in):: iter
    !local variables
    if(.not. self%initialized) then
        write(*,'(a)') 'ERROR: typ_cost_object_atombased is not initialized!'
        stop
    endif
    associate(ann_arr=>self%ann_arr,opt_ann=>self%opt_ann)
    associate(atoms_train=>self%atoms_train,atoms_valid=>self%atoms_valid)
    associate(symfunc_train=>self%symfunc_train,symfunc_valid=>self%symfunc_valid)
    call convert_opt_x_ann_arr(opt_ann,ann_arr)
    call ann_evaluate(parini,iter,ann_arr,symfunc_train,atoms_train,"train")
    call ann_evaluate(parini,iter,ann_arr,symfunc_valid,atoms_valid,"valid")
    end associate
    end associate
    end associate
end subroutine ann_evaluate_all_atombased
!*****************************************************************************************
end module mod_train_atombased
!*****************************************************************************************
