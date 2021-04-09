!> @file
!! Define operations to handle and define input files
!! @author
!!    Copyright (C) 2014-2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS

module f_input_file
  use dictionaries
  use yaml_strings, only: operator(.eqv.),f_strcpy
  use f_utils, only: f_zero
  implicit none

  private

  !> Error ids for this module.
  integer, parameter :: ERR_UNDEF=-1
  integer :: INPUT_VAR_NOT_IN_LIST = ERR_UNDEF
  integer :: INPUT_VAR_NOT_IN_RANGE = ERR_UNDEF
  integer :: INPUT_VAR_ILLEGAL = ERR_UNDEF

  character(len = *), parameter, public :: ATTRS = "_attributes"
  character(len = *), parameter :: PROF_KEY = "PROFILE_FROM"
  character(len = *), parameter :: USER_KEY = "USER_DEFINED"

  !> Reserved keywords
  character(len = *), parameter :: COMMENT = "COMMENT"         !< Short description of the key (documentation)
  character(len = *), parameter :: DESCRIPTION = "DESCRIPTION" !< Description of the key (doc)
  character(len = *), parameter :: RANGE = "RANGE"         !< Specify the range of permitted values
  character(len = *), parameter :: EXCLUSIVE = "EXCLUSIVE" !< Exclude some values
  character(len = *), parameter :: DEFAULT = "default"     !< Default value
  character(len = *), parameter :: COND = "CONDITION" !< Indicate a condition for the use (WHEN, ...)
  character(len = *), parameter :: WHEN = "WHEN"             !< Used in a condition
  character(len = *), parameter :: WHEN_NOT = "WHEN_NOT"
  character(len = *), parameter :: MASTER_KEY = "MASTER_KEY"
  character(len = *), parameter :: IMPORT_KEY = "import"

  !internal variable of the module, only used to dump errors.
  type(dictionary), pointer :: failed_exclusive

  !> for internal flib usage
  public :: input_file_errors
  public :: input_value_is_default,input_keys_get_profile,input_file_complete,input_file_dump
  public :: input_file_minimal

contains

  !> Callback routine when an error occurs
  subroutine abort_excl()
    use yaml_output
    use dictionaries
    implicit none
    call yaml_mapping_open("Allowed values for the variable")
    call yaml_dict_dump(failed_exclusive)
    call yaml_mapping_close()
    call f_err_severe()
  end subroutine abort_excl

  !> Callback routine for illegal input variables
  subroutine warn_illegal()
    implicit none

  end subroutine warn_illegal

  !> define the errors of the module
  subroutine input_file_errors()
    use dictionaries, only: f_err_define
    implicit none

    if (INPUT_VAR_NOT_IN_LIST == ERR_UNDEF) then
       call f_err_define(err_name='INPUT_VAR_NOT_IN_LIST',&
            err_msg='given value not in allowed list.',&
            err_action='choose a value from the list below.',&
            err_id=INPUT_VAR_NOT_IN_LIST,callback=abort_excl)
    end if
    if (INPUT_VAR_NOT_IN_RANGE == ERR_UNDEF) then
       call f_err_define(err_name='INPUT_VAR_NOT_IN_RANGE',&
            err_msg='given value not in allowed range.',&
            err_action='adjust the given value.',&
            err_id=INPUT_VAR_NOT_IN_RANGE)
    end if
    if (INPUT_VAR_ILLEGAL == ERR_UNDEF) then
       call f_err_define(err_name='INPUT_VAR_ILLEGAL',&
            err_msg='provided variable is not allowed in this context.',&
            err_action='correct or remove the input variable.',&
            err_id=INPUT_VAR_ILLEGAL,callback=warn_illegal)
    end if

  end subroutine input_file_errors


  !> Check and complete input file
  subroutine input_file_complete(inputdef,dict,imports,nocheck,verbose)
    use dynamic_memory
    use yaml_output
    implicit none
    !> Dictionary of input definitions
    type(dictionary), pointer :: inputdef
    !> Dictionary of the input files
    type(dictionary), pointer :: dict
    !> List of the keys which should not be checked
    type(dictionary), pointer, optional :: nocheck
    !> Dictionary of the preloaded input parameters associated to the importing
    !! for the input file
    type(dictionary), pointer, optional :: imports
    !> variable controlling the verbosity of the output, in particular
    !! when some errors have been raised in the input variables filling
    logical, intent(in), optional :: verbose

    !local variables
    logical :: localcheck,verb
    type(dictionary), pointer :: dict_tmp,iter,dict_tmp2

    call f_routine(id='input_file_complete')

    !if present imports, the user dictionary has to be saved for overriding
    if (present(imports) .and. (IMPORT_KEY .in. dict)) then
       if (associated(imports)) then
          dict_tmp => dict
          nullify(dict)
          !distinguish now the presence of a list or a scalar in the import keyword
          dict_tmp2 => dict_tmp//IMPORT_KEY
          if (dict_len(dict_tmp2) > 0) then
             !iterate on the list elements
             iter => dict_iter(dict_tmp2)
             do while(associated(iter))
                call dict_update(dict,imports//dict_value(iter))
                iter => dict_next(iter)
             end do
          else if (dict_size(dict_tmp2) > 0 ) then
             call f_err_throw(err_id = INPUT_VAR_ILLEGAL, &
                  err_msg ='The "'//IMPORT_KEY//'" key only accepts scalars or lists')
          else
             !Only one entry
             call dict_copy(src=imports//dict_value(dict_tmp2),dest=dict)
          end if
          !Remove COMMENT and DESCRIPTION entries
          if (COMMENT .in. dict) call dict_remove(dict,COMMENT)
          if (DESCRIPTION .in. dict) call dict_remove(dict,DESCRIPTION)
          !then override with the input defaults
          call dict_update(dict,dict_tmp)
          call dict_free(dict_tmp)
       end if
    end if
    verb=.false.
    if (present(verbose)) verb=verbose

    localcheck=.true.
    dict_tmp => dict_iter(inputdef)
    do while (associated(dict_tmp))
       if (present(nocheck)) localcheck = dict_key(dict_tmp) .notin. nocheck
       call input_keys_fill(inputdef,dict,trim(dict_key(dict_tmp)),localcheck,verb)
       dict_tmp => dict_next(dict_tmp)
    end do

    call f_release_routine()

  end subroutine input_file_complete
  

  subroutine input_keys_set(inputdef,userDef, dict, file, key)
    use dictionaries
    use yaml_output
    use dynamic_memory
    implicit none
    logical, intent(out) :: userDef
    !> dictionaries of the input definitions
    type(dictionary), pointer :: inputdef
    !> user input file
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file, key

    integer :: i
    type(dictionary), pointer :: ref,iter
    character(len = max_field_length) :: val, profile_
    double precision, dimension(2) :: rg
!!$    integer :: skeys
!!$    logical :: found
!!$    character(len = max_field_length), dimension(:), allocatable :: keys

    !    call f_routine(id='input_keys_set')

    ref => inputdef // file // key

    !profile_(1:max_field_length) = " "
    !if (trim(profile_) == "") profile_(1:max_field_length) = DEFAULT
    call f_strcpy(src=DEFAULT,dest=profile_)

    userDef = key .in. dict !(has_key(dict, key))
    if (userDef) then
       ! Key should be present only for some unmet conditions.
       if (.not.set_(dict, ref)) then
          !to see if the f_release_routine has to be controlled automatically
          ! call f_release_routine() !to be called before raising the error
          call f_err_throw(err_id = INPUT_VAR_ILLEGAL, &
               err_msg = trim(file) // "/" // trim(key) // " is not allowed in this context.")
          return
       end if
       val = dict // key

       ! There is already a value in dict.
       if (has_key(ref, val)) then
          ! The key was asking for a profile, we copy it.
          call dict_copy(dict // key, ref // val)
          profile_ = val
       else
          ! The value is user-defined, we validate it with info from dict from.
          if (has_key(ref, RANGE)) then
!!$             rg(1) = ref // RANGE // 0
!!$             rg(2) = ref // RANGE // 1
             rg = ref // RANGE
             call validate(dict // key, key, rg)
!!$          else if (has_key(ref, EXCLUSIVE)) then
!!$             failed_exclusive => ref // EXCLUSIVE
!!$             allocate(keys(dict_size(failed_exclusive)))
!!$             keys = dict_keys(failed_exclusive)
!!$             found = .false.
!!$             skeys = size(keys)
!!$             do i = 1, skeys, 1
!!$                found = trim(val) .eqv. trim(keys(i))
!!$                if (found) exit
!!$             end do
!!$             deallocate(keys)
!!$             if (.not. found) then
!!$                call f_err_throw(err_id = INPUT_VAR_NOT_IN_LIST, &
!!$                     err_msg = trim(key) // " = '" // trim(val) //&
!!$                     "' is not allowed, see above the allowed values.")
!!$                nullify(failed_exclusive)
!!$                return
!!$             end if
          else if (EXCLUSIVE .in. ref) then
             failed_exclusive => ref // EXCLUSIVE
             if (val .notin. failed_exclusive) then
                call f_err_throw(err_id = INPUT_VAR_NOT_IN_LIST, &
                     err_msg = trim(key) // " = '" // trim(val) //&
                     "' is not allowed, see above the allowed values.")
                nullify(failed_exclusive)
                return
             end if
          end if
       end if
    else

       ! Key should be present only for some unmet conditions.
       if (.not.set_(dict, ref)) then
!!$          call f_err_throw(err_id = INPUT_VAR_ILLEGAL, &
!!$               & err_msg = trim(file) // "/" // trim(key) // " has to be presentd with a master key.")
          !          call f_release_routine()
          !print *,trim(key)'XXXXXXXXXXXXx'
          return
       end if

       ! Hard-coded profile from key.
       if (PROF_KEY .in. ref) then
          val = ref // PROF_KEY !this retrieve the value of the driver key
          !if might be a profile
          if (val .in. dict) then
             profile_ = dict // val
          end if
       end if

       ! There is no value in dict, we take it from ref.
       !first check if the provided value is among the profiles of ref
       if (profile_ .notin. ref) then
          !it still might be one of the values of the profiles of ref
          nullify(iter)
          do while(iterating(iter,on=inputdef // file // val))
             if (dict_value(iter) .eqv. profile_) then
                profile_=dict_key(iter)
                exit
             end if
          end do
       end if
       if ( profile_ .notin. ref) profile_ = DEFAULT
       !still search if the chosen profile correspons to the value of another profile
       val = dict_value(ref // profile_)
       if (val .in. ref) then
          call dict_copy(dict // key, ref // val)
       else
          call dict_copy(dict // key, ref // profile_)
       end if
    end if

    ! Copy the comment.
    if (has_key(ref, COMMENT)) &
         call dict_copy(dict // (trim(key) // ATTRS) // COMMENT, ref // COMMENT)
    ! Save the source.
    if (userDef) &
         call set(dict // (trim(key) // ATTRS) // USER_KEY, .true.)
    if (profile_ /= DEFAULT) &
         call set(dict // (trim(key) // ATTRS) // PROF_KEY, profile_)

    !    call f_release_routine()
  contains

    function set_(dict, ref)
      implicit none
      type(dictionary), pointer :: dict, ref
      logical :: set_
      !local variables
      logical :: l1
      type(dictionary), pointer :: tmp,tmp0,tmp_not,iter
      character(max_field_length) :: mkey, val_master, val_when

!!$      set_ = .true.
!!$      if (has_key(ref, COND)) then
!!$         mkey = ref // COND // MASTER_KEY
!!$         if (.not. has_key(dict, mkey)) then
!!$            set_ = .false.
!!$            return
!!$         end if
!!$         val_master = dict // mkey
!!$         set_ = .false.
!!$         tmp => ref // COND // WHEN
!!$         do j = 0, dict_len(tmp) - 1, 1
!!$            val_when = tmp // j
!!$            set_ = set_ .or. (trim(val_master) .eqv. trim(val_when))
!!$         end do
!!$      end if

      set_ = COND .notin. ref
      if (set_) return !there are no conditions on the reference variable
      tmp0 => ref // COND
      mkey = tmp0 // MASTER_KEY
      set_ = mkey .in. dict
      if (.not. set_) return !the variable is not present, not coherent
      val_master = dict // mkey
      tmp = tmp0 .get. WHEN
      tmp_not = tmp0 .get. WHEN_NOT
      !call yaml_map('val_master',val_master)
      !call yaml_map('when',tmp)
      !call yaml_map('whennot',tmp_not)
      !call yaml_map('intmp',[(val_master .in. tmp),(val_master .notin. tmp_not)])
      l1=(val_master .in. tmp) .or. .not. associated(tmp)
      set_ = l1 .and. (val_master .notin. tmp_not)
      !call yaml_map('set_',set_)
      if (set_) return !still check if the value is coherent with the profile
      tmp0 => inputdef // file // mkey
      nullify(iter)
      do while(iterating(iter,on=tmp) .and. .not. set_)
         call f_zero(val_when)
         !call yaml_map('val',dict_value(iter))
         !call yaml_map('tmp0',tmp0)
         val_when = tmp0 .get. dict_value(iter)
         set_ = trim(val_master) .eqv. trim(val_when)
      end do
      if (.not. set_) return
      nullify(iter)
      do while(iterating(iter,on=tmp_not) .and. set_)
         call f_zero(val_when)
         !call yaml_map('valnot',dict_value(iter))
         val_when = tmp0 .get. dict_value(iter)
         set_ = .not. (trim(val_master) .eqv. trim(val_when))
      end do
    end function set_

    recursive subroutine validate(dict, key, rg)
      implicit none
      type(dictionary), pointer :: dict
      character(len = *), intent(in) :: key
      double precision, dimension(2), intent(in) :: rg

      character(len = max_field_length) :: val
      character(max_field_length), dimension(:), allocatable :: keys
      double precision :: var
      integer :: dlen, skeys

      if (associated(dict%child)) then
         if (dict_len(dict) >= 1) then
            ! List case.
            dlen = dict_len(dict)
            do i = 0, dlen - 1, 1
               call validate(dict // i, key, rg)
            end do
         else
            ! Dictionary case
            allocate(keys(dict_size(dict)))
            keys = dict_keys(dict)
            skeys = size(keys)
            do i = 1, skeys, 1
               call validate(dict // keys(i), key, rg)
            end do
            deallocate(keys)
         end if
      else
         var = dict
         if (var < rg(1) .or. var > rg(2)) then
            val = dict
            call f_err_throw(err_id = INPUT_VAR_NOT_IN_RANGE, &
                 & err_msg = trim(key) // " = '" // trim(val) // &
                 & "' not in range.")
            return
         end if
      end if
    end subroutine validate
  END SUBROUTINE input_keys_set

  subroutine input_keys_fill(inputdef,dict, file,check,verbose)
    use dynamic_memory
    use yaml_output
    implicit none
    !Arguments
    !> dictionaries of the input definitions
    type(dictionary), pointer :: inputdef
    !> user input file
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file
    logical, intent(in) :: check,verbose
    !Local variables
    !integer :: i
    logical :: user, hasUserDef,docheck
    type(dictionary), pointer :: ref,ref_iter,errs
    !character(len=max_field_length), dimension(:), allocatable :: keys

    docheck=check

    !    call f_routine(id='input_keys_fill')
    if (docheck) call input_keys_control(inputdef,dict,file)

    ref => inputdef // file

    ref_iter => dict_iter(inputdef // file)
    hasUserDef = .false.
!!$    call yaml_map('inputdef_now',inputdef // file)
    do while(associated(ref_iter))
       if (trim(dict_key(ref_iter)) /= DESCRIPTION) then
!!$          call yaml_map('ref_iter',dict_key(ref_iter))
!!$          call yaml_map('ref_iter_val',dict_value(ref_iter))
!!$          call yaml_map('ref_iter_data',dict_value(ref_iter))
!!$          call dump_dict_impl(ref_iter)
          !open a try-catch section to understand where the error is, if any
          call f_err_open_try()
          call input_keys_set(inputdef,user, dict // file, file, dict_key(ref_iter))
          call f_err_close_try(exceptions=errs)
          if (associated(errs)) then
             if (verbose) call yaml_map('List of error found while parsing input variables',errs)
             call dict_free(errs)
             call f_err_throw('Error(s) found in input_keys_fill for the field "'//trim(file)//&
                  '" and the key "'//trim(dict_key(ref_iter))//'", see details in the above message(s),'//&
                  ' before the input file printout.',&
                  err_id=INPUT_VAR_ILLEGAL)
          end if
          hasUserDef = (hasUserDef .or. user)
       end if
       ref_iter=> dict_next(ref_iter)
    end do

    !    call f_release_routine()
  END SUBROUTINE input_keys_fill


  !> Control if all the keys which are defined in a given field are associated with a true input variable
  subroutine input_keys_control(inputdef,dict,file)
    use dictionaries
    use yaml_output, only: yaml_map,yaml_warning
    use yaml_strings, only: yaml_toa
    implicit none
    !> dictionaries of the input definitions
    type(dictionary), pointer :: inputdef
    !> user input file
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: file
    !local variables
    type(dictionary), pointer :: dict_tmp,ref,dict_err,dict_it

    ref=> inputdef // file
    !parse all the keys of the dictionary
    dict_tmp=>dict_iter(dict//file)
    do while(associated(dict_tmp))
       if (.not. (dict_key(dict_tmp) .in. ref) .and. &
            & index(dict_key(dict_tmp), ATTRS) == 0) then
    !      call yaml_map('Allowed keys',dict_keys(ref))
          dict_it=>dict_iter(ref)
          call dict_init(dict_err)
          do while(associated(dict_it))
             if (trim(dict_key(dict_it)) /= DESCRIPTION) &
                  call add(dict_err,dict_key(dict_it))
             dict_it=>dict_next(dict_it)
          end do
          !dict_err=>list_new(.item. dict_keys(ref))
          call yaml_warning('Input file, section "'//file//&
            '"; invalid key "'//trim(dict_key(dict_tmp))//'".')
          call yaml_map('Allowed keys',dict_err)
          call dict_free(dict_err)
          call f_err_throw('An invalid key ('//trim(dict_key(dict_tmp))&
              //') has been found in section "'&
                //file//'". Check above the allowed keys.' ,&
            err_id=INPUT_VAR_ILLEGAL,callback=f_err_severe)
       end if
       dict_tmp=> dict_next(dict_tmp)
    end do
  end subroutine input_keys_control


  function input_keys_get_profile(dict, key, user_defined)
    use dictionaries
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: key
    logical, intent(out) :: user_defined
    character(len = max_field_length) :: input_keys_get_profile
    !local variables
    type(dictionary), pointer :: elem

    user_defined = .false.
    call f_strcpy(src=DEFAULT,dest=input_keys_get_profile)
    !input_keys_get_source(1:max_field_length) = DEFAULT
    elem = dict .get. (trim(key) // ATTRS)
    user_defined = elem .get. USER_KEY
    input_keys_get_profile = elem .get. PROF_KEY

!!$    if (has_key(dict, trim(key) // ATTRS)) then
!!$       if (has_key(dict // (trim(key) // ATTRS), USER_KEY)) &
!!$            & user_defined = dict // (trim(key) // ATTRS) // USER_KEY
!!$       if (has_key(dict // (trim(key) // ATTRS), PROF_KEY)) &
!!$            & input_keys_get_source = dict // (trim(key) // ATTRS) // PROF_KEY
!!$    end if
  end function input_keys_get_profile

  !> return true if the key has not been defined by the user and it is a default variable
  function input_value_is_default(dict, key) result(yes)
    implicit none
    type(dictionary), pointer :: dict
    character(len = *), intent(in) :: key
    logical :: yes
    !local variables
    logical :: user_defined
    character(len = max_field_length) :: prof

    prof = input_keys_get_profile(dict,key,user_defined)
    yes = (trim(prof) == DEFAULT) .and. .not. user_defined

  end function input_value_is_default

  recursive subroutine input_variable_dump(dict,userOnly)
    use yaml_output
    use dictionaries
    use f_utils, only: f_zero
    implicit none
    type(dictionary), pointer :: dict
    logical, intent(in) :: userOnly !< dump only the input variables which are user-defined

    logical :: flow, userDef
    integer :: i, dlen, skeys
    type(dictionary), pointer :: parent, attr, iter
    character(max_field_length) :: descr, tag, prof, output
    character(max_field_length), dimension(:), allocatable :: keys

    if (index(dict%data%key, ATTRS) > 0) return

    call f_zero(descr)! = " "
    call f_zero(tag)! = " "
    call f_zero(prof)! = " "
    userDef = .true.
    parent => dict%parent
    if (associated(parent)) then
       attr = parent .get. (trim(dict%data%key)//ATTRS)
       descr = attr .get. COMMENT
       prof = attr .get. PROF_KEY
       userDef = attr .get. USER_KEY
    end if

    if (dict_len(dict) > 0) then
       ! List case.
       if (userOnly .and. .not.userDef .and. trim(dict%data%key) /= "") return

       flow = (.not.associated(dict%child%child))
       if (.not.flow .and. trim(descr) /= "") then
          call yaml_sequence_open(trim(dict%data%key), tag = tag, advance = "no")
          call yaml_comment(trim(descr), tabbing = 50)
       else
          call yaml_sequence_open(trim(dict%data%key), tag = tag, flow=flow)
       end if
       dlen = dict_len(dict)
       do i = 0, dlen - 1, 1
          call yaml_sequence("", advance = "no")
          call input_variable_dump(dict // i,userOnly)
       end do
       if (flow .and. trim(descr) /= "") then
          call yaml_sequence_close(advance = "no")
          call yaml_comment(trim(descr), tabbing = 50)
       else
          call yaml_sequence_close()
       end if
    else if (dict_size(dict) > 0) then
       ! Dictionary case
       if (userOnly .and. .not.userDef) return

       if (len_trim(dict%data%key) > 0) &
            & call yaml_mapping_open(trim(dict%data%key),flow=.false.)
       iter => dict_next(dict)
       allocate(keys(dict_size(dict)))
       keys = dict_keys(dict)
       skeys = size(keys)
       do i = 1, skeys, 1
          call input_variable_dump(dict // keys(i),userOnly)
       end do
       deallocate(keys)
       if (len_trim(dict%data%key) > 0) call yaml_mapping_close()
    else if (associated(dict)) then
       ! Leaf case.
       if (dict%data%item >= 0) then
          ! List entry
          call yaml_sequence(trim(dict%data%value))
       else
          ! Dictionary entry
          if (userOnly .and. .not.userDef) return

          if (userOnly .and. trim(prof) /= "") then
             output = prof
          else
             output = dict%data%value
          end if
          if (trim(descr) /= "") then
             call yaml_map(trim(dict%data%key), trim(output), tag = tag, advance = "no")
             call yaml_comment(trim(descr), tabbing = 50)
          else
             call yaml_map(trim(dict%data%key), trim(output), tag = tag)
          end if
       end if
    end if

  end subroutine input_variable_dump


  !> This routine is used to create a minimal dictionary which can be used at the place
  !! of the one provided as an indication on the understood variables in inputdef
  subroutine input_file_minimal(inputdef,dict,minimal,nested,as_is)
    use dynamic_memory
    use dictionaries
    use yaml_output
    implicit none
    !Arguments
    type(dictionary), pointer :: inputdef             !< Dictionary of the input definitions
    type(dictionary), pointer, intent(in) :: dict     !< User input file
    type(dictionary), pointer, intent(in) :: as_is    !< Add keys in as_is and not in inputdef
    type(dictionary), pointer, intent(in) :: nested   !< Add other keys (extracted by subcategories of nested elements)
    type(dictionary), pointer, intent(out) :: minimal !< List of keys in dict that remain which require special treatments
                                                      !! which are also in inputdef (with different values) or in as_is
    !Local variables
    type(dictionary), pointer :: dict_tmp,min_cat,dict_tmp0
    character(len=max_field_length) :: category,category0
    logical :: cat_found

    call f_routine(id='input_file_minimal')

    nullify(minimal)

    !recursively search into the reference input variables

    dict_tmp => dict_iter(inputdef)!parameters)

    do while(associated(dict_tmp))
       !for any of the keys of parameters look at the corresponding value of the dictionary
       category=dict_key(dict_tmp)
       !call yaml_map('dict category',parameters//category)
       !print *,'category',trim(category),has_key(dict,category)
       !call yaml_map('dict category',dict_tmp)
       if (has_key(dict,category)) then
          call minimal_category(dict_tmp,dict//category,min_cat)
          if (associated(min_cat)) then
             if (.not. associated(minimal)) call dict_init(minimal)
             call set(minimal//category,min_cat)
          end if
       end if
       !stop
       dict_tmp => dict_next(dict_tmp)
    end do

    !then add other information to the minimal dictionary which is associated
    !to specific system parameters
    !in this way the default values can be extracted by the subcategories
    !of the nested elements
    dict_tmp0 => dict_iter(nested)
    do while(associated(dict_tmp0))
       category0=dict_value(dict_tmp0)
       !! basis set
       if (category0 .in. dict) then
          dict_tmp => dict_iter(dict//category0)
          do while(associated(dict_tmp))
             category=dict_key(dict_tmp)
             !Pb with stack (Cray - ftn 05/2014), solved with the cat_found temporary variable
             cat_found = category .in. inputdef//category0
             if (.not. cat_found .and. index(category,ATTRS) == 0 ) then
                !verify that no parameters correspond to default values
                call minimal_category(inputdef//category0,dict_tmp,min_cat)
                if (associated(min_cat)) then
                   if (.not. associated(minimal)) call dict_init(minimal)
                   call set(minimal//category0//category,min_cat)
                end if
             end if
             dict_tmp => dict_next(dict_tmp)
          end do
       end if
       dict_tmp0 => dict_next(dict_tmp0)
    end do
    !fragment dictionary has to be copied as-is
    !other variables have to follow the same treatment
    dict_tmp => dict_iter(as_is)
    do while(associated(dict_tmp))
       category=dict_value(dict_tmp)
       if (category .in. dict) then
          if (.not. associated(minimal)) call dict_init(minimal)
          call dict_copy(minimal//category,dict//category)
       end if
       dict_tmp => dict_next(dict_tmp)
    end do

    call f_release_routine()

  contains

    subroutine minimal_category(vars,input,minim)
      implicit none
      type(dictionary), pointer :: vars,input,minim
      !local variables
      logical :: profile_found
      character(len=max_field_length) :: def_var,var_prof,prof_var,scalar_tmp
      type(dictionary), pointer :: defvar,var
      nullify(minim)

      var=>dict_iter(vars)
      !        call yaml_map('var dict',var)

      do while(associated(var))
         def_var=dict_key(var)
         !           print *,'here2 ',trim(def_var),has_key(input,def_var)

         !search if the input data have values among the profiles
         if (has_key(input,def_var)) then
            profile_found=.false.
            !see if the dictionary has the PROF_KEY in its possibilities
!!$            call f_zero(prof_var)
!!$            call f_zero(profile_value)
!!$            profile_value=var .get. PROF_KEY
!!$            if (len_trim(profile_value) == 0) prof_var=input .get. profile_value
            prof_var(1:len(prof_var))=' '
            if (has_key(var,PROF_KEY)) then
               if (has_key(input,dict_value(var//PROF_KEY))) &
                    prof_var=dict_value(input//dict_value(var//PROF_KEY))
            end if

            defvar => dict_iter(var)
            !              call yaml_map('var dict inside',defvar)
            check_profile: do while(associated(defvar))
               !exclude keys for definition of the variable
               var_prof=dict_key(defvar)
               !              call yaml_map('key',var_prof)
!!$               if (trim(var_prof) /= COMMENT .and. trim(var_prof) /= COND .and.&
!!$                    trim(var_prof) /= RANGE .and. trim(var_prof) /= PROF_KEY .and. &
!!$                    trim(var_prof) /= EXCLUSIVE) then
               select case(trim(var_prof))
               case(COMMENT,COND,RANGE,PROF_KEY,EXCLUSIVE)
                  !do nothing in these cases
               case default
                  !check if some profile meets desired values
                  !call yaml_map('defvar',defvar)
                  !call yaml_map('input',input//def_var)
                  !call yaml_map('result',defvar == input//def_var)
                  !call yaml_map('var_prof',var_prof)
                  !call yaml_map('test',trim(var_prof) /= DEFAULT .and. var_prof /= prof_var)
                  !print *,'key',def_var
                  profile_found= (defvar == input//def_var)
                  !if it has not been found verify first that the list has not to be compacted
                  !in one sense
                  if (.not. profile_found) then
                     call list_scalar(input//def_var,scalar_tmp)
                     if (len_trim(scalar_tmp) /= 0) profile_found= dict_value(defvar) .eqv. scalar_tmp
                  end if
                  !or in another
                  if (.not. profile_found) then
                     call list_scalar(defvar,scalar_tmp)
                     if (len_trim(scalar_tmp) /= 0) profile_found= dict_value(input//def_var) .eqv. scalar_tmp
                  end if
                  if (profile_found) then
                     if (trim(var_prof) /= DEFAULT .and. var_prof /= prof_var) then
                        if (.not. associated(minim)) call dict_init(minim)
                        call set(minim//def_var,var_prof)
                     end if
                     exit check_profile
                  end if
               end select
!!$               end if
               defvar => dict_next(defvar)
            end do check_profile
            !the key has not been found among the profiles, therefore it should be entered as is
            if (.not. profile_found .and. len_trim(dict_value(input//def_var))/=0) then
               if (.not. associated(minim)) call dict_init(minim)
               call list_scalar(input//def_var,var_prof)
!!$               defvar => dict_iter(input//def_var)
!!$               call f_zero(var_prof)
!!$               !var_prof=repeat(' ',len(var_prof))
!!$               !if the dictionary is not a list there is no need to start
!!$               if (dict_len(input // def_var)==0) nullify(defvar)
!!$               compact_list: do while(associated(defvar))
!!$                  !if scalar, retrieve the value, otherwise exit
!!$                  if (dict_size(defvar) == 0 .and. dict_len(defvar)==0) then
!!$                     prof_var=defvar
!!$                  else
!!$                     !var_prof=repeat(' ',len(var_prof))
!!$                     call f_zero(var_prof)
!!$                     exit compact_list
!!$                  end if
!!$                  !check if all the values of the list are equal to the first one
!!$                  if (len_trim(var_prof) == 0) then
!!$                     var_prof=prof_var
!!$                  else
!!$                     !check if all the values are OK, otherwise exit at first failure
!!$                     if (var_prof /= prof_var) then
!!$                        !var_prof=repeat(' ',len(var_prof))
!!$                        call f_zero(var_prof)
!!$                        exit compact_list
!!$                     end if
!!$                  end if
!!$                  defvar => dict_next(defvar)
!!$               end do compact_list
               !if the dictionary is not a one-level list or if it is a list with different values
               ! copy it as-is
               if (len_trim(var_prof) == 0) then
                  call dict_copy(minim//def_var,input//def_var)
               else !otherwise put the scalar value associated
                  !but first we should check if it is associated to a profile
                  call set(minim//def_var,var_prof)
               end if
            end if
         end if
         var => dict_next(var)
      end do
    end subroutine minimal_category

  end subroutine input_file_minimal

  !>clean the list items if the dictionary is a list with all the items identical
  subroutine list_scalar(list,scalar)
    implicit none
    type(dictionary), pointer :: list
    character(len=max_field_length), intent(out) :: scalar !<empty if failed
    !local variables
    character(len=max_field_length) :: probe
    type(dictionary), pointer :: iter

    !the result is first cleaned
    call f_zero(scalar)
    !if the dictionary is not a list there is no need to start
    if (dict_len(list)==0) return
    iter => dict_iter(list)
    do while(associated(iter))
       !if scalar, retrieve the value, otherwise exit
       if (dict_size(iter) == 0 .and. dict_len(iter)==0) then
          probe=iter
       else
          call f_zero(scalar)
          exit
       end if
       !assign the first value of the list
       if (len_trim(scalar) == 0) then
          call f_strcpy(src=probe,dest=scalar)
       else
          !check if all the values are OK, otherwise exit at first failure
          if (scalar /= probe) then
             call f_zero(scalar)
             exit
          end if
       end if
       iter => dict_next(iter)
    end do

  end subroutine list_scalar


  !> Dump the dictionary of the input variables.
  !! Should dump only the keys relative to the input variables and
  !! print out warnings for the ignored keys
  subroutine input_file_dump(dict, userOnly,nodump_list,msg)
    use yaml_output
    implicit none
    type(dictionary), pointer :: dict   !< Dictionary to dump
    logical, intent(in), optional :: userOnly
    !> Message to be written as comment in the beginning
    !! do not write anything if present and empty
    character(len=*), intent(in), optional :: msg
    type(dictionary), pointer, optional :: nodump_list !<list containing keys not to dump

    !local variables
    integer, parameter :: natoms_dump=500
    logical :: userOnly_,todump

    type(dictionary), pointer :: iter

    userOnly_ = .false.
    if (present(userOnly)) userOnly_ = userOnly

    if (present(msg)) then
       if (len_trim(msg) /= 0) call yaml_comment(msg, hfill = "-")
    else
       call yaml_comment("Input parameters", hfill = "-")
    end if

!!$    iter => dict_iter(dict)
!!$    do while(associated(iter))
    nullify(iter)
    do while(iterating(iter,on=dict))
       todump=.true.
       if (present(nodump_list)) todump = dict_key(iter) .notin. nodump_list
       if (todump) call input_variable_dump(iter,userOnly_)
       !iter => dict_next(iter)
    end do

  end subroutine input_file_dump


end module f_input_file
