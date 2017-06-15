!> @file
!! Manage signals or hooks to call registered methods at
!! given entries in the code. Also defines a declarative
!! API to export objects and related method to language
!! bindings.
!! @author
!!    Copyright (C) 2016 BigDFT group
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS
!!
!! @defgroup FLIB_OBJECT  Define introspectable API for object method and signaling (flib)
!! @ingroup FLIB
!! @brief flib definition of introspectable object with associated methods and signals.
!! @details
!! module_f_objects provides routines to declare objects and their associated
!! methods and signals, so they can be called from language bindings.

!> Module used to define, emit and connect to signals.
!! @ingroup FLIB_OBJECT
!! @details
!! Objects are reference by their id, which is a string to identify
!! the class of a object. Use f_object_new() to register a new object class.
!! The constructor and destructor methods are used by language bindings to
!! create new objects of this class on the fly. Associated methods and functions
!! can be added to a class using f_object_add_method() and f_object_add_function().
module module_f_objects
  use dictionaries
  use f_precisions, only: f_address

  implicit none

  private

  integer, private, save :: ERROR_OBJECT
  type(dictionary), pointer :: class_library => null()

  integer, parameter :: MAX_ARGS_IMPLEMENTED = 5

  type signal_ctx
     character(len = max_field_length) :: obj_id, id
     type(dictionary), pointer :: sig
  end type signal_ctx

  public :: signal_ctx
  public :: f_object_new_, f_object_add_method, f_object_get_method
  public :: f_object_finalize
  public :: f_object_add_signal, f_object_signal_prepare
  public :: f_object_signal_add_arg_, f_object_signal_add_str
  public :: f_object_signal_emit, f_object_signal_connect

  integer, parameter :: MAX_KERNEL_ARGS = 7
  type kernel_ctx
     integer(f_address) :: callback
     integer :: callback_n_args
     
     integer :: n_args
     integer(f_address), dimension(MAX_KERNEL_ARGS) :: args
     integer :: n_strs
     integer, dimension(MAX_KERNEL_ARGS) :: strs
  end type kernel_ctx

  public :: kernel_ctx
  public :: f_object_kernel_new_
  public :: f_object_kernel_add_arg_, f_object_kernel_add_str

  public :: f_object_has_signal

contains

  subroutine ensure_init()
    implicit none

    if (associated(class_library)) return

    call dict_init(class_library)
    call f_err_define('ERROR_OBJECT', 'object class error', &
         ERROR_OBJECT, err_action='Check the code.')
  end subroutine ensure_init

  subroutine f_object_finalize()
    implicit none
    if (associated(class_library)) call dict_free(class_library)
  end subroutine f_object_finalize

  !> @private
  !! Don't use this method directy because it uses address instead of
  !! method names, prefer f_object_new().
  !! @internal
  subroutine f_object_new_(obj_id, constructor_add, destructor_add)
    implicit none
    character(len = *), intent(in) :: obj_id
    integer(f_address), optional, intent(in) :: constructor_add, destructor_add
    
    call ensure_init()
    if (f_err_raise(obj_id .in. class_library, &
         & "object '" // obj_id // "' already exists.", err_id = ERROR_OBJECT)) return

    call set(class_library // obj_id, "")
    if (present(constructor_add)) &
         & call f_object_add_method(obj_id, "constructor", constructor_add, 1)
    if (present(destructor_add)) &
         & call f_object_add_method(obj_id, "destructor", destructor_add, 0)
  end subroutine f_object_new_
  
  subroutine f_object_add_method(obj_id, id, method_add, n_args, isfunc)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    integer(f_address), intent(in) :: method_add
    integer, intent(in) :: n_args
    logical, optional, intent(in) :: isfunc

    call ensure_init()
    if (f_err_raise(.not. (obj_id .in. class_library) .and. obj_id /= "class", &
         & "object '" // obj_id // "' not defined.", err_id = ERROR_OBJECT)) return
    
    call set(class_library // obj_id // "methods" // id // "address", method_add)
    call set(class_library // obj_id // "methods" // id // "n_args", n_args)
    call set(class_library // obj_id // "methods" // id // "function", .false.)
    if (present(isfunc)) &
         & call set(class_library // obj_id // "methods" // id // "function", isfunc)
  end subroutine f_object_add_method

  subroutine f_object_get_method(obj_id, method_id, n_args, isfunc, callback)
    implicit none
    character(len = *), intent(in) :: obj_id, method_id
    integer, intent(out) :: n_args, isfunc
    integer(f_address), intent(out) :: callback

    logical :: isfunc_
    type(dictionary), pointer :: method,tmp
    character(len = max_field_length) :: obj_key, meth_key

    n_args = 0
    callback = 0
    write(obj_key, "(A)") obj_id(1:min(max_field_length, len(obj_id)))
    write(meth_key, "(A)") method_id(1:min(max_field_length, len(method_id)))

    call ensure_init()
    if ((obj_key .notin. class_library) .and. trim(obj_key) /= "class") return
    tmp =>class_library // obj_key // "methods" !for ftn internal compiler error
    if (meth_key .notin. tmp) return

    method => class_library // obj_key // "methods" // meth_key
    n_args = method // "n_args"
    callback = method // "address"
    isfunc = 0
    isfunc_ = method // "function"
    if (isfunc_) isfunc = 1
  end subroutine f_object_get_method

  subroutine f_object_add_signal(obj_id, id, n_args)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    integer, intent(in) :: n_args

    call ensure_init()
    if (f_err_raise(.not. (obj_id .in. class_library) .and. obj_id /= "class", &
         & "object '" // obj_id // "' not defined.", err_id = ERROR_OBJECT)) return
    if (f_err_raise(n_args < 0 .or. n_args > MAX_ARGS_IMPLEMENTED, &
         & "signal '" // id // "' defined with too many arguments.", &
         & err_id = ERROR_OBJECT)) return
    
    call set(class_library // obj_id // "signals" // id // "n_args", n_args)
  end subroutine f_object_add_signal

  function ensure_signal(obj_id, id) result(sig)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    type(dictionary), pointer :: tmp, sig

    call ensure_init()

    nullify(sig)
    if (f_err_raise(.not. (obj_id .in. class_library) .and. obj_id /= "class", &
         & "object '" // obj_id // "' not defined.", err_id = ERROR_OBJECT)) return
    tmp =>  class_library // obj_id // "signals" 
    if (f_err_raise(.not. (id .in. tmp), &
         & "signal '" // obj_id // "::" // id // "' not defined.", err_id = ERROR_OBJECT)) return
    sig => tmp // id
  end function ensure_signal

  function f_object_signal_prepare(obj_id, id, ctx, reset) result(emit)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    type(signal_ctx), intent(out) :: ctx
    logical, intent(in), optional :: reset

    logical :: emit, reset_

    emit = .false.

    ctx%sig => ensure_signal(obj_id, id)
    if (.not. associated(ctx%sig)) return
    write(ctx%obj_id, "(A)") obj_id
    write(ctx%id, "(A)") id

    reset_ = .true.
    if (present(reset)) reset_ = reset
    if (dict_len(ctx%sig // "arguments") > 0 .and. reset_) &
         & call dict_remove(ctx%sig, "arguments")
    if (dict_len(ctx%sig // "strings") > 0 .and. reset_) &
         & call dict_remove(ctx%sig, "strings")
    emit = (dict_len(ctx%sig // "hooks") > 0)
  end function f_object_signal_prepare

  subroutine f_object_signal_add_arg_(ctx, arg_add)
    implicit none
    type(signal_ctx), intent(in) :: ctx
    integer(f_address), intent(in) :: arg_add

    if (f_err_raise(.not. associated(ctx%sig), &
         & trim(ctx%obj_id) // "::" // trim(ctx%id) // " not defined", &
         & err_id = ERROR_OBJECT)) return

    call add(ctx%sig // "arguments", arg_add)
  end subroutine f_object_signal_add_arg_

  subroutine f_object_signal_add_str(ctx, arg)
    use f_precisions, only: f_loc
    implicit none
    type(signal_ctx), intent(in) :: ctx
    character(len = *), intent(inout) :: arg

    if (f_err_raise(.not. associated(ctx%sig), &
         & trim(ctx%obj_id) // "::" // trim(ctx%id) // " not defined", &
         & err_id = ERROR_OBJECT)) return

    call add(ctx%sig // "arguments", f_loc(arg))
    call add(ctx%sig // "strings", len(arg))
  end subroutine f_object_signal_add_str

  subroutine f_object_signal_emit(ctx)
    implicit none
    type(signal_ctx), intent(in) :: ctx
    
    type(dictionary), pointer :: iter
    integer(f_address) :: callback
    integer(f_address), dimension(MAX_ARGS_IMPLEMENTED) :: args
    integer, dimension(MAX_ARGS_IMPLEMENTED) :: lens
    integer :: n_args_signal, n_args_kernel, n_args
    integer :: n_strs_signal, n_strs_kernel, n_strs

    if (f_err_raise(.not. associated(ctx%sig), &
         & trim(ctx%obj_id) // "::" // trim(ctx%id) // " not defined", &
         & err_id = ERROR_OBJECT)) return

    n_args_signal = ctx%sig // "n_args"

    if (f_err_raise(dict_len(ctx%sig // "arguments") /= n_args_signal, &
         & "not enough packed arguments for signal '" // ctx%id // "'.", &
         & err_id = ERROR_OBJECT)) return

    if (n_args_signal > 0) args(1:n_args_signal) = ctx%sig // "arguments"
    n_strs_signal = dict_len(ctx%sig // "strings")
    if (n_strs_signal > 0) lens(1:n_strs_signal) = ctx%sig // "strings"

    iter => dict_iter(ctx%sig // "hooks")
    do while (associated(iter))
       ! Add the kernel arguments here.
       n_args_kernel = dict_len(iter // "arguments")
       n_args = n_args_signal + max(n_args_kernel, 0)
       if (n_args_kernel > 0) args(n_args_signal + 1:n_args) = iter // "arguments"
       
       n_strs_kernel = dict_len(iter // "strings")
       n_strs = n_strs_signal + max(n_strs_kernel, 0)
       if (n_strs_kernel > 0) lens(n_strs_signal + 1:n_strs) = iter // "strings"

       callback = iter // "address"
       
       select case(n_args)
       case (0)
          call callable_void(callback)
       case (1)
          select case(n_strs)
          case (0)
             call callable_arg(callback, args(1))
          case (1)
             call callable_str(callback, args(1), lens(1))
          end select
       case (2)
          select case(n_strs)
          case (0)
             call callable_arg_arg(callback, args(1), args(2))
          case (1)
             call callable_arg_str(callback, args(1), args(2), lens(1))
          case (2)
             call callable_str_str(callback, args(1), args(2), lens(1), lens(2))
          end select
       case (3)
          select case(n_strs)
          case (0)
             call callable_arg_arg_arg(callback, &
                  & args(1), args(2), args(3))
          case (1)
             call callable_arg_arg_str(callback, &
                  & args(1), args(2), args(3), lens(1))
          case (2)
             call callable_arg_str_str(callback, &
                  & args(1), args(2), args(3), lens(1), lens(2))
          case (3)
             call callable_str_str_str(callback, &
                  & args(1), args(2), args(3), lens(1), lens(2), lens(3))
          end select
       case (4)
          select case(n_strs)
          case (0)
             call callable_arg_arg_arg_arg(callback, &
                  & args(1), args(2), args(3), args(4))
          case (1)
             call callable_arg_arg_arg_str(callback, &
                  & args(1), args(2), args(3), args(4), lens(1))
          case (2)
             call callable_arg_arg_str_str(callback, &
                  & args(1), args(2), args(3), args(4), lens(1), lens(2))
          case (3)
             call callable_arg_str_str_str(callback, &
                  & args(1), args(2), args(3), args(4), lens(1), lens(2), lens(3))
          case (4)
             call callable_str_str_str_str(callback, &
                  & args(1), args(2), args(3), args(4), &
                  & lens(1), lens(2), lens(3), lens(4))
          end select
       case (5)
          select case(n_strs)
          case (0)
             call callable_arg_arg_arg_arg_arg(callback, &
                  & args(1), args(2), args(3), args(4), args(5))
          case (1)
             call callable_arg_arg_arg_arg_str(callback, &
                  & args(1), args(2), args(3), args(4), args(5), lens(1))
          case (2)
             call callable_arg_arg_arg_str_str(callback, &
                  & args(1), args(2), args(3), args(4), args(5), lens(1), lens(2))
          case (3)
             call callable_arg_arg_str_str_str(callback, &
                  & args(1), args(2), args(3), args(4), args(5), &
                  & lens(1), lens(2), lens(3))
          case (4)
             call callable_arg_str_str_str_str(callback, &
                  & args(1), args(2), args(3), args(4), args(5), &
                  & lens(1), lens(2), lens(3), lens(4))
          case (5)
             call callable_str_str_str_str_str(callback, &
                  & args(1), args(2), args(3), args(4), args(5), &
                  & lens(1), lens(2), lens(3), lens(4), lens(5))
          end select
       end select
       iter => dict_next(iter)
    end do
  end subroutine f_object_signal_emit

  subroutine f_object_signal_connect(obj_id, id, kernel, sid)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    type(kernel_ctx), intent(in) :: kernel
    integer, intent(out), optional :: sid

    type(dictionary), pointer :: sig, hook
    integer :: n_args_signal, i, sid_

    if (present(sid)) sid = -1

    sig => ensure_signal(obj_id, id)
    if (.not. associated(sig)) return

    n_args_signal = class_library // obj_id // "signals" // id // "n_args"
    if (f_err_raise(n_args_signal + kernel%n_args /= kernel%callback_n_args, &
         & "kernel don't have the right number of arguments for signal '" // &
         & obj_id // "::" // id // "'.", err_id = ERROR_OBJECT)) return

    ! Get the last hook, to retrieve its id.
    hook => class_library // obj_id // "signals" // id // "hooks"
    if (dict_len(hook) > 0) then
       hook => hook // (dict_len(hook) - 1)
       sid_ = hook // "id"
       sid_ = sid_ + 1
    else
       sid_ = 0
    end if

    call dict_init(hook)
    call set(hook // "id", sid_)
    call set(hook // "address", kernel%callback)
    do i = 1, kernel%n_args
       call add(hook // "arguments", kernel%args(i))
    end do
    do i = 1, kernel%n_strs
       call add(hook // "strings", kernel%strs(i))
    end do

    call add(class_library // obj_id // "signals" // id // "hooks", hook)

    if (present(sid)) sid = sid_
  end subroutine f_object_signal_connect

  function f_object_has_signal(obj_id, id)
    implicit none
    character(len = *), intent(in) :: obj_id, id
    logical :: f_object_has_signal

    call ensure_init()

    f_object_has_signal = (obj_id .in. class_library)
    if (f_object_has_signal) &
         & f_object_has_signal = (id .in. class_library // obj_id // "signals")
  end function f_object_has_signal

  function f_object_kernel_new_(callback_add, n_args) result(ctx)
    implicit none
    integer(f_address), intent(in) :: callback_add
    integer, intent(in) :: n_args
    type(kernel_ctx) :: ctx

    if (f_err_raise(n_args > MAX_ARGS_IMPLEMENTED, "Too many arguments for kernel", &
         & err_id = ERROR_OBJECT)) return

    ctx%callback = callback_add
    ctx%callback_n_args = n_args
    ctx%n_args = 0
    ctx%n_strs = 0
  end function f_object_kernel_new_

  subroutine f_object_kernel_add_arg_(ctx, arg_add)
    implicit none
    type(kernel_ctx), intent(inout) :: ctx
    integer(f_address), intent(in) :: arg_add

    if (ctx%n_args == MAX_KERNEL_ARGS) then
       call f_err_throw("Too many arguments for kernel", err_id = ERROR_OBJECT)
       return
    end if
    ctx%n_args = ctx%n_args + 1
    ctx%args(ctx%n_args) = arg_add
  end subroutine f_object_kernel_add_arg_

  subroutine f_object_kernel_add_str(ctx, arg)
    use f_precisions, only: f_loc
    implicit none
    type(kernel_ctx), intent(inout) :: ctx
    character(len = *), intent(in) :: arg

    if (ctx%n_args == MAX_KERNEL_ARGS) then
       call f_err_throw("Too many arguments for kernel", err_id = ERROR_OBJECT)
       return
    end if
    ctx%n_args = ctx%n_args + 1
    ctx%args(ctx%n_args) = f_loc(arg)
    ctx%n_strs = ctx%n_strs + 1
    ctx%strs(ctx%n_args) = len(arg)
  end subroutine f_object_kernel_add_str
end module module_f_objects



!> @relates module_f_objects
!! Define a new object with given routines for constructor
!! and destructors.
subroutine f_object_new(id, constructor, destructor)
  use f_precisions
  use module_f_objects, only: wrapper_new => f_object_new_
  character(len = *), intent(in) :: id
  external :: constructor, destructor
  
  call wrapper_new(id, f_loc(constructor), f_loc(destructor))
end subroutine f_object_new

!> @relates module_f_objects
!! Define the name of the method with \c id. The number of arguments
!! the method is using is also mandatory to allow runtime checks.
subroutine f_object_add_method(obj_id, id, method, n_args)
  use f_precisions
  use module_f_objects, only: wrapper_add => f_object_add_method
  character(len = *), intent(in) :: obj_id !< Object class identifier
  character(len = *), intent(in) :: id     !< Method name
  integer, intent(in) :: n_args            !< Number of argumnets
  external :: method                       !< Subroutine corresponding to the method
  
  call wrapper_add(obj_id, id, f_loc(method), n_args)
end subroutine f_object_add_method

!> @relates module_f_objects
subroutine f_object_add_function(obj_id, id, method, n_args)
  use f_precisions
  use module_f_objects, only: wrapper_add => f_object_add_method
  character(len = *), intent(in) :: obj_id, id
  integer, intent(in) :: n_args
  external :: method
  
  call wrapper_add(obj_id, id, f_loc(method), n_args, .true.)
end subroutine f_object_add_function

!> @relates module_f_objects
subroutine f_object_get_method(obj_id, method_id, n_args, isfunc, callback)
  use f_precisions
  use module_f_objects, only: wrapper_get => f_object_get_method
  character(len = *), intent(in) :: obj_id, method_id
  integer, intent(out) :: n_args, isfunc
  integer(f_address), intent(out) :: callback
  
  call wrapper_get(obj_id, method_id, n_args, isfunc, callback)
end subroutine f_object_get_method

!> @relates module_f_objects
!! Add a new argument to the stack before calling f_object_signal_emit().
!! Use "class" as obj_id to add an argument to a generic signal.
subroutine f_object_signal_add_arg(ctx, arg)
  use f_precisions
  use module_f_objects, only: wrapper_add => f_object_signal_add_arg_, signal_ctx
  type(signal_ctx), intent(in) :: ctx !< signal context.
  external :: arg !< A variable, whatever kind

  call wrapper_add(ctx, f_loc(arg))
end subroutine f_object_signal_add_arg

!> @relates module_f_objects
subroutine f_object_kernel_new(ctx, hook, n_args)
  use f_precisions
  use module_f_objects, only: wrapper_new => f_object_kernel_new_, kernel_ctx
  type(kernel_ctx), intent(out) :: ctx
  integer, intent(in) :: n_args
  external :: hook

  ctx = wrapper_new(f_loc(hook), n_args)
end subroutine f_object_kernel_new

!> @relates module_f_objects
subroutine f_object_kernel_add_arg(ctx, arg)
  use f_precisions
  use module_f_objects, only: wrapper_add => f_object_kernel_add_arg_, kernel_ctx
  type(kernel_ctx), intent(inout) :: ctx !< kernel context.
  external :: arg !< A variable, whatever kind

  call wrapper_add(ctx, f_loc(arg))
end subroutine f_object_kernel_add_arg

subroutine f_object_signal_connect_bind(obj, sig, ctx, id)
  use module_f_objects, only: f_object_signal_connect, kernel_ctx
  character(len = *), intent(in) :: obj, sig
  type(kernel_ctx), intent(out) :: ctx
  integer, intent(out) :: id

  call f_object_signal_connect(obj, sig, ctx, id)
end subroutine f_object_signal_connect_bind
