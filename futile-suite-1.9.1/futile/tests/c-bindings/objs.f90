module my_objects

  type my_object
     character(len = 10) :: label
     integer, dimension(:), pointer :: data
     double precision :: mean
  end type my_object

contains

  subroutine my_object_alloc(pobj, add)
    use f_precisions
    type(my_object), pointer :: pobj
    integer(f_address), intent(out) :: add

    allocate(pobj)
    call my_object_nullify(pobj)
    add = f_loc(pobj)
  end subroutine my_object_alloc

  subroutine my_object_dealloc(pobj)
    type(my_object), pointer :: pobj

    call my_object_free(pobj)
    deallocate(pobj)
  end subroutine my_object_dealloc

  subroutine my_object_nullify(obj)
    type(my_object), intent(out) :: obj
    write(obj%label, "(A)") ""
    nullify(obj%data)
    obj%mean = 0.
  end subroutine my_object_nullify

  subroutine my_object_set_data(obj, label, data, ln)
    use dynamic_memory
    type(my_object), intent(inout) :: obj
    character(len = *), intent(in) :: label
    integer, intent(in) :: ln
    integer, dimension(ln), intent(in) :: data

    if (associated(obj%data)) call f_free_ptr(obj%data)
    if (ln > 0) then
       obj%data = f_malloc_ptr(size(data), id = "data")
       obj%data = data
       obj%mean = real(sum(obj%data)) / real(ln)
    else
       call my_object_nullify(obj)
    end if
    write(obj%label, "(A)") label
  end subroutine my_object_set_data

  subroutine my_object_get_data(obj, data)
    use dynamic_memory
    use f_python
    type(my_object), intent(in) :: obj
    type(ndarray), intent(out) :: data

    data = toNdArray_ptr(obj%data)
  end subroutine my_object_get_data

  subroutine my_object_get_size(obj, ln)
    type(my_object), intent(in) :: obj
    integer, intent(out) :: ln

    ln = -1
    if (associated(obj%data)) ln = size(obj%data)
  end subroutine my_object_get_size

  subroutine my_object_get_mean(obj, mean)
    type(my_object), intent(in) :: obj
    double precision, intent(out) :: mean

    mean = obj%mean
  end subroutine my_object_get_mean

  subroutine my_object_serialize(obj)
    use yaml_output
    use yaml_strings
    type(my_object), intent(in) :: obj

    integer :: i

    call yaml_sequence_open(trim(obj%label))
    do i = 1, size(obj%data), 1
       call yaml_sequence(advance = "NO")
       call yaml_scalar(yaml_toa(obj%data(i)))
    end do
    call yaml_sequence_close()
  end subroutine my_object_serialize

  subroutine my_object_free(obj)
    use dynamic_memory
    type(my_object), intent(inout) :: obj
    if (associated(obj%data)) call f_free_ptr(obj%data)
    call my_object_nullify(obj)
  end subroutine my_object_free

end module my_objects


program test
  use my_objects
  use f_precisions
  
  type(my_object) :: obj
  integer :: ierr

  call f_lib_initialize()

  call f_object_new("my_object", my_object_alloc, my_object_dealloc)
  call f_object_add_method("my_object", "set_data", my_object_set_data, 3)
  call f_object_add_method("my_object", "serialize", my_object_serialize, 0)
  call f_object_add_method("my_object", "get_size", my_object_get_size, 1)
  call f_object_add_method("my_object", "get_mean", my_object_get_mean, 1)
  call f_object_add_method("my_object", "get_data", my_object_get_data, 1)

  call f_object_add_method("class", "version", version, 0)

  call my_object_nullify(obj)
  call my_object_set_data(obj, "fortran", (/ 1, 2, 3 /), 3)

  call f_python_initialize(0, 1, 0, 1)

  call f_python_execute('print " nproc: %d" % futile.nproc', ierr)
  call f_python_execute("futile.version()", ierr)

  call f_python_execute("import numpy", ierr)

  call f_python_add_object("my_object", "obj", obj)
  !call f_python_execute('obj = futile.FObject("my_object", %ld)' % f_loc(obj))

  call f_python_execute("obj.serialize()", ierr)
  call f_python_execute('obj.set_data("python", (4,5,6,7), 4)', ierr)
  call f_python_execute("obj.serialize()", ierr)
  call f_python_execute("print ' get_size: %d' % obj.get_size(futile.SCALAR_I4)", ierr)
  call f_python_execute("print ' get_mean: %g' % obj.get_mean(0.)", ierr)
  call f_python_execute("(data,) = obj.get_data(futile.ARRAY)", ierr)
  call f_python_execute("print ' get_data: %s' % data", ierr)

  call f_python_execute("print ' modifying data:',; data[2] = 42; print 'yes'", ierr)
  call f_python_execute("obj.serialize()", ierr)

  call f_python_execute('obj.set_data("numpy", numpy.array((123, 456), dtype = numpy.int32), 2)', ierr)
  call f_python_execute("obj.serialize()", ierr)

  call f_python_execute('obj2 = futile.FObject("my_object")', ierr)
  call f_python_execute('obj2.set_data("python new", (42, ), 1)', ierr)
  call f_python_execute("obj2.serialize()", ierr)

  call f_python_finalize()

  call my_object_free(obj)

  call f_lib_finalize()

contains

  subroutine version()
    use f_lib_package
    use yaml_output
    
    call yaml_map("Futile version", package_version)
  end subroutine version
end program test
