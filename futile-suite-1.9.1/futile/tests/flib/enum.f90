program test_enum
  use yaml_output
  use f_regtests
  implicit none

  call f_lib_initialize()

  call yaml_sequence_open("Tests")
  call run(creator)
  call run(equality)
  call run(difference)
  call run(convert_i)
  call run(convert_s)
  call run(attributes)
  call yaml_sequence_close()
  
  call f_lib_finalize()

contains

  subroutine creator(label)
    use f_enums
    implicit none
    character(len = *), intent(out) :: label

    type(f_enumerator) :: e

    label = "Creator"
    call nullify_f_enum(e)
    call verify(e == f_enumerator_null(), "equals to null enum")
  end subroutine creator

  subroutine equality(label)
    use f_enums
    implicit none
    character(len = *), intent(out) :: label

    type(f_enumerator) :: e, o, d, v
    integer, parameter :: ival = 123
    integer, parameter :: oval = 999
    character(len = *), parameter :: name = 'TEST_VALUE'
    character(len = *), parameter :: other = 'OTHER_VALUE'

    label = "Equality"

    e = f_enumerator(name, ival, null())
    o = f_enumerator(name, ival, null())
    d = f_enumerator(name, oval, null())
    v = f_enumerator(other, ival, null())

    call verify(e == ival, "equality with an integer")
    call verify(ival == e, "equality with an integer")
    call verify(.not. e == oval, "inexact equality with an integer")
    call verify(.not. oval == e, "inexact equality with an integer")
    call verify(e == name, "equality with a string")
    call verify(.not. e == other, "inexact equality with a string")
    call verify(e == e, "equality with itself")
    call verify(e == o, "equality with an other similar enum")
    call verify(.not. e == d, "inexact equality with a different enum")
    call verify(.not. e == v, "inexact equality with a different enum")
  end subroutine equality

  subroutine difference(label)
    use f_enums
    implicit none
    character(len = *), intent(out) :: label

    type(f_enumerator) :: e, o, d, v
    integer, parameter :: ival = 123
    integer, parameter :: oval = 999
    character(len = *), parameter :: name = 'TEST_VALUE'
    character(len = *), parameter :: other = 'OTHER_VALUE'

    label = "Inequality"

    e = f_enumerator(name, ival, null())
    o = f_enumerator(name, ival, null())
    d = f_enumerator(name, oval, null())
    v = f_enumerator(other, ival, null())

    call verify(e /= oval, "inequality with an integer")
    call verify(.not. e /= ival, "inexact inequality with an integer")
    call verify(e /= other, "inequality with a string")
    call verify(.not. e /= name, "inexact inequality with a string")
    call verify(e /= d, "inequality with a different enum")
    call verify(e /= v, "inequality with an other different enum")
    call verify(.not. e /= e, "inexact inequality with itself")
    call verify(.not. e /= o, "inexact inequality with a similar enum")
  end subroutine difference

  subroutine convert_i(label)
    use f_enums
    implicit none
    character(len = *), intent(out) :: label

    type(f_enumerator) :: e
    integer, parameter :: ival = 123
    character(len = *), parameter :: name = 'TEST_VALUE'

    label = "Conversion to an integer"

    e = f_enumerator(name, ival, null())
    call compare(toi(e), ival)
  end subroutine convert_i
  
  subroutine convert_s(label)
    use f_enums
    implicit none
    character(len = *), intent(out) :: label

    type(f_enumerator) :: e
    integer, parameter :: ival = 123
    character(len = *), parameter :: name = 'TEST_VALUE'

    label = "Conversion to a string"

    e = f_enumerator(name, ival, null())
    call compare(toa(e), name)
  end subroutine convert_s
  
  subroutine attributes(label)
    use f_enums
    implicit none
    character(len = *), intent(out) :: label

    type(f_enumerator) :: e
    type(f_enumerator), target :: family, o
    integer, parameter :: ival = 123
    integer, parameter :: fval = 999
    character(len = *), parameter :: name = 'TEST_VALUE'
    character(len = *), parameter :: surname = 'GENERIC_VALUE'

    label = "Attributes"

    e = f_enumerator(name, ival, null())
    o = f_enumerator(name, ival, null())
    family = f_enumerator(surname, fval, null())
    call verify(.not. (e .hasattr. family), "no family")
    call f_enum_attr(e, family)
    call verify(e .hasattr. family, "with a family")
    call verify(.not. (e .hasattr. o), "not with another family")
    call verify(e .hasattr. fval, "with a family number")
    call verify(.not. (e .hasattr. ival), "not with another family number")
    call verify(e .hasattr. surname, "with a family name")
    call verify(.not. (e .hasattr. name), "not with another family name")
    call f_enum_attr(e, o)
    call verify(e .hasattr. o, "with a second family")
    call verify(family .hasattr. o, "second family by transitivity")
  end subroutine attributes
  
end program test_enum
