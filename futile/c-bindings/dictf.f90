subroutine bind_dict_move_to_key(dict, exists, key)
  use dictionaries, only: dictionary, operator(//), has_key
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists
  character(len = *), intent(in) :: key

  exists = has_key(dict, key(1:len(key)))
  if (exists) dict => dict // key(1:len(key))
END SUBROUTINE bind_dict_move_to_key


subroutine bind_dict_move_to_item(dict, exists, id)
  use dictionaries, only: dictionary, operator(//), dict_len
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists
  integer, intent(in) :: id

  exists = (id < dict_len(dict) .and. id >= 0)
  if (exists) dict => dict // id
END SUBROUTINE bind_dict_move_to_item


subroutine bind_dict_insert(dict, key)
  use dictionaries, only: dictionary, operator(//)
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: key

  ! This is a patch for Intel, to be corrected properly later.
  dict => dict // key(1:len(key))
END SUBROUTINE bind_dict_insert


subroutine bind_dict_append(dict)
  use dictionaries, only: dictionary, operator(//), dict_len
  implicit none
  type(dictionary), pointer :: dict

  dict => dict // dict_len(dict)
END SUBROUTINE bind_dict_append


subroutine bind_dict_set(dict, val)
  use dictionaries, only: dictionary, set
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: val

  ! This is a patch for Intel, to be corrected properly later.
  call set(dict, val(1:len(val)))
END SUBROUTINE bind_dict_set


subroutine bind_dict_dump(dict, unit)
  use dictionaries, only: dictionary
  use yaml_output, only: yaml_dict_dump
  implicit none
  type(dictionary), pointer :: dict
  integer, intent(in) :: unit

  if (unit < 0) then
     call yaml_dict_dump(dict)
  else
     call yaml_dict_dump(dict, unit = unit)
  end if
END SUBROUTINE bind_dict_dump


subroutine bind_dict_dump_to_file(dict, path)
  use dictionaries, only: dictionary
  use yaml_output, only: yaml_dict_dump, yaml_set_stream, yaml_close_stream, yaml_get_default_stream
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: path

  integer :: unit

  call yaml_set_stream(filename = trim(path), &
       & record_length = 92, setdefault = .true., tabbing = 0)
  call yaml_get_default_stream(unit = unit)
  call yaml_dict_dump(dict, unit = unit)
  call yaml_close_stream(unit = unit)
END SUBROUTINE bind_dict_dump_to_file


subroutine bind_dict_parse(dict, buf)
  use dictionaries, only: dictionary, operator(//), dict_len,operator(.pop.),dict_free
  use yaml_parse, only: yaml_parse_from_string
  implicit none
  type(dictionary), pointer :: dict
  character(len = *), intent(in) :: buf
  type(dictionary), pointer :: dict_load

  nullify(dict_load)
  call yaml_parse_from_string(dict_load, buf)
  if (dict_len(dict_load) == 1) then
     dict => dict_load .pop. 0
  end if
  call dict_free(dict_load)
END SUBROUTINE bind_dict_parse


subroutine bind_dict_pop(dict, exists, key)
  use dictionaries, only: dictionary, has_key, dict_remove, dict_init
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists
  character(len = *), intent(in) :: key

  exists = (has_key(dict, key(1:len(key))))
  if (exists) then
     call dict_remove(dict, key(1:len(key)))
     if (.not. associated(dict)) call dict_init(dict)
  end if
END SUBROUTINE bind_dict_pop


subroutine bind_dict_value(dict, buf)
  use dictionaries, only: dictionary, max_field_length, dict_value
  implicit none
  type(dictionary), pointer :: dict
  character(len = max_field_length), intent(out) :: buf
  
  buf = dict_value(dict)
END SUBROUTINE bind_dict_value


subroutine bind_dict_key(dict, buf)
  use dictionaries, only: dictionary, max_field_length, dict_key
  implicit none
  type(dictionary), pointer :: dict
  character(len = max_field_length), intent(out) :: buf
  
  buf = dict_key(dict)
END SUBROUTINE bind_dict_key


subroutine bind_dict_iter(dict, exists)
  use dictionaries, only: dictionary, dict_iter
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists

  type(dictionary), pointer :: start

  start => dict_iter(dict)
  exists = associated(start)
  if (exists) dict => start
END SUBROUTINE bind_dict_iter


subroutine bind_dict_next(dict, exists)
  use dictionaries, only: dictionary, dict_next
  implicit none
  type(dictionary), pointer :: dict
  logical, intent(out) :: exists

  type(dictionary), pointer :: next

  next => dict_next(dict)
  exists = associated(next)
  if (exists) dict => next
END SUBROUTINE bind_dict_next


subroutine bind_dict_len(dict, ln)
  use dictionaries, only: dictionary, dict_len
  implicit none
  type(dictionary), pointer :: dict
  integer, intent(out) :: ln

  ln = dict_len(dict)
END SUBROUTINE bind_dict_len


subroutine bind_dict_size(dict, ln)
  use dictionaries, only: dictionary, dict_size
  implicit none
  type(dictionary), pointer :: dict
  integer, intent(out) :: ln

  ln = dict_size(dict)
END SUBROUTINE bind_dict_size


subroutine bind_dict_copy(dict, ref)
  use dictionaries, only: dictionary, dict_copy
  implicit none
  type(dictionary), pointer :: dict, ref

  call dict_copy(dict, ref)
END SUBROUTINE bind_dict_copy


subroutine bind_dict_update(dict, ref)
  use dictionaries, only: dictionary, dict_update
  implicit none
  type(dictionary), pointer :: dict, ref

  call dict_update(dict, ref)
END SUBROUTINE bind_dict_update


subroutine bind_dict_init(dict)
  use dictionaries, only: dictionary, dict_init
  implicit none
  type(dictionary), pointer :: dict

  call dict_init(dict)
END SUBROUTINE bind_dict_init

subroutine bind_dict_free(dict)
  use dictionaries_base, only: dictionary, dict_free
  implicit none
  type(dictionary), pointer :: dict

  call dict_free(dict)
end subroutine bind_dict_free
