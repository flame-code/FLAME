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

