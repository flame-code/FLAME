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

subroutine bind_dict_set_double(dict,key,keylen,val)
  use f_precisions
  use dictionaries
  type(dictionary), pointer :: dict
  integer, intent(in) :: keylen
  character(len=keylen), intent(in) :: key
  real(f_double), intent(in) :: val
  call set(dict // key, val)
end subroutine bind_dict_set_double

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

subroutine bind_dict_set_double_array(dict,key,keylen,array,arraylen)
  use f_precisions
  use dictionaries
  type(dictionary), pointer :: dict
  integer, intent(in) :: keylen,arraylen
  character(len=keylen), intent(in) :: key
  real(f_double), dimension(arraylen), intent(in) :: array

  call set(dict // key, array)
end subroutine bind_dict_set_double_array

subroutine bind_dict_set_double_matrix(dict,key,keylen,matrix,lx,ly)
  use f_precisions
  use dictionaries
  type(dictionary), pointer :: dict
  integer, intent(in) :: keylen,lx,ly
  character(len=keylen), intent(in) :: key
  real(f_double), dimension(lx,ly), intent(in) :: matrix
  !local variables
  integer :: i
  do i=1,ly
     call set(dict // key // (i-1), matrix(:,i))
  end do
end subroutine bind_dict_set_double_matrix

subroutine bind_dict_set_string(dict,key,val)
  use dictionaries
  use yaml_strings, only: convert_f_char_ptr
  type(dictionary), pointer :: dict
  character, dimension(*), intent(in) :: key,val
  !local variables
  character(len=max_field_length) :: key_,val_
  call convert_f_char_ptr(src=key,dest=key_)
  call convert_f_char_ptr(src=val,dest=val_)
  call set(dict // trim(key_), trim(val_))
end subroutine bind_dict_set_string

subroutine bind_dict_set_dict(dict,key,keylen,val)
  use dictionaries
  type(dictionary), pointer :: dict,val
  integer, intent(in) :: keylen
  character(len=keylen), intent(in) :: key

  call set(dict // key, val)
end subroutine bind_dict_set_dict

subroutine bind_dict_add_dict(dict,val)
  use dictionaries
  type(dictionary), pointer :: dict,val

  call add(dict, val)
end subroutine bind_dict_add_dict

subroutine bind_dict_add_double(dict,val)
  use dictionaries
  use f_precisions
  type(dictionary), pointer :: dict
  real(f_double), intent(in) :: val

  call add(dict, val)
end subroutine bind_dict_add_double

subroutine bind_dict_add_int(dict,val)
  use dictionaries
  use f_precisions
  type(dictionary), pointer :: dict
  integer(f_integer), intent(in) :: val

  call add(dict, val)
end subroutine bind_dict_add_int

subroutine bind_dict_add_char(dict,val)
  use dictionaries
  use f_precisions
  use yaml_strings
  type(dictionary), pointer :: dict
  character, dimension(*), intent(in) :: val

  character(len = max_field_length) :: buf
  
  call convert_f_char_ptr(val, buf)
  call add(dict, buf)
end subroutine bind_dict_add_char


subroutine bind_dict_get_double_array(dict,key,keylen,array,arraylen,istat)
  use f_precisions
  use dictionaries
  type(dictionary), pointer :: dict
  integer, intent(in) :: keylen,arraylen
  character(len=keylen), intent(in) :: key
  integer, intent(out) :: istat
  real(f_double), dimension(arraylen), intent(out) :: array
  
  istat=1
  if (key .in. dict) then
     istat=0
     array=dict//key
  end if

end subroutine bind_dict_get_double_array

subroutine bind_dict_get_dict(dict,key,keylen,val,istat)
  use f_precisions
  use dictionaries
  type(dictionary), pointer :: dict,val
  integer, intent(in) :: keylen
  character(len=keylen), intent(in) :: key
  integer, intent(out) :: istat
  
  istat=1
  if (key .in. dict) then
     istat=0
     val=>dict//key
  end if

end subroutine bind_dict_get_dict

subroutine bind_iter_null(iter)
  use dictionaries
  type(dictionary), pointer :: iter
  nullify(iter)
end subroutine bind_iter_null

subroutine bind_iterate(iter,dict,istat)
  use dictionaries
  type(dictionary), pointer :: iter,dict
  istat=1
  if (associated(iter)) then
     iter=>dict_next(iter)
  else
     iter => dict_iter(dict)
  end if
  if(associated(iter)) istat=0

  !if (iterate(iter,on=dict)) istat=0
end subroutine bind_iterate
