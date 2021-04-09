  c_allocation=.false.
  if(len_trim(m%info) > 0) then
     info => yaml_load(m%info)
     val=' '
     val=info .get. INFO_TYPE_KEY
     call dict_free(info)
     select case(trim(val))
     case(INFO_SHARED_TYPE)
        c_allocation=.true.
     end select
  end if
  if (c_allocation) then
     p = smpi_shared_malloc(product(m%shape(1:(m%rank-1)))*(m%shape(m%rank)+padding)*f_sizeof(d),&
          trim(m%array_id)//char(0), int(m%rank,f_long))
     call c_f_pointer(p, array, int(m%shape(1:m%rank)))
     !call f_map_ptr(m%lbounds,m%ubounds,array_tmp,array)
     !here we might add a fortran call to the function for the reshaping of the pointer
     ierror=0
     include 'allocate-inc.f90'
     return
  end if
