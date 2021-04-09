  alignment=0
  if(len_trim(m%info) > 0) then
    info => yaml_load(m%info)
    alignment = dict_get(info,INFO_ALIGNMENT_KEY,0)
    call dict_free(info)
  end if
  if (alignment /=0) then
   !  print *, "size", product(m%shape(1:(m%rank-1)))*(m%shape(m%rank)+padding)*f_sizeof(d),&
   !          "align", alignment, "div",product(m%shape(1:(m%rank-1)))*(m%shape(m%rank)+padding)*f_sizeof(d)/alignment
    alignment = f_sizeof(d)
    aligned_size = product(m%shape(1:(m%rank-1)))* &
                          (m%shape(m%rank)+padding)*f_sizeof(d)
    ! Pad the size being allocated to be a multiple of the alignment
    aligned_size = CEILING((aligned_size-1.0)/alignment) * alignment
    p= aligned_alloc(alignment,aligned_size)
    call c_f_pointer(p, array, int(m%shape(1:m%rank)))
    ierror=0
    include 'allocate-inc.f90'
    return 
  end if
