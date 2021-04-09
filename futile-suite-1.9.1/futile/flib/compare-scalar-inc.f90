    character(len = *), intent(in), optional :: label

    if (val /= expected) then
       if (present(label)) then
          call f_err_throw(label // ": expected " // yaml_toa(expected) // ", value " // yaml_toa(val))
       else
          call f_err_throw("expected " // yaml_toa(expected) // ", value " // yaml_toa(val))
       end if
    end if
