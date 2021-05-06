    character(len = *), intent(in), optional :: label

    if (any(val /= expected)) then
       if (present(label)) then
          call f_err_throw(label // ": maximum difference of " // yaml_toa(maxval(val - expected)))
       else
          call f_err_throw("maximum difference of " // yaml_toa(maxval(val - expected)))
       end if
    end if
