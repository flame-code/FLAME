subroutine bind_yaml_swap_stream( &
    new_unit, &
    old_unit, &
    ierr)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: new_unit
  integer, intent(out) :: old_unit
  integer, intent(out) :: ierr

  call yaml_swap_stream( &
    new_unit, &
    old_unit, &
    ierr)
end subroutine bind_yaml_swap_stream

subroutine bind_yaml_output_errors( &
    )
  use yaml_strings
  use dictionaries
  use exception_callbacks, only: f_err_set_all_errors_callback => f_err_set_all_errors_callback, &
    f_err_set_last_error_callback => f_err_set_last_error_callback
  use f_precisions
  use yaml_output
  implicit none

  call yaml_output_errors( &
)
end subroutine bind_yaml_output_errors

subroutine bind_yaml_set_default_stream( &
    unit, &
    ierr)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: unit
  integer, intent(out) :: ierr

  call yaml_set_default_stream( &
    unit, &
    ierr)
end subroutine bind_yaml_set_default_stream

subroutine bind_yaml_get_default_stream( &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(out) :: unit

  call yaml_get_default_stream( &
    unit)
end subroutine bind_yaml_get_default_stream

subroutine bind_yaml_stream_connected( &
    filename_len, &
    filename, &
    unit, &
    istat)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: filename_len
  character(len = filename_len), intent(in) :: filename
  integer, intent(out) :: unit
  integer, optional, intent(out) :: istat

  call yaml_stream_connected( &
    filename, &
    unit, &
    istat)
end subroutine bind_yaml_stream_connected

subroutine bind_yaml_set_stream( &
    unit, &
    filename_len, &
    filename, &
    istat, &
    tabbing, &
    record_length, &
    position_len, &
    position, &
    setdefault)
  use dictionaries
  use yaml_strings
  use f_precisions
  use f_utils, only: f_utils_recl => f_utils_recl, &
    f_get_free_unit => f_get_free_unit, &
    f_open_file => f_open_file
  use yaml_output
  implicit none
  integer, optional, intent(in) :: unit
  integer, intent(in) :: filename_len
  character(len = filename_len), optional, intent(in) :: filename
  integer, optional, intent(out) :: istat
  integer, optional, intent(in) :: tabbing
  integer, optional, intent(in) :: record_length
  integer, intent(in) :: position_len
  character(len = position_len), optional, intent(in) :: position
  logical, optional, intent(in) :: setdefault

  call yaml_set_stream( &
    unit, &
    filename, &
    istat, &
    tabbing, &
    record_length, &
    position, &
    setdefault)
end subroutine bind_yaml_set_stream

subroutine bind_yaml_stream_attributes( &
    unit, &
    stream_unit, &
    icursor, &
    flowrite, &
    itab_active, &
    iflowlevel, &
    ilevel, &
    ilast, &
    indent, &
    indent_previous, &
    record_length)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, optional, intent(in) :: unit
  integer, optional, intent(in) :: stream_unit
  integer, optional, intent(out) :: icursor
  logical, optional, intent(out) :: flowrite
  integer, optional, intent(out) :: itab_active
  integer, optional, intent(out) :: iflowlevel
  integer, optional, intent(out) :: ilevel
  integer, optional, intent(out) :: ilast
  integer, optional, intent(out) :: indent
  integer, optional, intent(out) :: indent_previous
  integer, optional, intent(out) :: record_length

  call yaml_stream_attributes( &
    unit, &
    stream_unit, &
    icursor, &
    flowrite, &
    itab_active, &
    iflowlevel, &
    ilevel, &
    ilast, &
    indent, &
    indent_previous, &
    record_length)
end subroutine bind_yaml_stream_attributes

subroutine bind_yaml_new_document( &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, optional, intent(in) :: unit

  call yaml_new_document( &
    unit)
end subroutine bind_yaml_new_document

subroutine bind_yaml_flush_document( &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, optional, intent(in) :: unit

  call yaml_flush_document( &
    unit)
end subroutine bind_yaml_flush_document

subroutine bind_yaml_close_stream( &
    unit, &
    istat)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, optional, intent(in) :: unit
  integer, optional, intent(out) :: istat

  call yaml_close_stream( &
    unit, &
    istat)
end subroutine bind_yaml_close_stream

subroutine bind_yaml_close_all_streams( &
    )
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none

  call yaml_close_all_streams( &
)
end subroutine bind_yaml_close_all_streams

subroutine bind_yaml_dict_inspect( &
    dict)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  type(dictionary), pointer :: dict

  call yaml_dict_inspect( &
    dict)
end subroutine bind_yaml_dict_inspect

subroutine bind_dump_progress_bar( &
    bar, &
    step, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use f_utils
  use yaml_output
  implicit none
  type(f_progress_bar), intent(inout) :: bar
  integer, optional, intent(in) :: step
  integer, optional, intent(in) :: unit

  call dump_progress_bar( &
    bar, &
    step, &
    unit)
end subroutine bind_dump_progress_bar

subroutine bind_yaml_cite( &
    paper_len, &
    paper, &
    unit)
  use f_bibliography
  use dictionaries
  use yaml_strings
  use f_precisions
  use f_utils
  use yaml_output
  implicit none
  integer, intent(in) :: paper_len
  character(len = paper_len), intent(in) :: paper
  integer, optional, intent(in) :: unit

  call yaml_cite( &
    paper, &
    unit)
end subroutine bind_yaml_cite

subroutine bind_yaml_bib_dump( &
    citations, &
    unit)
  use f_bibliography
  use dictionaries
  use yaml_strings
  use f_precisions
  use f_utils
  use yaml_output
  implicit none
  type(dictionary), pointer :: citations
  integer, optional, intent(in) :: unit

  call yaml_bib_dump( &
    citations, &
    unit)
end subroutine bind_yaml_bib_dump

subroutine bind_yaml_scalar( &
    message_len, &
    message, &
    advance_len, &
    advance, &
    unit, &
    hfill_len, &
    hfill)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: message_len
  character(len = message_len), intent(in) :: message
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: hfill_len
  character(len = hfill_len), optional, intent(in) :: hfill

  call yaml_scalar( &
    message, &
    advance, &
    unit, &
    hfill)
end subroutine bind_yaml_scalar

subroutine bind_yaml_mapping_open( &
    mapname_len, &
    mapname, &
    label_len, &
    label, &
    tag_len, &
    tag, &
    flow, &
    tabbing, &
    advance_len, &
    advance, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), optional, intent(in) :: mapname
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: tag_len
  character(len = tag_len), optional, intent(in) :: tag
  logical, optional, intent(in) :: flow
  integer, optional, intent(in) :: tabbing
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit

  call yaml_mapping_open( &
    mapname, &
    label, &
    tag, &
    flow, &
    tabbing, &
    advance, &
    unit)
end subroutine bind_yaml_mapping_open

subroutine bind_yaml_mapping_close( &
    advance_len, &
    advance, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit

  call yaml_mapping_close( &
    advance, &
    unit)
end subroutine bind_yaml_mapping_close

subroutine bind_yaml_sequence_open( &
    mapname_len, &
    mapname, &
    label_len, &
    label, &
    tag_len, &
    tag, &
    flow, &
    tabbing, &
    advance_len, &
    advance, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), optional, intent(in) :: mapname
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: tag_len
  character(len = tag_len), optional, intent(in) :: tag
  logical, optional, intent(in) :: flow
  integer, optional, intent(in) :: tabbing
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit

  call yaml_sequence_open( &
    mapname, &
    label, &
    tag, &
    flow, &
    tabbing, &
    advance, &
    unit)
end subroutine bind_yaml_sequence_open

subroutine bind_yaml_sequence_close( &
    advance_len, &
    advance, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit

  call yaml_sequence_close( &
    advance, &
    unit)
end subroutine bind_yaml_sequence_close

subroutine bind_yaml_newline( &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, optional, intent(in) :: unit

  call yaml_newline( &
    unit)
end subroutine bind_yaml_newline

subroutine bind_yaml_sequence( &
    seqvalue_len, &
    seqvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    padding)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: seqvalue_len
  character(len = seqvalue_len), optional, intent(in) :: seqvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, optional, intent(in) :: padding

  call yaml_sequence( &
    seqvalue, &
    label, &
    advance, &
    unit, &
    padding)
end subroutine bind_yaml_sequence

subroutine bind_yaml_dict_dump( &
    dict, &
    unit, &
    flow, &
    verbatim)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  type(dictionary), pointer, intent(in) :: dict
  integer, optional, intent(in) :: unit
  logical, optional, intent(in) :: flow
  logical, optional, intent(in) :: verbatim

  call yaml_dict_dump( &
    dict, &
    unit, &
    flow, &
    verbatim)
end subroutine bind_yaml_dict_dump

subroutine bind_yaml_dict_dump_all( &
    dict, &
    unit, &
    flow, &
    verbatim)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  type(dictionary), pointer, intent(in) :: dict
  integer, optional, intent(in) :: unit
  logical, optional, intent(in) :: flow
  logical, optional, intent(in) :: verbatim

  call yaml_dict_dump_all( &
    dict, &
    unit, &
    flow, &
    verbatim)
end subroutine bind_yaml_dict_dump_all

subroutine bind_yaml_release_document( &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, optional, intent(in) :: unit

  call yaml_release_document( &
    unit)
end subroutine bind_yaml_release_document

subroutine bind_yaml_map( &
    mapname_len, &
    mapname, &
    mapvalue_len, &
    mapvalue, &
    label_len, &
    label, &
    tag_len, &
    tag, &
    advance_len, &
    advance, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_len
  character(len = mapvalue_len), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: tag_len
  character(len = tag_len), optional, intent(in) :: tag
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    tag, &
    advance, &
    unit)
end subroutine bind_yaml_map

subroutine bind_yaml_map_dict( &
    mapname_len, &
    mapname, &
    mapvalue, &
    label_len, &
    label, &
    unit, &
    flow)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  type(dictionary), pointer, intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, optional, intent(in) :: unit
  logical, optional, intent(in) :: flow

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    unit, &
    flow)
end subroutine bind_yaml_map_dict

subroutine bind_yaml_map_enum( &
    mapname_len, &
    mapname, &
    mapvalue, &
    label_len, &
    label, &
    unit, &
    flow)
  use f_enums
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  type(f_enumerator), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, optional, intent(in) :: unit
  logical, optional, intent(in) :: flow

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    unit, &
    flow)
end subroutine bind_yaml_map_enum

subroutine bind_yaml_map_li( &
    mapname_len, &
    mapname, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer(kind = f_long), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_li

subroutine bind_yaml_map_i( &
    mapname_len, &
    mapname, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer(kind = f_integer), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_i

subroutine bind_yaml_map_f( &
    mapname_len, &
    mapname, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  real(kind = f_simple), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_f

subroutine bind_yaml_map_d( &
    mapname_len, &
    mapname, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  real(kind = f_double), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_d

subroutine bind_yaml_map_l( &
    mapname_len, &
    mapname, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  logical, intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_l

subroutine bind_yaml_map_dv( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  real(kind = f_double), dimension(mapvalue_dim_0), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_dv

subroutine bind_yaml_map_rv( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  real(kind = f_simple), dimension(mapvalue_dim_0), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_rv

subroutine bind_yaml_map_cv( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue_len, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer, intent(in) :: mapvalue_len
  character(len = mapvalue_len), dimension(mapvalue_dim_0), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_cv

subroutine bind_yaml_map_lv( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  logical, dimension(mapvalue_dim_0), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_lv

subroutine bind_yaml_map_liv( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer(kind = f_long), dimension(mapvalue_dim_0), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_liv

subroutine bind_yaml_map_iv( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer(kind = f_integer), dimension(mapvalue_dim_0), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_iv

subroutine bind_yaml_map_dm( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue_dim_1, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer, intent(in) :: mapvalue_dim_1
  real(kind = f_double), dimension(mapvalue_dim_0, mapvalue_dim_1), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_dm

subroutine bind_yaml_map_rm( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue_dim_1, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer, intent(in) :: mapvalue_dim_1
  real(kind = f_simple), dimension(mapvalue_dim_0, mapvalue_dim_1), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_rm

subroutine bind_yaml_map_im( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue_dim_1, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer, intent(in) :: mapvalue_dim_1
  integer(kind = f_integer), dimension(mapvalue_dim_0, mapvalue_dim_1), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_im

subroutine bind_yaml_map_lm( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue_dim_1, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer, intent(in) :: mapvalue_dim_1
  logical, dimension(mapvalue_dim_0, mapvalue_dim_1), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_lm

subroutine bind_yaml_map_dt( &
    mapname_len, &
    mapname, &
    mapvalue_dim_0, &
    mapvalue_dim_1, &
    mapvalue_dim_2, &
    mapvalue, &
    label_len, &
    label, &
    advance_len, &
    advance, &
    unit, &
    fmt_len, &
    fmt)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: mapname_len
  character(len = mapname_len), intent(in) :: mapname
  integer, intent(in) :: mapvalue_dim_0
  integer, intent(in) :: mapvalue_dim_1
  integer, intent(in) :: mapvalue_dim_2
  real(kind = f_double), dimension(mapvalue_dim_0, mapvalue_dim_1, mapvalue_dim_2), intent(in) :: mapvalue
  integer, intent(in) :: label_len
  character(len = label_len), optional, intent(in) :: label
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: fmt_len
  character(len = fmt_len), optional, intent(in) :: fmt

  call yaml_map( &
    mapname, &
    mapvalue, &
    label, &
    advance, &
    unit, &
    fmt)
end subroutine bind_yaml_map_dt

subroutine bind_yaml_warning_str( &
    message, &
    level, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  type(f_string), intent(in) :: message
  integer, optional, intent(in) :: level
  integer, optional, intent(in) :: unit

  call yaml_warning( &
    message, &
    level, &
    unit)
end subroutine bind_yaml_warning_str

subroutine bind_yaml_warning_c( &
    message_len, &
    message, &
    level, &
    unit)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: message_len
  character(len = message_len), intent(in) :: message
  integer, optional, intent(in) :: level
  integer, optional, intent(in) :: unit

  call yaml_warning( &
    message, &
    level, &
    unit)
end subroutine bind_yaml_warning_c

subroutine bind_yaml_comment_str( &
    message, &
    advance_len, &
    advance, &
    unit, &
    hfill_len, &
    hfill, &
    tabbing)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  type(f_string), intent(in) :: message
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: hfill_len
  character(len = hfill_len), optional, intent(in) :: hfill
  integer, optional, intent(in) :: tabbing

  call yaml_comment( &
    message, &
    advance, &
    unit, &
    hfill, &
    tabbing)
end subroutine bind_yaml_comment_str

subroutine bind_yaml_comment_c( &
    message_len, &
    message, &
    advance_len, &
    advance, &
    unit, &
    hfill_len, &
    hfill, &
    tabbing)
  use dictionaries
  use yaml_strings
  use f_precisions
  use yaml_output
  implicit none
  integer, intent(in) :: message_len
  character(len = message_len), intent(in) :: message
  integer, intent(in) :: advance_len
  character(len = advance_len), optional, intent(in) :: advance
  integer, optional, intent(in) :: unit
  integer, intent(in) :: hfill_len
  character(len = hfill_len), optional, intent(in) :: hfill
  integer, optional, intent(in) :: tabbing

  call yaml_comment( &
    message, &
    advance, &
    unit, &
    hfill, &
    tabbing)
end subroutine bind_yaml_comment_c

