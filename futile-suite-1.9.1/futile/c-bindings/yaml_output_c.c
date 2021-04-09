#include "yaml_output.h"
#include <config.h>
#include <string.h>

void FC_FUNC_(bind_yaml_swap_stream, BIND_YAML_SWAP_STREAM)(const int*,
  int*,
  int*);
void yaml_swap_stream(int new_unit,
  int old_unit,
  int ierr)
{
  FC_FUNC_(bind_yaml_swap_stream, BIND_YAML_SWAP_STREAM)
    (&new_unit, &old_unit, &ierr);
}

void FC_FUNC_(bind_yaml_output_errors, BIND_YAML_OUTPUT_ERRORS)(void);
void yaml_output_errors(void)
{
  FC_FUNC_(bind_yaml_output_errors, BIND_YAML_OUTPUT_ERRORS)
    ();
}

void FC_FUNC_(bind_yaml_set_default_stream, BIND_YAML_SET_DEFAULT_STREAM)(const int*,
  int*);
void yaml_set_default_stream(int unit,
  int ierr)
{
  FC_FUNC_(bind_yaml_set_default_stream, BIND_YAML_SET_DEFAULT_STREAM)
    (&unit, &ierr);
}

void FC_FUNC_(bind_yaml_get_default_stream, BIND_YAML_GET_DEFAULT_STREAM)(int*);
void yaml_get_default_stream(int unit)
{
  FC_FUNC_(bind_yaml_get_default_stream, BIND_YAML_GET_DEFAULT_STREAM)
    (&unit);
}

void FC_FUNC_(bind_yaml_stream_connected, BIND_YAML_STREAM_CONNECTED)(const size_t*,
  const char*,
  int*,
  int*,
  size_t);
void yaml_stream_connected(const char* filename,
  int unit,
  int (*istat))
{
  size_t filename_len = filename ? strlen(filename) : 0;
  FC_FUNC_(bind_yaml_stream_connected, BIND_YAML_STREAM_CONNECTED)
    (&filename_len, filename, &unit, istat, filename_len);
}

void FC_FUNC_(bind_yaml_set_stream, BIND_YAML_SET_STREAM)(const int*,
  const size_t*,
  const char*,
  int*,
  const int*,
  const int*,
  const size_t*,
  const char*,
  const bool*,
  size_t,
  size_t);
void yaml_set_stream(int (*unit),
  const char (*filename),
  int (*istat),
  int (*tabbing),
  int (*record_length),
  const char (*position),
  bool (*setdefault))
{
  size_t filename_len = filename ? strlen(filename) : 0;
  size_t position_len = position ? strlen(position) : 0;
  FC_FUNC_(bind_yaml_set_stream, BIND_YAML_SET_STREAM)
    (unit, &filename_len, filename, istat, tabbing, record_length, &position_len, position, setdefault, filename_len, position_len);
}

void FC_FUNC_(bind_yaml_stream_attributes, BIND_YAML_STREAM_ATTRIBUTES)(const int*,
  const int*,
  int*,
  bool*,
  int*,
  int*,
  int*,
  int*,
  int*,
  int*,
  int*);
void yaml_stream_attributes(int (*unit),
  int (*stream_unit),
  int (*icursor),
  bool (*flowrite),
  int (*itab_active),
  int (*iflowlevel),
  int (*ilevel),
  int (*ilast),
  int (*indent),
  int (*indent_previous),
  int (*record_length))
{
  FC_FUNC_(bind_yaml_stream_attributes, BIND_YAML_STREAM_ATTRIBUTES)
    (unit, stream_unit, icursor, flowrite, itab_active, iflowlevel, ilevel, ilast, indent, indent_previous, record_length);
}

void FC_FUNC_(bind_yaml_new_document, BIND_YAML_NEW_DOCUMENT)(const int*);
void yaml_new_document(int (*unit))
{
  FC_FUNC_(bind_yaml_new_document, BIND_YAML_NEW_DOCUMENT)
    (unit);
}

void FC_FUNC_(bind_yaml_flush_document, BIND_YAML_FLUSH_DOCUMENT)(const int*);
void yaml_flush_document(int (*unit))
{
  FC_FUNC_(bind_yaml_flush_document, BIND_YAML_FLUSH_DOCUMENT)
    (unit);
}

void FC_FUNC_(bind_yaml_close_stream, BIND_YAML_CLOSE_STREAM)(const int*,
  int*);
void yaml_close_stream(int (*unit),
  int (*istat))
{
  FC_FUNC_(bind_yaml_close_stream, BIND_YAML_CLOSE_STREAM)
    (unit, istat);
}

void FC_FUNC_(bind_yaml_close_all_streams, BIND_YAML_CLOSE_ALL_STREAMS)(void);
void yaml_close_all_streams(void)
{
  FC_FUNC_(bind_yaml_close_all_streams, BIND_YAML_CLOSE_ALL_STREAMS)
    ();
}

void FC_FUNC_(bind_yaml_dict_inspect, BIND_YAML_DICT_INSPECT)(f90_dictionary_pointer*);
void yaml_dict_inspect(f90_dictionary_pointer* dict)
{
  FC_FUNC_(bind_yaml_dict_inspect, BIND_YAML_DICT_INSPECT)
    (dict);
}

void FC_FUNC_(bind_dump_progress_bar, BIND_DUMP_PROGRESS_BAR)(f90_f_progress_bar*,
  const int*,
  const int*);
void dump_progress_bar(f90_f_progress_bar* bar,
  int (*step),
  int (*unit))
{
  FC_FUNC_(bind_dump_progress_bar, BIND_DUMP_PROGRESS_BAR)
    (bar, step, unit);
}

void FC_FUNC_(bind_yaml_cite, BIND_YAML_CITE)(const size_t*,
  const char*,
  const int*,
  size_t);
void yaml_cite(const char* paper,
  int (*unit))
{
  size_t paper_len = paper ? strlen(paper) : 0;
  FC_FUNC_(bind_yaml_cite, BIND_YAML_CITE)
    (&paper_len, paper, unit, paper_len);
}

void FC_FUNC_(bind_yaml_bib_dump, BIND_YAML_BIB_DUMP)(f90_dictionary_pointer*,
  const int*);
void yaml_bib_dump(f90_dictionary_pointer* citations,
  int (*unit))
{
  FC_FUNC_(bind_yaml_bib_dump, BIND_YAML_BIB_DUMP)
    (citations, unit);
}

void FC_FUNC_(bind_yaml_scalar, BIND_YAML_SCALAR)(const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t);
void yaml_scalar(const char* message,
  const char (*advance),
  int (*unit),
  const char (*hfill))
{
  size_t message_len = message ? strlen(message) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t hfill_len = hfill ? strlen(hfill) : 0;
  FC_FUNC_(bind_yaml_scalar, BIND_YAML_SCALAR)
    (&message_len, message, &advance_len, advance, unit, &hfill_len, hfill, message_len, advance_len, hfill_len);
}

void FC_FUNC_(bind_yaml_mapping_open, BIND_YAML_MAPPING_OPEN)(const size_t*,
  const char*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const bool*,
  const int*,
  const size_t*,
  const char*,
  const int*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_mapping_open(const char (*mapname),
  const char (*label),
  const char (*tag),
  bool (*flow),
  int (*tabbing),
  const char (*advance),
  int (*unit))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t tag_len = tag ? strlen(tag) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  FC_FUNC_(bind_yaml_mapping_open, BIND_YAML_MAPPING_OPEN)
    (&mapname_len, mapname, &label_len, label, &tag_len, tag, flow, tabbing, &advance_len, advance, unit, mapname_len, label_len, tag_len, advance_len);
}

void FC_FUNC_(bind_yaml_mapping_close, BIND_YAML_MAPPING_CLOSE)(const size_t*,
  const char*,
  const int*,
  size_t);
void yaml_mapping_close(const char (*advance),
  int (*unit))
{
  size_t advance_len = advance ? strlen(advance) : 0;
  FC_FUNC_(bind_yaml_mapping_close, BIND_YAML_MAPPING_CLOSE)
    (&advance_len, advance, unit, advance_len);
}

void FC_FUNC_(bind_yaml_sequence_open, BIND_YAML_SEQUENCE_OPEN)(const size_t*,
  const char*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const bool*,
  const int*,
  const size_t*,
  const char*,
  const int*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_sequence_open(const char (*mapname),
  const char (*label),
  const char (*tag),
  bool (*flow),
  int (*tabbing),
  const char (*advance),
  int (*unit))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t tag_len = tag ? strlen(tag) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  FC_FUNC_(bind_yaml_sequence_open, BIND_YAML_SEQUENCE_OPEN)
    (&mapname_len, mapname, &label_len, label, &tag_len, tag, flow, tabbing, &advance_len, advance, unit, mapname_len, label_len, tag_len, advance_len);
}

void FC_FUNC_(bind_yaml_sequence_close, BIND_YAML_SEQUENCE_CLOSE)(const size_t*,
  const char*,
  const int*,
  size_t);
void yaml_sequence_close(const char (*advance),
  int (*unit))
{
  size_t advance_len = advance ? strlen(advance) : 0;
  FC_FUNC_(bind_yaml_sequence_close, BIND_YAML_SEQUENCE_CLOSE)
    (&advance_len, advance, unit, advance_len);
}

void FC_FUNC_(bind_yaml_newline, BIND_YAML_NEWLINE)(const int*);
void yaml_newline(int (*unit))
{
  FC_FUNC_(bind_yaml_newline, BIND_YAML_NEWLINE)
    (unit);
}

void FC_FUNC_(bind_yaml_sequence, BIND_YAML_SEQUENCE)(const size_t*,
  const char*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const int*,
  size_t,
  size_t,
  size_t);
void yaml_sequence(const char (*seqvalue),
  const char (*label),
  const char (*advance),
  int (*unit),
  int (*padding))
{
  size_t seqvalue_len = seqvalue ? strlen(seqvalue) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  FC_FUNC_(bind_yaml_sequence, BIND_YAML_SEQUENCE)
    (&seqvalue_len, seqvalue, &label_len, label, &advance_len, advance, unit, padding, seqvalue_len, label_len, advance_len);
}

void FC_FUNC_(bind_yaml_dict_dump, BIND_YAML_DICT_DUMP)(const f90_dictionary_pointer*,
  const int*,
  const bool*,
  const bool*);
void yaml_dict_dump(const f90_dictionary_pointer* dict,
  int (*unit),
  bool (*flow),
  bool (*verbatim))
{
  FC_FUNC_(bind_yaml_dict_dump, BIND_YAML_DICT_DUMP)
    (dict, unit, flow, verbatim);
}

void FC_FUNC_(bind_yaml_dict_dump_all, BIND_YAML_DICT_DUMP_ALL)(const f90_dictionary_pointer*,
  const int*,
  const bool*,
  const bool*);
void yaml_dict_dump_all(const f90_dictionary_pointer* dict,
  int (*unit),
  bool (*flow),
  bool (*verbatim))
{
  FC_FUNC_(bind_yaml_dict_dump_all, BIND_YAML_DICT_DUMP_ALL)
    (dict, unit, flow, verbatim);
}

void FC_FUNC_(bind_yaml_release_document, BIND_YAML_RELEASE_DOCUMENT)(const int*);
void yaml_release_document(int (*unit))
{
  FC_FUNC_(bind_yaml_release_document, BIND_YAML_RELEASE_DOCUMENT)
    (unit);
}

void FC_FUNC_(bind_yaml_map, BIND_YAML_MAP)(const size_t*,
  const char*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  size_t,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map(const char* mapname,
  const char* mapvalue,
  const char (*label),
  const char (*tag),
  const char (*advance),
  int (*unit))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t mapvalue_len = mapvalue ? strlen(mapvalue) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t tag_len = tag ? strlen(tag) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  FC_FUNC_(bind_yaml_map, BIND_YAML_MAP)
    (&mapname_len, mapname, &mapvalue_len, mapvalue, &label_len, label, &tag_len, tag, &advance_len, advance, unit, mapname_len, mapvalue_len, label_len, tag_len, advance_len);
}

void FC_FUNC_(bind_yaml_map_dict, BIND_YAML_MAP_DICT)(const size_t*,
  const char*,
  const f90_dictionary_pointer*,
  const size_t*,
  const char*,
  const int*,
  const bool*,
  size_t,
  size_t);
void yaml_map_dict(const char* mapname,
  const f90_dictionary_pointer* mapvalue,
  const char (*label),
  int (*unit),
  bool (*flow))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  FC_FUNC_(bind_yaml_map_dict, BIND_YAML_MAP_DICT)
    (&mapname_len, mapname, mapvalue, &label_len, label, unit, flow, mapname_len, label_len);
}

void FC_FUNC_(bind_yaml_map_enum, BIND_YAML_MAP_ENUM)(const size_t*,
  const char*,
  const f90_f_enumerator*,
  const size_t*,
  const char*,
  const int*,
  const bool*,
  size_t,
  size_t);
void yaml_map_enum(const char* mapname,
  const f90_f_enumerator* mapvalue,
  const char (*label),
  int (*unit),
  bool (*flow))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  FC_FUNC_(bind_yaml_map_enum, BIND_YAML_MAP_ENUM)
    (&mapname_len, mapname, mapvalue, &label_len, label, unit, flow, mapname_len, label_len);
}

void FC_FUNC_(bind_yaml_map_li, BIND_YAML_MAP_LI)(const size_t*,
  const char*,
  const long*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_li(const char* mapname,
  long mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_li, BIND_YAML_MAP_LI)
    (&mapname_len, mapname, &mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_i, BIND_YAML_MAP_I)(const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_i(const char* mapname,
  int mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_i, BIND_YAML_MAP_I)
    (&mapname_len, mapname, &mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_f, BIND_YAML_MAP_F)(const size_t*,
  const char*,
  const float*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_f(const char* mapname,
  float mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_f, BIND_YAML_MAP_F)
    (&mapname_len, mapname, &mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_d, BIND_YAML_MAP_D)(const size_t*,
  const char*,
  const double*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_d(const char* mapname,
  double mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_d, BIND_YAML_MAP_D)
    (&mapname_len, mapname, &mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_l, BIND_YAML_MAP_L)(const size_t*,
  const char*,
  const bool*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_l(const char* mapname,
  bool mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_l, BIND_YAML_MAP_L)
    (&mapname_len, mapname, &mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_dv, BIND_YAML_MAP_DV)(const size_t*,
  const char*,
  const size_t*,
  const double*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_dv(const char* mapname,
  size_t mapvalue_dim_0,
  const double* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_dv, BIND_YAML_MAP_DV)
    (&mapname_len, mapname, &mapvalue_dim_0, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_rv, BIND_YAML_MAP_RV)(const size_t*,
  const char*,
  const size_t*,
  const float*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_rv(const char* mapname,
  size_t mapvalue_dim_0,
  const float* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_rv, BIND_YAML_MAP_RV)
    (&mapname_len, mapname, &mapvalue_dim_0, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_cv, BIND_YAML_MAP_CV)(const size_t*,
  const char*,
  const size_t*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_cv(const char* mapname,
  size_t mapvalue_dim_0,
  const char* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t mapvalue_len = mapvalue ? strlen(mapvalue) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_cv, BIND_YAML_MAP_CV)
    (&mapname_len, mapname, &mapvalue_dim_0, &mapvalue_len, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, mapvalue_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_lv, BIND_YAML_MAP_LV)(const size_t*,
  const char*,
  const size_t*,
  const bool*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_lv(const char* mapname,
  size_t mapvalue_dim_0,
  const bool* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_lv, BIND_YAML_MAP_LV)
    (&mapname_len, mapname, &mapvalue_dim_0, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_liv, BIND_YAML_MAP_LIV)(const size_t*,
  const char*,
  const size_t*,
  const long*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_liv(const char* mapname,
  size_t mapvalue_dim_0,
  const long* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_liv, BIND_YAML_MAP_LIV)
    (&mapname_len, mapname, &mapvalue_dim_0, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_iv, BIND_YAML_MAP_IV)(const size_t*,
  const char*,
  const size_t*,
  const int*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_iv(const char* mapname,
  size_t mapvalue_dim_0,
  const int* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_iv, BIND_YAML_MAP_IV)
    (&mapname_len, mapname, &mapvalue_dim_0, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_dm, BIND_YAML_MAP_DM)(const size_t*,
  const char*,
  const size_t*,
  const size_t*,
  const double*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_dm(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const double* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_dm, BIND_YAML_MAP_DM)
    (&mapname_len, mapname, &mapvalue_dim_0, &mapvalue_dim_1, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_rm, BIND_YAML_MAP_RM)(const size_t*,
  const char*,
  const size_t*,
  const size_t*,
  const float*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_rm(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const float* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_rm, BIND_YAML_MAP_RM)
    (&mapname_len, mapname, &mapvalue_dim_0, &mapvalue_dim_1, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_im, BIND_YAML_MAP_IM)(const size_t*,
  const char*,
  const size_t*,
  const size_t*,
  const int*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_im(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const int* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_im, BIND_YAML_MAP_IM)
    (&mapname_len, mapname, &mapvalue_dim_0, &mapvalue_dim_1, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_lm, BIND_YAML_MAP_LM)(const size_t*,
  const char*,
  const size_t*,
  const size_t*,
  const bool*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_lm(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const bool* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_lm, BIND_YAML_MAP_LM)
    (&mapname_len, mapname, &mapvalue_dim_0, &mapvalue_dim_1, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_map_dt, BIND_YAML_MAP_DT)(const size_t*,
  const char*,
  const size_t*,
  const size_t*,
  const size_t*,
  const double*,
  const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  size_t,
  size_t,
  size_t,
  size_t);
void yaml_map_dt(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  size_t mapvalue_dim_2,
  const double* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt))
{
  size_t mapname_len = mapname ? strlen(mapname) : 0;
  size_t label_len = label ? strlen(label) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t fmt_len = fmt ? strlen(fmt) : 0;
  FC_FUNC_(bind_yaml_map_dt, BIND_YAML_MAP_DT)
    (&mapname_len, mapname, &mapvalue_dim_0, &mapvalue_dim_1, &mapvalue_dim_2, mapvalue, &label_len, label, &advance_len, advance, unit, &fmt_len, fmt, mapname_len, label_len, advance_len, fmt_len);
}

void FC_FUNC_(bind_yaml_warning_str, BIND_YAML_WARNING_STR)(const f90_f_string*,
  const int*,
  const int*);
void yaml_warning_str(const f90_f_string* message,
  int (*level),
  int (*unit))
{
  FC_FUNC_(bind_yaml_warning_str, BIND_YAML_WARNING_STR)
    (message, level, unit);
}

void FC_FUNC_(bind_yaml_warning_c, BIND_YAML_WARNING_C)(const size_t*,
  const char*,
  const int*,
  const int*,
  size_t);
void yaml_warning_c(const char* message,
  int (*level),
  int (*unit))
{
  size_t message_len = message ? strlen(message) : 0;
  FC_FUNC_(bind_yaml_warning_c, BIND_YAML_WARNING_C)
    (&message_len, message, level, unit, message_len);
}

void FC_FUNC_(bind_yaml_comment_str, BIND_YAML_COMMENT_STR)(const f90_f_string*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  const int*,
  size_t,
  size_t);
void yaml_comment_str(const f90_f_string* message,
  const char (*advance),
  int (*unit),
  const char (*hfill),
  int (*tabbing))
{
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t hfill_len = hfill ? strlen(hfill) : 0;
  FC_FUNC_(bind_yaml_comment_str, BIND_YAML_COMMENT_STR)
    (message, &advance_len, advance, unit, &hfill_len, hfill, tabbing, advance_len, hfill_len);
}

void FC_FUNC_(bind_yaml_comment_c, BIND_YAML_COMMENT_C)(const size_t*,
  const char*,
  const size_t*,
  const char*,
  const int*,
  const size_t*,
  const char*,
  const int*,
  size_t,
  size_t,
  size_t);
void yaml_comment_c(const char* message,
  const char (*advance),
  int (*unit),
  const char (*hfill),
  int (*tabbing))
{
  size_t message_len = message ? strlen(message) : 0;
  size_t advance_len = advance ? strlen(advance) : 0;
  size_t hfill_len = hfill ? strlen(hfill) : 0;
  FC_FUNC_(bind_yaml_comment_c, BIND_YAML_COMMENT_C)
    (&message_len, message, &advance_len, advance, unit, &hfill_len, hfill, tabbing, message_len, advance_len, hfill_len);
}

