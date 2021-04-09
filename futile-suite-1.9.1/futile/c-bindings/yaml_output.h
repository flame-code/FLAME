#ifndef YAML_OUTPUT_H
#define YAML_OUTPUT_H

#include "futile_cst.h"
#include "f_progress_bar.h"
#include "f_string.h"
#include "dict.h"
#include "enum.h"


void yaml_swap_stream(int new_unit,
  int old_unit,
  int ierr);

void yaml_output_errors(void);

void yaml_set_default_stream(int unit,
  int ierr);

void yaml_get_default_stream(int unit);

void yaml_stream_connected(const char* filename,
  int unit,
  int (*istat));

void yaml_set_stream(int (*unit),
  const char (*filename),
  int (*istat),
  int (*tabbing),
  int (*record_length),
  const char (*position),
  bool (*setdefault));

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
  int (*record_length));

void yaml_new_document(int (*unit));

void yaml_flush_document(int (*unit));

void yaml_close_stream(int (*unit),
  int (*istat));

void yaml_close_all_streams(void);

void yaml_dict_inspect(f90_dictionary_pointer* dict);

void dump_progress_bar(f90_f_progress_bar* bar,
  int (*step),
  int (*unit));

void yaml_cite(const char* paper,
  int (*unit));

void yaml_bib_dump(f90_dictionary_pointer* citations,
  int (*unit));

void yaml_scalar(const char* message,
  const char (*advance),
  int (*unit),
  const char (*hfill));

void yaml_mapping_open(const char (*mapname),
  const char (*label),
  const char (*tag),
  bool (*flow),
  int (*tabbing),
  const char (*advance),
  int (*unit));

void yaml_mapping_close(const char (*advance),
  int (*unit));

void yaml_sequence_open(const char (*mapname),
  const char (*label),
  const char (*tag),
  bool (*flow),
  int (*tabbing),
  const char (*advance),
  int (*unit));

void yaml_sequence_close(const char (*advance),
  int (*unit));

void yaml_newline(int (*unit));

void yaml_sequence(const char (*seqvalue),
  const char (*label),
  const char (*advance),
  int (*unit),
  int (*padding));

void yaml_dict_dump(const f90_dictionary_pointer* dict,
  int (*unit),
  bool (*flow),
  bool (*verbatim));

void yaml_dict_dump_all(const f90_dictionary_pointer* dict,
  int (*unit),
  bool (*flow),
  bool (*verbatim));

void yaml_release_document(int (*unit));

void yaml_map(const char* mapname,
  const char* mapvalue,
  const char (*label),
  const char (*tag),
  const char (*advance),
  int (*unit));

void yaml_map_dict(const char* mapname,
  const f90_dictionary_pointer* mapvalue,
  const char (*label),
  int (*unit),
  bool (*flow));

void yaml_map_enum(const char* mapname,
  const f90_f_enumerator* mapvalue,
  const char (*label),
  int (*unit),
  bool (*flow));

void yaml_map_li(const char* mapname,
  long mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_i(const char* mapname,
  int mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_f(const char* mapname,
  float mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_d(const char* mapname,
  double mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_l(const char* mapname,
  bool mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_dv(const char* mapname,
  size_t mapvalue_dim_0,
  const double* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_rv(const char* mapname,
  size_t mapvalue_dim_0,
  const float* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_cv(const char* mapname,
  size_t mapvalue_dim_0,
  const char* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_lv(const char* mapname,
  size_t mapvalue_dim_0,
  const bool* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_liv(const char* mapname,
  size_t mapvalue_dim_0,
  const long* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_iv(const char* mapname,
  size_t mapvalue_dim_0,
  const int* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_dm(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const double* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_rm(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const float* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_im(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const int* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_lm(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  const bool* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_map_dt(const char* mapname,
  size_t mapvalue_dim_0,
  size_t mapvalue_dim_1,
  size_t mapvalue_dim_2,
  const double* mapvalue,
  const char (*label),
  const char (*advance),
  int (*unit),
  const char (*fmt));

void yaml_warning_str(const f90_f_string* message,
  int (*level),
  int (*unit));

void yaml_warning_c(const char* message,
  int (*level),
  int (*unit));

void yaml_comment_str(const f90_f_string* message,
  const char (*advance),
  int (*unit),
  const char (*hfill),
  int (*tabbing));

void yaml_comment_c(const char* message,
  const char (*advance),
  int (*unit),
  const char (*hfill),
  int (*tabbing));

#endif
