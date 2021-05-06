#ifndef DICT_H
#define DICT_H
#include "futile_cst.h"
#include <stdbool.h>

F_DEFINE_TYPE(dictionary);

#define max_field_length 256

typedef struct {
  const f90_dictionary_pointer* parent;
  f90_dictionary_pointer dict;
  char key[max_field_length + 1];
  char value[max_field_length + 1];
} f90_dictionary_iterator;

f90_dictionary_pointer* dict_new(void);

void dict_set_string(f90_dictionary_pointer* dict, const char* key ,const char* value);
void dict_set_double_array(f90_dictionary_pointer* dict, const char* key ,const double * array, size_t len);
void dict_set_double(f90_dictionary_pointer* dict, const char* key ,double value);
void dict_set_double_matrix(f90_dictionary_pointer* dict, const char* key , const double * array, size_t lenx, size_t leny);

void dict_set_dict(f90_dictionary_pointer* dict, const char* key , f90_dictionary_pointer* value);

void dict_add_dict(f90_dictionary_pointer* dict, f90_dictionary_pointer* value);

void dict_add_int(f90_dictionary_pointer* dict, int value);
void dict_add_double(f90_dictionary_pointer* dict, double value);
void dict_add_string(f90_dictionary_pointer* dict, const char *value);

bool dict_get_double_array(f90_dictionary_pointer* dict, const char* key , double * array, size_t len);

bool dict_get_dict(f90_dictionary_pointer* dict, const char* key , f90_dictionary_pointer* value);

void dict_iter_new(f90_dictionary_iterator *iter, const f90_dictionary_pointer *parent);

bool iterate(f90_dictionary_iterator *iter);

void dict_init(f90_dictionary_pointer *dict);

void dict_free(f90_dictionary_pointer *dict);

void futile_dicts_initialize();
void futile_dicts_finalize();

#endif
