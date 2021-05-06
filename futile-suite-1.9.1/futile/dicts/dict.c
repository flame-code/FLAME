#include <config.h>
#include "dict.h"
#include "dict-fapi.h"

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

#define VALID(I) (I==0) ? true : false

#define RSTRIP(buf) {int i; for (i = max_field_length-1; i >= 0 && buf[i] == ' '; i--) buf[i] = '\0'; buf[max_field_length] = '\0';}


void dict_init(f90_dictionary_pointer *dict)
{
  FC_FUNC_(bind_dict_init, BIND_DICT_INIT)(dict);
}

void dict_free(f90_dictionary_pointer *dict)
{
  FC_FUNC_(bind_dict_free, BIND_DICT_FREE)(dict);
}

void dict_set_string(f90_dictionary_pointer* dict, const char* key ,const char* value)
{
  FC_FUNC_(bind_dict_set_string, BIND_DICT_SET_STRING)(dict, key, value);
}

void dict_set_double(f90_dictionary_pointer* dict, const char* key ,double value)
{
  int klen = strlen(key);
  FC_FUNC_(bind_dict_set_double, BIND_DICT_SET_DOUBLE)(dict, key, &klen, &value, klen);
}

void dict_set_double_array(f90_dictionary_pointer* dict, const char* key , const double * array, size_t len)
{
  int flen = len;
  int klen = strlen(key);
  FC_FUNC_(bind_dict_set_double_array, BIND_DICT_SET_DOUBLE_ARRAY)(dict, key, &klen, array, &flen, klen);
}

void dict_set_double_matrix(f90_dictionary_pointer* dict, const char* key , const double * array, size_t lenx, size_t leny)
{
  int flenx = lenx;
  int fleny = leny;
  int klen = strlen(key);
  FC_FUNC_(bind_dict_set_double_matrix, BIND_DICT_SET_DOUBLE_MATRIX)(dict, key, &klen, array, &flenx, &fleny, klen);
}


void dict_set_dict(f90_dictionary_pointer* dict, const char* key , f90_dictionary_pointer* value)
{
  int klen = strlen(key);
  FC_FUNC_(bind_dict_set_dict,BIND_DICT_SET_DICT)(dict, key, &klen, value, klen);
}

void dict_add_dict(f90_dictionary_pointer* dict, f90_dictionary_pointer* value)
{
  FC_FUNC_(bind_dict_add_dict,BIND_DICT_ADD_DICT)(dict,value);
}

void dict_add_double(f90_dictionary_pointer* dict, double value)
{
  FC_FUNC_(bind_dict_add_double,BIND_DICT_ADD_DOUBLE)(dict,&value);
}

void dict_add_int(f90_dictionary_pointer* dict, int value)
{
  FC_FUNC_(bind_dict_add_int,BIND_DICT_ADD_INT)(dict,&value);
}

void dict_add_string(f90_dictionary_pointer* dict, const char *value)
{
  FC_FUNC_(bind_dict_add_char,BIND_DICT_ADD_CHAR)(dict,value);
}


bool dict_get_double_array(f90_dictionary_pointer* dict, const char* key , double * array, size_t len)
{
  int istat;
  int flen = len;
  int klen = strlen(key);
  FC_FUNC_(bind_dict_get_double_array, BIND_DICT_GET_DOUBLE_ARRAY)(dict, key, &klen, array, &flen, &istat, klen);
  return VALID(istat);
}

bool dict_get_dict(f90_dictionary_pointer* dict, const char* key , f90_dictionary_pointer* value)
{
  int istat;
  int klen = strlen(key);
  FC_FUNC_(bind_dict_get_dict,BIND_DICT_GET_DICT)(dict, key, &klen, value, &istat, klen);
  return VALID(istat);
}

void dict_iter_new(f90_dictionary_iterator *iter, const f90_dictionary_pointer *parent)
{
  iter->parent = parent;
  FC_FUNC_(bind_iter_null,BIND_ITER_NULL)(&iter->dict);
}

bool iterate(f90_dictionary_iterator* iter)
{
  int istat;

  FC_FUNC_(bind_iterate,BIND_ITERATE)(&iter->dict,iter->parent,&istat);
  if (!VALID(istat))
    return false;

  FC_FUNC_(bind_dict_key,BIND_DICT_KEY)(&iter->dict,iter->key,max_field_length);
  RSTRIP(iter->key);
  FC_FUNC_(bind_dict_value,BIND_DICT_VALUE)(&iter->dict, iter->value,max_field_length);
  RSTRIP(iter->value);

  return true;
}

void futile_dicts_initialize()
{
  FC_FUNC_(f_dicts_initialize, F_DICTS_INITIALIZE)();
}
void futile_dicts_finalize()
{
  FC_FUNC_(f_dicts_finalize, F_DICTS_FINALIZE)();
}
