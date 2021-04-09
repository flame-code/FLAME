#include "dict.h"

#include <string.h>
#include <stdio.h>

#define ASSERT(T) {if (!(T)) {fprintf(stdout, "%s: failed\n", #T); return FALSE;}}
#define ALLEQ(T,U,N) {int it; for(it=0;it < N;it=it+1) { ASSERT(T[it]==U[it]);}}

gboolean test_f90_dict()
{
  f90_dictionary_pointer dict,dict_tmp;
  f90_dictionary_iterator it;
  const double arr[] = {1.2, 2.3, 3.4, 4.5};
  double arr_res[4];
  int i;
  double val;

  dict_init(&dict);
  /*fill a dict*/
  dict_set_double_array(&dict,"array",arr,3);
  /*try to retrieve it with the wrong key*/
  ASSERT(!dict_get_double_array(&dict,"arraywrong",arr_res,3));
  /*now add an extra element*/
  ASSERT(dict_get_dict(&dict,"array",&dict_tmp));
  dict_add_double(&dict_tmp,arr[3]);
  /*try to retrieve it with the good key*/
  ASSERT(dict_get_double_array(&dict,"array",arr_res,4));
  /*compare element per element*/
  ALLEQ(arr,arr_res,4);
  ASSERT(dict_get_dict(&dict,"array",&dict_tmp));
  for (i = 0, dict_iter_new(&it, &dict_tmp); iterate(&it); i++)
    ASSERT(arr[i] == atof(it.value));
  dict_free(&dict);

  return TRUE;
}

#define RUN(T) {if (!T) return 1; else fprintf(stdout, "%s: OK\n", #T);}

int main(int argc, char **argv)
{
#ifdef GLIB_MAJOR_VERSION
  g_type_init();
#endif 
  RUN(test_f90_dict());
  return 0;
}
