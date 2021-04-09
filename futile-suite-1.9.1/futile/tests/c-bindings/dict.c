#include "misc.h"
#include "tree.h"

#include <string.h>
#include <stdio.h>

#define ASSERT(T) {if (!(T)) {fprintf(stdout, "%s: failed\n", #T); return FALSE;}}
#define ALLEQ(T,U,N) {int it; for(it=0;it < N;it=it+1) { ASSERT(T[it]==U[it]);}}

gboolean test_dict_new()
{
  FutileTree *dict;
  FutileTreeIter iter;

  dict = futile_tree_new(NULL);
  ASSERT(dict != (FutileTree*)0);
  futile_tree_unref(dict);

  dict = futile_tree_new(&iter);
  ASSERT(F_TYPE(iter.pointer) == F_TYPE(dict->root));
  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_add()
{
  FutileTree *dict;
  FutileTreeIter iter;
  gchar *str;

  dict = futile_tree_new(&iter);
  
  futile_tree_iter_add(&iter, "0.");
  ASSERT(futile_tree_iter_len(&iter) == 1);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "__list__"));
  g_free(str);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_set()
{
  FutileTree *dict;
  FutileTreeIter iter;
  gchar *str;

  dict = futile_tree_new(&iter);
  
  futile_tree_iter_set(&iter, "nitermax", "45");
  ASSERT(futile_tree_iter_len(&iter) == 1);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "__dict__"));
  g_free(str);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_pop()
{
  FutileTree *dict;
  FutileTreeIter iter;
  gboolean valid;

  dict = futile_tree_new(&iter);
  
  futile_tree_iter_set(&iter, "nitermax", "45");
  futile_tree_iter_set(&iter, "ncong", "4");
  ASSERT(futile_tree_iter_len(&iter) == 2);

  valid = futile_tree_iter_pop(&iter, "nitermax2");
  ASSERT(valid == FALSE);
  valid = futile_tree_iter_pop(&iter, "nitermax");
  ASSERT(valid == TRUE);
  ASSERT(futile_tree_iter_len(&iter) == 1);
  valid = futile_tree_iter_pop(&iter, "ncong");
  ASSERT(valid == TRUE);
  ASSERT(futile_tree_iter_len(&iter) == 0);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_insert()
{
  FutileTree *dict;
  FutileTreeIter root, iter;
  gchar *str;

  dict = futile_tree_new(&root);

  futile_tree_iter_insert(&root, "nitermax", &iter);
  ASSERT(futile_tree_iter_len(&root) == 1);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, "nitermax"));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  futile_tree_iter_set(&iter, NULL, "45");
  ASSERT(futile_tree_iter_len(&root) == 1);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "45"));
  g_free(str);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_append()
{
  FutileTree *dict;
  FutileTreeIter root, iter;
  gchar *str;

  dict = futile_tree_new(&root);

  futile_tree_iter_append(&root, &iter);
  ASSERT(futile_tree_iter_len(&root) == 1);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  futile_tree_iter_set(&iter, NULL, "0.");
  ASSERT(futile_tree_iter_len(&root) == 1);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "0."));
  g_free(str);
  futile_tree_iter_append(&root, &iter);
  ASSERT(futile_tree_iter_len(&root) == 2);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_at_key()
{
  FutileTree *dict;
  FutileTreeIter root, iter;
  gchar *str;
  gboolean exists;

  dict = futile_tree_new(&root);

  futile_tree_iter_set(&root, "idsx", "6");
  futile_tree_iter_set(&root, "nitermax", "45");
  ASSERT(futile_tree_iter_len(&root) == 2);

  exists = futile_tree_iter_at_key(&root, "ncong", &iter);
  ASSERT(exists == FALSE);
  exists = futile_tree_iter_at_key(&root, "nitermax", &iter);
  ASSERT(exists == TRUE);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, "nitermax"));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "45"));
  g_free(str);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_at_item()
{
  FutileTree *dict;
  FutileTreeIter root, iter;
  gchar *str;
  gboolean exists;

  dict = futile_tree_new(&root);

  futile_tree_iter_append(&root, &iter);
  futile_tree_iter_set(&iter, NULL, "0.");
  futile_tree_iter_append(&root, &iter);
  futile_tree_iter_set(&iter, NULL, "1.");
  futile_tree_iter_append(&root, &iter);
  futile_tree_iter_set(&iter, NULL, "2.");
  ASSERT(futile_tree_iter_len(&root) == 3);

  exists = futile_tree_iter_at_item(&root, 5, &iter);
  ASSERT(exists == FALSE);
  exists = futile_tree_iter_at_item(&root, 1, &iter);
  ASSERT(exists == TRUE);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "1."));
  g_free(str);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_loop()
{
  FutileTree *dict;
  FutileTreeIter root, iter;
  gchar *str;
  gboolean exists;
  const gchar *arr[] = {"0.", "1.", "2.", NULL};

  dict = futile_tree_new(&root);

  futile_tree_iter_set_array(&root, NULL, arr);
  ASSERT(futile_tree_iter_len(&root) == 3);

  exists = futile_tree_iter_first(&root, &iter);
  ASSERT(exists == TRUE);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "0."));
  g_free(str);
  exists = futile_tree_iter_next(&iter);
  ASSERT(exists == TRUE);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "1."));
  g_free(str);
  exists = futile_tree_iter_next(&iter);
  ASSERT(exists == TRUE);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "2."));
  g_free(str);
  exists = futile_tree_iter_next(&iter);
  ASSERT(exists == FALSE);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_set_array()
{
  FutileTree *dict;
  FutileTreeIter iter, lst;
  const gchar *arr[] = {"0.", "1.", "2.", NULL};
  gchar *str;

  dict = futile_tree_new(&iter);
  
  futile_tree_iter_set_array(&iter, "coord", arr);

  ASSERT(futile_tree_iter_len(&iter) == 1);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "__dict__"));
  g_free(str);
  futile_tree_iter_at_key(&iter, "coord", &lst);
  ASSERT(futile_tree_iter_len(&lst) == 3);
  str = futile_tree_iter_key(&lst);
  ASSERT(!strcmp(str, "coord"));
  g_free(str);
  str = futile_tree_iter_value(&lst);
  ASSERT(!strcmp(str, "__list__"));
  g_free(str);
  futile_tree_iter_at_item(&lst, 1, &iter);
  ASSERT(futile_tree_iter_len(&iter) == 0);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "1."));
  g_free(str);

  futile_tree_unref(dict);

  return TRUE;
}

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


gboolean test_dict_set_dict()
{
  FutileTree *dict, *leaf;
  FutileTreeIter iter, lst;
  gchar *str;
  const gchar *arr[] = {"0.", "1.", "2.", NULL};
  gboolean exists;

  leaf = futile_tree_new(&iter);
  futile_tree_iter_set_array(&iter, "Zn", arr);

  dict = futile_tree_new(&iter);
  futile_tree_iter_append(&iter, &lst);
  futile_tree_iter_set(&lst, "units", "bohr");
  futile_tree_iter_append(&iter, &lst);
  futile_tree_iter_set_tree(&lst, NULL, leaf);

  futile_tree_unref(leaf);

  ASSERT(futile_tree_iter_len(&iter) == 2);
  str = futile_tree_iter_key(&iter);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "__list__"));
  g_free(str);
  exists = futile_tree_iter_at_item(&iter, 1, &lst);
  ASSERT(exists == TRUE);
  ASSERT(futile_tree_iter_len(&lst) == 1);
  str = futile_tree_iter_key(&lst);
  ASSERT(!strcmp(str, ""));
  g_free(str);
  str = futile_tree_iter_value(&lst);
  ASSERT(!strcmp(str, "__dict__"));
  g_free(str);
  exists = futile_tree_iter_at_key(&lst, "Zn", &iter);
  ASSERT(exists == TRUE);
  str = futile_tree_iter_value(&iter);
  ASSERT(!strcmp(str, "__list__"));
  g_free(str);
  ASSERT(futile_tree_iter_len(&iter) == 3);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_from_yaml()
{
  const gchar *yaml = "---\n"
"cell: [10., 20., 30.]\n"
"positions:\n"
"- Si: [0., 0., 0.]\n"
"- Si: [1., 2., 3.]\n"
"\n";
  FutileTree *dict;
  FutileTreeIter iter, lst;
  gboolean exists;
  
  dict = futile_tree_new_from_yaml(yaml, &iter);

  ASSERT(futile_tree_iter_len(&iter) == 2);
  exists = futile_tree_iter_at_key(&iter, "positions", &lst);
  ASSERT(exists == TRUE);
  ASSERT(futile_tree_iter_len(&lst) == 2);

  futile_tree_unref(dict);

  return TRUE;
}

gboolean test_dict_by_hand()
{
  FutileTree *dict;
  FutileTreeIter iter, lst, at, coord;
  gboolean exists;
  
  dict = futile_tree_new(&iter);

  futile_tree_iter_insert(&iter, "cell", &lst);
  futile_tree_iter_add(&lst, "10.");
  futile_tree_iter_add(&lst, "20.");
  futile_tree_iter_add(&lst, "30.");

  futile_tree_iter_insert(&iter, "positions", &lst);

  futile_tree_iter_append(&lst, &at);
  futile_tree_iter_insert(&at, "Si", &coord);
  futile_tree_iter_add(&coord, "0.");
  futile_tree_iter_add(&coord, "0.");
  futile_tree_iter_add(&coord, "0.");

  futile_tree_iter_append(&lst, &at);
  futile_tree_iter_insert(&at, "Si", &coord);
  futile_tree_iter_add(&coord, "1.");
  futile_tree_iter_add(&coord, "2.");
  futile_tree_iter_add(&coord, "3.");

  ASSERT(futile_tree_iter_len(&iter) == 2);
  exists = futile_tree_iter_at_key(&iter, "positions", &lst);
  ASSERT(exists == TRUE);
  ASSERT(futile_tree_iter_len(&lst) == 2);

  futile_tree_unref(dict);

  return TRUE;
}

#define RUN(T) {if (!T) return 1; else fprintf(stdout, "%s: OK\n", #T);}

int main(int argc, char **argv)
{
#ifdef GLIB_MAJOR_VERSION
  g_type_init();
#endif
  
  RUN(test_dict_new());
  RUN(test_dict_add());
  RUN(test_dict_set());
  RUN(test_dict_pop());
  RUN(test_dict_insert());
  RUN(test_dict_append());
  RUN(test_dict_at_key());
  RUN(test_dict_at_item());
  RUN(test_dict_loop());
  RUN(test_dict_set_array());
  RUN(test_dict_set_dict());
  RUN(test_dict_by_hand());

  RUN(test_f90_dict());

  futile_initialize();
  RUN(test_dict_from_yaml());
  futile_finalize();


  return 0;
}
