#include <string.h>

#include <config.h>

#include "tree.h"
#include "dict-fapi.h"

/******************************/
/* FutileTree data structure */
/******************************/
#ifdef GLIB_MAJOR_VERSION
G_DEFINE_TYPE(FutileTree, futile_tree, G_TYPE_OBJECT)

static void futile_tree_dispose(GObject *tree);
static void futile_tree_finalize(GObject *tree);

static void futile_tree_class_init(FutileTreeClass *klass)
{
  /* Connect the overloading methods. */
  G_OBJECT_CLASS(klass)->dispose      = futile_tree_dispose;
  G_OBJECT_CLASS(klass)->finalize     = futile_tree_finalize;
  /* G_OBJECT_CLASS(klass)->set_property = visu_data_set_property; */
  /* G_OBJECT_CLASS(klass)->get_property = visu_data_get_property; */
}
#endif

static void futile_tree_init(FutileTree *obj)
{
#ifdef GLIB_MAJOR_VERSION
  memset((void*)((char*)obj + sizeof(GObject)), 0, sizeof(FutileTree) - sizeof(GObject));
#else
  memset(obj, 0, sizeof(FutileTree));
  G_OBJECT(obj)->ref_count = 1;
#endif

  /* g_message("New dict %p.", (gpointer)obj); */
}
#ifdef GLIB_MAJOR_VERSION
static void futile_tree_dispose(GObject *obj)
{
  FutileTree *tree = FUTILE_TREE(obj);

  if (tree->dispose_has_run)
    return;
  tree->dispose_has_run = TRUE;

  /* Chain up to the parent class */
  G_OBJECT_CLASS(futile_tree_parent_class)->dispose(obj);
}
#endif
static void futile_tree_finalize(GObject *obj)
{
  FutileTree *tree = FUTILE_TREE(obj);

  /* g_message("Killing %p (%p).", (gpointer)obj, tree->root); */
  if (F_TYPE(tree->root))
    FC_FUNC_(bind_dict_free, BIND_DICT_FREE)(&tree->root);

#ifdef GLIB_MAJOR_VERSION
  G_OBJECT_CLASS(futile_tree_parent_class)->finalize(obj);
#endif
}
void futile_tree_unref(FutileTree *tree)
{
  g_object_unref(G_OBJECT(tree));
#ifdef GLIB_MAJOR_VERSION
#else
  if (G_OBJECT(tree)->ref_count <= 0)
    {
      futile_tree_finalize(G_OBJECT(tree));
      g_free(tree);
    }
#endif
}

/**
 * futile_tree_new:
 * @root: (allow-none) (out caller-allocates):
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
FutileTree *futile_tree_new(FutileTreeIter *root)
{
  FutileTree *tree;

#ifdef GLIB_MAJOR_VERSION
  tree = FUTILE_TREE(g_object_new(FUTILE_TREE_TYPE, NULL));
#else
  tree = g_malloc(sizeof(FutileTree));
  futile_tree_init(tree);
#endif
  FC_FUNC_(bind_dict_init, BIND_DICT_INIT)(&tree->root);

  if (root)
    {
      root->tree = tree;
      root->pointer = tree->root;
    }

  return tree;
}
/**
 * futile_tree_new_from_yaml:
 * @buf: 
 * @root: (allow-none) (out caller-allocates):
 *
 * Pouet.
 *
 * Returns: (transfer full):
 **/
FutileTree *futile_tree_new_from_yaml(const gchar *buf, FutileTreeIter *root)
{
  FutileTree *tree;

#ifdef GLIB_MAJOR_VERSION
  tree = FUTILE_TREE(g_object_new(FUTILE_TREE_TYPE, NULL));
#else
  tree = g_malloc(sizeof(FutileTree));
  futile_tree_init(tree);
#endif
  FC_FUNC_(bind_dict_parse, BIND_DICT_PARSE)(&tree->root, buf, strlen(buf));

  if (root)
    {
      root->tree = tree;
      root->pointer = tree->root;
    }

  return tree;
}
FutileTree *futile_tree_new_from_fortran(f90_dictionary_pointer dictf)
{
  FutileTree *tree;

#ifdef GLIB_MAJOR_VERSION
  tree = FUTILE_TREE(g_object_new(FUTILE_TREE_TYPE, NULL));
#else
  tree = g_malloc(sizeof(FutileTree));
  futile_tree_init(tree);
#endif
  tree->root = dictf;

  return tree;
}
/**
 * futile_tree_move_to_key:
 * @tree: 
 * @iter: (out caller-allocates) (allow-none):
 * @key:
 *
 * 
 *
 * Returns: TRUE, if @key exists.
 **/
gboolean futile_tree_iter_at_key(const FutileTreeIter *at, const gchar *key,
                                 FutileTreeIter *leaf)
{
  int exists;
  FutileTreeIter iter;

  iter = *at;
  FC_FUNC_(bind_dict_move_to_key, BIND_DICT_MOVE_TO_KEY)(&iter.pointer, &exists, key, strlen(key));
  if (leaf)
    *leaf = iter;
  return (exists) ? TRUE : FALSE;
}
/**
 * futile_tree_move_to_item:
 * @dict: 
 * @iter: (out caller-allocates) (allow-none):
 * @id:
 *
 * 
 *
 * Returns: TRUE, if @key exists.
 **/
gboolean futile_tree_iter_at_item(const FutileTreeIter *at, guint id,
                                  FutileTreeIter *leaf)
{
  int exists;
  FutileTreeIter iter;

  iter = *at;
  FC_FUNC_(bind_dict_move_to_item, BIND_DICT_MOVE_TO_ITEM)(&iter.pointer, &exists, (int*)&id);
  if (leaf)
    *leaf = iter;
  return (exists) ? TRUE : FALSE;
}
gboolean futile_tree_iter_first(const FutileTreeIter *at, FutileTreeIter *first)
{
  int exists;
  FutileTreeIter iter;

  iter = *at;
  FC_FUNC_(bind_dict_iter, BIND_DICT_ITER)(&iter.pointer, &exists);
  if (first)
    *first = iter;
  return (exists) ? TRUE : FALSE;
}
/**
 * futile_tree_next:
 * @dict: 
 * @iter: (out caller-allocates) (allow-none):
 *
 * 
 *
 * Returns: 
 **/
gboolean futile_tree_iter_next(FutileTreeIter *iter)
{
  int exists;

  FC_FUNC_(bind_dict_next, BIND_DICT_NEXT)(&iter->pointer, &exists);
  return (exists) ? TRUE : FALSE;
}
/**
 * futile_tree_insert:
 * @at: 
 * @key: 
 * @leaf: (out caller-allocates) (allow-none):
 *
 * Pouet.
 **/
void futile_tree_iter_insert(const FutileTreeIter *at, const gchar *key,
                             FutileTreeIter *leaf)
{
  FutileTreeIter iter;

  iter = *at;
  FC_FUNC_(bind_dict_insert, BIND_DICT_INSERT)(&iter.pointer, key, strlen(key));
  if (leaf)
    *leaf = iter;
}
/**
 * futile_tree_append:
 * @at: 
 * @leaf: (out caller-allocates) (allow-none):
 *
 * Pouet.
 **/
void futile_tree_iter_append(const FutileTreeIter *at, FutileTreeIter *leaf)
{
  FutileTreeIter iter;

  iter = *at;
  FC_FUNC_(bind_dict_append, BIND_DICT_APPEND)(&iter.pointer);
  if (leaf)
    *leaf = iter;
}
void futile_tree_iter_add(FutileTreeIter *iter, const gchar *value)
{
  f90_dictionary_pointer root;
  
  root = iter->pointer;
  FC_FUNC_(bind_dict_append, BIND_DICT_APPEND)(&iter->pointer);
  FC_FUNC_(bind_dict_set, BIND_DICT_SET)(&iter->pointer, value, strlen(value));
  iter->pointer = root;
}
/**
 * futile_tree_set:
 * @dict: 
 * @id: (allow-none): 
 * @value: 
 *
 * Pouet.
 **/
void futile_tree_iter_set(FutileTreeIter *iter, const gchar *id, const gchar *value)
{
  f90_dictionary_pointer root;
  
  root = iter->pointer;
  if (id)
    FC_FUNC_(bind_dict_insert, BIND_DICT_INSERT)(&iter->pointer, id, strlen(id));
  FC_FUNC_(bind_dict_set, BIND_DICT_SET)(&iter->pointer, value, strlen(value));
  iter->pointer = root;
}
/**
 * futile_tree_set_array:
 * @dict: 
 * @id: (allow-none):
 * @value: (array zero-terminated=1):
 *
 * 
 **/
void futile_tree_iter_set_array(FutileTreeIter *iter, const gchar *id, const gchar **value)
{
  guint i;
  f90_dictionary_pointer root, key;

  root = iter->pointer;
  if (id)
    FC_FUNC_(bind_dict_insert, BIND_DICT_INSERT)(&iter->pointer, id, strlen(id));
  key = iter->pointer;
  for (i = 0; value[i]; i++)
    {
      iter->pointer = key;
      FC_FUNC_(bind_dict_append, BIND_DICT_APPEND)(&iter->pointer);
      FC_FUNC_(bind_dict_set, BIND_DICT_SET)(&iter->pointer, value[i], strlen(value[i]));
    }
  iter->pointer = root;
}
/**
 * futile_tree_set_dict:
 * @dict: 
 * @id: (allow-none):
 * @value:
 *
 * 
 **/
void  futile_tree_iter_set_tree(FutileTreeIter *iter, const gchar *id,
                                const FutileTree *tree)
{
  f90_dictionary_pointer root;
  
  root = iter->pointer;
  if (id)
    FC_FUNC_(bind_dict_insert, BIND_DICT_INSERT)(&iter->pointer, id, strlen(id));
  FC_FUNC_(bind_dict_update, BIND_DICT_UPDATE)(&iter->pointer, &tree->root);
  iter->pointer = root;
}
gboolean futile_tree_iter_pop(FutileTreeIter *iter, const gchar *key)
{
  int exists;
  gboolean reset;

  reset = (F_TYPE(iter->pointer) == F_TYPE(iter->tree->root));
  FC_FUNC_(bind_dict_pop, BIND_DICT_POP)(&iter->pointer, &exists, key, strlen(key));
  if (reset)
    iter->tree->root = iter->pointer;
  return (exists) ? TRUE : FALSE;
}
/**
 * futile_tree_value:
 * @dict: 
 *
 * 
 *
 * Returns: (transfer full):
 **/
gchar* futile_tree_iter_value(const FutileTreeIter *iter)
{
  char buf[max_field_length + 1];
  int i;
  guint ln;
  gchar *out;
  
  buf[max_field_length] = ' ';
  FC_FUNC_(bind_dict_value, BIND_DICT_VALUE)(&iter->pointer, buf, max_field_length);
  for (i = max_field_length; i >= 0 && buf[i] == ' '; i--)
    buf[i] = '\0';
  ln = max_field_length - i;
  out = g_malloc(sizeof(gchar) * (ln + 1));
  memcpy(out, buf, sizeof(gchar) * (ln + 1));
  return out;
}
/**
 * futile_tree_key:
 * @dict: 
 *
 * 
 *
 * Returns: (transfer full):
 **/
gchar* futile_tree_iter_key(const FutileTreeIter *iter)
{
  char buf[max_field_length + 1];
  int i;
  guint ln;
  gchar *out;
  
  buf[max_field_length] = ' ';
  FC_FUNC_(bind_dict_key, BIND_DICT_KEY)(&iter->pointer, buf, max_field_length);
  for (i = max_field_length; i >= 0 && buf[i] == ' '; i--)
    buf[i] = '\0';
  ln = max_field_length - i;
  out = g_malloc(sizeof(gchar) * (ln + 1));
  memcpy(out, buf, sizeof(gchar) * (ln + 1));
  return out;
}
guint futile_tree_iter_len(const FutileTreeIter *iter)
{
  int ln;

  FC_FUNC_(bind_dict_len, BIND_DICT_LEN)(&iter->pointer, &ln);
  if (ln > 0)
    return (guint)ln;
  FC_FUNC_(bind_dict_size, BIND_DICT_SIZE)(&iter->pointer, &ln);
  if (ln > 0)
    return (guint)ln;
  return 0;
}
void futile_tree_dump(FutileTree *tree, gint unit)
{
  FC_FUNC_(bind_dict_dump, BIND_DICT_DUMP)(&tree->root, &unit);
}
void futile_tree_dump_to_file(FutileTree *tree, const gchar *filename)
{
  FC_FUNC_(bind_dict_dump_to_file, BIND_DICT_DUMP_TO_FILE)(&tree->root, filename, strlen(filename));
}
/*********************************/
