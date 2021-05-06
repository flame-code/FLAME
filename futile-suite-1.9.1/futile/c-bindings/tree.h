#ifndef TREE_H
#define TREE_H

#include "futile_cst.h"

#include "dict.h"

/******************************/
/* FutileTree data structure */
/******************************/
#ifdef GLIB_MAJOR_VERSION
#define FUTILE_TREE_TYPE    (futile_tree_get_type())
#define FUTILE_TREE(obj)                                               \
  (G_TYPE_CHECK_INSTANCE_CAST(obj, FUTILE_TREE_TYPE, FutileTree))
#define FUTILE_TREE_CLASS(klass)                                       \
  (G_TYPE_CHECK_CLASS_CAST(klass, FUTILE_TREE_TYPE, FutileTreeClass))
#define FUTILE_TREE_GET_CLASS(obj)                                     \
  (G_TYPE_INSTANCE_GET_CLASS(obj, FUTILE_TREE_TYPE, FutileTreeClass))
#define FUTILE_IS_CLASS_TREE(klass)                    \
  (G_TYPE_CHECK_CLASS_TYPE(klass, FUTILE_TREE_TYPE))
#define FUTILE_IS_TYPE_TREE(obj)                       \
  (G_TYPE_CHECK_INSTANCE_TYPE(obj, FUTILE_TREE_TYPE))

typedef struct _FutileTreeClass FutileTreeClass;
struct _FutileTreeClass
{
  GObjectClass parent;
};
GType futile_tree_get_type(void);
#else
#define FUTILE_TREE_TYPE    (999)
#define FUTILE_TREE(obj)    ((FutileTree*)obj)
#endif

typedef struct _FutileTree FutileTree;
struct _FutileTree
{
  /* Object management. */
  GObject parent;
  gboolean dispose_has_run;

  /* Private. */
  f90_dictionary_pointer root;
};
typedef struct _FutileTreeIter FutileTreeIter;
struct _FutileTreeIter
{
  FutileTree *tree;
  f90_dictionary_pointer pointer;
};

FutileTree *futile_tree_new      (FutileTreeIter *root);
FutileTree *futile_tree_new_from_yaml(const gchar *buf, FutileTreeIter *root);
void  futile_tree_unref          (FutileTree *tree);


void  futile_tree_iter_insert    (const FutileTreeIter *at, const gchar *key,
                                  FutileTreeIter *leaf);
void  futile_tree_iter_append    (const FutileTreeIter *at, FutileTreeIter *leaf);
gboolean futile_tree_iter_at_key (const FutileTreeIter *at, const gchar *key,
                                  FutileTreeIter *leaf);
gboolean futile_tree_iter_at_item(const FutileTreeIter *at, guint id,
                                  FutileTreeIter *leaf);
gboolean futile_tree_iter_first  (const FutileTreeIter *at, FutileTreeIter *first);
gboolean futile_tree_iter_next   (FutileTreeIter *iter);

void  futile_tree_iter_add       (FutileTreeIter *iter, const gchar *value);
void  futile_tree_iter_set       (FutileTreeIter *iter,
                                  const gchar *id, const gchar *value);
void  futile_tree_iter_set_array (FutileTreeIter *iter,
                                  const gchar *id, const gchar **value);
void  futile_tree_iter_set_tree  (FutileTreeIter *iter, const gchar *id,
                                  const FutileTree *tree);

gboolean futile_tree_iter_pop    (FutileTreeIter *iter, const gchar *key);

gchar* futile_tree_iter_key      (const FutileTreeIter *iter);
gchar* futile_tree_iter_value    (const FutileTreeIter *iter);

guint futile_tree_iter_len       (const FutileTreeIter *iter);

void  futile_tree_dump           (FutileTree *tree, gint unit);
void  futile_tree_dump_to_file   (FutileTree *tree, const gchar *filename);
/*********************************/

#endif
