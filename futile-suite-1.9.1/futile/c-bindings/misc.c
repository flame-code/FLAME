#include <config.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>

#include "misc.h"

void FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)(void);
void FC_FUNC_(f_lib_finalize, F_LIB_FINALIZE)(void);
void FC_FUNC_(f_object_get_method, F_OBJECT_GET_METHOD)(const char *obj_id, const char *meth_id, int *n_args, int *isfunc, void **callback, int ln_obj_id, int ln_meth_id);
void FC_FUNC_(f_python_ndarray_init, F_PYTHON_NDARRAY_INIT)(void);
void FC_FUNC_(f_python_ndarray_get, F_PYTHON_NDARRAY_GET)(void *ndarray, void *data, int *ndims, int *shapes, char *type, int ln);

void futile_initialize()
{
  FC_FUNC_(f_lib_initialize, F_LIB_INITIALIZE)();
}
void futile_finalize()
{
  FC_FUNC_(f_lib_finalize, F_LIB_FINALIZE)();
}

gboolean futile_object_get_method(FutileMethod *meth,
                                  const char *obj_id, const char *meth_id)
{
  void *callback;

  if (!meth_id || !meth_id[0])
    return FALSE;

  if (obj_id && obj_id[0])
    FC_FUNC_(f_object_get_method, F_OBJECT_GET_METHOD)(obj_id, meth_id, (int*)&meth->n_args, &meth->isfunc, &callback, strlen(obj_id), strlen(meth_id));
  else
    FC_FUNC_(f_object_get_method, F_OBJECT_GET_METHOD)("class", meth_id, (int*)&meth->n_args, &meth->isfunc, &callback, strlen("class"), strlen(meth_id));
  if (!callback)
    return FALSE;
  
  meth->callback = (FutileMethodFortranFunc)callback;

  memset(meth->args, '\0', sizeof(meth->args));
  meth->n_strs = 0;
  memset(meth->strlens, '\0', sizeof(meth->strlens));
  return TRUE;
}

static FutileArg* _add_arg(FutileMethod *meth, void *arg)
{
  unsigned int i;

  for (i = 0; i < FUTILE_METHOD_ARG_MAX; i++)
    if (!meth->args[i].arg)
      {
        meth->args[i].arg = arg;
        meth->args[i].type = FUTILE_POINTER;
        meth->args[i].size = 1;
        return meth->args + i;
      }
  return NULL;
}

void futile_object_method_add_arg(FutileMethod *meth, void *arg)
{
  _add_arg(meth, arg);
}
void futile_object_method_add_arg_str(FutileMethod *meth, char *arg)
{
  _add_arg(meth, arg);
  meth->strlens[meth->n_strs++] = strlen(arg);
}
void futile_object_method_add_arg_full(FutileMethod *meth, void *arg,
                                       FutileNumeric type, size_t size,
                                       void *buffer, FutileDestroyFunc freeBuffer)
{
  FutileArg *storage;

  storage = _add_arg(meth, arg);
  if (!storage)
    return;

  storage->type = type;
  storage->size = size;
  storage->data = buffer;
  storage->freefunc = freeBuffer;
}

void* futile_object_method_get_arg_arr(FutileMethod *meth, unsigned int i)
{
  if (i >= FUTILE_METHOD_ARG_MAX)
    return NULL;

  return meth->args[i].arg;
}

void futile_object_method_clean(FutileMethod *meth)
{
  int i;

  for (i = 0; i < meth->n_args; i++)
    if (meth->args[i].freefunc && meth->args[i].data)
      meth->args[i].freefunc(meth->args[i].data);
}

void futile_object_method_execute(FutileMethod *meth)
{
  switch (meth->n_args)
    {
    case (0):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback();
          break;
        case (1):
          if (meth->callback) meth->callback(meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (1):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0].arg);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (2):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (3):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->args[2].arg);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->args[2].arg, meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->args[2].arg, meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    case (4):
      switch (meth->n_strs)
        {
        case (0):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->args[2].arg, meth->args[3].arg);
          break;
        case (1):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->args[2].arg, meth->args[3].arg, meth->strlens[0]);
          break;
        case (2):
          if (meth->callback) meth->callback(meth->args[0].arg, meth->args[1].arg, meth->args[2].arg, meth->args[3].arg, meth->strlens[0], meth->strlens[1]);
          break;
        default:
          break;
        }
      break;
    default:
      fprintf(stderr, "Implement %d arguments.\n", meth->n_args);
      break;
    }
}

void* futile_object_ndarray_new(void **paddress)
{
  FutileMethod meth;
  void *address;

  if (!futile_object_get_method(&meth, NULL, "ndarray_new"))
    {
      FC_FUNC_(f_python_ndarray_init, F_PYTHON_NDARRAY_INIT)();
      if (!futile_object_get_method(&meth, NULL, "ndarray_new"))
        return NULL;
    }

  meth.n_args += 1;
  *paddress = malloc(sizeof(void*) * FUTILE_F_POINTER_SIZE);
  futile_object_method_add_arg(&meth, &paddress);
  futile_object_method_add_arg(&meth, &address);
  futile_object_method_execute(&meth);
  return address;
}

void futile_object_ndarray_free(void *paddress)
{
  FutileMethod meth;

  if (!futile_object_get_method(&meth, NULL, "ndarray_free"))
    return;

  meth.n_args += 1;
  futile_object_method_add_arg(&meth, &paddress);
  futile_object_method_execute(&meth);
  free(paddress);
}

void* futile_object_ndarray_get(void *ndarray,
                                int *ndims, int shapes[7], FutileNumeric *type)
{
  char t[2];
  void *data;

  FC_FUNC_(f_python_ndarray_get, F_PYTHON_NDARRAY_GET)(ndarray, &data, ndims, shapes, t, 2);
  if (t[0] == 'i' && t[1] == '4')
    *type = FUTILE_INTEGER_4;
  else if (t[0] == 'f' && t[1] == '8')
    *type = FUTILE_REAL_8;
  else
    return NULL;
  return data;
}
