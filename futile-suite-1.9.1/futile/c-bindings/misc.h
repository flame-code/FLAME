#ifndef MISC_H
#define MISC_H

#include "futile_cst.h"

void futile_initialize();
void futile_finalize();

#define FUTILE_F_POINTER_SIZE 64

#define FUTILE_METHOD_ARG_MAX 32

typedef enum _FutileNumeric
  {
    FUTILE_POINTER,
    FUTILE_INTEGER_4,
    FUTILE_REAL_8
  } FutileNumeric;

typedef union
{
  int ival;
  double dval;
} FutileNumericValue;

typedef void (*FutileDestroyFunc)(void *data);
typedef void (*FutileMethodFortranFunc)();

typedef void (*FObjectCallable)(void);
typedef struct _FObjectKernel
{
  char kernel[256];
} FObjectKernel;

typedef struct _FutileArg FutileArg;
struct _FutileArg
{
  FutileNumeric type;
  size_t size;

  void *arg;

  void *data;
  FutileDestroyFunc freefunc;
};

typedef struct _FutileMethod FutileMethod;
struct _FutileMethod
{
  unsigned int n_args;
  unsigned int n_strs;
  int isfunc;

  FutileMethodFortranFunc callback;

  FutileArg args[FUTILE_METHOD_ARG_MAX];
  int strlens[FUTILE_METHOD_ARG_MAX];
};

void futileC_object_kernel_new(FObjectKernel *kernel, FObjectCallable callback, unsigned int nargs);
void futileC_object_kernel_add_arg(FObjectKernel *kernel, void *arg);
int futileC_object_signal_connect(const char *obj, const char *sig, FObjectKernel *kernel);

gboolean futile_object_get_method(FutileMethod *meth,
                                  const char *obj_id, const char *meth_id);
void futile_object_method_add_arg(FutileMethod *meth, void *arg);
void futile_object_method_add_arg_str(FutileMethod *meth, char *arg);
void futile_object_method_add_arg_full(FutileMethod *meth, void *arg,
                                       FutileNumeric type, size_t size,
                                       void *buffer, FutileDestroyFunc freeBuffer);
void futile_object_method_execute(FutileMethod *meth);
void* futile_object_method_get_arg_arr(FutileMethod *meth, unsigned int i);
void futile_object_method_clean(FutileMethod *meth);

void* futile_object_ndarray_new(void **paddress);
void futile_object_ndarray_free(void *paddress);
void* futile_object_ndarray_get(void *ndarray,
                                int *ndims, int shapes[7], FutileNumeric *type);
#endif
