#include <Python.h>
#include <structmember.h>

#include <misc.h>

#include "config.h"
#ifdef HAVE_PYTHON_NUMPY
#include <numpy/ndarrayobject.h>
#endif

typedef enum
  {
    F_PY_SCALAR_I4,
    F_PY_SCALAR_R8,
    F_PY_ARRAY
  } FPyTypeFlags;

typedef struct {
  PyObject_HEAD
  FPyTypeFlags flags;
} FPyArgFlags;

static PyMemberDef f_py_argflags_members[] = {
  {NULL}  /* Sentinel */
};

static PyMethodDef f_py_argflags_methods[] = {
  {NULL}  /* Sentinel */
};

static void f_py_argflags_dealloc(FPyArgFlags *self)
{
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* f_py_argflags_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  FPyArgFlags *self;

  self = (FPyArgFlags*)type->tp_alloc(type, 0);

  return (PyObject*)self;
}

static int f_py_argflags_init(FPyArgFlags *self, PyObject *args, PyObject *kwds)
{
  int type;

  static char *kwlist[] = {"type", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "i", kwlist, &type))
    return -1;

  if (type < F_PY_SCALAR_I4 || type > F_PY_ARRAY)
    {
      PyErr_SetString(PyExc_TypeError, "unknown flag");
      return -1;
    }

  self->flags = type;

  return 0;
}

static PyTypeObject FPyArgFlagsType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "futile.ArgFlags",         /*tp_name*/
  sizeof(FPyArgFlags),       /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)f_py_argflags_dealloc, /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  "argument flags",          /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  f_py_argflags_methods,     /* tp_methods */
  f_py_argflags_members,     /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)f_py_argflags_init,/* tp_init */
  0,                         /* tp_alloc */
  f_py_argflags_new,           /* tp_new */
};

static int f_py_argflags_check(PyObject *obj, FPyTypeFlags type)
{
  return (obj->ob_type == &FPyArgFlagsType &&
          ((FPyArgFlags*)obj)->flags == type);
}

static void* f_py_tuple_to_arr(PyObject *tuple, FutileNumeric *type, Py_ssize_t *size)
{
  Py_ssize_t i;
  PyObject *first;
  int *iarr;
  double *darr;

  *type = FUTILE_INTEGER_4;

  *size = PyTuple_GET_SIZE(tuple);
  if (*size == 0)
    return NULL;

  first = PyTuple_GET_ITEM(tuple, 0);
  if (PyInt_Check(first))
    {
      *type = FUTILE_INTEGER_4;
      iarr = malloc(sizeof(int) * *size); // should be long here.
      for (i = 0; i < *size; i++)
        iarr[i] = (int)PyInt_AS_LONG(PyTuple_GET_ITEM(tuple, i));
      return iarr;
    }
  else if (PyFloat_Check(first))
    {
      *type = FUTILE_REAL_8;
      darr = malloc(sizeof(double) * *size);
      for (i = 0; i < *size; i++)
        darr[i] = (double)PyFloat_AS_DOUBLE(PyTuple_GET_ITEM(tuple, i));
      return darr;
    }
  else
    return NULL;
}
#ifdef HAVE_PYTHON_NUMPY
static void* f_py_array_to_arr(PyObject *array, FutileNumeric *type, Py_ssize_t *size, PyObject **farr)
{
  PyArray_Descr *descr;
  npy_intp *indices;
  void *arr;

  *type = FUTILE_INTEGER_4;
  *size = 0;
  *farr = PyArray_FromAny(array, NULL, 1, 7, NPY_INOUT_FARRAY, NULL);
  if (*farr == NULL)
    return NULL;

  indices = malloc(sizeof(npy_intp) * PyArray_NDIM(*farr));
  memset(indices, '\0', sizeof(npy_intp) * PyArray_NDIM(*farr));
  arr = PyArray_GetPtr((PyArrayObject*)*farr, indices);
  free(indices);

  *size = PyArray_Size(*farr);
  descr = PyArray_DESCR(*farr);
  switch (descr->kind)
    {
    case 'i':
      if (descr->elsize != 4)
        {
          Py_DECREF(*farr);
          return NULL;
        }
      *type = FUTILE_INTEGER_4;
      break;
    case 'f':
      if (descr->elsize != 8)
        {
          Py_DECREF(*farr);
          return NULL;
        }
      *type = FUTILE_REAL_8;
      break;
    default:
      PyErr_SetString(PyExc_TypeError, "unsupported type");
      Py_DECREF(*farr);
      return NULL;
    }

  return arr;
}
static PyObject* f_py_arr_to_array(void *data, int ndims, int shapes[7], int type)
{
  npy_intp dims[7];
  int type_num;
  int i;
  
  switch (type)
    {
    case (FUTILE_INTEGER_4):
      type_num = NPY_INT32;
      break;
    case (FUTILE_REAL_8):
      type_num = NPY_DOUBLE;
      break;
    default:
      PyErr_SetString(PyExc_TypeError, "unsupported type");
      return NULL;
    }
  for (i = 0; i < ndims; i++)
    dims[i] = shapes[ndims - i - 1];
  return PyArray_SimpleNewFromData(ndims, dims, type_num, data);
}
#endif

typedef struct {
  PyObject_HEAD
  PyObject *obj_id;
  void *address;
  void *paddress;
} FPyObject;

static int f_py_object_check(PyObject *o);

typedef struct {
  PyObject_HEAD
  FPyObject *self;
  FutileMethod meth;
} FPyMethod;


static PyMemberDef f_py_method_members[] = {
  {NULL}  /* Sentinel */
};

static PyMethodDef f_py_method_methods[] = {
  {NULL}  /* Sentinel */
};

static void f_py_method_dealloc(FPyMethod *self)
{
  Py_XDECREF(self->self);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* f_py_method_new_from_meth(PyTypeObject *type,
                                           const FutileMethod *meth,
                                           FPyObject *obj)
{
  FPyMethod *self;

  self = (FPyMethod*)type->tp_alloc(type, 0);
  self->meth = *meth;
  Py_XINCREF(obj);
  self->self = obj;
  if (obj != NULL)
    self->meth.n_args += 1;

  return (PyObject*)self;
}

static PyObject* f_py_method_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  FPyMethod *self;

  self = (FPyMethod*)type->tp_alloc(type, 0);

  return (PyObject*)self;
}

static int f_py_method_init(FPyMethod *self, PyObject *args, PyObject *kwds)
{
  const char *meth_id;
  FutileMethod meth;

  static char *kwlist[] = {"meth_id", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &meth_id))
    return -1;

  if (!futile_object_get_method(&meth, NULL, meth_id))
    {
      PyErr_SetString(PyExc_TypeError, "no method");
      return -1;
    }

  self->meth = meth;

  return 0;
}

#ifdef HAVE_PYTHON_NUMPY
static void pyFree(void *obj)
{
  Py_DECREF(obj);
}
#endif

static PyObject* f_py_method_call(FPyMethod *self, PyObject *args, PyObject *kwds)
{
  PyObject *iterator = PyObject_GetIter(args);
  PyObject *item, *ret;
#ifdef HAVE_PYTHON_NUMPY
  PyObject *farr;
#endif

  int i, j, i_args, ndims, shapes[7];
  Py_ssize_t size;
  FutileNumeric type;
  void *arr, *parr;

  if (self->self && self->self->address)
    futile_object_method_add_arg(&self->meth, self->self->address);

  if (iterator == NULL)
    return NULL;

  i_args = (self->self && self->self->address) ? 1 : 0;
  for (i = i_args; i < self->meth.n_args; i += 1)
    {
      item = PyIter_Next(iterator);
      if (item == NULL)
        {
          PyErr_SetString(PyExc_TypeError, "missing arguments");
          Py_DECREF(iterator);
          return NULL;
        }

      /* fprintf(stderr, "%s\n", PyString_AsString(PyObject_Str(PyObject_Type(item)))); */
      if (PyInt_Check(item) || f_py_argflags_check(item, F_PY_SCALAR_I4))
        {
          arr = malloc(sizeof(int));
          if (PyInt_Check(item))
            *(int*)arr = (int)PyInt_AS_LONG(item);
          futile_object_method_add_arg_full(&self->meth, arr, FUTILE_INTEGER_4,
                                            1, arr, free);
        }
      else if (PyFloat_Check(item) || f_py_argflags_check(item, F_PY_SCALAR_R8))
        {
          arr = malloc(sizeof(double));
          if (PyFloat_Check(item))
            *(double*)arr = PyFloat_AS_DOUBLE(item);
          futile_object_method_add_arg_full(&self->meth, arr, FUTILE_REAL_8,
                                            1, arr, free);
        }
      else if (f_py_argflags_check(item, F_PY_ARRAY))
        {
          arr = futile_object_ndarray_new(&parr);
          if (arr == NULL)
            {
              Py_DECREF(item);
              Py_DECREF(iterator);
              futile_object_method_clean(&self->meth);
              return NULL;
            }
          futile_object_method_add_arg_full(&self->meth, arr, FUTILE_POINTER,
                                            1, parr, futile_object_ndarray_free);
        }
      else if (PyString_Check(item))
        {
          futile_object_method_add_arg_str(&self->meth, PyString_AsString(item));
        }
      else if (PyTuple_Check(item))
        {
          arr = f_py_tuple_to_arr(item, &type, &size);
          futile_object_method_add_arg_full(&self->meth, arr,
                                            type, size, arr, free);
        }
#ifdef HAVE_PYTHON_NUMPY
      else if (PyArray_Check(item))
        {
          arr = f_py_array_to_arr(item, &type, &size, &farr);
          futile_object_method_add_arg_full(&self->meth, arr,
                                            type, size, farr, pyFree);
        }
#endif
      else if (f_py_object_check(item))
        {
          futile_object_method_add_arg(&self->meth, ((FPyObject*)item)->address);
        }
      else
        {
          PyErr_SetString(PyExc_TypeError, "waiting for type futile.FObject");
          Py_DECREF(item);
          Py_DECREF(iterator);
          futile_object_method_clean(&self->meth);
          return NULL;
        }
      /* release reference when done */
      Py_DECREF(item);
    }
  Py_DECREF(iterator);

  futile_object_method_execute(&self->meth);

  /* Convert back values. */
  ret = PyTuple_New(self->meth.n_args);
  j = 0;
  iterator = PyObject_GetIter(args);
  for (i = i_args; i < self->meth.n_args; i += 1)
    {
      item = PyIter_Next(iterator);

      if (PyInt_Check(item) || f_py_argflags_check(item, F_PY_SCALAR_I4))
        {
          arr = futile_object_method_get_arg_arr(&self->meth, i);
          if (f_py_argflags_check(item, F_PY_SCALAR_I4) ||
              PyInt_AS_LONG(item) != (long)*(int*)arr)
            PyTuple_SET_ITEM(ret, j++, PyInt_FromLong((long)*(int*)arr));
        }
      else if (PyFloat_Check(item) || f_py_argflags_check(item, F_PY_SCALAR_R8))
        {
          arr = futile_object_method_get_arg_arr(&self->meth, i);
          if (f_py_argflags_check(item, F_PY_SCALAR_R8) ||
              PyFloat_AS_DOUBLE(item) != *(double*)arr)
            PyTuple_SET_ITEM(ret, j++, PyFloat_FromDouble(*(double*)arr));
        }
      else if (f_py_argflags_check(item, F_PY_ARRAY))
        {
          parr = futile_object_method_get_arg_arr(&self->meth, i);
          arr = futile_object_ndarray_get(parr, &ndims, shapes, &type);
#ifdef HAVE_PYTHON_NUMPY
          farr = f_py_arr_to_array(arr, ndims, shapes, type);
          if (farr != NULL)
            PyTuple_SET_ITEM(ret, j++, farr);
#endif
        }
      Py_DECREF(item);
    }
  Py_DECREF(iterator);

  futile_object_method_clean(&self->meth);

  if (j > 0)
    {
      _PyTuple_Resize(&ret, j);
      return ret;
    }
  else
    {
      Py_DECREF(ret);
      Py_RETURN_NONE;
    }
}

static PyTypeObject FPyMethodType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "futile.FMethod",          /*tp_name*/
  sizeof(FPyMethod),         /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)f_py_method_dealloc, /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  (ternaryfunc)f_py_method_call, /*tp_call*/
  0,                         /*tp_str*/
  0,                         /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  "Fortran methods",         /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  f_py_method_methods,       /* tp_methods */
  f_py_method_members,       /* tp_members */
  0,                         /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)f_py_method_init,/* tp_init */
  0,                         /* tp_alloc */
  f_py_method_new,           /* tp_new */
};

static PyMemberDef f_py_object_members[] = {
  {NULL}  /* Sentinel */
};

static PyObject* f_py_object_get_obj_id(FPyObject *self, void *closure)
{
  Py_INCREF(self->obj_id);
  return self->obj_id;
}

static int f_py_object_set_obj_id(FPyObject *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the obj_id attribute");
    return -1;
  }
  
  if (!PyString_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The obj_id attribute value must be a string");
    return -1;
  }
      
  Py_XDECREF(self->obj_id);
  Py_INCREF(value);
  self->obj_id = value;    

  return 0;
}

static PyObject* f_py_object_get_address(FPyObject *self, void *closure)
{
  return PyLong_FromVoidPtr(self->address);
}

static int f_py_object_set_address(FPyObject *self, PyObject *value, void *closure)
{
  if (value == NULL) {
    PyErr_SetString(PyExc_TypeError, "Cannot delete the address attribute");
    return -1;
  }
  
  if (!PyLong_Check(value) && !PyInt_Check(value)) {
    PyErr_SetString(PyExc_TypeError, 
                    "The address attribute value must be an address");
    return -1;
  }
      
  self->address = PyLong_AsVoidPtr(value);

  return 0;
}
static PyGetSetDef f_py_object_getseters[] = {
  {"obj_id", 
   (getter)f_py_object_get_obj_id, (setter)f_py_object_set_obj_id,
   "object id",
   NULL},
  {"address", 
   (getter)f_py_object_get_address, (setter)f_py_object_set_address,
   "Fortran object address",
   NULL},
  {NULL}  /* Sentinel */
};

static PyObject* f_py_object_getattro(FPyObject *self, PyObject *attr_name)
{
  FutileMethod meth;
  const char *obj_id, *meth_id;

  obj_id = PyString_AsString(self->obj_id);
  meth_id = PyString_AsString(attr_name);

  /* Look for a dynamic method as exported by Fortran. */
  if (futile_object_get_method(&meth, obj_id, meth_id))
    return f_py_method_new_from_meth(&FPyMethodType, &meth, self);
  else
    return PyObject_GenericGetAttr((PyObject*)self, attr_name);
}

static PyMethodDef f_py_object_methods[] = {
  {NULL}  /* Sentinel */
};

static void f_py_object_dealloc(FPyObject *self)
{
  FutileMethod meth;

  if (self->paddress &&
      futile_object_get_method(&meth, PyString_AsString(self->obj_id), "destructor"))
    {
      meth.n_args += 1;
      futile_object_method_add_arg(&meth, &self->paddress);
      futile_object_method_execute(&meth);
      free(self->paddress);
    }
    
  Py_XDECREF(self->obj_id);
  self->ob_type->tp_free((PyObject*)self);
}

static PyObject* f_py_object_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
  FPyObject *self;

  self = (FPyObject*)type->tp_alloc(type, 0);
  if (self != NULL) {
    self->obj_id = PyString_FromString("");
    if (self->obj_id == NULL)
      {
        Py_DECREF(self);
        return NULL;
      }
        
    self->address = (void*)0;
  }

  return (PyObject*)self;
}

static int f_py_object_init(FPyObject *self, PyObject *args, PyObject *kwds)
{
  PyObject *obj_id = NULL, *address = NULL;
  FutileMethod meth;

  static char *kwlist[] = {"obj_id", "address", NULL};

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O|O", kwlist, &obj_id, &address))
    return -1; 

  if (f_py_object_set_obj_id(self, obj_id, (void*)0) < 0)
    return -1;
  if (address && f_py_object_set_address(self, address, (void*)0) < 0)
    return -1;
  else if (!address &&
           futile_object_get_method(&meth, PyString_AsString(self->obj_id), "constructor"))
    {
      meth.n_args += 1;
      self->paddress = malloc(sizeof(void) * FUTILE_F_POINTER_SIZE);
      futile_object_method_add_arg(&meth, &self->paddress);
      futile_object_method_add_arg(&meth, &self->address);
      futile_object_method_execute(&meth);
    }

  return 0;
}

static PyTypeObject FPyObjectType = {
  PyObject_HEAD_INIT(NULL)
  0,                         /*ob_size*/
  "futile.FObject",          /*tp_name*/
  sizeof(FPyObject),         /*tp_basicsize*/
  0,                         /*tp_itemsize*/
  (destructor)f_py_object_dealloc, /*tp_dealloc*/
  0,                         /*tp_print*/
  0,                         /*tp_getattr*/
  0,                         /*tp_setattr*/
  0,                         /*tp_compare*/
  0,                         /*tp_repr*/
  0,                         /*tp_as_number*/
  0,                         /*tp_as_sequence*/
  0,                         /*tp_as_mapping*/
  0,                         /*tp_hash */
  0,                         /*tp_call*/
  0,                         /*tp_str*/
  (getattrofunc)f_py_object_getattro, /*tp_getattro*/
  0,                         /*tp_setattro*/
  0,                         /*tp_as_buffer*/
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, /*tp_flags*/
  "Fortran objects",         /* tp_doc */
  0,                         /* tp_traverse */
  0,                         /* tp_clear */
  0,                         /* tp_richcompare */
  0,                         /* tp_weaklistoffset */
  0,                         /* tp_iter */
  0,                         /* tp_iternext */
  f_py_object_methods,       /* tp_methods */
  f_py_object_members,       /* tp_members */
  f_py_object_getseters,     /* tp_getset */
  0,                         /* tp_base */
  0,                         /* tp_dict */
  0,                         /* tp_descr_get */
  0,                         /* tp_descr_set */
  0,                         /* tp_dictoffset */
  (initproc)f_py_object_init,/* tp_init */
  0,                         /* tp_alloc */
  f_py_object_new,           /* tp_new */
};

static int f_py_object_check(PyObject *o)
{
  return PyObject_IsInstance(o, (PyObject*)&FPyObjectType);
}

static PyMethodDef FutilePyMethods[] = {
  {NULL, NULL, 0, NULL}        /* Sentinel */
};

static PyObject* f_py_futile_getattro(PyObject *self, PyObject *attr_name)
{
  FutileMethod meth;
  const char *meth_id;

  if (strcmp(PyModule_GetName(self), "futile"))
    return PyObject_GenericGetAttr(self, attr_name);

  meth_id = PyString_AsString(attr_name);

  /* Look for a dynamic method as exported by Fortran. */
  if (futile_object_get_method(&meth, NULL, meth_id))
    return f_py_method_new_from_meth(&FPyMethodType, &meth, NULL);
  else
    return PyObject_GenericGetAttr(self, attr_name);
}


PyMODINIT_FUNC
initfutile(void)
{
  PyObject *futile;
  PyTypeObject *futileType;
  FPyArgFlags *flag;

  if (PyType_Ready(&FPyMethodType) < 0)
    return;
  if (PyType_Ready(&FPyObjectType) < 0)
    return;
  if (PyType_Ready(&FPyArgFlagsType) < 0)
    return;

  futile = Py_InitModule("futile", FutilePyMethods);
  if (futile == NULL)
    return;

#ifdef HAVE_PYTHON_NUMPY
  import_array();
#endif

  Py_INCREF(&FPyMethodType);
  PyModule_AddObject(futile, "FMethod", (PyObject*)&FPyMethodType);
  Py_INCREF(&FPyObjectType);
  PyModule_AddObject(futile, "FObject", (PyObject*)&FPyObjectType);
  Py_INCREF(&FPyArgFlagsType);
  PyModule_AddObject(futile, "ArgFlags", (PyObject*)&FPyArgFlagsType);

  flag = PyObject_New(FPyArgFlags, &FPyArgFlagsType);
  flag->flags = F_PY_SCALAR_I4;
  PyModule_AddObject(futile, "SCALAR_I4", (PyObject*)flag);
  flag = PyObject_New(FPyArgFlags, &FPyArgFlagsType);
  flag->flags = F_PY_SCALAR_R8;
  PyModule_AddObject(futile, "SCALAR_R8", (PyObject*)flag);
  flag = PyObject_New(FPyArgFlags, &FPyArgFlagsType);
  flag->flags = F_PY_ARRAY;
  PyModule_AddObject(futile, "ARRAY", (PyObject*)flag);

  futileType = futile->ob_type;
  futileType->tp_getattro = (getattrofunc)f_py_futile_getattro;
}
