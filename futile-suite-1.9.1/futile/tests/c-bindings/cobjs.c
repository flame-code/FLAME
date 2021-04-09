#include "config.h"
#include <stdlib.h>
#include <stdio.h>

#include "misc.h"
#include "futileC_objects.h"


void FC_FUNC_(obj_print, OBJ_PRINT)(void *obj);
void FC_FUNC_(obj_set, OBJ_SET)(void *obj, int *val);

void onInit(void *obj, int *val0)
{
  FC_FUNC_(obj_set, OBJ_SET)(obj, val0);
  free(val0);
}

void onSet(void *obj)
{
  FC_FUNC_(obj_print, OBJ_PRINT)(obj);
}

void FC_FUNC_(connect_in_c, CONNECT_IN_C)()
{
  FObjectKernel kernel;
  int *val0;

  val0 = malloc(sizeof(int));
  *val0 = 123;
  futileC_object_kernel_new(&kernel, (FObjectCallable)onInit, 2);
  futileC_object_kernel_add_arg(&kernel, val0);
  futileC_object_signal_connect("my_object", "init", &kernel);

  futileC_object_kernel_new(&kernel, (FObjectCallable)onSet, 1);
  futileC_object_signal_connect("my_object", "set", &kernel);
}
