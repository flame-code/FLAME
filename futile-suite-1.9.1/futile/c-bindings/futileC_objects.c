/**
 * \file
 *    C bindings for the f_object API.
 */

/**
 * \author
 *    Copyright (C) 2016 BigDFT group
 *    This file is distributed under the terms of the
 *    GNU General Public License, see ~/COPYING file
 *    or http://www.gnu.org/copyleft/gpl.txt .
 *    For the list of contributors, see ~/AUTHORS
 */

#include "config.h"
#include "misc.h"

#include <string.h>

typedef void* FAddress;

void FC_FUNC_(f_object_kernel_new,
              F_OBJECT_KERNEL_NEW)(FObjectKernel *kernel, FAddress callback, int *nargs);
void FC_FUNC_(f_object_kernel_add_arg,
              F_OBJECT_KERNEL_ADD_ARG)(FObjectKernel *kernel, void *arg);
void FC_FUNC_(f_object_signal_connect_bind,
              F_OBJECT_SIGNAL_CONNECT_BIND)(const char *obj, const char *sig, FObjectKernel *kernel, int *id, int obj_len, int sig_len);

void futileC_object_kernel_new(FObjectKernel *kernel, FObjectCallable callback, unsigned int nargs)
{
  FC_FUNC_(f_object_kernel_new,
           F_OBJECT_KERNEL_NEW)(kernel, (FAddress)callback, (int*)&nargs);
}

void futileC_object_kernel_add_arg(FObjectKernel *kernel, void *arg)
{
  FC_FUNC_(f_object_kernel_add_arg,
           F_OBJECT_KERNEL_ADD_ARG)(kernel, arg);
}

int futileC_object_signal_connect(const char *obj, const char *sig, FObjectKernel *kernel)
{
  int id;
  FC_FUNC_(f_object_signal_connect_bind,
           F_OBJECT_SIGNAL_CONNECT_BIND)(obj, sig, kernel, &id, strlen(obj), strlen(sig));
  return id;
}
