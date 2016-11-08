/**
 * \file
 *    Manual on how to use plugins.
 */

/**
 * \author
 *    Copyright (C) 2016 BigDFT group
 *    This file is distributed under the terms of the
 *    GNU General Public License, see ~/COPYING file
 *    or http://www.gnu.org/copyleft/gpl.txt .
 *    For the list of contributors, see ~/AUTHORS
 * \example test_plugin.f90
 *    Manual on how to use plugins.
 */

#include <stdlib.h>
#include <dlfcn.h>
#include <string.h>

#include <stdio.h>

#include "config.h"
#include "utils.h"

#define STR(A) str(A)
#define str(s) #s
#define FC_INIT FC_FUNC(init, INIT)

static unsigned int n_args = 0;
static void *args[FFUNC_MAX_ARGS] = { NULL };
static unsigned int strs[FFUNC_MAX_ARGS] = { 0 };

void FC_FUNC_(plugin_reset_arg, PLUGIN_RESET_ARG)()
{
  n_args = 0;
}
void FC_FUNC_(plugin_add_arg, PLUGIN_ADD_ARG)(void *arg)
{
  if (n_args >= FFUNC_MAX_ARGS)
    return;
  
  args[n_args++] = arg;
}
void FC_FUNC_(plugin_add_str, PLUGIN_ADD_STR)(void *str, unsigned int ln)
{
  if (n_args >= FFUNC_MAX_ARGS)
    return;

  args[n_args] = str;
  strs[n_args++] = ln;
}

void FC_FUNC_(plugin_load, PLUGIN_LOAD)(const char *name, int *ierr, unsigned int ln)
{
  char *libname;
  unsigned int i;
  void *handle;
  void *func;

  libname = malloc(sizeof(char) * (ln + 9));
  memcpy(libname, "./lib", sizeof(char) * 5);
  memcpy(libname + 5, name, sizeof(char) * ln);
  for (i = ln + 4; i > 4 && libname[i] == ' '; i--);
  memcpy(libname + ++i, ".so", sizeof(char) * 3);
  libname[i + 3] = '\0';

  handle = dlopen(libname, RTLD_NOW);
  free(libname);
  if (!handle)
    {
      *ierr = 1;
      return;
    }

  /* Try first Fortran init version. */
  func = dlsym(handle, STR(FC_INIT));
  if (!func)
    {
      /* Try then C init version. */
      func = dlsym(handle, "init");
      if (!func)
        {
          *ierr = 2;
          return;
        }
    }

  *ierr = 0;
  f_func_call(func, n_args, args, strs);
}

void FC_FUNC_(plugin_error, PLUGIN_ERROR)(char errMess[256])
{
  unsigned int i;
  char *mess;

  memset(errMess, ' ', sizeof(char) * 256);

  mess = dlerror();
  for (i = 0; mess[i] && i < 256; i++)
    errMess[i] = mess[i];
}
