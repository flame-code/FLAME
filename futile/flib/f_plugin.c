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

#define STR(A) str(A)
#define str(s) #s
#define FC_INIT FC_FUNC(init, INIT)

typedef void (*PluginInitFunc)(void);

void FC_FUNC_(plugin_load, PLUGIN_LOAD)(const char *name, int *ierr, unsigned int ln)
{
  char *libname;
  unsigned int i;
  void *handle;
  void *func;
  PluginInitFunc init;

  libname = malloc(sizeof(char) * (ln + 7));
  memcpy(libname, "./lib", sizeof(char) * 5);
  memcpy(libname + 5, name, sizeof(char) * ln);
  for (i = ln + 5; i > 4 && libname[i] == ' '; i--);
  memcpy(libname + i, ".so", sizeof(char) * 3);
  libname[i + 3] = '\0';

  handle = dlopen(libname, RTLD_LAZY);
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
  *(void **) (&init) = func;
  init();

  *ierr = 0;
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
