/*
!> @file
!!  Routines for the overriding of the malloc function to be used for debug cases
!! @author
!!    Copyright (C) 2007-2013 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
*/
#include <config.h>
#define _GNU_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <dlfcn.h>

#define FC_STR(str) #str
#define FC_QUOTE(str) FC_STR(str)
#define FC_SYMBOL(sym,SYM) FC_QUOTE(FC_FUNC(sym,SYM))

/* #define RTLD_NEXT RTLD_DEFAULT (this has been put for windows) */

static void* (*real_malloc)(size_t)=NULL;

size_t get_padding(void);

static void mtrace_init(void)
{
    real_malloc = dlsym(RTLD_NEXT, "malloc");
    if (NULL == real_malloc) {
        fprintf(stderr, "Error in `dlsym`: %s\n", dlerror());
    }
}

void *malloc(size_t size)
{
    if(real_malloc==NULL) {
        mtrace_init();
    }

    void *p = NULL;
    size_t repad=0;
    char *repad_env=getenv("FUTILE_REPAD_SIZE_INTERNAL");
    /*printf("Got environment again:%s \n",repad_env);*/
    if (NULL != repad_env){
      int nvars = sscanf(repad_env,"%zu",&repad);
      if (nvars != 1) {
        fprintf(stderr, "Error in `scanning`: %d\n", nvars);
      }
      /*fprintf(stderr, "getval(%d) ", (int) repad);*/
    }

    size_t actual_size = size + repad;
    if (repad != 0) {
      fprintf(stderr, "malloc(%d) ", (int) size);
      fprintf(stderr, " [padding(%d)] = ", (int) repad);
    }
    p = real_malloc(actual_size);
    if (repad != 0 ) fprintf(stderr, "%p\n", p);
    return p;
}
