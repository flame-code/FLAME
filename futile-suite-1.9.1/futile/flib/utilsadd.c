/*
!> @file
!!  Routines to extract address values from different objects
!!  needed to be named differently to circumvent g95 limitation
!! @author
!!    Copyright (C) 2007-2011 BigDFT group 
!!    This file is distributed under the terms of the
!!    GNU General Public License, see ~/COPYING file
!!    or http://www.gnu.org/copyleft/gpl.txt .
!!    For the list of contributors, see ~/AUTHORS 
*/

#include <config.h>

#define _GNU_SOURCE

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "utils.h"

void FC_FUNC_(call_external_c, CALL_EXTERNAL_C)(void *callback(),void *address())
{
  //  *address=0;
  address = callback;
  printf("\n test address = %p; \n", (void*)callback);
  callback();
  printf("\n test address2 = %p , %ld; \n", (void*)callback,sizeof(callback));
  address();
  printf("\n test address3 = %p , %lld; \n", (void*)callback,(long long int)callback);
  return;
}

void FC_FUNC_(c_memcopy, C_MEMCOPY)(void *to, long long int * fromadd, long long int *ln)
{
  char *from = (char*)*fromadd;
  size_t nbytes = sizeof(char) * *ln;
  memcpy(to, from, nbytes);
}

void FC_FUNC_(callable_str, CALLABLE_STR)(FFunc_str *func, void **a, int *ln)
{
  if (func && a)
    (*func)(*a, *ln);
}

void FC_FUNC_(callable_arg_arg, CALLABLE_ARG_ARG)(FFunc_arg_arg *func, void **a, void **a2)
{
  if (func && a && a2)
    (*func)(*a, *a2);
}

void FC_FUNC_(callable_arg_str, CALLABLE_ARG_STR)(FFunc_arg_str *func, void **a, void **a2, int *ln)
{
  if (func && a && a2)
    (*func)(*a, *a2, *ln);
}

void FC_FUNC_(callable_str_str, CALLABLE_STR_STR)(FFunc_str_str *func, void **a, void **a2, int *ln, int *ln2)
{
  if (func && a && a2)
    (*func)(*a, *a2, *ln, *ln2);
}

void FC_FUNC_(callable_arg_arg_arg, CALLABLE_ARG_ARG_ARG)(FFunc_arg_arg_arg *func, void **a, void **a2, void **a3)
{
  if (func && a && a2 && a3)
    (*func)(*a, *a2, *a3);
}

void FC_FUNC_(callable_arg_arg_str, CALLABLE_ARG_ARG_STR)(FFunc_arg_arg_str *func, void **a, void **a2, void **a3, int *ln)
{
  if (func && a && a2 && a3)
    (*func)(*a, *a2, *a3, *ln);
}

void FC_FUNC_(callable_arg_str_str, CALLABLE_ARG_STR_STR)(FFunc_arg_str_str *func, void **a, void **a2, void **a3, int *ln, int *ln2)
{
  if (func && a && a2 && a3)
    (*func)(*a, *a2, *a3, *ln, *ln2);
}

void FC_FUNC_(callable_str_str_str, CALLABLE_STR_STR_STR)(FFunc_str_str_str *func, void **a, void **a2, void **a3, int *ln, int *ln2, int *ln3)
{
  if (func && a && a2 && a3)
    (*func)(*a, *a2, *a3, *ln, *ln2, *ln3);
}

void FC_FUNC_(callable_arg_arg_arg_arg, CALLABLE_ARG_ARG_ARG_ARG)(FFunc_arg_arg_arg_arg *func, void **a, void **a2, void **a3, void **a4)
{
  if (func && a && a2 && a3 && a4)
    (*func)(*a, *a2, *a3, *a4);
}

void FC_FUNC_(callable_arg_arg_arg_str, CALLABLE_ARG_ARG_ARG_STR)(FFunc_arg_arg_arg_str *func, void **a, void **a2, void **a3, void **a4, int *ln)
{
  if (func && a && a2 && a3 && a4)
    (*func)(*a, *a2, *a3, *a4, *ln);
}

void FC_FUNC_(callable_arg_arg_str_str, CALLABLE_ARG_ARG_STR_STR)(FFunc_arg_arg_str_str *func, void **a, void **a2, void **a3, void **a4, int *ln, int *ln2)
{
  if (func && a && a2 && a3 && a4)
    (*func)(*a, *a2, *a3, *a4, *ln, *ln2);
}

void FC_FUNC_(callable_arg_str_str_str, CALLABLE_ARG_STR_STR_STR)(FFunc_arg_str_str_str *func, void **a, void **a2, void **a3, void **a4, int *ln, int *ln2, int *ln3)
{
  if (func && a && a2 && a3 && a4)
    (*func)(*a, *a2, *a3, *a4, *ln, *ln2, *ln3);
}

void FC_FUNC_(callable_str_str_str_str, CALLABLE_STR_STR_STR_STR)(FFunc_str_str_str_str *func, void **a, void **a2, void **a3, void **a4, int *ln, int *ln2, int *ln3, int *ln4)
{
  if (func && a && a2 && a3 && a4)
    (*func)(*a, *a2, *a3, *a4, *ln, *ln2, *ln3, *ln4);
}


void FC_FUNC_(callable_arg_arg_arg_arg_arg, CALLABLE_ARG_ARG_ARG_ARG_ARG)(FFunc_arg_arg_arg_arg_arg *func, void **a, void **a2, void **a3, void **a4, void **a5)
{
  if (func && a && a2 && a3 && a4 && a5)
    (*func)(*a, *a2, *a3, *a4, *a5);
}

void FC_FUNC_(callable_arg_arg_arg_arg_str, CALLABLE_ARG_ARG_ARG_ARG_STR)(FFunc_arg_arg_arg_arg_str *func, void **a, void **a2, void **a3, void **a4, void **a5, int *ln)
{
  if (func && a && a2 && a3 && a4 && a5)
    (*func)(*a, *a2, *a3, *a4, *a5, *ln);
}

void FC_FUNC_(callable_arg_arg_arg_str_str, CALLABLE_ARG_ARG_ARG_STR_STR)(FFunc_arg_arg_arg_str_str *func, void **a, void **a2, void **a3, void **a4, void **a5, int *ln, int *ln2)
{
  if (func && a && a2 && a3 && a4 && a5)
    (*func)(*a, *a2, *a3, *a4, *a5, *ln, *ln2);
}

void FC_FUNC_(callable_arg_arg_str_str_str, CALLABLE_ARG_ARG_STR_STR_STR)(FFunc_arg_arg_str_str_str *func, void **a, void **a2, void **a3, void **a4, void **a5, int *ln, int *ln2, int *ln3)
{
  if (func && a && a2 && a3 && a4 && a5)
    (*func)(*a, *a2, *a3, *a4, *a5, *ln, *ln2, *ln3);
}

void FC_FUNC_(callable_arg_str_str_str_str, CALLABLE_ARG_STR_STR_STR_STR)(FFunc_arg_str_str_str_str *func, void **a, void **a2, void **a3, void **a4, void **a5, int *ln, int *ln2, int *ln3, int *ln4)
{
  if (func && a && a2 && a3 && a4 && a5)
    (*func)(*a, *a2, *a3, *a4, *a5, *ln, *ln2, *ln3, *ln4);
}

void FC_FUNC_(callable_str_str_str_str_str, CALLABLE_STR_STR_STR_STR_STR)(FFunc_str_str_str_str_str *func, void **a, void **a2, void **a3, void **a4, void **a5, int *ln, int *ln2, int *ln3, int *ln4, int *ln5)
{
  if (func && a && a2 && a3 && a4 && a5)
    (*func)(*a, *a2, *a3, *a4, *a5, *ln, *ln2, *ln3, *ln4, *ln5);
}

void f_func_call(void *func, unsigned int n_args, void* args[FFUNC_MAX_ARGS],
                 const unsigned int strs[FFUNC_MAX_ARGS])
{
  FFunc_void func_0;
  FFunc_str func_1;
  FFunc_str_str func_2;
  FFunc_str_str_str func_3;
  FFunc_str_str_str_str func_4;
  FFunc_str_str_str_str_str func_5;

  switch (n_args)
    {
    case 0:
      *(&func_0) = func;
      func_0();
      return;
    case 1:
      *(&func_1) = func;
      func_1(args[0], strs[0]);
      return;
    case 2:
      *(&func_2) = func;
      func_2(args[0], args[1], strs[0], strs[1]);
      return;
    case 3:
      *(&func_3) = func;
      func_3(args[0], args[1], args[2], strs[0], strs[1], strs[2]);
      return;
    case 4:
      *(&func_4) = func;
      func_4(args[0], args[1], args[2], args[3], strs[0], strs[1], strs[2], strs[3]);
      return;
    case 5:
      *(&func_5) = func;
      func_5(args[0], args[1], args[2], args[3], args[4],
             strs[0], strs[1], strs[2], strs[3], strs[4]);
      return;
    default:
      return;
    }
}

//Symbol duplications for fortran interfaces

void FC_FUNC(geti1, GETI1)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(geti2, GETI2)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(geti3, GETI3)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(geti4, GETI4)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getc1, GETC1)(int *len,void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getl1, GETL1)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getl2, GETL2)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getl3, GETL3)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getr1, GETR1)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getr2, GETR2)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getr3, GETR3)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getr4, GETR4)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp1, GETDP1)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp2, GETDP2)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp3, GETDP3)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp4, GETDP4)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp5, GETDP5)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp6, GETDP6)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp7, GETDP7)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getz2, GETZ2)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getz3, GETZ3)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp1ptr, GETDP1PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp2ptr, GETDP2PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp3ptr, GETDP3PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp4ptr, GETDP4PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp5ptr, GETDP5PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getdp6ptr, GETDP6PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(geti1ptr, GETI1PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(geti2ptr, GETI2PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(geti3ptr, GETI3PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(geti4ptr, GETI4PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getl1ptr, GETL1PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getl2ptr, GETL2PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getl3ptr, GETL3PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getc1ptr, GETC1ptr)(int *len,void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}

void FC_FUNC(getz1ptr, GETZ1PTR)(void *ptr,long long int *address)
{
  *address=0;
  *address = (long long int)ptr;
  return;
}
