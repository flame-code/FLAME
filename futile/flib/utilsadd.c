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

void FC_FUNC_(c_memcopy, C_MEMCOPY)(void *to, long long int * fromadd, int *ln)
{
  char *from = (char*)*fromadd;
  memcpy(to, from, sizeof(char) * *ln);
}

void FC_FUNC_(call_external_c_fromadd, CALL_EXTERNAL_C_FROMADD)(long long int * add)
{
  void * ext;
  //  long long int *ext_data;

  ext=(void*) *add;
  //*ext_data=0;

  //  *address=0;
  //callback1= (void*) *add;
  //printf("\n test NEW address = %p; \n", (void*) *add);
  //callback1();
  //*addredss();
  FC_FUNC_(call_external_f,CALL_EXTERNAL_F)(ext);//,ext_data);

  // printf("\n test NEW address = %lld; \n", *add);
  //  printf("\n test NEW address3 = %p , %lld; \n", (void*)callback,*address);
  return;
}

typedef void (*f_subroutine_data)(void *dataadd);
void FC_FUNC_(call_external_c_fromadd_data, CALL_EXTERNAL_C_FROMADD_DATA)(f_subroutine_data *add, void **dataadd)
{
  if (add && dataadd)
    (*add)(*dataadd);
}

typedef void (*f_subroutine_data_str)(void *dataadd, int ln);
void FC_FUNC_(call_external_c_fromadd_data_str, CALL_EXTERNAL_C_FROMADD_DATA_STR)(f_subroutine_data_str *add, void **dataadd, int *ln)
{
  if (add && dataadd)
    (*add)(*dataadd, *ln);
}

typedef void (*f_subroutine_data_data)(void *dataadd, void *dataadd2);
void FC_FUNC_(call_external_c_fromadd_data_data, CALL_EXTERNAL_C_FROMADD_DATA_DATA)(f_subroutine_data_data *add, void **dataadd, void **dataadd2)
{
  if (add && dataadd && dataadd2)
    (*add)(*dataadd, *dataadd2);
}

typedef void (*f_subroutine_data_data_str)(void *dataadd, void *dataadd2, int ln);
void FC_FUNC_(call_external_c_fromadd_data_data_str, CALL_EXTERNAL_C_FROMADD_DATA_DATA_STR)(f_subroutine_data_data_str *add, void **dataadd, void **dataadd2, int *ln)
{
  if (add && dataadd && dataadd2)
    (*add)(*dataadd, *dataadd2, *ln);
}

typedef void (*f_subroutine_data_data_str_str)(void *dataadd, void *dataadd2, int ln, int ln2);
void FC_FUNC_(call_external_c_fromadd_data_data_str_str, CALL_EXTERNAL_C_FROMADD_DATA_DATA_STR_STR)(f_subroutine_data_data_str_str *add, void **dataadd, void **dataadd2, int *ln, int *ln2)
{
  if (add && dataadd && dataadd2)
    (*add)(*dataadd, *dataadd2, *ln, *ln2);
}

typedef void (*f_subroutine_data_data_data)(void *dataadd, void *dataadd2, void *dataadd3);
void FC_FUNC_(call_external_c_fromadd_data_data_data, CALL_EXTERNAL_C_FROMADD_DATA_DATA_DATA)(f_subroutine_data_data_data *add, void **dataadd, void **dataadd2, void **dataadd3)
{
  if (add && dataadd && dataadd2 && dataadd3)
    (*add)(*dataadd, *dataadd2, *dataadd3);
}

typedef void (*f_subroutine_data_data_data_str)(void *dataadd, void *dataadd2, void **dataadd3, int ln);
void FC_FUNC_(call_external_c_fromadd_data_data_data_str, CALL_EXTERNAL_C_FROMADD_DATA_DATA_DATA_STR)(f_subroutine_data_data_data_str *add, void **dataadd, void **dataadd2, void **dataadd3, int *ln)
{
  if (add && dataadd && dataadd2 && dataadd3)
    (*add)(*dataadd, *dataadd2, *dataadd3, *ln);
}

typedef void (*f_subroutine_data_data_data_str_str)(void *dataadd, void *dataadd2, void **dataadd3, int ln, int ln2);
void FC_FUNC_(call_external_c_fromadd_data_data_data_str_str, CALL_EXTERNAL_C_FROMADD_DATA_DATA_DATA_STR_STR)(f_subroutine_data_data_data_str_str *add, void **dataadd, void **dataadd2, void **dataadd3, int *ln, int *ln2)
{
  if (add && dataadd && dataadd2 && dataadd3)
    (*add)(*dataadd, *dataadd2, *dataadd3, *ln, *ln2);
}

typedef void (*f_subroutine_data_data_data_str_str_str)(void *dataadd, void *dataadd2, void **dataadd3, int ln, int ln2, int ln3);
void FC_FUNC_(call_external_c_fromadd_data_data_data_str_str_str, CALL_EXTERNAL_C_FROMADD_DATA_DATA_DATA_STR_STR_STR)(f_subroutine_data_data_data_str_str_str *add, void **dataadd, void **dataadd2, void **dataadd3, int *ln, int *ln2, int *ln3)
{
  if (add && dataadd && dataadd2 && dataadd3)
    (*add)(*dataadd, *dataadd2, *dataadd3, *ln, *ln2, *ln3);
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
