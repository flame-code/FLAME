/**
 * \file
 *    Common callable definitions.
 */

/**
 * \author
 *    Copyright (C) 2016 BigDFT group
 *    This file is distributed under the terms of the
 *    GNU General Public License, see ~/COPYING file
 *    or http://www.gnu.org/copyleft/gpl.txt .
 *    For the list of contributors, see ~/AUTHORS
 */

#define FFUNC_MAX_ARGS 3

#include <addresses.h>

typedef void (*FFunc_str)(void *a, int ln);

typedef void (*FFunc_arg_arg)(void *a, void *a2);
typedef void (*FFunc_arg_str)(void *a, void *a2, int ln);
typedef void (*FFunc_str_str)(void *a, void *a2, int ln, int ln2);

typedef void (*FFunc_arg_arg_arg)(void *a, void *a2, void *a3);
typedef void (*FFunc_arg_arg_str)(void *a, void *a2, void *a3, int ln);
typedef void (*FFunc_arg_str_str)(void *a, void *a2, void *a3, int ln, int ln2);
typedef void (*FFunc_str_str_str)(void *a, void *a2, void *a3,
                                  int ln, int ln2, int ln3);

typedef void (*FFunc_arg_arg_arg_arg)(void *a, void *a2, void *a3, void *a4);
typedef void (*FFunc_arg_arg_arg_str)(void *a, void *a2, void *a3, void *a4, int ln);
typedef void (*FFunc_arg_arg_str_str)(void *a, void *a2, void *a3, void *a4,
                                      int ln, int ln2);
typedef void (*FFunc_arg_str_str_str)(void *a, void *a2, void *a3, void *a4,
                                      int ln, int ln2, int ln3);
typedef void (*FFunc_str_str_str_str)(void *a, void *a2, void *a3, void *a4,
                                      int ln, int ln2, int ln3, int ln4);

typedef void (*FFunc_arg_arg_arg_arg_arg)(void *a, void *a2, void *a3, void *a4, void *a5);
typedef void (*FFunc_arg_arg_arg_arg_str)(void *a, void *a2, void *a3, void *a4, void *a5,
                                          int ln);
typedef void (*FFunc_arg_arg_arg_str_str)(void *a, void *a2, void *a3, void *a4, void *a5,
                                          int ln, int ln2);
typedef void (*FFunc_arg_arg_str_str_str)(void *a, void *a2, void *a3, void *a4, void *a5,
                                          int ln, int ln2, int ln3);
typedef void (*FFunc_arg_str_str_str_str)(void *a, void *a2, void *a3, void *a4, void *a5,
                                          int ln, int ln2, int ln3, int ln4);
typedef void (*FFunc_str_str_str_str_str)(void *a, void *a2, void *a3, void *a4, void *a5,
                                          int ln, int ln2, int ln3, int ln4, int ln5);

void f_func_call(void *func, unsigned int n_args, void* args[FFUNC_MAX_ARGS],
                 const unsigned int strs[FFUNC_MAX_ARGS]);
