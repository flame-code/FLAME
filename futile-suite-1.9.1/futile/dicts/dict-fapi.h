/** @file
 * Bindings for the BigDFT package
 * @author
 * Copyright (C) 2013-2015 BigDFT group (DC)
 * This file is distributed under the terms of the
 * GNU General Public License, see ~/COPYING file
 * or http://www.gnu.org/copyleft/gpl.txt .
 * For the list of contributors, see ~/AUTHORS
**/
#ifndef BINDINGS_API_H
#define BINDINGS_API_H

#include "dict.h"

#undef hz

/* bind_dict_append futile/c-bindings/dictf.f90:36 */
/* Fortran header:
subroutine bind_dict_append(dict)
use dictionaries, only: dictionary, operator(//), dict_len
implicit none
type(dictionary), pointer :: dict
*/
void FC_FUNC_(bind_dict_append, BIND_DICT_APPEND)(f90_dictionary_pointer *dict);
/* bind_dict_dump futile/c-bindings/dictf.f90:56 */
/* Fortran header:
subroutine bind_dict_dump(dict, unit)
use dictionaries, only: dictionary
use yaml_output, only: yaml_dict_dump
implicit none
type(dictionary), pointer :: dict
integer, intent(in) :: unit
*/
void FC_FUNC_(bind_dict_dump, BIND_DICT_DUMP)(f90_dictionary_pointer *dict, 
                                              const int *unit);
/* bind_dict_dump_to_file futile/c-bindings/dictf.f90:71 */
/* Fortran header:
subroutine bind_dict_dump_to_file(dict, path)
use dictionaries, only: dictionary
use yaml_output, only: yaml_dict_dump, yaml_set_stream, yaml_close_stream, yaml_get_default_stream
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: path

integer :: unit
*/
void FC_FUNC_(bind_dict_dump_to_file, BIND_DICT_DUMP_TO_FILE)(f90_dictionary_pointer *dict, 
                                                              const char *path, 
                                                              int str_ln_1);
/* bind_dict_free futile/c-bindings/dictf.f90:214 */
/* Fortran header:
subroutine bind_dict_free(dict)
use dictionaries_base, only: dict_free
implicit none
type(dictionary), pointer :: dict
*/
void FC_FUNC_(bind_dict_free, BIND_DICT_FREE)(f90_dictionary_pointer *dict);
/* bind_dict_init futile/c-bindings/dictf.f90:206 */
/* Fortran header:
subroutine bind_dict_init(dict)
use dictionaries, only: dictionary, dict_init
implicit none
type(dictionary), pointer :: dict
*/
void FC_FUNC_(bind_dict_init, BIND_DICT_INIT)(f90_dictionary_pointer *dict);
/* bind_dict_insert futile/c-bindings/dictf.f90:25 */
/* Fortran header:
subroutine bind_dict_insert(dict, key)
use dictionaries, only: dictionary, operator(//)
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: key
*/
void FC_FUNC_(bind_dict_insert, BIND_DICT_INSERT)(f90_dictionary_pointer *dict, 
                                                  const char *key, 
                                                  int str_ln_1);
/* bind_dict_iter futile/c-bindings/dictf.f90:140 */
/* Fortran header:
subroutine bind_dict_iter(dict, exists)
use dictionaries, only: dictionary, dict_iter
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists

type(dictionary), pointer :: start
*/
void FC_FUNC_(bind_dict_iter, BIND_DICT_ITER)(f90_dictionary_pointer *dict, 
                                              int *exists);
/* bind_dict_key futile/c-bindings/dictf.f90:130 */
/* Fortran header:
subroutine bind_dict_key(dict, buf)
use dictionaries, only: dictionary, max_field_length, dict_key
implicit none
type(dictionary), pointer :: dict
character(len = max_field_length), intent(out) :: buf
*/
void FC_FUNC_(bind_dict_key, BIND_DICT_KEY)(const f90_dictionary_pointer *dict, 
                                            char *buf, 
                                            int str_ln_1);
/* bind_dict_len futile/c-bindings/dictf.f90:168 */
/* Fortran header:
subroutine bind_dict_len(dict, ln)
use dictionaries, only: dictionary, dict_len
implicit none
type(dictionary), pointer :: dict
integer, intent(out) :: ln
*/
void FC_FUNC_(bind_dict_len, BIND_DICT_LEN)(const f90_dictionary_pointer *dict, 
                                            int *ln);
/* bind_dict_move_to_item futile/c-bindings/dictf.f90:13 */
/* Fortran header:
subroutine bind_dict_move_to_item(dict, exists, id)
use dictionaries, only: dictionary, operator(//), dict_len
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists
integer, intent(in) :: id
*/
void FC_FUNC_(bind_dict_move_to_item, BIND_DICT_MOVE_TO_ITEM)(f90_dictionary_pointer *dict, 
                                                              int *exists, 
                                                              const int *id);
/* bind_dict_move_to_key futile/c-bindings/dictf.f90:1 */
/* Fortran header:
subroutine bind_dict_move_to_key(dict, exists, key)
use dictionaries, only: dictionary, operator(//), has_key
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists
character(len = *), intent(in) :: key
*/
void FC_FUNC_(bind_dict_move_to_key, BIND_DICT_MOVE_TO_KEY)(f90_dictionary_pointer *dict, 
                                                            int *exists, 
                                                            const char *key, 
                                                            int str_ln_1);
/* bind_dict_next futile/c-bindings/dictf.f90:154 */
/* Fortran header:
subroutine bind_dict_next(dict, exists)
use dictionaries, only: dictionary, dict_next
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists

type(dictionary), pointer :: next
*/
void FC_FUNC_(bind_dict_next, BIND_DICT_NEXT)(f90_dictionary_pointer *dict, 
                                              int *exists);
/* bind_dict_parse futile/c-bindings/dictf.f90:88 */
/* Fortran header:
subroutine bind_dict_parse(dict, buf)
use dictionaries, only: dictionary, operator(//), dict_len,operator(.pop.),dict_free
use yaml_parse, only: yaml_parse_from_string
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: buf
type(dictionary), pointer :: dict_load
*/
void FC_FUNC_(bind_dict_parse, BIND_DICT_PARSE)(f90_dictionary_pointer *dict, 
                                                const char *buf, 
                                                int str_ln_1);
/* bind_dict_pop futile/c-bindings/dictf.f90:105 */
/* Fortran header:
subroutine bind_dict_pop(dict, exists, key)
use dictionaries, only: dictionary, has_key, dict_remove, dict_init
implicit none
type(dictionary), pointer :: dict
logical, intent(out) :: exists
character(len = *), intent(in) :: key
*/
void FC_FUNC_(bind_dict_pop, BIND_DICT_POP)(f90_dictionary_pointer *dict, 
                                            int *exists, 
                                            const char *key, 
                                            int str_ln_1);
/* bind_dict_set futile/c-bindings/dictf.f90:45 */
/* Fortran header:
subroutine bind_dict_set(dict, val)
use dictionaries, only: dictionary, set
implicit none
type(dictionary), pointer :: dict
character(len = *), intent(in) :: val
*/
void FC_FUNC_(bind_dict_set, BIND_DICT_SET)(f90_dictionary_pointer *dict, 
                                            const char *val, 
                                            int str_ln_1);
/* bind_dict_size futile/c-bindings/dictf.f90:178 */
/* Fortran header:
subroutine bind_dict_size(dict, ln)
use dictionaries, only: dictionary, dict_size
implicit none
type(dictionary), pointer :: dict
integer, intent(out) :: ln
*/
void FC_FUNC_(bind_dict_size, BIND_DICT_SIZE)(const f90_dictionary_pointer *dict, 
                                              int *ln);
/* bind_dict_update futile/c-bindings/dictf.f90:197 */
/* Fortran header:
subroutine bind_dict_update(dict, ref)
use dictionaries, only: dictionary, dict_update
implicit none
type(dictionary), pointer :: dict, ref
*/
void FC_FUNC_(bind_dict_update, BIND_DICT_UPDATE)(f90_dictionary_pointer *dict, 
                                                  const f90_dictionary_pointer *ref);
/* bind_dict_value futile/c-bindings/dictf.f90:120 */
/* Fortran header:
subroutine bind_dict_value(dict, buf)
use dictionaries, only: dictionary, max_field_length, dict_value
implicit none
type(dictionary), pointer :: dict
character(len = max_field_length), intent(out) :: buf
*/
void FC_FUNC_(bind_dict_value, BIND_DICT_VALUE)(const f90_dictionary_pointer *dict, 
                                                char *buf, 
                                                int str_ln_1);
/* bind_dict_set_dict ./c-bindings/dictf.f90:233 */
/* Fortran header:
subroutine bind_dict_set_dict(dict,key,keylen,val)
use dictionaries
type(dictionary), pointer :: dict,val
integer, intent(in) :: keylen
character(len=keylen), intent(in) :: key
*/
void FC_FUNC_(bind_dict_set_dict, BIND_DICT_SET_DICT)(f90_dictionary_pointer *dict, 
                                                      const char *key, 
                                                      const int *keylen, 
                                                      f90_dictionary_pointer *val,
						      int str_ln_1);
/* bind_dict_set_double_array ./c-bindings/dictf.f90:222 */
/* Fortran header:
subroutine bind_dict_set_double_array(dict,key,keylen,array,arraylen)
use f_precisions
use dictionaries
type(dictionary), pointer :: dict
integer, intent(in) :: keylen,arraylen
character(len=keylen), intent(in) :: key
real(f_double), dimension(arraylen), intent(in) :: array
*/
void FC_FUNC_(bind_dict_set_double_array, BIND_DICT_SET_DOUBLE_ARRAY)(f90_dictionary_pointer *dict,
								      const char *key, 
                                                                      const int *keylen, 
                                                                      const double *array, 
                                                                      const int *arraylen,
								      int str_ln_1);

void FC_FUNC_(bind_dict_set_double_matrix, BIND_DICT_SET_DOUBLE_MATRIX)(f90_dictionary_pointer *dict,
								      const char *key, 
                                                                      const int *keylen, 
                                                                      const double *array, 
                                                                      const int *arraylenx,
								      const int *arrayleny,
								      int str_ln_1);

void FC_FUNC_(bind_dict_set_double, BIND_DICT_SET_DOUBLE)(f90_dictionary_pointer *dict,
							  const char *key, 
							  const int *keylen, 
							  const double *array,
							  int str_ln_1);

void FC_FUNC_(bind_dict_set_string, BIND_DICT_SET_STRING)(f90_dictionary_pointer *dict,
							  const char* key, 
							  const char* value);

void FC_FUNC_(bind_dict_get_double_array, BIND_DICT_GET_DOUBLE_ARRAY)(f90_dictionary_pointer *dict,
								      const char *key, 
                                                                      const int *keylen, 
                                                                      double *array, 
                                                                      const int *arraylen,
								      int *istat,
								      int str_ln_1);

void FC_FUNC_(bind_dict_get_dict, BIND_DICT_GET_DICT)(f90_dictionary_pointer *dict, 
                                                      const char *key, 
                                                      const int *keylen, 
                                                      f90_dictionary_pointer *val,
						      int *istat,
						      int str_ln_1);

void FC_FUNC_(bind_dict_add_dict, BIND_DICT_ADD_DICT)(f90_dictionary_pointer *dict, 
                                                      f90_dictionary_pointer *val);


void FC_FUNC_(bind_dict_add_double, BIND_DICT_ADD_DOUBLE)(f90_dictionary_pointer *dict, 
							  const double *val);

void FC_FUNC_(bind_dict_add_int, BIND_DICT_ADD_INT)(f90_dictionary_pointer *dict, 
                                                    const int *val);

void FC_FUNC_(bind_dict_add_char, BIND_DICT_ADD_CHAR)(f90_dictionary_pointer *dict, 
                                                      const char *val);

 /*subroutine bind_dict_add_dict(dict,val)
  use dictionaries
  type(dictionary), pointer :: dict,val*/


#endif
