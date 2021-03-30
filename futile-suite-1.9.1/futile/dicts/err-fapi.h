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

void FC_FUNC_(bind_err_define_with_action, BIND_ERR_DEFINE_WITH_ACTION)(int*, const char*, const int*, const char*, const int*, const char *, const int*, int, int, int);
void FC_FUNC_(bind_err_define, BIND_ERR_DEFINE)(int*, const char*, const int*, const char*, const int*, int, int);

void FC_FUNC_(bind_err_open_try, BIND_ERR_OPEN_TRY)(void);
void FC_FUNC_(bind_err_close_try, BIND_ERR_CLOSE_TRY)(void);

void FC_FUNC_(bind_err_throw_by_id, BIND_ERR_THROW_BY_ID)(const int*);
void FC_FUNC_(bind_err_throw_by_id_with_mess, BIND_ERR_THROW_BY_ID_WITH_MESS)(const char*, const int*, const int*, int);

void FC_FUNC_(bind_err_throw_by_name, BIND_ERR_THROW_BY_NAME)(const char*, const int*, int);
void FC_FUNC_(bind_err_throw_by_name_with_mess, BIND_ERR_THROW_BY_NAME_WITH_MESS)(const char*, const int*, const char*, const int*, int, int);
void FC_FUNC_(bind_err_check, BIND_ERR_CHECK)(int*);
void FC_FUNC_(bind_err_check_by_id, BIND_ERR_CHECK_BY_ID)(int*, const int*);
void FC_FUNC_(bind_err_check_by_name, BIND_ERR_CHECK_BY_NAME)(int*, const char*, const int*, int);

#endif
