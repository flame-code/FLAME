#ifndef ERR_H
#define ERR_H
#include "futile_cst.h"
#include <limits.h>
#include <stdbool.h>

#define ERR_NOT_DEFINE UINT_MAX
#define ERR_ANY (UINT_MAX - 1)

typedef unsigned int ErrId;

ErrId err_define(const char *name, const char *message, const char *action);

void err_open_try(void);
void err_close_try(void);

void err_throw_by_id(const char *message, ErrId id);
void err_throw_by_name(const char *message, const char *name);

bool err_check(ErrId id);
bool err_check_by_name(const char *name);

#endif
