#include <config.h>
#include "err.h"
#include "err-fapi.h"

#include <string.h>

ErrId err_define(const char *name, const char *message, const char *action)
{
  int nlen, mlen;
  int id;

  if (!name || !message)
    return ERR_NOT_DEFINE;

  nlen = strlen(name);
  mlen = strlen(message);
  if (action)
    {
      int len = strlen(action);
      FC_FUNC_(bind_err_define_with_action, BIND_ERR_DEFINE_WITH_ACTION)(&id, name, &nlen, message, &mlen, action, &len, nlen, mlen, len);
    }
  else
    FC_FUNC_(bind_err_define, BIND_ERR_DEFINE)(&id, name, &nlen, message, &mlen, nlen, mlen);

  return (ErrId)id;
}

void err_open_try(void)
{
  FC_FUNC_(bind_err_open_try, BIND_ERR_OPEN_TRY)();
}

void err_close_try(void)
{
  FC_FUNC_(bind_err_close_try, BIND_ERR_CLOSE_TRY)();
}

void err_throw_by_id(const char *message, ErrId id)
{
  int iid = (int)id;
  if (message)
    {
      int len = strlen(message);
      FC_FUNC_(bind_err_throw_by_id_with_mess, BIND_ERR_THROW_BY_ID_WITH_MESS)(message, &len, &iid, len);
    }
  else
    FC_FUNC_(bind_err_throw_by_id, BIND_ERR_THROW_BY_ID)(&iid);
}

void err_throw_by_name(const char *message, const char *name)
{
  int nlen;

  if (!name)
    return;
  
  nlen = strlen(name);
  if (message)
    {
      int len = strlen(message);
      FC_FUNC_(bind_err_throw_by_name_with_mess, BIND_ERR_THROW_BY_NAME_WITH_MESS)(message, &len, name, &nlen, len, nlen);
    }
  else
    FC_FUNC_(bind_err_throw_by_name, BIND_ERR_THROW_BY_NAME)(name, &nlen, nlen);
}

bool err_check(ErrId id)
{
  int ret;
  
  if (id == ERR_ANY)
    FC_FUNC_(bind_err_check, BIND_ERR_CHECK)(&ret);
  else
    {
      int iid = (int)id;
      FC_FUNC_(bind_err_check_by_id, BIND_ERR_CHECK_BY_ID)(&ret, &iid);
    }
  return ret == 1 ? true : false;
}

bool err_check_by_name(const char *name)
{
  int ret;
  int nlen;

  if (!name)
    return false;

  nlen = strlen(name);
  FC_FUNC_(bind_err_check_by_name, BIND_ERR_CHECK_BY_NAME)(&ret, name, &nlen, nlen);
  return ret == 1 ? true : false;
}
