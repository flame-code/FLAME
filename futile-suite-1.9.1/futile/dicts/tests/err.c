#include "dict.h"
#include "err.h"

#include <string.h>
#include <stdio.h>
#include <stdbool.h>

#define ASSERT(T) {if (!(T)) {fprintf(stdout, "%s: failed\n", #T); return false;}}
#define ALLEQ(T,U,N) {int it; for(it=0;it < N;it=it+1) { ASSERT(T[it]==U[it]);}}

bool test_f90_err_define()
{
  ErrId id;
  
  id = err_define("ERROR_TEST", "an error for test.", NULL);
  ASSERT(id != ERR_NOT_DEFINE);

  return true;
}

bool test_f90_err_by_id()
{
  ErrId id;

  id = err_define("ERROR_TEST1", "an error for test.", NULL);
  ASSERT(id != ERR_NOT_DEFINE);

  err_open_try();
  err_throw_by_id("something is wrong!", id);

  ASSERT(err_check(ERR_ANY));
  ASSERT(err_check(id));
  ASSERT(err_check_by_name("ERROR_TEST1"));

  err_close_try();
  ASSERT(!err_check(ERR_ANY));

  return true;
}

bool test_f90_err_by_name()
{
  ErrId id;

  id = err_define("ERROR_TEST2", "an error for test.", NULL);
  ASSERT(id != ERR_NOT_DEFINE);

  err_open_try();
  err_throw_by_name("something is wrong!", "ERROR_TEST2");

  ASSERT(err_check(ERR_ANY));
  ASSERT(err_check(id));
  ASSERT(err_check_by_name("ERROR_TEST2"));

  err_close_try();
  ASSERT(!err_check(ERR_ANY));

  return true;
}

#define RUN(T) {if (!T) return 1; else fprintf(stdout, "%s: OK\n", #T);}

int main(int argc, char **argv)
{
#ifdef GLIB_MAJOR_VERSION
  g_type_init();
#endif
  futile_dicts_initialize();
  RUN(test_f90_err_define());
  RUN(test_f90_err_by_id());
  RUN(test_f90_err_by_name());
  futile_dicts_finalize();
  return 0;
}
