/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef FC_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* Define to a macro mangling the given C identifier (in lower and upper
   case), which must not contain underscores, for linking with Fortran. */
#define FC_FUNC(name,NAME) name ## _

/* As FC_FUNC, but for C identifiers containing underscores. */
#define FC_FUNC_(name,NAME) name ## _

/* Define to 1 if you have the `aligned_alloc' function. */
#define HAVE_ALIGNED_ALLOC 1

/* Define to 1 if you have the `clock_gettime' function. */
#define HAVE_CLOCK_GETTIME 1

/* Flush(6) can be used safely in fortran */
#define HAVE_FC_FLUSH 1

/* call flush(6) can be used safely in fortran */
/* #undef HAVE_FC_FLUSH_SUB */

/* get_command_argument() can be used safely in Fortran */
#define HAVE_FC_GET_COMMAND_ARGUMENT 1

/* If set, we can call glib.h */
/* #undef HAVE_GLIB */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define to 1 if you have the `dl' library (-ldl). */
#define HAVE_LIBDL 1

/* Define to 1 if you have the `rt' library (-lrt). */
#define HAVE_LIBRT 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* use MPI2 capabilities. */
#define HAVE_MPI2 1

/* use MPI_INIT_THREAD */
#define HAVE_MPI_INIT_THREAD 1

/* Define to 1 if you have the <numpy/ndarrayobject.h> header file. */
/* #undef HAVE_NUMPY_NDARRAYOBJECT_H */

/* if set, we can call Python.h */
/* #undef HAVE_PYTHON */

/* if set, we can call numpy/ndarrayobject.h */
/* #undef HAVE_PYTHON_NUMPY */

/* use SIMGRID allocators. */
/* #undef HAVE_SIMGRID_SHARED_ALLOCS */

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if you have the `strndup' function. */
#define HAVE_STRNDUP 1

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <time.h> header file. */
#define HAVE_TIME_H 1

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Name of package */
#define PACKAGE "futile"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "Damien.Caliste@cea.fr"

/* Define to the full name of this package. */
#define PACKAGE_NAME "Futile"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "Futile 1.8"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "futile"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.8"

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* Version number of package */
#define VERSION "1.8"
