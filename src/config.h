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

/* compile the code with debugging options */
/* #undef HAVE_DEBUG */

/* If set, we can call glib.h */
/* #undef HAVE_GLIB */

/* Define to 1 if you have the <inttypes.h> header file. */
/* #undef HAVE_INTTYPES_H */

/* libarchive is linkable. */
/* #undef HAVE_LIB_ARCHIVE */

/* Define to 1 if you have the <memory.h> header file. */
/* #undef HAVE_MEMORY_H */

/* use MPI2 capabilities. */
#define HAVE_MPI2 1

/* use MPI3 capabilities (like MPI_IALLREDUCE and MPI_IALLTOALLV). */
/* #undef HAVE_MPI3 */

/* use MPI_INIT_THREAD */
#define HAVE_MPI_INIT_THREAD 1

/* If set, we can call Python.h */
/* #undef HAVE_PYTHON */

/* Define to 1 if you have the <stdint.h> header file. */
/* #undef HAVE_STDINT_H */

/* Define to 1 if you have the <stdlib.h> header file. */
/* #undef HAVE_STDLIB_H */

/* Define to 1 if you have the <strings.h> header file. */
/* #undef HAVE_STRINGS_H */

/* Define to 1 if you have the <string.h> header file. */
/* #undef HAVE_STRING_H */

/* Define to 1 if you have the `strndup' function. */
#define HAVE_STRNDUP 1

/* Define to 1 if you have the <sys/stat.h> header file. */
/* #undef HAVE_SYS_STAT_H */

/* Define to 1 if you have the <sys/types.h> header file. */
/* #undef HAVE_SYS_TYPES_H */

/* Define to 1 if you have the <unistd.h> header file. */
/* #undef HAVE_UNISTD_H */

/* Name of package */
#define PACKAGE "FLAME"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "alirezagh76@gmail.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "FLAME - TO BE FILLED LATER"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "FLAME - VERSION TO BE FILLED LATER"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "FLAME"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "VERSION TO BE FILLED LATER"

/* Define to 1 if you have the ANSI C header files. */
/* #undef STDC_HEADERS */

/* Version number of package */
#define VERSION "VERSION TO BE FILLED LATER"
