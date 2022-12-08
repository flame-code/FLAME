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

/* use MPI2 capabilities. */
#define HAVE_MPI2 1

/* use MPI_INIT_THREAD */
#define HAVE_MPI_INIT_THREAD 1

/* Name of package */
#define PACKAGE "FLAME"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "alirezagh76@gmail.com"

/* Define to the full name of this package. */
#define PACKAGE_NAME "FLAME - Fully Loaded library for Atomistic Modelling Environment"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "FLAME - Fully Loaded library for Atomistic Modelling Environment 1.0"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "FLAME"

/* Define to the home page for this package. */
#define PACKAGE_URL ""

/* Define to the version of this package. */
#define PACKAGE_VERSION "1.0"

/* Version number of package */
#define VERSION "1.0"
