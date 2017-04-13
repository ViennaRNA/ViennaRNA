dnl Process this file with autoconf to produce a configure script.

dnl Every other copy of the package version number gets its value from here
AC_INIT([ViennaRNA],[2.3.5],[rna@tbi.univie.ac.at],[ViennaRNA], [http://www.tbi.univie.ac.at/RNA])
AC_PREREQ([2.60])
AC_CONFIG_AUX_DIR(config)

AC_SUBST(VERSION)

AC_CONFIG_MACRO_DIR([m4])

AC_USE_SYSTEM_EXTENSIONS

dnl change the default installation prefix
dnl AC_PREFIX_DEFAULT(/usr/local/RNA)

AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([1.11])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config.h])

dnl configure options
AM_WITH_DMALLOC

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP

if test "$GCC" = yes; then
   AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
   AC_DEFINE(UNUSED,, only for gcc)
fi

AC_PROG_INSTALL
AM_PROG_AR
AC_PROG_MAKE_SET
AC_PROG_GREP
AC_PROG_MKDIR_P
AM_PROG_CC_C_O
AC_PROG_SED
AC_PROG_LN_S
AC_PROG_CC_C99

AC_PROG_LIBTOOL
LT_INIT([pic-only, static])
AC_SUBST(LIBTOOL_DEPS)


dnl Checks for libraries.
AC_CHECK_LIB([m], [exp])

dnl Checks for header files.
AC_HEADER_STDC
AC_HEADER_STDBOOL
AC_CHECK_HEADERS([malloc.h float.h limits.h stdlib.h string.h strings.h unistd.h math.h stdarg.h])

dnl Checks for funtions
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor strdup strstr strchr strrchr strstr strtol strtoul pow rint sqrt erand48 memset memmove erand48 asprintf vasprintf])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE([])

AC_RNA_INIT

AC_OUTPUT

AC_RNA_NOTICE
