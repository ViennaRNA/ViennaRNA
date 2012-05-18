dnl Process this file with autoconf to produce a configure script.

dnl Every other copy of the package version number gets its value from here
AC_INIT([ViennaRNA],[2.0.7],[rna@tbi.univie.ac.at],[ViennaRNA], [http://www.tbi.univie.ac.a/~ivo/RNA])
AC_PREREQ([2.59])

AC_SUBST(VERSION)

m4_include([m4/package.m4])

dnl change the default installation prefix
dnl AC_PREFIX_DEFAULT(/usr/local/RNA)

AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([1.11])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AM_CONFIG_HEADER([config.h])

dnl configure options
AM_WITH_DMALLOC

AC_ARG_VAR(PERLPREFIX, [Prefix for installation of the Perl module. 
By default the perl module is installed wherever Perl thinks it belongs, ignoring the $prefix used by configure])


ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_CANONICAL_HOST

dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CC_FPIC

AC_ARG_VAR(FPIC,[compiler flag for generating PIC code that can be used in a shared library])

if test "$GCC" = yes; then
   AC_DEFINE(UNUSED, __attribute__ ((unused)), avoid warning about unused variables)
else
   AC_DEFINE(UNUSED,, only for gcc)
fi

AC_PROG_INSTALL
AC_PROG_RANLIB
AC_PROG_MAKE_SET
AC_PROG_GREP
AM_PROG_CC_C_O


dnl Check for OpenMP and prepare flags
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)

dnl Checks for libraries.
AC_CHECK_LIB([m], [exp])

dnl Checks for header files.
AC_HEADER_STDC
AC_HEADER_STDBOOL
AC_CHECK_HEADERS([malloc.h float.h limits.h stdlib.h string.h strings.h unistd.h unistd.h math.h])

dnl Checks for funtions
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor strdup strstr strchr strrchr strstr strtol strtoul pow rint sqrt erand48 memset memmove erand48])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE([])

AC_RNA_INIT

AC_OUTPUT

AC_RNA_NOTICE
