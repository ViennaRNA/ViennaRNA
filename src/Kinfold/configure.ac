#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)

AC_INIT([kinfold], [1.3], [rna@tbi.univie.ac.at], [Kinfold])
AC_CONFIG_SRCDIR([cache_util.h])

AM_SILENT_RULES([yes])
dnl Every other copy of the package version number gets its value from here
AM_INIT_AUTOMAKE([1.11])

dnl configure options
AM_WITH_DMALLOC

dnl Checks for programs.
AC_PROG_CC
dnl AC_PROG_MAKE_SET

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS(config.h)

AC_SUBST(VERSION)

ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_CANONICAL_HOST

dnl Checks for library functions.
AC_CHECK_FUNCS([strdup memset strchr])

PKG_PROG_PKG_CONFIG

## Search for the RNAlib2 and set variables for Makefiles.
##
## NOTE, we only check for RNAlib2 if the corresponding
## variables are not set in the environment. Also note,
## linking against API v3.0 of the ViennaRNA Package is
## default, so we check old API variables first. If they
## are not set, we use pkg-config to determine whether
## we find API3, or API2 within the system. If all cases
## fail, the configure script fails as well
AS_IF([test "x$VRNA2_CFLAGS" != "x" && test "x$VRNA2_LIBS" != "x"],
    [
        HAVE_LIBRNA_API2=1
        AC_MSG_NOTICE([Using VRNA API 2 as provided in environment variables])
    ],[
        PKG_CHECK_MODULES([VRNA],
                          [RNAlib2 >= 2.2.0],
                          [],
                          [PKG_CHECK_MODULES([VRNA2], [RNAlib2 >= 2.0])])
    ])

AS_IF([test "x$VRNA_CFLAGS" != "x" && test "x$VRNA_LIBS" != "x"],
    [
        AC_DEFINE([HAVE_LIBRNA_API3], [1], [Link against RNAlib v3.0 API])
        HAVE_LIBRNA_API3=1
        AC_SUBST([VRNA_CFLAGS])
        AC_SUBST([VRNA_LIBS])
    ],[
    AS_IF([test "x$VRNA2_CFLAGS" != "x" && test "x$VRNA2_LIBS" != "x"],
        [ AC_DEFINE([HAVE_LIBRNA_API2], [1], [Link against RNAlib v2.0 API])
          HAVE_LIBRNA_API2=1
          AC_SUBST([VRNA2_CFLAGS])
          AC_SUBST([VRNA2_LIBS])
        ])
])

AM_CONDITIONAL([WITH_LIBRNA_API3], [test "x${HAVE_LIBRNA_API3}" = "x1"])
AM_CONDITIONAL([WITH_LIBRNA_API2], [test "x${HAVE_LIBRNA_API2}" = "x1"])

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(stdlib.h string.h strings.h unistd.h)

dnl Checks for libraries.
dnl Replace `main' with a function in -lm:
AC_CHECK_LIB(m, exp)

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

AC_CONFIG_FILES([Example/Makefile Makefile])

AC_OUTPUT
