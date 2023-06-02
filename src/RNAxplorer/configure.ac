AC_INIT([RNAxplorer],[0.9.0],[ronny@tbi.univie.ac.at],[RNAxplorer],[https://github.com/ViennaRNA/RNAxplorer])

dnl minimally required ViennaRNA version
VRNA_Version=2.4.12
AC_PREREQ([2.59])
AC_CONFIG_AUX_DIR(config)

AC_SUBST(VERSION)

AC_CONFIG_MACRO_DIR([m4])

AM_SILENT_RULES([yes])

AM_INIT_AUTOMAKE(-Wall -Werror foreign)


#Create debug conditional
# prevent setting xFLAGS to default of -g -O2
if test x"$CFLAGS" = x""; then
  AC_SUBST(CFLAGS, [ ])
fi
if test x"$CPPFLAGS" = x""; then
  AC_SUBST(CPPFLAGS, [ ])
fi

AC_ARG_ENABLE([debug],
              AC_HELP_STRING([--enable-debug],
              [compile with debug symbols @<:@default=no@:>@]),
              [want_debug="$enableval"], [want_debug=no])

if test "x$want_debug" = "xyes"; then
  CFLAGS="$CFLAGS -O0 -ggdb"
  CPPFLAGS="$CPPFLAGS -O0 -ggdb"
  AC_DEFINE([DEBUG], 1, [Define for debugging])
else
  CFLAGS="$CFLAGS -O2"
  CPPFLAGS="$CPPFLAGS -O2"
fi

#add lapack, (needed for multidimensional sampling)
AX_LAPACK([
      AC_CHECK_HEADERS([openblas/lapacke.h lapacke/lapacke.h lapacke.h],[rx_lapacke_header=yes;break;],)
      AS_IF([test "x$rx_lapacke_header" = "xyes"],[
        LDFLAGS="$LDFLAGS -llapack"
      ],[
        AC_MSG_ERROR([could not find lapacke.h])
      ])
      ],[AC_MSG_ERROR([could not find required version of Lapack library])])

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
AM_PROG_CC_C_O

#Output this variables to the makefiles

AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

#Check headers
AC_CHECK_HEADERS(malloc.h strings.h unistd.h)
AC_C_INLINE   #check if C compiler uses INLINE

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config.h])

AM_PROG_AR

AC_PROG_LIBTOOL
LT_INIT([pic-only, static])
AC_SUBST(LIBTOOL_DEPS)

dnl Check for OpenMP and prepare flags
AC_OPENMP
AC_SUBST(OPENMP_CFLAGS)
CPPFLAGS="$CPPFLAGS $OPENMP_CFLAGS"
LDFLAGS="$LDFLAGS $OPENMP_CFLAGS"

dnl Vienna RNA package library path support, if not installed in standard directory/path
VRNA_OK=true

AC_ARG_WITH([vrna],
    AC_HELP_STRING(
        [--with-vrna=PREFIX],
        [Installation directory of the Vienna RNA library]
      )
)

dnl require rna lib linking
ENABLE_LIBRNA=yes


dnl ------------------------------
dnl Vienna RNA lib configuration
dnl
dnl check whether pkg-config is available
AC_CHECK_PROG([PKG_CONFIG_AVAILABLE],pkg-config,[true],[false])

VRNALIB_OK=false

dnl if pkg-config is available, use it to configure the Vienna RNA library

AS_IF([test "$PKG_CONFIG_AVAILABLE" = "true"],[
    dnl if with_vrna is set, then use it to set the path for pkg_config
    if test "$with_vrna" != "" ; then
        export PKG_CONFIG_PATH="$with_vrna/lib/pkgconfig:$PKG_CONFIG_PATH"
    fi
    #Search for the RNAlib2 and set variables for Makefiles if found.  _CFLAGS and _LIBS can be overriden with the modue Name from PKG_CHECK_MODULES
    PKG_CHECK_MODULES([VRNA], [RNAlib2 >= $VRNA_Version],
                  [ AC_SUBST([VRNA_CFLAGS])
                    AC_SUBST([VRNA_LIBS])
                    VRNALIB_OK=true
                  ],
                  [VRNALIB_OK=false])
])

dnl if pkg-config is not available or something else went wrong with pkg-config,
dnl then guess the configuration
dnl
AS_IF([test "$VRNALIB_OK" != "true" -a "$PKG_CONFIG_AVAILABLE" != "true" ],
[
  AC_MSG_NOTICE([Guess the Vienna RNA library configuration:])
  if test "$with_vrna" = "" ; then
    if test -e "/usr/local/include/ViennaRNA" ; then
      vrna_prefix="/usr/local"
    else
      vrna_prefix="/usr"
    fi
  else
    vrna_prefix="$with_vrna"
  fi
  VRNA_LIBS="-L$vrna_prefix/lib -lRNA"
  dnl if libgomp is available, then add it to VRNA_LIBS; otherwise,
  dnl assume (optimistically) that libRNA works without libgomp
  dnl AC_CHECK_LIB([gomp],[main],[VRNA_LIBS="$VRNA_LIBS -lgomp"])
  VRNA_CFLAGS="-I$vrna_prefix/include"
  AC_MSG_NOTICE([Guessing of ViennaRNA flags may be avoided by installing pkg-config.])
  AC_CHECK_LIB([RNA],[vrna_fold],[VRNALIB_OK=true],[VRNALIB_OK=false])
])

AS_IF([test "$VRNALIB_OK" != "true"],
[
  AC_MSG_ERROR([
                ============================================================
                   The Vienna RNA library (>= $VRNA_Version) could not be found! 
                   Please install it to a standard location or 
                   specify its path with --with-vrna=<path-to-viennaRNA>!
                ============================================================
                ])
])

# add vrna to libs
##LIBS="$VRNA_LIBS $LIBS"




dnl Checks for libraries.
AC_CHECK_LIB([m], [exp])

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([malloc.h float.h limits.h stdlib.h string.h strings.h unistd.h unistd.h math.h])
dnl ----------------------------------------
dnl for Vienna RNA
dnl     defines HAVE_STRDUP if strdup is available
AC_CHECK_FUNCS([strdup])

AC_RXP_INIT

AC_OUTPUT([packaging/RNAxplorer.spec packaging/PKGBUILD])

AC_RXP_NOTICE

