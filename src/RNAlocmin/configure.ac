AC_INIT(RNAlocmin,2.1.0)

AC_CONFIG_MACRO_DIR([.])

#Output this variables to the makefiles

AC_SUBST(VERSION)
ISODATE=`date +%Y-%m-%d`
AC_SUBST(ISODATE)

AC_USE_SYSTEM_EXTENSIONS

AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([-Wall -Werror foreign])

#Check headers
AC_CONFIG_HEADERS(config.h) 

AC_CANONICAL_HOST

#check for compilers c++ and c
AC_PROG_CXX
AC_PROG_CXXCPP
AC_PROG_CC
AC_PROG_CPP

AX_CXX_COMPILE_STDCXX([11])

AC_CHECK_FUNCS([strchr strdup strtol])
AC_CHECK_HEADERS([limits.h])
AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_FUNC_STRTOD
AC_TYPE_SIZE_T

# Search for the RNAlib2 and set variables for Makefiles if found.
# _CFLAGS and _LIBS can be overriden with the module name from
# PKG_CHECK_MODULES
PKG_CHECK_MODULES([VRNA], [RNAlib2],
                  [ AC_SUBST([VRNA_CFLAGS])
                    AC_SUBST([VRNA_LIBS])
                  ],
                  [])

AC_OUTPUT(Makefile)

