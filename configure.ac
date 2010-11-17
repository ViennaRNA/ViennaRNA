dnl Process this file with autoconf to produce a configure script.
AC_INIT([ViennaRNA],[2.0.0],[rna@tbi.univie.ac.at],[ViennaRNA])
dnl Every other copy of the package version number gets its value from here
AM_INIT_AUTOMAKE
dnl change the default installation prefix
dnl AC_PREFIX_DEFAULT(/usr/local/RNA)

# libsvm version number (there should be a libsvm-${SVM_VERSION}.tar.gz in the ViennaRNA directory)
SVM_VERSION=2.91
AC_SUBST([SVM_SOURCE_ARCHIVE], [libsvm-${SVM_VERSION}.tar.gz])


dnl configure options
AM_WITH_DMALLOC
AC_ARG_WITH(cluster, [  --with-cluster          build AnalyseSeqs and AnalyseDists])
AC_ARG_WITH(perl, [  --without-perl          don't build Perl module])
AC_ARG_WITH(forester, [  --without-forester          don't build RNAforester program])
AC_ARG_WITH(kinfold, [  --without-kinfold          don't build Kinfold program])
AC_ARG_WITH(svm, [  --without-svm          don't use svm classifiers])
AC_ARG_ENABLE(floatpf, [  --enable-floatpf        use float instead of double for partition function], [])

AC_ARG_VAR(PERLPREFIX, [Prefix for installation of the Perl module. 
By default the perl module is installed wherever Perl thinks it belongs, ignoring the $prefix used by configure])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AM_CONFIG_HEADER([config.h])

AC_SUBST(VERSION)

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
AM_PROG_CC_C_O

dnl check for auxilary programs to compile

AC_PATH_PROG(PerlCmd, perl)
if test -n "$PerlCmd"; then
  if $PerlCmd -e 'require 5.004'; then :
  else
     AC_MSG_RESULT(You need Perl version 5.004 or higher for the Perl module)
     with_perl="no"
   fi
fi

if test -z "$PerlCmd"; then
    AC_MSG_RESULT(No suitable Perl found -- will not build Perl module)
    AC_MSG_RESULT(You may set the PerlCmd environment variable to point to
    a suitable perl binary)
    with_perl="no"
fi

if test "x${with_perl}" != xno; then
  if test -f $srcdir/Perl/Makefile.am; then
    AC_CONFIG_SUBDIRS([Perl])
  else
    AC_MSG_WARN([Perl subdirectory does not exist])
    with_perl="no"
  fi
fi

if test "$with_kinfold" != "no"; then
   if test -f $srcdir/Kinfold/Makefile.am; then
     AC_CONFIG_SUBDIRS([Kinfold])
   else
     AC_MSG_WARN([Kinfold subdirectory does not exist])
     with_kinfold="no"
   fi
fi

if test "$with_forester" != "no"; then
   if test -f $srcdir/RNAforester/Makefile.am; then
     AC_CONFIG_SUBDIRS([RNAforester])
   else
     AC_MSG_WARN([RNAforester subdirectory does not exist])
     with_forester="no"
  fi
fi

if test "$with_cluster" != "no"; then
   if test -f $srcdir/Cluster/Makefile.am; then
     AC_CONFIG_SUBDIRS([Cluster])
   else
     AC_MSG_WARN([Cluster subdirectory does not exist])
     with_cluster="no"
   fi
fi


dnl check for necessity of compiling provided SVM library
if test "x${with_svm}" != xno; then
  AC_MSG_NOTICE([performing libSVM specific checks])
  # Checks for header files.
  AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h strings.h unistd.h])
  # Checks for typedefs, structures, and compiler characteristics.
  AC_HEADER_STDBOOL
  AC_C_INLINE
  AC_TYPE_SIZE_T
  # Checks for library functions.
  AC_FUNC_MALLOC
  AC_FUNC_REALLOC
  AC_FUNC_STRTOD
  AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strdup strndup strrchr strstr strtol strtoul])
  AC_MSG_CHECKING([whether provided libSVM source is present...])
  UNCOMPRESS_LIBSVM=0
  if test -f $srcdir/libsvm-${SVM_VERSION}/svm.cpp; then
    AC_SUBST([LIBSVM_DIR], [$srcdir/libsvm-${SVM_VERSION}/svm.cpp])
  else
    UNCOMPRESS_LIBSVM=1
  fi
  if test -f $srcdir/libsvm-${SVM_VERSION}/svm.h; then
    AC_SUBST([LIBSVM_DIR],[$srcdir/libsvm-${SVM_VERSION}])
  else
    UNCOMPRESS_LIBSVM=1
  fi
  if test "x${UNCOMPRESS_LIBSVM}" = x1; then
    if test -f $srcdir/libsvm-${SVM_VERSION}.tar.gz; then
      AC_MSG_RESULT(yes)
      AC_PATH_PROG([TAR], [tar], [])
      if test "x${TAR}" = x; then
        AC_MSG_NOTICE([==> tar command not found!])
        AC_MSG_ERROR([the tar command is mandatory for unpacking libsvm sources])
      fi
      AC_MSG_NOTICE([unpacking libSVM source...])
      ${TAR} xfz $srcdir/libsvm-${SVM_VERSION}.tar.gz
      AC_SUBST([LIBSVM_DIR],[$srcdir/libsvm-${SVM_VERSION}])
      AC_DEFINE([USE_SVM], [1], [use support vector machine])
    else
      AC_MSG_RESULT(no)
      AC_MSG_WARN([no suitable libsvm package found, deactivating SVM capabilities])
      with_svm="no"
    fi
  else
    AC_DEFINE([USE_SVM], [1], [use support vector machine])
    AC_MSG_RESULT(yes)
  fi
  AC_SUBST([CXXLD],[${CXX}]) # this is rather a hack for RNALfold.c linking correctly
fi

AC_OPENMP

AM_CONDITIONAL(MAKE_PERL_EXT, test "$with_perl" != "no")

AM_CONDITIONAL(MAKE_CLUSTER, test "$with_cluster" = "yes")

AM_CONDITIONAL(MAKE_FORESTER, test "$with_forester" != "no")

AM_CONDITIONAL(MAKE_KINFOLD, test "$with_kinfold" != "no")

AM_CONDITIONAL(WITH_LIBSVM, test "$with_svm" != "no")



dnl Checks for libraries.
dnl AC_REPLACE_GNU_GETOPT

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS(malloc.h strings.h unistd.h math.h)
AC_CHECK_FUNCS([memset memmove])

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_TYPE_SIZE_T
AC_C_INLINE([])

dnl Checks for library functions.
AC_CHECK_FUNCS(strdup strstr strchr erand48)

AC_CONFIG_FILES([Makefile ViennaRNA.spec Utils/Makefile Progs/Makefile lib/Makefile man/Makefile Cluster/Makefile H/Makefile Perl/Makefile])
AC_CONFIG_FILES([man/cmdlopt.sh],[chmod +x man/cmdlopt.sh])
AC_OUTPUT

eval _bindir=$(eval echo $bindir)
eval _libdir=$(eval echo $libdir)
eval _includedir=${includedir}
eval _datadir=$datadir
eval _mandir=$mandir

AC_MSG_NOTICE(
[
Configure successful with the following options:

  Perl Extension:      ${with_perl:-yes}
  Analyse{Dists,Seqs}: ${with_cluster:-no}
  Kinfold:             ${with_kinfold:-yes}
  RNAforester:         ${with_forester:-yes}
  SVM:                 ${with_svm:-yes}

Files will be installed in the following directories:

  Executables:  $_bindir
  Libraries:    $_libdir
  Header files: $_includedir
  Extra Data:   $_datadir
  Man pages:    $_mandir
])
