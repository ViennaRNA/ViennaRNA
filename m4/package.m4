# This file is part of Autoconf.                       -*- Autoconf -*-

# ViennaRNA Package 2011 Ronny Lorenz
#

##-----------------##
## Private macros. ##
##-----------------##

AC_DEFUN([AC_RNA_TEST_FILE],[
AC_MSG_CHECKING([for $1])
if test -f $1 ; then
  AC_MSG_RESULT([yes])
  $2
else
  AC_MSG_RESULT([no])
  $3
fi
])

AC_DEFUN([AC_RNA_TEST_DIR],[
AC_MSG_CHECKING([whether directory \"$1\" exists])
if test -d $1 ; then
  AC_MSG_RESULT([yes])
  $2
else
  AC_MSG_RESULT([no])
  $3
fi
])

AC_DEFUN([AC_RNA_PACKAGE_IF_ENABLED],[
if test "x$with_$1" != xno; then
  $2
fi
])

#
# AC_RNA_WITH_PACKAGE(package-name, package-dir, package-description, default)
#
# This macro handles additional package inclusion
# Parameters:
#       package-name:       a lowercase name of the optional package
#                           which is used for $with_package_name variables
#                           and --with[out]-package-name options in
#                           configure script
#
#       package-dir:        directory path to the optional package sources
#
#       package-desciption: a very brief description used for the package
#                           specific help output in configure script
#
#       default:            activate the package as default
#                           Values: "yes" or "no"
#
# Example: AC_RNA_WITH_PACKAGE([foo], [path/to/foo], [the incredible Foo program], yes)
#

AC_DEFUN([AC_RNA_WITH_PACKAGE],[

# announce the option to include it in configure script
AC_ARG_WITH([$1],
            [ifelse($4, yes,
              [AS_HELP_STRING([--without-$1], [don't build $3])],
              [AS_HELP_STRING([--with-$1], [build $3])])])

# check if enabling the package makes sense at configure-time
# and deactivate it if not


AC_RNA_PACKAGE_IF_ENABLED([$1],[
  AC_RNA_TEST_DIR([$2],
    [with_$1=$with_$1],
    [with_$1="no"])
])

AC_RNA_PACKAGE_IF_ENABLED([$1],[
  for i in $5; do
    AC_RNA_TEST_FILE([$2/$i],
      [with_$1=$with_$1],
      [with_$1="no"])
  done
])

])


##----------------##
## Public macros. ##
##----------------##



AC_DEFUN([AC_RNA_INIT],[

SVM_VERSION=2.91

AC_RNA_WITH_PACKAGE([perl],
                    [Perl],
                    [Perl module],
                    [yes],
                    [Makefile.am])
AC_RNA_WITH_PACKAGE([kinfold],
                    [Kinfold],
                    [Kinfold program],
                    [yes],
                    [Makefile.am])
AC_RNA_WITH_PACKAGE([forester],
                    [RNAforester],
                    [RNAforester program],
                    [yes],
                    [Makefile.am])
AC_RNA_WITH_PACKAGE([cluster],
                    [Cluster],
                    [AnalyseSeqs and AnalyseDists],
                    [no],
                    [Makefile.am])
AC_RNA_WITH_PACKAGE([svm],
                    [libsvm-${SVM_VERSION}],
                    [svm classifiers],
                    [yes],
                    [svm.cpp svm.h])

AC_SUBST( [SVM_SOURCE_ARCHIVE],
          [libsvm-${SVM_VERSION}.tar.gz])

AC_PATH_PROG(PerlCmd, perl)
if test -n "$PerlCmd"; then
  if $PerlCmd -e 'require 5.004'; then :
  else
     AC_MSG_WARN([You need Perl version 5.004 or higher for the Perl module])
     with_perl="no"
   fi
fi

if test -z "$PerlCmd"; then
    AC_MSG_WARN([No suitable Perl found -- will not build Perl module])
    AC_MSG_WARN([You may set the PerlCmd environment variable to point to
    a suitable perl binary])
    with_perl="no"
fi

AC_RNA_PACKAGE_IF_ENABLED([perl],[
  AC_CONFIG_SUBDIRS([Perl])
  AC_CONFIG_FILES([Perl/Makefile Perl/Makefile.PL])
])

AC_RNA_PACKAGE_IF_ENABLED([kinfold],[
  AC_CONFIG_SUBDIRS([Kinfold])
])

AC_RNA_PACKAGE_IF_ENABLED([forester],[
  AC_CONFIG_SUBDIRS([RNAforester])
])

AC_RNA_PACKAGE_IF_ENABLED([cluster],[
  AC_CONFIG_SUBDIRS([Cluster])
  AC_CONFIG_FILES([Cluster/Makefile])
])

AC_RNA_PACKAGE_IF_ENABLED([svm],[
  AC_MSG_NOTICE([performing libSVM specific checks])
  # Checks for header files.
  AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h strings.h unistd.h])
  # Checks for typedefs, structures, and compiler characteristics.
  AC_HEADER_STDBOOL
  # Checks for library functions.
  AC_FUNC_MALLOC
  AC_FUNC_REALLOC
  AC_FUNC_STRTOD
  AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strdup strndup strrchr strstr strtol strtoul])
  AC_SUBST([CXXLD],[${CXX}]) # this is rather a hack for RNALfold.c linking correctly
])


AM_CONDITIONAL(MAKE_PERL_EXT, test "$with_perl" != "no")
AM_CONDITIONAL(MAKE_KINFOLD, test "$with_kinfold" != "no")
AM_CONDITIONAL(MAKE_FORESTER, test "$with_forester" != "no")
AM_CONDITIONAL(MAKE_CLUSTER, test "$with_cluster" = "yes")
AM_CONDITIONAL(WITH_LIBSVM, test "$with_svm" != "no")


AC_CONFIG_FILES([Makefile ViennaRNA.spec Utils/Makefile Progs/Makefile lib/Makefile man/Makefile H/Makefile])

AC_CONFIG_FILES([man/cmdlopt.sh],[chmod +x man/cmdlopt.sh])

])

AC_DEFUN([AC_RNA_NOTICE],[

# get directory paths

eval _bindir=$(eval echo $bindir)
eval _libdir=$(eval echo $libdir)
eval _includedir=${includedir}
eval _datadir=$datadir
eval _mandir=$mandir
eval _docdir=$docdir

# Notify the user

AC_MSG_NOTICE(
[
Configure successful with the following options:

  Perl Extension:       ${with_perl:-yes}
  Analyse{Dists,Seqs}:  ${with_cluster:-no}
  Kinfold:              ${with_kinfold:-yes}
  RNAforester:          ${with_forester:-yes}
  SVM:                  ${with_svm:-yes}
  Documentation:        ${with_documentation:-no}
    (HTML):             ${with_documentation_html:-no}
    (PDF):              ${with_documentation_pdf:-no}
-
Files will be installed in the following directories:

  Executables:    $_bindir
  Libraries:      $_libdir
  Header files:   $_includedir
  Extra Data:     $_datadir
  Man pages:      $_mandir
  Documentation:  $_docdir
])

])



