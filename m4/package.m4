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
# AC_RNA_ADD_PACKAGE( package-name,
#                     package-description,
#                     default-on,
#                     [action-if-not-default],
#                     [action-if-default],
#                     [files to check])
#
# This macro handles additional package inclusion
# Parameters:
#       package-name:       a lowercase name of the optional package
#                           which is used for $with_package_name variables
#                           and --with[out]-package-name options in
#                           configure script
#                           The package_name must consist of alphanumeric
#                           characters including the dash only! Every
#                           occuring dash will be replaced by a '-' char
#                           in the --with[out]-package-name option
#
#       package-desciption: a very brief description used for the package
#                           specific help output in configure script
#
#       default-on:         package build | installed by default
#                           Values: "yes" or "no"
#
# Example: AC_RNA_ADD_PACKAGE([foo], [the incredible Foo program], [yes], [with_foo=no], [with_foo=yes], [file1 file2])
#

AC_DEFUN([AC_RNA_ADD_PACKAGE],[

# announce the option to include it in configure script
AC_ARG_WITH(m4_translit([[$1]], [_], [-]),
            [ifelse($3, yes,
              [AS_HELP_STRING([--without-m4_translit([$1], [_], [-])], [don't build | install $2])],
              [AS_HELP_STRING([--with-m4_translit([$1], [_], [-])], [build | install $2])])],
            [$4],
            [$5])

# check if enabling the package makes sense at configure-time
# and deactivate it if not

AC_RNA_PACKAGE_IF_ENABLED([$1],[
  for i in $6; do
    AC_RNA_TEST_FILE([$i],
      [with_$1=$with_$1],
      [with_$1=no])
  done
])

])

AC_DEFUN([AC_RNA_DOCUMENTATION_INIT],[

AC_CHECK_PROGS(doxygen, [doxygen],no)
AC_CHECK_PROGS(pdflatex,[pdflatex],no)
AC_CHECK_PROGS(latex,[latex],no)
AC_CHECK_PROGS(makeindex,[makeindex],no)
AC_CHECK_PROGS(dot,[dot],no)
AC_CHECK_PROGS(egrep,[egrep],no)

])

##----------------##
## Public macros. ##
##----------------##

AC_DEFUN([AC_RNA_INIT],[

SVM_VERSION=2.91

dnl add packages to the configure process

AC_RNA_ADD_PACKAGE( [perl],
                    [Perl module],
                    [yes],
                    [with_perl=no],
                    [with_perl=yes],
                    [Perl/Makefile.am])
AC_RNA_ADD_PACKAGE( [kinfold],
                    [Kinfold program],
                    [yes],
                    [with_kinfold=no],
                    [with_kinfold=yes],
                    [Kinfold/Makefile.am])
AC_RNA_ADD_PACKAGE( [forester],
                    [RNAforester program],
                    [yes],
                    [with_forester=no],
                    [with_forester=yes],
                    [RNAforester/Makefile.am])
AC_RNA_ADD_PACKAGE( [cluster],
                    [AnalyseSeqs and AnalyseDists],
                    [no],
                    [with_cluster=yes],
                    [with_cluster=no],
                    [Cluster/Makefile.am])
AC_RNA_ADD_PACKAGE( [svm],
                    [svm classifiers],
                    [yes],
                    [libsvm-${SVM_VERSION}/svm.cpp libsvm-${SVM_VERSION}/svm.h],
                    [with_svm=no],
                    [with_svm=yes])
AC_RNA_ADD_PACKAGE( [doc_pdf],
                    [PDF RNAlib reference manual],
                    [yes],
                    [with_doc_pdf=no],
                    [with_doc_pdf=yes])
AC_RNA_ADD_PACKAGE( [doc_html],
                    [HTML RNAlib reference manual],
                    [yes],
                    [with_doc_html=no],
                    [with_doc_html=yes])
AC_RNA_ADD_PACKAGE( [doc],
                    [RNAlib reference manual],
                    [yes],
                    [ with_doc=no
                      with_doc_pdf=no
                      with_doc_html=no],
                    [with_doc=yes])

dnl do option specific things
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
  AC_SUBST([CXXLD],[${CXX}]) # this is rather a hack for RNALfold.c linking correctly
])


AM_CONDITIONAL(MAKE_PERL_EXT, test "$with_perl" != "no")
AM_CONDITIONAL(MAKE_KINFOLD, test "$with_kinfold" != "no")
AM_CONDITIONAL(MAKE_FORESTER, test "$with_forester" != "no")
AM_CONDITIONAL(MAKE_CLUSTER, test "$with_cluster" = "yes")
AM_CONDITIONAL(WITH_LIBSVM, test "$with_svm" != "no")

AC_RNA_DOCUMENTATION_INIT


AC_CONFIG_FILES([Makefile ViennaRNA.spec Utils/Makefile Progs/Makefile lib/Makefile man/Makefile H/Makefile])

AC_CONFIG_FILES([man/cmdlopt.sh],[chmod +x man/cmdlopt.sh])

])

AC_DEFUN([AC_RNA_NOTICE],[

# get directory paths

eval _bindir=$(eval printf "%s" $bindir)
eval _libdir=$(eval printf "%s" $libdir)
eval _includedir=$(eval printf "%s" $includedir)
eval _datadir=$(eval printf "%s" $datadir)
eval _mandir=$(eval printf "%s" $mandir)
eval _docdir=$(eval printf "%s" $docdir)
eval _htmldir=$(eval printf "%s" $htmldir)
eval _pdfdir=$(eval printf "%s" $pdfdir)

# Notify the user

AC_MSG_NOTICE(
[
Configure successful with the following options:

  Perl Extension:       ${with_perl:-yes}
  Analyse{Dists,Seqs}:  ${with_cluster:-no}
  Kinfold:              ${with_kinfold:-yes}
  RNAforester:          ${with_forester:-yes}
  SVM:                  ${with_svm:-yes}
  Documentation:        ${with_doc:-no}
    (HTML):             ${with_doc_html:-no}
    (PDF):              ${with_doc_pdf:-no}
-
Files will be installed in the following directories:

  Executables:    $_bindir
  Libraries:      $_libdir
  Header files:   $_includedir
  Extra Data:     $_datadir
  Man pages:      $_mandir
  Documentation:  $_docdir
    (HTML):       $_htmldir
    (PDF):        $_pdfdir
])

])



