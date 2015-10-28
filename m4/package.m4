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

AC_DEFUN([AC_RNA_PACKAGE_IF_ENABLED],[
if test "x$with_$1" != "xno"; then
  $2
fi
])

AC_DEFUN([AC_RNA_FEATURE_IF_ENABLED],[
if test "x$enable_$1" != "xno"; then
  $2
fi
])

#
# AC_RNA_ADD_PACKAGE( package-name,
#                     package-description,
#                     default-on,
#                     [action-if-not-default],
#                     [action-if-default],
#                     [files to check for])
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
            [ifelse([$3], [yes],
              [AS_HELP_STRING([--without-m4_translit([$1], [_], [-])], [do not build | install $2])],
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

AC_DEFUN([AC_RNA_ADD_FEATURE],[

# announce the feature for inclusion in the configure script

AC_ARG_ENABLE(m4_translit([[$1]], [_], [-]),
            [ifelse([$3], [yes],
              [AS_HELP_STRING([--disable-m4_translit([$1], [_], [-])], [disable feature: $2])],
              [AS_HELP_STRING([--enable-m4_translit([$1], [_], [-])], [enable feature: $2])])],
            [$4],
            [$5])


])

# AC_RNA_DOCUMENTATION_INIT(PROJECT_NAME, [config-file], [documentation-output-directory])
#
#
AC_DEFUN([AC_RNA_DOCUMENTATION_INIT],[

AC_PATH_PROG(doxygen, [doxygen],no)
AC_PATH_PROG(pdflatex,[pdflatex],no)
AC_PATH_PROG(latex,[latex],no)
AC_PATH_PROG(makeindex,[makeindex],no)
AC_PATH_PROG(dot,[dot],no)
AC_PATH_PROG(egrep,[egrep],no)
AC_PATH_PROG(perl,[perl],no)

DOXYGEN_PDFLATEX_WORKARROUND=yes

# check whether we are able to generate the doxygen documentation
AC_RNA_PACKAGE_IF_ENABLED([doc],[
  if test "x$doxygen" != xno;
  then
    # test for programs necessary in order to use doxygen

    if test "x$DOXYGEN_PDFLATEX_WORKARROUND" = xno;
    then
    # this is a workarround for older versions of doxygen as installed e.g. in fc12 where
    # pdflatex usage does not work

      if test "x$pdflatex" = xno;
      then
        if test "x$latex" = xno;
        then
          AC_MSG_WARN([neither latex or pdflatex exists on your system!])
          AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
          doxygen=no
        else
          _latex_cmd=$latex
        fi
      else
        _latex_cmd=$pdflatex
      fi
    else
      if test "x$latex" = xno;
      then
        AC_MSG_WARN([neither latex or pdflatex exists on your system!])
        AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
        doxygen=no
        _latex_cmd=
      else
        AC_MSG_WARN([due to a bug in older versions of doxygen, latex will be used for reference manual generation even if pdflatex is available])
        _latex_cmd=$latex
        pdflatex=no
      fi
    fi

    if test "x$makeindex" = xno;
    then
      AC_MSG_WARN([makeindex command not found on your system!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen=no
    fi

    if test "x$egrep" = xno;
    then
      AC_MSG_WARN([egrep command not found on your system!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen=no
    fi

    if test "x$dot" = xno;
    then
      AC_MSG_WARN([dot command not found on your system!])
      AC_MSG_WARN([deactivating graph output in reference manual!])
    fi

    if test "x$perl" = xno;
    then
      AC_MSG_WARN([perl command not found on your system!])
      AC_MSG_WARN([deactivating automatic (re)generation of reference manual!])
      doxygen=no
    fi

  fi
])


# setup everything in order to generate the doxygen configfile

AC_RNA_PACKAGE_IF_ENABLED([doc],[

  AC_SUBST([DOXYGEN_PROJECT_NAME], [$1-$PACKAGE_VERSION])
  AC_SUBST([DOXYGEN_SRCDIR], [$srcdir])
  AC_SUBST([DOXYGEN_DOCDIR], [ifelse([$3], [], [doc], [$3])])
  AC_SUBST([DOXYGEN_CONF], [ifelse([$2], [], [doxygen.conf], [$2])])


# prepare the config file for doxygen if we are able to generate a reference manual
  if test "x$doxygen" != xno;
  then

    AC_SUBST([DOXYGEN_CMD_LATEX], [$_latex_cmd])
    AC_SUBST([DOXYGEN_CMD_MAKEINDEX], [$makeindex])
    AC_SUBST([DOXYGEN_HAVE_DOT],[ifelse([$dot], [no], [NO], [YES])])
    AC_SUBST([DOXYGEN_WITH_PDFLATEX], [ifelse([$pdflatex],[no],[NO],[YES])])
    AC_SUBST([DOXYGEN_GENERATE_HTML], [ifelse([$with_doc_html], [no], [NO], [YES])])
    AC_SUBST([DOXYGEN_GENERATE_LATEX], [ifelse([$with_doc_pdf], [no], [NO], [YES])])

    AC_CONFIG_FILES([${DOXYGEN_DOCDIR}/${DOXYGEN_CONF}])

  else

# otherwise check if a generated reference manual already exists

    AC_RNA_PACKAGE_IF_ENABLED([doc_pdf],[
      AC_RNA_TEST_FILE( [$DOXYGEN_DOCDIR/$DOXYGEN_PROJECT_NAME.pdf],
                        [with_doc_pdf=yes],
                        [with_doc_pdf=no])])

    AC_RNA_PACKAGE_IF_ENABLED([doc_html],[
      AC_RNA_TEST_FILE( [$DOXYGEN_DOCDIR/html/index.html],
                        [with_doc_html=yes],
                        [with_doc_html=no])])

    if test "x$with_doc_pdf" = "x$with_doc_html";
    then
      if test "x$with_doc_pdf" = xno;
      then
        with_doc=no
      fi
    fi
  fi
])

AC_SUBST([REFERENCE_MANUAL_PDF_NAME], [ifelse([$with_doc_pdf],
                                              [no],
                                              [],
                                              [$DOXYGEN_PROJECT_NAME.pdf])])
AC_SUBST([REFERENCE_MANUAL_TAGFILE],  [ifelse([$doxygen],
                                              [no],
                                              [],
                                              [$DOXYGEN_PROJECT_NAME.tag])])


# setup variables used in Makefile.am

# Define ${htmldir} if the configure script was created with a version of
# autoconf older than 2.60
# Alternatively, if ${htmldir} is exactly '${docdir}', append a /html to
# separate html files from rest of doc.
# Otherwise, just append the PACKAGE_NAME to the htmldir
if test "x${htmldir}" = "x";
then
  AC_MSG_WARN([resetting htmldir])
  htmldir="${docdir}/html"
fi

if test "x${htmldir}" = 'x${docdir}';
then
  htmldir="${docdir}/html"
else
  htmldir=${htmldir}/${PACKAGE_NAME}
fi

AC_SUBST(htmldir, [${htmldir}])

#

AM_CONDITIONAL(WITH_REFERENCE_MANUAL, test "x$with_doc" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_BUILD, test "x$doxygen" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_PDF, test "x$with_doc_pdf" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_HTML, test "x$with_doc_html" != xno)

])

##----------------##
## Public macros. ##
##----------------##

AC_DEFUN([AC_RNA_INIT],[

SVM_VERSION=3.20
with_pf_float=no

dnl add packages to the configure process

AC_RNA_ADD_PACKAGE( [perl],
                    [Perl interface],
                    [yes],
                    [with_perl=no],
                    [with_perl=yes],
                    [${srcdir}/interfaces/Perl/Makefile.am])
AC_RNA_ADD_PACKAGE( [python],
                    [Python interface],
                    [yes],
                    [with_python=no],
                    [with_python=yes],
                    [${srcdir}/interfaces/Python/Makefile.am])
##AC_RNA_ADD_PACKAGE( [ruby],
##                    [Ruby interface],
##                    [no],
##                    [with_ruby=yes],
##                    [with_ruby=no],
##                    [${srcdir}/interfaces/Ruby/Makefile.am])
AC_RNA_ADD_PACKAGE( [kinfold],
                    [Kinfold program],
                    [yes],
                    [with_kinfold=no],
                    [with_kinfold=yes],
                    [${srcdir}/src/Kinfold/Makefile.am])
AC_RNA_ADD_PACKAGE( [forester],
                    [RNAforester program],
                    [yes],
                    [with_forester=no],
                    [with_forester=yes],
                    [${srcdir}/src/RNAforester/Makefile.am])
AC_RNA_ADD_PACKAGE( [cluster],
                    [AnalyseSeqs and AnalyseDists],
                    [no],
                    [with_cluster=yes],
                    [with_cluster=no],
                    [${srcdir}/src/Cluster/Makefile.am])
AC_RNA_ADD_PACKAGE( [svm],
                    [svm classifiers],
                    [yes],
                    [with_svm=no],
                    [with_svm=yes],
                    [${srcdir}/src/libsvm-${SVM_VERSION}/svm.cpp ${srcdir}/src/libsvm-${SVM_VERSION}/svm.h])
AC_RNA_ADD_PACKAGE( [json],
                    [json in/out support],
                    [yes],
                    [with_json=no],
                    [with_json=yes],
                    [${srcdir}/src/json/json.c ${srcdir}/src/json/json.h])
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
AC_RNA_ADD_PACKAGE( [check],
                    [Unit tests],
                    [yes],
                    [with_check=no],
                    [with_check=yes])
AC_RNA_ADD_PACKAGE( [gsl],
                    [GNU Scientific Library],
                    [yes],
                    [with_gsl=no],
                    [with_gsl=yes])
AC_RNA_ADD_FEATURE( [boustrophedon],
                    [Boustrophedon scheme for stochastic backtracking],
                    [yes],
                    [enable_boustrophedon=no],
                    [enable_boustrophedon=yes])
AC_RNA_ADD_FEATURE( [gen_hard_constraints],
                    [Generic hard constraints],
                    [no],
                    [enable_gen_hard_constraints=yes],
                    [enable_gen_hard_constraints=no])

## begin with initialization according to configure-time specific options

## The following test ensures the right type for FLT_OR_DBL in the SWIG RNAlib interface
AC_MSG_CHECKING([whether float precision is used for partition function arrays instead of double precision])
bla=`${GREP} "^#define LARGE_PF" ${srcdir}/src/ViennaRNA/dp_matrices.h`
if test "x$bla" = "x";
then
  with_pf_float=yes
fi
AC_MSG_RESULT([$with_pf_float])
AM_CONDITIONAL([WITH_LARGE_PF], [test "$with_pf_float" != "yes"])

# check prerequisites for Perl interface
if test "x$PerlCmd" = "x";
then
  AC_PATH_PROG(PerlCmd, perl)
fi
ifelse([$PerlCmd], [],[
  AC_MSG_WARN([No suitable Perl found -- will not build Perl module])
  AC_MSG_WARN([You may set the PerlCmd environment variable to point to
              a suitable perl binary])
  with_perl="no"
],[
  if $PerlCmd -e 'require 5.004'; then :
  else
    AC_MSG_WARN([You need Perl version 5.004 or higher for the Perl module])
    with_perl="no"
  fi
])


# prepare all files for perl interface
AC_RNA_PACKAGE_IF_ENABLED([perl],[
  AC_DEFINE([WITH_PERL_INTERFACE], [1], [Create the perl interface to RNAlib])
  AC_SUBST([PERL_INTERFACE], [Perl])
  AC_CONFIG_FILES([interfaces/Perl/Makefile interfaces/Perl/Makefile.PL])
])

# check prerequisites for Python interface
AC_RNA_PACKAGE_IF_ENABLED([python],[
if test "x$PythonCmd" = "x";
then
  AC_PATH_PROG(PythonCmd, python)
fi
ifelse([$PythonCmd], [], [
  AC_MSG_WARN([No suitable Python found -- will not build Python extension])
  AC_MSG_WARN([You may set the PythonCmd environment variable to point to
                a suitable python binary])
  with_python="no"
  enabled_but_failed_python="(No python executable found)"
],[
  version_test=`$PythonCmd ${srcdir}/interfaces/Python/version_test.py`
  if test "x$version_test" = "xok"; then :
  else
    AC_MSG_WARN([You need Python >= 2.5 and < 3.0 to build the Python extension])
    AC_MSG_WARN([You may set the PythonCmd environment variable to point to
                a suitable python binary])
    with_python="no"
    enabled_but_failed_python="(Python executable is of wrong version)"
  fi
])
])

# prepare all files for python interface
AC_RNA_PACKAGE_IF_ENABLED([python],[
  AC_DEFINE([WITH_PYTHON_INTERFACE], [1], [Create the python interface to RNAlib])
  AC_SUBST([PYTHON_INTERFACE], [Python])
  AC_CONFIG_FILES([interfaces/Python/Makefile interfaces/Python/setup.py])
])


AC_RNA_PACKAGE_IF_ENABLED([kinfold],[
  AC_CONFIG_SUBDIRS([src/Kinfold])
])

AC_RNA_PACKAGE_IF_ENABLED([forester],[
  AC_CONFIG_SUBDIRS([src/RNAforester])
])

AC_RNA_PACKAGE_IF_ENABLED([cluster],[
  AC_DEFINE([WITH_CLUSTER], [1], [Analyse{Dists,Seqs}])
  AC_SUBST([CLUSTER_DIR], [Cluster])
  AC_CONFIG_FILES([src/Cluster/Makefile])
])

AC_RNA_PACKAGE_IF_ENABLED([svm],[
  AC_SUBST([LIBSVM_DIR], [libsvm-${SVM_VERSION}])
  AC_SUBST([WITH_SVM], [USE_SVM])
  AC_SUBST([CXXLD],[${CXX}]) # this is rather a hack for RNALfold.c linking correctly
  AC_DEFINE([USE_SVM], [1], [Compute z-scores for RNALfold])
])

AC_RNA_PACKAGE_IF_ENABLED([json],[
  AC_SUBST([WITH_JSON], [WITH_JSON_SUPPORT])
  AC_DEFINE([WITH_JSON_SUPPORT], [1], [Add JSON support for input and output functions])
])

# check prerequisties for unit testing support
AC_RNA_PACKAGE_IF_ENABLED([check],[
PKG_CHECK_MODULES([CHECK], [check >= 0.9.4], [],[
  AC_MSG_WARN([check not found -- will not build Unit tests])
  with_check="no"
  enabled_but_failed_check="(check framework not found)"
])
])

AC_RNA_PACKAGE_IF_ENABLED([check],[
  AC_DEFINE([WITH_CHECK], [1], [Include Unit tests])
  AC_SUBST([CHECK_DIR], [tests])
  AC_CONFIG_FILES([tests/Makefile])
])

# check prerequisties for gsl support
AC_RNA_PACKAGE_IF_ENABLED([gsl],[
  AC_CHECK_LIB([m],[cos])
  AC_CHECK_LIB([gslcblas],[cblas_dgemm])
  AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

  if test "$ac_cv_lib_gsl_gsl_blas_dgemm" != yes; then
    AC_MSG_WARN("Can't find libgsl. Falling back to default implementation.")
    enabled_but_failed_gsl="(libgsl not found)"
    with_gsl="no"
  fi
])

AC_RNA_PACKAGE_IF_ENABLED([gsl],[
  AC_DEFINE([WITH_GSL], [1], [Use GNU Scientific Library])
])

## Add preprocessor define statement for Boustrophedon scheme in stochastic backtracking in part_func.c
AC_RNA_FEATURE_IF_ENABLED([boustrophedon],[
  AC_DEFINE([WITH_BOUSTROPHEDON], [1], [Use Boustrophedon scheme for stochastic backtracking])
])

## Add preprocessor define statement for generlaized hard constraints feature
AC_RNA_FEATURE_IF_ENABLED([gen_hard_constraints],[
  AC_DEFINE([WITH_GEN_HC], [1], [Provide generic hard constraints])
  GENERIC_HC_DEF=-DWITH_GEN_HC
])
AC_SUBST(GENERIC_HC_DEF)

## Add linker flag for OpenMP in pkg-config file
AC_RNA_FEATURE_IF_ENABLED([openmp],[
  LIBGOMPFLAG=-lgomp
])
AC_SUBST(LIBGOMPFLAG)

AM_CONDITIONAL(MAKE_KINFOLD, test "$with_kinfold" != "no")
AM_CONDITIONAL(MAKE_FORESTER, test "$with_forester" != "no")
AM_CONDITIONAL(MAKE_CLUSTER, test "$with_cluster" = "yes")
AM_CONDITIONAL(WITH_LIBSVM, test "$with_svm" != "no")
AM_CONDITIONAL(WITH_JSON, test "$with_json" != "no")
AM_CONDITIONAL(WITH_CHECK, test "$with_check" != "yes")
AM_CONDITIONAL(WITH_GSL, test "$with_gsl" != "no")


AC_RNA_DOCUMENTATION_INIT([RNAlib])

AC_CONFIG_FILES([misc/Makefile misc/ViennaRNA.spec misc/PKGBUILD])
AC_CONFIG_FILES([interfaces/Makefile])
AC_CONFIG_FILES([Makefile RNAlib2.pc src/Utils/Makefile src/bin/Makefile src/Makefile man/Makefile src/ViennaRNA/Makefile doc/Makefile])

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
Configured successful with the following options:

RNAlib Scripting Interfaces:
  Perl Interface:           ${with_perl:-yes}       $enabled_but_failed_perl
  Python Interface:         ${with_python:-yes}     $enabled_but_failed_python

Extra Programs:
  Analyse{Dists,Seqs}:      ${with_cluster:-no}
  Kinfold:                  ${with_kinfold:-yes}
  RNAforester:              ${with_forester:-yes}

Other Options:
  SVM:                      ${with_svm:-yes}
  JSON:                     ${with_json:-yes}
  GSL:                      ${with_gsl:-yes}        $enabled_but_failed_gsl
  Boustrophedon:            ${enable_boustrophedon:-yes}
  Generic Hard Constraints: ${enable_gen_hard_constraints:-no}
  OpenMP:                   ${enable_openmp:-yes}

Documentation:              ${with_doc:-no}
    (HTML):                 ${with_doc_html:-no}
    (PDF):                  ${with_doc_pdf:-no}

Unit Tests:
  check:                    ${with_check:-yes}      $enabled_but_failed_check

-
Files will be installed in the following directories:

  Executables:    $_bindir
  Libraries:      $_libdir
  Header files:   $_includedir
  Extra Data:     $_datadir
  Man pages:      $_mandir
  Documentation:  $_docdir
    (HTML):       $(eval printf "%s" $_htmldir)
    (PDF):        $(eval printf "%s" $_pdfdir)
])

])



