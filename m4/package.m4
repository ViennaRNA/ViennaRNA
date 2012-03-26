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

    AC_CONFIG_FILES([${DOXYGEN_CONF}])

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
                                              [$DOXYGEN_DOCDIR/$DOXYGEN_PROJECT_NAME.pdf])])
AC_SUBST([REFERENCE_MANUAL_TAGFILE],  [ifelse([$doxygen],
                                              [no],
                                              [],
                                              [$DOXYGEN_DOCDIR/$DOXYGEN_PROJECT_NAME.tag])])


# setup variables used in Makefile.am
AM_CONDITIONAL(WITH_REFERENCE_MANUAL, test "x$with_doc" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_BUILD, test "x$doxygen" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_PDF, test "x$with_doc_pdf" != xno)
AM_CONDITIONAL(WITH_REFERENCE_MANUAL_HTML, test "x$with_doc_html" != xno)

])

##----------------##
## Public macros. ##
##----------------##

AC_DEFUN([AC_RNA_INIT],[

SVM_VERSION=2.91
with_pf_float=no

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
                    [with_svm=no],
                    [with_svm=yes],
                    [libsvm-${SVM_VERSION}/svm.cpp libsvm-${SVM_VERSION}/svm.h])
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

## begin with initialization according to configure-time specific options

AC_PATH_PROG(PerlCmd, perl)
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


AC_RNA_PACKAGE_IF_ENABLED([perl],[
  ## The following test ensures the right type for FLT_OR_DBL in the SWIG RNAlib interface
  AC_MSG_CHECKING([whether float precision is used for partition function arrays instead of double precision])
  bla=`${GREP} "^#define LARGE_PF" H/data_structures.h`
  if test "x$bla" = "x";
  then
    with_pf_float=yes
  fi
  AC_MSG_RESULT([$with_pf_float])
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
  AC_SUBST([LIBSVM_DIR], [libsvm-${SVM_VERSION}])
  AC_SUBST([CXXLD],[${CXX}]) # this is rather a hack for RNALfold.c linking correctly
  AC_DEFINE([USE_SVM], [1], [Compute z-scores for RNALfold])
])

AM_CONDITIONAL(WITH_LARGE_PF, test "$with_pf_float" != "yes")
AM_CONDITIONAL(MAKE_PERL_EXT, test "$with_perl" != "no")
AM_CONDITIONAL(MAKE_KINFOLD, test "$with_kinfold" != "no")
AM_CONDITIONAL(MAKE_FORESTER, test "$with_forester" != "no")
AM_CONDITIONAL(MAKE_CLUSTER, test "$with_cluster" = "yes")
AM_CONDITIONAL(WITH_LIBSVM, test "$with_svm" != "no")

AC_RNA_DOCUMENTATION_INIT([RNAlib])


AC_CONFIG_FILES([Makefile ViennaRNA.spec RNAlib2.pc Utils/Makefile Progs/Makefile lib/Makefile man/Makefile H/Makefile])

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
    (HTML):       $(eval printf "%s" $_htmldir)/html
    (PDF):        $(eval printf "%s" $_pdfdir)
])

])



