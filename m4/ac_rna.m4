# This file is part of Autoconf.                       -*- Autoconf -*-

# ViennaRNA Package 2011 Ronny Lorenz
#


##----------------##
## Public macros. ##
##----------------##

AC_DEFUN([AC_RNA_INIT],[

AX_COMPILER_VENDOR
AC_CANONICAL_HOST

AC_ARG_VAR([PERL],[The perl interpreter])
AC_PATH_PROGS([PERL],[$PERL perl],[no])
AS_IF([test "$PERL" == "no"],[
  AC_MSG_ERROR([Perl is required to install and run the ViennaRNA Package])
])

##--------------------------##
## Additional Compile Flags ##
##--------------------------##

AC_LANG_PUSH([C])
AX_CHECK_COMPILE_FLAG([-fno-strict-aliasing], [
  AX_APPEND_FLAG(["-fno-strict-aliasing"], [RNA_CFLAGS])
],[],[],[])
AC_LANG_POP([C])

AC_LANG_PUSH([C++])
AX_CHECK_COMPILE_FLAG([-fno-strict-aliasing], [
  AX_APPEND_FLAG(["-fno-strict-aliasing"], [RNA_CXXFLAGS])
],[],[],[])
AC_LANG_POP([C++])

AX_CHECK_LINK_FLAG([-fno-strict-aliasing], [
  AX_APPEND_FLAG(["-fno-strict-aliasing"], [RNA_LDFLAGS])
], [],[], [])

##------------------##
## Enable Features  ##
##------------------##

RNA_ENABLE_C11
RNA_ENABLE_OSX
RNA_ENABLE_LTO
RNA_ENABLE_SVM
RNA_ENABLE_JSON
RNA_ENABLE_GSL
RNA_ENABLE_OPENMP
RNA_ENABLE_PTHREADS
RNA_ENABLE_BOUSTROPHEDON
RNA_ENABLE_NR_SAMPLE_HASH
RNA_ENABLE_FLOATPF
RNA_ENABLE_DEPRECATION_WARNINGS
RNA_ENABLE_COLORED_TTY
RNA_ENABLE_STATIC_BIN
RNA_ENABLE_SIMD
RNA_ENABLE_VECTORIZE
RNA_ENABLE_MPFR

## Set post conditions for feature
## settings

RNA_FEATURE_POST

##--------------------##
## Enable scripting   ##
## language interface ##
##--------------------##

RNA_ENABLE_SWIG_INTERFACES

##------------------##
## Enable Reference ##
## Manual           ##
##------------------##

RNA_ENABLE_DOXYGEN_REFMAN([RNAlib])

##--------------------##
## Enable Tutorial    ##
##--------------------##
RNA_ENABLE_TUTORIAL([tutorial])

##--------------------##
## Enable Unit tests  ##
##--------------------##
RNA_ENABLE_UNIT_TESTS

##----------------------------------##
## Check general build dependencies ##
##----------------------------------##
RNA_CHECK_BUILD_REQUIREMENTS

##--------------------##
## Enable Subpackages ##
##--------------------##

RNA_ENABLE_PKG_KINFOLD
RNA_ENABLE_PKG_FORESTER
RNA_ENABLE_PKG_CLUSTER
RNA_ENABLE_PKG_KINWALKER
RNA_ENABLE_PKG_RNALOCMIN

##--------------------##
## Prepare Version    ##
## Macros             ##
##--------------------##
VRNA_VERSION_MAJOR=vrna_version_major
VRNA_VERSION_MINOR=vrna_version_minor
VRNA_VERSION_PATCH=vrna_version_patch

AC_SUBST(VRNA_VERSION_MAJOR)
AC_SUBST(VRNA_VERSION_MINOR)
AC_SUBST(VRNA_VERSION_PATCH)

##------------------##
## Prepare files    ##
##------------------##

AC_CONFIG_FILES([src/ViennaRNA/vrna_config.h])
AC_CONFIG_FILES([misc/Makefile])
AC_CONFIG_FILES([interfaces/Makefile])
AC_CONFIG_FILES([Makefile RNAlib2.pc])
AC_CONFIG_FILES([src/Utils/Makefile src/bin/Makefile src/Makefile src/ViennaRNA/Makefile])
AC_CONFIG_FILES([src/ViennaRNA/static/Makefile])
AC_CONFIG_FILES([man/Makefile doc/Makefile RNA-Tutorial/Makefile])
AC_CONFIG_FILES([man/cmdlopt.sh],[chmod +x man/cmdlopt.sh])
AC_CONFIG_FILES([examples/Makefile])
AC_CONFIG_FILES([packaging/viennarna.spec packaging/PKGBUILD])
AC_CONFIG_FILES([packaging/win_installer_archlinux_i686.nsi packaging/win_installer_archlinux_x86_64.nsi])
AC_CONFIG_FILES([packaging/win_installer_fedora_i686.nsi packaging/win_installer_fedora_x86_64.nsi])

])

AC_DEFUN([AC_RNA_NOTICE],[

# get directory paths

eval _bindir=$(eval printf "%s" $bindir)
eval _libdir=$(eval printf "%s" $libdir)
eval _includedir=$(eval printf "%s" $includedir)
eval _datadir=$(eval printf "%s" $datadir)
eval _mandir=$(eval printf "%s" $mandir)

AS_IF([test $with_perl = "yes"],[
  eval _perl_arch_dir=$(eval printf "%s" "$prefix" ${PERL_ARCH_RELATIVE_INSTALL_DIR})
  eval _perl_lib_dir=$(eval printf "%s" "$prefix" ${PERL_LIB_RELATIVE_INSTALL_DIR})
  ], [
  _perl_arch_dir=""
  _perl_lib_dir=""
  _perl_install="Not to be installed"
])
AS_IF([test $with_python = "yes"],[
  eval _python2_arch_dir=$(eval printf "%s" ${py2execdir})
  eval _python2_lib_dir=$(eval printf "%s" ${python2dir})
  ],[
    _python2_arch_dir=""
    _python2_lib_dir=""
    _python2_install="Not to be installed $python2_enabled_but_failed"
])
AS_IF([test $with_python3 = "yes"],[
  eval _python3_arch_dir=$(eval printf "%s" ${py3execdir})
  eval _python3_lib_dir=$(eval printf "%s" ${python3dir})
  ],[
    _python3_arch_dir=""
    _python3_lib_dir=""
    _python3_install="Not to be installed $python3_enabled_but_failed"
])
AS_IF([test "x$enable_universal_binary" != "xno"],[
  _osx_arch=$osx_arch
  ],[
  osx_arch="no"
])

AS_IF([test "x$with_doc" != "xno"],[
  eval _docdir=$(eval printf "%s" $docdir)
  AS_IF([test "x$with_doc_html" != "xno"],[
    eval _htmldir=$(eval printf "%s" $htmldir)],[
    _htmldir=""
    ])
  AS_IF([test "x$with_doc_pdf" != "xno"],[
    eval _pdfdir=$(eval printf "%s" $pdfdir)
    eval _pdfdir2=$(eval printf "%s" $_pdfdir)],[
    _pdfdir=""
  ])
  ],[
  _docdir="Not to be installed"
])

AS_IF([test "x$with_tutorial" != "xno"],[
  eval _pdfdir=$(eval printf "%s" $pdfdir)],[
  _pdfdir=""
])


AS_IF([test "x$ac_rna_warning" != "x"],[
  ac_rna_warning="
==================================================
Warning:
$ac_rna_warning
==================================================
"
])


m4_map_args([ AC_RNA_COLOR_RESULT_PACKAGE],
            [kinfold],
            [forester],
            [cluster],
            [rnalocmin],
            [kinwalker],
            [svm],
            [gsl],
            [json],
            [perl],
            [python],
            [python3],
            [doc_pdf],
            [doc_html],
            [tutorial_pdf],
            [tutorial_html],
            [check])

m4_map_args([ AC_RNA_COLOR_RESULT_FEATURE],
            [mpfr],
            [boustrophedon],
            [NRhash],
            [c11],
            [tty_colors],
            [floatpf],
            [warn_deprecated],
            [vectorize],
            [simd],
            [lto],
            [pthreads],
            [openmp],
            [unittests],
            [check_perl],
            [check_python],
            [check_python3],
            [macosx_installer],
            [macosx_sdk])

m4_map_args([ AC_RNA_COLOR_RESULT_SIMPLE],
            [osx_arch],
            [_bindir],
            [_libdir],
            [_includedir],
            [_mandir],
            [_datadir],
            [_docdir],
            [_htmldir],
            [_pdfdir2],
            [_perl_install],
            [_perl_arch_dir],
            [_perl_lib_dir],
            [_python2_install],
            [_python2_arch_dir],
            [_python2_lib_dir],
            [_python3_install],
            [_python3_arch_dir],
            [_python3_lib_dir])

# Notify the user

AC_RNA_STRING_APPEND_FORMAT_BOLD([ac_rna_name_string], [ViennaRNA Package ${PACKAGE_VERSION}])
AC_RNA_STRING_APPEND_FORMAT_BOLD([ac_rna_final_msg], [You can run 'make', 'make check', and 'make install' now!])

AC_MSG_NOTICE([

======================================
    $ac_rna_name_string
======================================

Successfully configured with the following options:

Sub Packages
------------
  * Kinfold                   : ${result_kinfold}
  * RNAforester               : ${result_forester}
  * Analyse{Dists,Seqs}       : ${result_cluster}
  * RNAlocmin                 : ${result_rnalocmin}
  * Kinwalker                 : ${result_kinwalker}

Extra Libraries
---------------
  * Support Vector Machine    : ${result_svm}
  * GNU Scientific Library    : ${result_gsl}
  * GNU MPFR                  : ${result_mpfr}
  * JSON                      : ${result_json}

Features
--------
  * Boustrophedon             : ${result_boustrophedon}
  * Use hash for NR Sampling  : ${result_NRhash}
  * C11 features              : ${result_c11}
  * TTY colors                : ${result_tty_colors}
  * Float Precision(PF}       : ${result_floatpf}
  * Deprecation Warnings      : ${result_warn_deprecated}

Optimizations
-------------
  * Auto Vectorization        : ${result_vectorize}
  * Explicit SIMD Extension   : ${result_simd} ${simd_failed}
  * Link Time Optimization    : ${result_lto}
  * POSIX Threads             : ${result_pthreads}
  * OpenMP                    : ${result_openmp}

Scripting Language Interfaces
-----------------------------
  * Perl 5                    : ${result_perl}
  * Python 2                  : ${result_python}
  * Python 3                  : ${result_python3}

Documentation
-------------
  * Reference Manual (PDF)    : ${result_doc_pdf} ${doc_pdf_failed}
  * Reference Manual (HTML)   : ${result_doc_html} ${doc_html_failed}
  * Tutorial (PDF)            : ${result_tutorial_pdf} ${tutorial_pdf_failed}
  * Tutorial (HTML)           : ${result_tutorial_html} ${tutorial_html_failed}

Unit Tests
----------
  * Executable Programs       : ${result_unittests}
  * C-Library                 : ${result_check}
  * Perl 5 Interface          : ${result_check_perl}
  * Python 2 Interface        : ${result_check_python}
  * Python 3 Interface        : ${result_check_python3}

MacOS X
-------
  * Universal Binary          : ${result_osx_arch}
  * Installer                 : ${result_macosx_installer}
  * SDK                       : ${result_macosx_sdk}

Install Directories
-------------------
  * Executables               : $result__bindir
  * Libraries                 : $result__libdir
  * Header files              : $result__includedir
  * Extra Data                : $result__datadir
  * Man pages                 : $result__mandir
  * Documentation             : $result__docdir
      (HTML)                  : $result__htmldir
      (PDF)                   : $result__pdfdir2
  * Perl5 Interface           : $result__perl_install
      (binaries)              : $result__perl_arch_dir
      (scripts)               : $result__perl_lib_dir
  * Python2 Interface         : $result__python2_install
      (binaries)              : $result__python2_arch_dir
      (scripts)               : $result__python2_lib_dir
  * Python3 Interface         : $result__python3_install
      (binaries)              : $result__python3_arch_dir
      (scripts)               : $result__python3_lib_dir
$ac_rna_warning
$ac_rna_final_msg])
])



