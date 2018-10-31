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
RNA_ENABLE_FLOATPF
RNA_ENABLE_DEPRECATION_WARNINGS
RNA_ENABLE_COLORED_TTY
RNA_ENABLE_STATIC_BIN
RNA_ENABLE_SSE
RNA_ENABLE_VECTORIZE

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
])

AS_IF([test "x$with_doc" != "xno"],[
  eval _docdir=$(eval printf "%s" $docdir)
  AS_IF([test "x$with_doc_html" != "xno"],[
    eval _htmldir=$(eval printf "%s" $htmldir)],[
    _htmldir=""
    ])
  AS_IF([test "x$with_doc_pdf" != "xno"],[
    eval _pdfdir=$(eval printf "%s" $pdfdir)],[
    _pdfdir=""
  ])
  ],[
  _docdir="Not to be installed"
])

AS_IF([test "x$with_tutorial" != "xno"],[
  eval _pdfdir=$(eval printf "%s" $pdfdir)],[
  _pdfdir=""
])

# Notify the user

AC_MSG_NOTICE([

======================================
    ViennaRNA Package ${PACKAGE_VERSION}
======================================

Sub Packages
------------
  * Kinfold                   : ${with_kinfold:-no}
  * RNAforester               : ${with_forester:-no}
  * Analyse{Dists,Seqs}       : ${with_cluster:-no}
  * RNAlocmin                 : ${with_rnalocmin:-no}
  * Kinwalker                 : ${with_kinwalker:-no}

Extra Libraries
---------------
  * Support Vector Machine    : ${with_svm:-no}
  * GNU Scientific Library    : ${with_gsl:-no}
  * JSON                      : ${with_json:-no}

Features
--------
  * Boustrophedon             : ${enable_boustrophedon:-no}
  * C11 features              : ${enable_c11:-no}
  * TTY colors                : ${enable_tty_colors:-no}
  * Float Precision(PF}       : ${enable_floatpf:-no}
  * Deprecation Warnings      : ${enable_warn_deprecated:-no}

Optimizations
-------------
  * Auto Vectorization        : ${enable_vectorize:-no}
  * Explicit SIMD Extension   : ${enable_sse:-no} ${simd_failed}
  * Link Time Optimization    : ${enable_lto:-no}
  * POSIX Threads             : ${enable_pthreads:-no}
  * OpenMP                    : ${enable_openmp:-no}

Scripting Language Interfaces
-----------------------------
  * Perl 5                    : ${with_perl:-no}
  * Python 2                  : ${with_python:-no}
  * Python 3                  : ${with_python3:-no}

Documentation
-------------
  * Reference Manual (PDF)    : ${with_doc_pdf:-no}
  * Reference Manual (HTML)   : ${with_doc_html:-no}
  * Tutorial (PDF)            : ${with_tutorial_pdf:-no}
  * Tutorial (HTML)           : ${with_tutorial_html:-no}

Unit Tests
----------
  * Executable Programs       : ${enable_unittests:-no}
  * C-Library                 : ${with_check:-no}
  * Perl 5 Interface          : ${enable_check_perl:-no}
  * Python 2 Interface        : ${enable_check_python:-no}
  * Python 3 Interface        : ${enable_check_python3:-no}

MacOS X
-------
  * Universal Binary          : ${osx_arch:-no}
  * Installer                 : ${enable_macosx_installer:-no}
  * SDK                       : ${enable_macosx_sdk:-no}

Install Directories
-------------------
  * Executables               : $_bindir
  * Libraries                 : $_libdir
  * Header files              : $_includedir
  * Extra Data                : $_datadir
  * Man pages                 : $_mandir
  * Documentation             : $_docdir
      (HTML)                  : $(eval printf "%s" $_htmldir)
      (PDF)                   : $(eval printf "%s" $_pdfdir)
  * Perl5 Interface           : $_perl_install
      (binaries)              : $_perl_arch_dir
      (scripts)               : $_perl_lib_dir
  * Python2 Interface         : $_python2_install
      (binaries)              : $_python2_arch_dir
      (scripts)               : $_python2_lib_dir
  * Python3 Interface         : $_python3_install
      (binaries)              : $_python3_arch_dir
      (scripts)               : $_python3_lib_dir

You can run 'make', 'make check' and 'make install' now!
])
])



