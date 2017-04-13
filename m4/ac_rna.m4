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
RNA_ENABLE_BOUSTROPHEDON
RNA_ENABLE_FLOATPF
RNA_ENABLE_DEPRECATION_WARNINGS
RNA_ENABLE_COLORED_TTY
RNA_ENABLE_STATIC_BIN
RNA_ENABLE_SSE

## Set post conditions for feature
## settings

RNA_FEATURE_POST

##--------------------##
## Enable Subpackages ##
##--------------------##

RNA_ENABLE_PKG_KINFOLD
RNA_ENABLE_PKG_FORESTER
RNA_ENABLE_PKG_CLUSTER
RNA_ENABLE_PKG_KINWALKER

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
## Enable Unit tests  ##
##--------------------##
RNA_ENABLE_UNIT_TESTS


##------------------##
## Prepare files    ##
##------------------##

AC_CONFIG_FILES([src/ViennaRNA/vrna_config.h])
AC_CONFIG_FILES([misc/Makefile])
AC_CONFIG_FILES([interfaces/Makefile])
AC_CONFIG_FILES([Makefile RNAlib2.pc src/Utils/Makefile src/bin/Makefile src/Makefile man/Makefile src/ViennaRNA/Makefile doc/Makefile])
AC_CONFIG_FILES([man/cmdlopt.sh],[chmod +x man/cmdlopt.sh])
AC_CONFIG_FILES([packaging/viennarna.spec packaging/PKGBUILD])
AC_CONFIG_FILES([packaging/win_installer_archlinux_i686.nsi packaging/win_installer_archlinux_x86_64.nsi packaging/win_installer_fedora_i686.nsi packaging/win_installer_fedora_x86_64.nsi])

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
    _python2_install="Not to be installed"
])
AS_IF([test $with_python3 = "yes"],[
  eval _python3_arch_dir=$(eval printf "%s" ${py3execdir})
  eval _python3_lib_dir=$(eval printf "%s" ${python3dir})
  ],[
    _python3_arch_dir=""
    _python3_lib_dir=""
    _python3_install="Not to be installed"
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

## collect enabled swig interfaces
RNA_GET_SWIG_INTERFACES(_swig_languages_enabled)

## collect enabled unit tests
RNA_GET_UNIT_TESTS(_unit_tests_enabled)

## collect enabled subpackages
RNA_GET_SUBPACKAGES(_packages_enabled)

## collect additional options
RNA_GET_FEATURE(_features_enabled)

## collect doxygen reference manual settings
RNA_GET_DOXYGEN_REFMAN(_refman_enabled)

## collect MacOSX config (if any)
RNA_GET_MACOSX_CONFIG(_macosx_enabled)
AS_IF([test "x$_macosx_enabled" != "x"], [
  AC_RNA_APPEND_VAR_COMMA(_features_enabled, [$_macosx_enabled])
])

# Notify the user

AC_MSG_NOTICE([

################################################
#         ViennaRNA Package ${PACKAGE_VERSION}             ##
#                                             ##
# configured successfully with the following  ##
# options:                                    ##
################################################


  * Extra Programs:

    ( ${_packages_enabled} )

  * Other Options:

    ( ${_features_enabled} )

  * RNAlib Scripting Language Interfaces:

    ( ${_swig_languages_enabled} )

  * RNAlib Documentation:

    ( ${_refman_enabled} )

  * Unit Tests will be performed for:

    ( ${_unit_tests_enabled} )


##############################################
# Files will be installed in the following  ##
# directories:                              ##
##############################################

  Executables:        $_bindir
  Libraries:          $_libdir
  Header files:       $_includedir
  Extra Data:         $_datadir
  Man pages:          $_mandir
  Documentation:      $_docdir
    (HTML):           $(eval printf "%s" $_htmldir)
    (PDF):            $(eval printf "%s" $_pdfdir)
  Perl5 Interface:    $_perl_install
    (binaries):       $_perl_arch_dir
    (scripts):        $_perl_lib_dir
  Python2 Interface:  $_python2_install
    (binaries):       $_python2_arch_dir
    (scripts):        $_python2_lib_dir
  Python3 Interface:  $_python3_install
    (binaries):       $_python3_arch_dir
    (scripts):        $_python3_lib_dir
])

])



