
AC_DEFUN([RNA_ENABLE_SWIG_INTERFACES],[

  AX_REQUIRE_DEFINED([AX_PKG_SWIG])

  RNA_ADD_PACKAGE([swig],
                  [SWIG scripting language interfaces],
                  [yes],
                  [with_swig=no],
                  [with_swig=yes],
                  [${srcdir}/interfaces/Makefile.am])

  AS_IF([test "x$with_swig" != "xno"],[
    AX_PKG_SWIG(2.0.0, [has_swig="yes"], [has_swig="no"])
  ])

  RNA_ENABLE_SWIG_PERL
  RNA_ENABLE_SWIG_PYTHON
  RNA_ENABLE_SWIG_PYTHON3

##RNA_ADD_PACKAGE( [ruby],
##                    [Ruby interface],
##                    [no],
##                    [with_ruby=yes],
##                    [with_ruby=no],
##                    [${srcdir}/interfaces/Ruby/Makefile.am])

])
AC_DEFUN([RNA_ENABLE_SWIG_PERL],[

  RNA_ADD_PACKAGE([perl],
                  [Perl interface],
                  [yes],
                  [with_perl=no],
                  [with_perl=yes],
                  [${srcdir}/interfaces/Perl/Makefile.am])


  ## check for perl requirements
  AS_IF([test "x$with_perl" != "xno"],[
    ## if swig is not available, check whether we already have swig generated sources
    if test "x$has_swig" != "xyes"
    then
      AC_RNA_TEST_FILE([${srcdir}/interfaces/Perl/RNA_wrap.c],[],[
        with_perl="no"
      ])
      AC_RNA_TEST_FILE([${srcdir}/interfaces/Perl/RNA.pm],[],[
        with_perl="no"
      ])
    fi
  ])

  RNA_PACKAGE_IF_ENABLED([perl],[
    AX_PERL_EXT
    if test "x$PERL" = "x"; then
      AC_MSG_ERROR([Perl is required to build.])
      [enable_perl_status="Perl is required to build."]
    fi
    AX_PERL_EXT_FLAGS([PERLXS_CFLAGS], [PERLXS_LDFLAGS])
    AC_SUBST([PERLXS_CFLAGS])
    AC_SUBST([PERLXS_LDFLAGS])
  ])

  # prepare all files for perl interface
  RNA_PACKAGE_IF_ENABLED([perl],[
    # Compose the correct installation path for perl modules
    #
    # here we actually have to account for INSTALLDIRS env variable, which can be
    #
    # site    = where the local systems administrator installs packages to
    # vendor  = where system packages are installed to, or
    # core    = where perl core packages are installed
    #
    # The default selection is 'site', but upon packaging for a specific distribution
    # we might want the user to set this to 'vendor'
    #
    AS_IF([ test "x$INSTALLDIRS" == "xvendor" ],[
      PERL_ARCH_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_VENDORARCH} | sed "s,${PERL_EXT_VENDORPREFIX},,"`
      PERL_LIB_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_VENDORLIB} | sed "s,${PERL_EXT_VENDORPREFIX},,"`
      ],[
      PERL_ARCH_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_SITEARCH} | sed "s,${PERL_EXT_SITEPREFIX},,"`
      PERL_LIB_RELATIVE_INSTALL_DIR=`echo ${PERL_EXT_SITELIB} | sed "s,${PERL_EXT_SITEPREFIX},,"`
    ])
    AC_SUBST(PERL_ARCH_RELATIVE_INSTALL_DIR)
    AC_SUBST(PERL_LIB_RELATIVE_INSTALL_DIR)

    AC_DEFINE([WITH_PERL_INTERFACE], [1], [Create the perl interface to RNAlib])
    AC_SUBST([PERL_INTERFACE], [Perl])
    AC_CONFIG_FILES([interfaces/Perl/Makefile])
  ])

])

AC_DEFUN([RNA_ENABLE_SWIG_PYTHON],[

  RNA_ADD_PACKAGE([python],
                  [Python interface],
                  [yes],
                  [with_python=no],
                  [with_python=yes],
                  [${srcdir}/interfaces/Python/Makefile.am])


  ## check for python requirements
  AS_IF([test "x$with_python" != "xno"],[
    ## if swig is not available, check whether we already have swig generated sources
    if test "x$has_swig" != "xyes"
    then
      AC_RNA_TEST_FILE([${srcdir}/interfaces/Python/RNA_wrap.c],[],[
        with_python="no"
      ])
      AC_RNA_TEST_FILE([${srcdir}/interfaces/Python/RNA.py],[],[
        with_python="no"
      ])
    fi
  ])

  AS_IF([test "x$with_python" != "xno"],[
    AX_PYTHON_DEVEL([< '3.0.0'])
    AM_PATH_PYTHON
    AX_SWIG_PYTHON
##    pythondir=$PYTHON_SITE_PKG
##    pyexecdir=$PYTHON_SITE_PKG_EXEC

    AC_SUBST(PYTHONDIR,$pythondir)
    AC_SUBST(PKGPYTHONDIR,$pkgpythondir)
    AC_SUBST(PYEXECDIR,$pyexecdir)
    AC_SUBST(PKGPYEXECDIR,$pkgpyexecdir)

    AC_DEFINE([WITH_PYTHON_INTERFACE], [1], [Create the python interface to RNAlib])
    AC_SUBST([PYTHON_INTERFACE], [Python])
    AC_CONFIG_FILES([interfaces/Python/Makefile])
  ])
])

AC_DEFUN([RNA_ENABLE_SWIG_PYTHON3],[

  RNA_ADD_PACKAGE([python3],
                  [Python3 interface],
                  [yes],
                  [with_python3=no],
                  [with_python3=yes],
                  [${srcdir}/interfaces/Python3/Makefile.am])


  ## check for python requirements
  AS_IF([test "x$with_python3" != "xno"],[
    ## if swig is not available, check whether we already have swig generated sources
    if test "x$has_swig" != "xyes"
    then
      AC_RNA_TEST_FILE([${srcdir}/interfaces/Python3/RNA_wrap.c],[],[
        with_python3="no"
      ])
      AC_RNA_TEST_FILE([${srcdir}/interfaces/Python3/RNA.py],[],[
        with_python3="no"
      ])
    fi
  ])

  AS_IF([test "x$with_python3" != "xno"],[
    # (AM_PATH_PYTHON) cannot be used for multiple Python version at once
    if test -z "$PYTHON3" ; then
      AC_PATH_PROG([PYTHON3], [python3])
    fi
    AC_ARG_VAR(PYTHON3, [Path to Python3 interpreter])

    AC_PATH_PROG([PYTHON3_CONFIG], [python3-config], [no])
                [if test "$PYTHON3_CONFIG" = "no"]
                [then]
                    [echo "The python3-config program was not found in the search path. Please ensure"]
                    [echo "that it is installed and its directory is included in the search path."]
                    [echo "Then run configure again before attempting to build OpenSCAP."]
                    [exit 1]
                [fi]

    AC_MSG_CHECKING([for Python3 include path])
    PYTHON3_INCLUDES=`python3-config --includes 2> /dev/null`
    AC_SUBST(PYTHON3_INCLUDES)
    AC_MSG_RESULT([$PYTHON3_INCLUDES])

    AC_MSG_CHECKING([for Python3 compile flags])
    PYTHON3_CFLAGS=`python3-config --cflags 2> /dev/null`
    AC_SUBST(PYTHON3_CFLAGS)
    AC_MSG_RESULT([$PYTHON3_CFLAGS])

    AC_MSG_CHECKING([for Python3 link flags])
    PYTHON3_LIBS=`python3-config --libs 2> /dev/null`
    AC_SUBST(PYTHON3_LIBS)
    AC_MSG_RESULT([$PYTHON3_LIBS])

    AC_MSG_CHECKING([for directory to install Python3 scripts in])
    if test -z "$PYTHON3_DIR" ; then
      # the string concatenation below is just a trick to prevent substitution
      PYTHON3_DIR=`$PYTHON3 -c "import distutils.sysconfig; \
            print(distutils.sysconfig.get_python_lib(0,0,prefix='$' '{prefix}'))"`
    fi
    AC_SUBST(python3dir, $PYTHON3_DIR)
    AC_MSG_RESULT([$PYTHON3_DIR])
    AC_ARG_VAR(PYTHON3_DIR, [Directory to install python3 scripts in])

    AC_MSG_CHECKING([for directory to install architecture dependent python3 things in])
    if test -z "$PYTHON3_EXECDIR" ; then
      PYTHON3_EXECDIR=`$PYTHON3 -c "import distutils.sysconfig; \
            print(distutils.sysconfig.get_python_lib(1,0,prefix='$' '{exec_prefix}'))"`
    fi
    AC_SUBST(py3execdir, $PYTHON3_EXECDIR)
    AC_MSG_RESULT([$PYTHON3_EXECDIR])
    AC_ARG_VAR(PYTHON3_EXECDIR, [Directory to install architecture dependent python3 things in])

    AC_MSG_CHECKING([for Python3 module extension])
    dnl Usually ".so", but for example, Mac OS X uses ".dylib".
    PYTHON3_SO=`$PYTHON3 -c "import distutils.sysconfig; \
            print(distutils.sysconfig.get_config_vars('SO')[[0]])"`
    AC_SUBST(PYTHON3_SO)
    AC_MSG_RESULT([$PYTHON3_SO])

    AC_MSG_CHECKING([for Python3 tag for cached compiled scripts])
    PYTHON3_CACHE_TAG=`$PYTHON3 -c "import imp; \
            print(imp.get_tag())"`
    AC_SUBST(PYTHON3_CACHE_TAG)
    AC_MSG_RESULT([$PYTHON3_CACHE_TAG])

    AC_MSG_CHECKING([for Python3 extension of cached and optimized bytecode])
    PYTHON3_CACHE_OPT1_EXT=`$PYTHON3 -c "import imp,sys; \
            print('%s.pyo'%imp.get_tag()) if sys.version_info.minor<5 \
            else print('{1}{2}'.format(*imp.util.cache_from_source('',optimization=1).rpartition(imp.get_tag())))"`
    AC_SUBST(PYTHON3_CACHE_OPT1_EXT)
    AC_MSG_RESULT([$PYTHON3_CACHE_OPT1_EXT])

    AC_DEFINE([WITH_PYTHON3_INTERFACE], [1], [Create the Python3 interface to RNAlib])
    AC_SUBST([PYTHON3_INTERFACE], [Python3])
    AC_CONFIG_FILES([interfaces/Python3/Makefile])
  ])
])
