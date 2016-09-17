
AC_DEFUN([RNA_GET_SWIG_INTERFACES],[
  ## collect the scripting language interfaces we gonna build
  _swig_packages=0
  AS_IF([test "x$with_perl" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [Perl 5])
    _swig_packages=1
  ])
  AS_IF([test "x$with_python" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [Python 2])
    _swig_packages=1
  ])
  AS_IF([test "x$with_python3" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [Python 3])
    _swig_packages=1
  ])
  AS_IF([test "x$with_swig" != "xyes" || test "$_swig_packages" -eq "0"],[
    AC_RNA_APPEND_VAR_COMMA($1, [None])
  ])
])

AC_DEFUN([RNA_ENABLE_SWIG_INTERFACES],[

  AX_REQUIRE_DEFINED([AX_PKG_SWIG])

  RNA_ADD_PACKAGE([swig],
                  [SWIG scripting language interfaces],
                  [yes],[],[],
                  [${srcdir}/interfaces/Makefile.am])

  AS_IF([test "x$with_swig" != "xno"],[
    wants_swig="yes"
    AX_PKG_SWIG(2.0.0, [has_swig="yes"], [has_swig="no"])
  ],[
    wants_swig="no"
  ])

  AM_CONDITIONAL(HAS_SWIG, test "x$has_swig" != "xno")
  
  RNA_ENABLE_SWIG_PERL
  RNA_ENABLE_SWIG_PYTHON
  RNA_ENABLE_SWIG_PYTHON3

])


AC_DEFUN([RNA_ENABLE_SWIG_PERL],[

  RNA_ADD_PACKAGE([perl],
                  [Perl interface],
                  [yes],[],[],
                  [${srcdir}/interfaces/Perl/Makefile.am])


  ## check for perl requirements
  AS_IF([test "x$with_perl" != "xno"],[
    AS_IF([test "x$wants_swig" = "xno"],[
      with_perl="no"
    ], [
      ## if swig is not available, check whether we already have swig generated sources
      if test "x$has_swig" != "xyes"
      then
        AC_RNA_TEST_FILE([${srcdir}/interfaces/Perl/RNA_wrap.cpp],[],[
          with_perl="no"
        ])
        AC_RNA_TEST_FILE([${srcdir}/interfaces/Perl/RNA.pm],[],[
          with_perl="no"
        ])
      fi
    ])
  ])

  RNA_PACKAGE_IF_ENABLED([perl],[
    AX_PERL_EXT
    if test "x$PERL" = "x"; then
      AC_MSG_ERROR([Perl is required to build.])
      [enable_perl_status="Perl is required to build."]
    fi
    AX_PERL_EXT_FLAGS([PERLXS_CFLAGS], [PERLXS_LDFLAGS])
    AX_PERL_EXT_LINK_CHECK([with_perl])
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
    AC_CONFIG_FILES([interfaces/Perl/Makefile interfaces/Perl/version.i])
  ])

])

AC_DEFUN([RNA_ENABLE_SWIG_PYTHON],[

  AX_REQUIRE_DEFINED([AX_PYTHON2_DEVEL])

  RNA_ADD_PACKAGE([python],
                  [Python interface],
                  [yes],[],[],
                  [${srcdir}/interfaces/Python/Makefile.am])


  ## check for python requirements
  AS_IF([test "x$with_python" != "xno"],[
    AS_IF([test "x$wants_swig" = "xno"],[
      with_python="no"
    ],[
      ## if swig is not available, check whether we already have swig generated sources
      if test "x$has_swig" != "xyes"
      then
        AC_RNA_TEST_FILE([${srcdir}/interfaces/Python/RNA_wrap.cpp],[],[
          with_python="no"
        ])
        AC_RNA_TEST_FILE([${srcdir}/interfaces/Python/RNA.py],[],[
          with_python="no"
        ])
      fi
    ])
  ])

  AS_IF([test "x$with_python" != "xno"],[

    ## check for python2 config
    AX_PYTHON2_DEVEL

    if test "x$python2_enabled_but_failed" != "x"
    then
      with_python="no"
    else
      AC_SUBST(PYTHON2DIR,$python2dir)
      AC_SUBST(PKGPYTHON2DIR,$pkgpython2dir)
      AC_SUBST(PYEXEC2DIR,$py2execdir)
      AC_SUBST(PKGPYEXEC2DIR,$pkgpy2execdir)

      AC_DEFINE([WITH_PYTHON2_INTERFACE], [1], [Create the python2 interface to RNAlib])
      AC_SUBST([PYTHON2_INTERFACE], [Python])
      AC_CONFIG_FILES([interfaces/Python/Makefile interfaces/Python/version.i])
    fi
  ])
])

AC_DEFUN([RNA_ENABLE_SWIG_PYTHON3],[

  AX_REQUIRE_DEFINED([AX_PYTHON3_DEVEL])

  RNA_ADD_PACKAGE([python3],
                  [Python3 interface],
                  [yes],[],[],
                  [${srcdir}/interfaces/Python3/Makefile.am])


  ## check for python requirements
  AS_IF([test "x$with_python3" != "xno"],[
    AS_IF([test "x$wants_swig" = "xno"],[
      with_python3="no"
    ],[
      ## if swig is not available, check whether we already have swig generated sources
      if test "x$has_swig" != "xyes"
      then
        AC_RNA_TEST_FILE([${srcdir}/interfaces/Python3/RNA_wrap.cpp],[],[
          with_python3="no"
        ])
        AC_RNA_TEST_FILE([${srcdir}/interfaces/Python3/RNA.py],[],[
          with_python3="no"
        ])
      fi
    ])
  ])

  AS_IF([test "x$with_python3" != "xno"],[

    ## check for python3 config
    AX_PYTHON3_DEVEL

    if test "x$python3_enabled_but_failed" != "x"
    then
      with_python3="no"
    else
      AC_DEFINE([WITH_PYTHON3_INTERFACE], [1], [Create the Python3 interface to RNAlib])
      AC_SUBST([PYTHON3_INTERFACE], [Python3])
    fi

    AC_CONFIG_FILES([interfaces/Python3/Makefile interfaces/Python3/version.i])
  ])
])
