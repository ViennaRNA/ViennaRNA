
AC_DEFUN([RNA_GET_UNIT_TESTS],[
  ## collect the unit tests we gonna perform
  _unittests_active=0
  AS_IF([test "x$with_check" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [C-Library])
    _unittests_active=1
  ])
  AS_IF([test "x$enable_check_perl" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [Perl 5])
    _unittests_active=1
  ])
  AS_IF([test "x$enable_check_python" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [Python 2])
    _unittests_active=1
  ])
  AS_IF([test "x$enable_check_python3" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [Python 3])
    _unittests_active=1
  ])
  AS_IF([test "x$enable_unittests" = "xno" || test "$_unittests_active" -eq "0"],[
    AC_RNA_APPEND_VAR_COMMA($1, [None])
  ])
])

#
# Unit tests with libcheck
#

AC_DEFUN([RNA_ENABLE_UNIT_TESTS],[

  RNA_ADD_FEATURE([unittests],
                  [Unit tests],
                  [yes])

  RNA_ADD_PACKAGE([check],
                  [C-library Unit tests],
                  [yes])

  RNA_ADD_FEATURE([check_executables],
                  [Unit tests for executables],
                  [yes])

  RNA_ADD_FEATURE([check_perl],
                  [Perl interface Unit tests],
                  [yes])

  RNA_ADD_FEATURE([check_python],
                  [Python interface Unit tests],
                  [yes])

  RNA_ADD_FEATURE([check_python3],
                  [Python3 interface Unit tests],
                  [yes])

  RNA_FEATURE_IF_ENABLED([unittests],[

    # check prerequisties for unit testing support
    RNA_PACKAGE_IF_ENABLED([check],[
      PKG_CHECK_MODULES([CHECK], [check >= 0.9.4], [],[
        AC_MSG_WARN([check not found -- will not build C-library Unit tests])
        with_check="no"
        enabled_but_failed_check="(check framework not found)"
      ])
    ])

    RNA_PACKAGE_IF_ENABLED([check],[
      AC_DEFINE([WITH_CHECK], [1], [Include C-library Unit tests])
    ])

    RNA_FEATURE_IF_ENABLED([check_executables],[
      AC_ARG_VAR([DIFF],[the 'diff' program to use for test output comparison])
      AC_PATH_PROG([DIFF],[diff],[])
      if test "x$DIFF" = "x"; then
        AC_MSG_WARN([diff not found -- deactivating check for executables!])
        AC_MSG_NOTICE([==> Set DIFF environment variable if present in non-standard path!])
        enable_check_executables="no"
      fi
    ])

    RNA_FEATURE_IF_ENABLED([check_perl],[
      ## switch off if perl interface is not going to be build
      if test "x$with_perl" != "xyes" ; then
        enable_check_perl="no"
      fi
    ])

    RNA_FEATURE_IF_ENABLED([check_python],[
      ## switch off if python interface is not going to be build
      if test "x$with_python" != "xyes" ; then
        enable_check_python="no"
      fi
    ])

    RNA_FEATURE_IF_ENABLED([check_python3],[
      ## switch off if python3 interface is not going to be build
      if test "x$with_python3" != "xyes" ; then
        enable_check_python3="no"
      fi
    ])

  ])

  AS_IF([test "x$enable_unittests" = "xno"],[
    with_check="no"
    enable_check_perl="no"
    enable_check_python="no"
    enable_check_python3="no"
  ])

  AM_CONDITIONAL(WITH_UNIT_TESTS, test "x$enable_unittests" != "xno")
  AM_CONDITIONAL(WITH_CHECK, test "x$with_check" != "xno")
  AM_CONDITIONAL(WITH_EXECUTABLE_TESTS, test "x$with_check_executables" != "xno")
  AM_CONDITIONAL(WITH_PERL_TESTS, test "x$enable_check_perl" != "xno")
  AM_CONDITIONAL(WITH_PYTHON_TESTS, test "x$enable_check_python" != "xno")
  AM_CONDITIONAL(WITH_PYTHON3_TESTS, test "x$enable_check_python3" != "xno")

  AC_CONFIG_FILES([tests/Makefile tests/RNApath.py tests/RNApath.pm tests/test-env.sh])
])

