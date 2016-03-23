
#
# Unit tests with libcheck
#

AC_DEFUN([RNA_ENABLE_UNIT_TESTS],[

  RNA_ADD_FEATURE([unittests],
                  [Unit tests],
                  [yes],
                  [enable_unittests=no],
                  [enable_unittests=yes])

  RNA_ADD_PACKAGE([check],
                  [C-library Unit tests],
                  [yes],
                  [with_check=no],
                  [with_check=yes])

  RNA_ADD_FEATURE([check-perl],
                  [Perl interface Unit tests],
                  [yes],
                  [enable_check_perl=no],
                  [enable_check_perl=yes])

  RNA_ADD_FEATURE([check-python],
                  [Python interface Unit tests],
                  [yes],
                  [enable_check_python=no],
                  [enable_check_python=yes])

  RNA_ADD_FEATURE([check-python3],
                  [Python3 interface Unit tests],
                  [yes],
                  [enable_check_python3=no],
                  [enable_check_python3=yes])

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
      AC_SUBST([CHECK_DIR], [tests])
      AC_CONFIG_FILES([tests/Makefile])
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
  AM_CONDITIONAL(WITH_PERL_TESTS, test "x$enable_check_perl" != "xno")
  AM_CONDITIONAL(WITH_PYTHON_TESTS, test "x$enable_check_python" != "xno")
  AM_CONDITIONAL(WITH_PYTHON3_TESTS, test "x$enable_check_python3" != "xno")
])

