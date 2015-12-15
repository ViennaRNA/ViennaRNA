#
# Unit tests with libcheck
#

AC_DEFUN([RNA_ENABLE_UNIT_TESTS],[

  RNA_ADD_PACKAGE([check],
                  [Unit tests],
                  [yes],
                  [with_check=no],
                  [with_check=yes])

  # check prerequisties for unit testing support
  RNA_PACKAGE_IF_ENABLED([check],[
  PKG_CHECK_MODULES([CHECK], [check >= 0.9.4], [],[
    AC_MSG_WARN([check not found -- will not build Unit tests])
    with_check="no"
    enabled_but_failed_check="(check framework not found)"
  ])
  ])

  RNA_PACKAGE_IF_ENABLED([check],[
    AC_DEFINE([WITH_CHECK], [1], [Include Unit tests])
    AC_SUBST([CHECK_DIR], [tests])
    AC_CONFIG_FILES([tests/Makefile])
  ])

  AM_CONDITIONAL(WITH_CHECK, test "$with_check" != "yes")
])

