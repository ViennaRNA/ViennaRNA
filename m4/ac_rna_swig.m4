
AC_DEFUN([RNA_ENABLE_SWIG_INTERFACES],[
  RNA_ADD_PACKAGE([perl],
                  [Perl interface],
                  [yes],
                  [with_perl=no],
                  [with_perl=yes],
                  [${srcdir}/interfaces/Perl/Makefile.am])

  RNA_ADD_PACKAGE([python],
                  [Python interface],
                  [yes],
                  [with_python=no],
                  [with_python=yes],
                  [${srcdir}/interfaces/Python/Makefile.am])

##RNA_ADD_PACKAGE( [ruby],
##                    [Ruby interface],
##                    [no],
##                    [with_ruby=yes],
##                    [with_ruby=no],
##                    [${srcdir}/interfaces/Ruby/Makefile.am])

  # check prerequisites for Perl interface
  AC_ARG_VAR([PerlCmd], [Perl executable])
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
  RNA_PACKAGE_IF_ENABLED([perl],[
    AC_DEFINE([WITH_PERL_INTERFACE], [1], [Create the perl interface to RNAlib])
    AC_SUBST([PERL_INTERFACE], [Perl])
    AC_CONFIG_FILES([interfaces/Perl/Makefile interfaces/Perl/Makefile.PL])
  ])


  # check prerequisites for Python interface
  RNA_PACKAGE_IF_ENABLED([python],[
  AC_ARG_VAR([PythonCmd], [Python executable])
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
  RNA_PACKAGE_IF_ENABLED([python],[
    AC_DEFINE([WITH_PYTHON_INTERFACE], [1], [Create the python interface to RNAlib])
    AC_SUBST([PYTHON_INTERFACE], [Python])
    AC_CONFIG_FILES([interfaces/Python/Makefile interfaces/Python/setup.py])
  ])

])

