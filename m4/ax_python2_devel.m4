
AC_DEFUN([AX_PYTHON2_DEVEL],[

    # (AM_PATH_PYTHON) cannot be used for multiple Python version at once
    if test -z "$PYTHON2" ; then
      AC_PATH_PROGS([PYTHON2], [python2 python2.7 python2.6], [no])
    fi
    AC_ARG_VAR(PYTHON2, [Path to Python2 interpreter (e.g.: /usr/bin/python2)])

    if test "${PYTHON2}" != "no" ; then
      if test -z "${PYTHON2_CONFIG}" ; then
        AC_PATH_PROGS([PYTHON2_CONFIG], [python2-config python2.7-config python2.6-config], [no])
      fi
      AC_ARG_VAR(PYTHON2_CONFIG, [Path to Python2 config tool (e.g.: /usr/bin/python2-config)])

      if test "$PYTHON2_CONFIG" = "no"
      then
        AC_MSG_WARN([
The python2-config program was not found in the search path. Please ensure
that it is installed and its directory is included in the search path.
])
        python2_enabled_but_failed="python2-config missing"
      else

        AC_MSG_CHECKING([for Python2 include path])
        PYTHON2_INCLUDES=`${PYTHON2_CONFIG} --includes 2> /dev/null`
        AC_SUBST(PYTHON2_INCLUDES)
        AC_MSG_RESULT([$PYTHON2_INCLUDES])

        AC_MSG_CHECKING([for Python2 compile flags])
        PYTHON2_CFLAGS=`${PYTHON2_CONFIG} --cflags 2> /dev/null`
        AC_SUBST(PYTHON2_CFLAGS)
        AC_MSG_RESULT([$PYTHON2_CFLAGS])

        AC_MSG_CHECKING([for Python2 libs])
        PYTHON2_LIBS=`${PYTHON2_CONFIG} --libs 2> /dev/null`
        AC_SUBST(PYTHON2_LIBS)
        AC_MSG_RESULT([$PYTHON2_LIBS])

        AC_MSG_CHECKING([for Python2 ldflags])
        case "$host" in
          # Handle OSX Python extensions differently
          # see: http://blog.tim-smith.us/2015/09/python-extension-modules-os-x/
          #
          *-darwin* | *-macos10*)
            PYTHON2_LDFLAGS="-bundle -undefined dynamic_lookup"
            ;;
          *)
            PYTHON2_LDFLAGS=`${PYTHON2_CONFIG} --ldflags 2> /dev/null`
            ;;
        esac
        AC_SUBST(PYTHON2_LDFLAGS)
        AC_MSG_RESULT([$PYTHON2_LDFLAGS])

        AC_MSG_CHECKING([for directory to install Python2 scripts in])
        if test -z "$PYTHON2_DIR" ; then
          # the string concatenation below is just a trick to prevent substitution
          PYTHON2_DIR=`$PYTHON2 -c "import distutils.sysconfig; \
                print(distutils.sysconfig.get_python_lib(0,0,prefix='$' '{prefix}'))"`
        fi
        AC_SUBST(python2dir, $PYTHON2_DIR)
        AC_MSG_RESULT([$PYTHON2_DIR])
        AC_ARG_VAR(PYTHON2_DIR, [Directory to install python2 scripts in])

        AC_MSG_CHECKING([for directory to install architecture dependent python2 things in])
        if test -z "$PYTHON2_EXECDIR" ; then
          PYTHON2_EXECDIR=`$PYTHON2 -c "import distutils.sysconfig; \
                print(distutils.sysconfig.get_python_lib(1,0,prefix='$' '{exec_prefix}'))"`
        fi
        AC_SUBST(py2execdir, $PYTHON2_EXECDIR)
        AC_MSG_RESULT([$PYTHON2_EXECDIR])
        AC_ARG_VAR(PYTHON2_EXECDIR, [Directory to install architecture dependent python2 things in])

        AC_MSG_CHECKING([for Python2 module extension])
        dnl Usually ".so", but for example, Mac OS X uses ".dylib".
        PYTHON2_SO=`$PYTHON2 -c "import distutils.sysconfig; \
                print(distutils.sysconfig.get_config_vars('SO')[[0]])"`
        AC_SUBST(PYTHON2_SO)
        AC_MSG_RESULT([$PYTHON2_SO])

        #
        # final check to see if everything compiles alright
        #
        AC_MSG_CHECKING([consistency of all components of python development environment])
        # save current global flags
        ac_save_LIBS="$LIBS"
        ac_save_CFLAGS="$CFLAGS"
        ac_save_CPPFLAGS="$CPPFLAGS"
        LIBS="$ac_save_LIBS $PYTHON2_LIBS"
        CFLAGS="$ac_save_CFLAGS $PYTHON2_CFLAGS"
        CPPFLAGS="$ac_save_CPPFLAGS $PYTHON2_INCLUDES"
        AC_LANG_PUSH([C])
        AC_LINK_IFELSE([
                AC_LANG_PROGRAM([[#include <Python.h>]],
                                [[Py_Initialize();]])
                ],[pythonexists=yes],[pythonexists=no])
        AC_LANG_POP([C])
        # turn back to default flags
        CPPFLAGS="$ac_save_CPPFLAGS"
        LIBS="$ac_save_LIBS"
        CFLAGS="$ac_save_CFLAGS"

        AC_MSG_RESULT([$pythonexists])

        if test ! "x$pythonexists" = "xyes"; then
          AC_MSG_FAILURE([
  Could not link test program to Python2. Maybe the main Python2 library has been
  installed in some non-standard library path. If so, pass it to configure,
  via the LIBS environment variable.
  Example: ./configure LIBS="-L/usr/non-standard-path/python/lib"
  ============================================================================
   ERROR!
   You probably have to install the development version of the Python package
   for your distribution.  The exact name of this package varies among them.
  ============================================================================
          ])
            python2_enabled_but_failed="could not link to python2 library"
        fi

    fi
  else
    python2_enabled_but_failed="python2 executable missing"
  fi
])
