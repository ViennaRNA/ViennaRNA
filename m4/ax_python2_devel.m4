
AC_DEFUN([AX_PYTHON2_DEVEL],[

    # (AM_PATH_PYTHON) cannot be used for multiple Python version at once
    if test -z "$PYTHON2" ; then
      AC_PATH_PROGS([PYTHON2], [python2 python2.7 python2.6], [no])
    fi
    AC_ARG_VAR(PYTHON2, [Path to Python2 interpreter (e.g.: /usr/bin/python2)])

    if test "${PYTHON2}" != "no" ; then
      if $PYTHON2 -c 'import distutils.sysconfig' 2> /dev/null
      then
        AC_MSG_CHECKING([for Python2 include path])
        PYTHON2_INC=`$PYTHON2 -c "import distutils.sysconfig; \
            print(distutils.sysconfig.get_python_inc())"`
        AC_SUBST(PYTHON2_INC)
        AC_MSG_RESULT([$PYTHON2_INC])

        # check if the Python.h header is available in the include directory
        AC_MSG_CHECKING([for $PYTHON2_INC/Python.h])
        if test -f "$PYTHON2_INC/Python.h"
        then
          AC_MSG_RESULT(yes)

          AC_MSG_CHECKING([for Python2 ldflags])
          case "$host" in
            # Handle OSX Python extensions differently
            # see: http://blog.tim-smith.us/2015/09/python-extension-modules-os-x/
            #
            *-darwin* | *-macos10*)
              PYTHON2_LDFLAGS="-bundle -undefined dynamic_lookup"
              ;;
            *)
              ;;
          esac
          AC_SUBST(PYTHON2_LDFLAGS)
          AC_MSG_RESULT([$PYTHON2_LDFLAGS])

          AC_MSG_CHECKING([for Python2 extension linker])
          PYTHON2_LD=`$PYTHON2 -c "import distutils.sysconfig; \
                  print(distutils.sysconfig.get_config_var('LDSHARED'))"`
          AC_SUBST(PYTHON2_LD)
          AC_MSG_RESULT([$PYTHON2_LD])

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
        else
          AC_MSG_RESULT(no)
          AC_MSG_WARN([
**********************************************************************
Python.h not found!
You probably need to install python-dev, python-devel, or something similar.
**********************************************************************
])
          python2_enabled_but_failed="Python.h missing"
        fi
      else
        AC_MSG_WARN([
**********************************************************************
Failed to import distutils.sysconfig!
**********************************************************************
])
          python2_enabled_but_failed="Can't import distutils.sysconfig"
      fi
    else
      python2_enabled_but_failed="python2 executable missing"
    fi
])
