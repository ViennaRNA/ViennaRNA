
AC_DEFUN([AX_PYTHON3_DEVEL],[

    # (AM_PATH_PYTHON) cannot be used for multiple Python version at once
    if test -z "$PYTHON3" ; then
      AC_PATH_PROGS([PYTHON3], [python3 python39 python3.9 python38 python3.8 python37 python3.7 python36 python3.6 python35 python3.5 python34 python3.4], [no])
    fi
    AC_ARG_VAR(PYTHON3, [Path to Python3 interpreter (e.g.: /usr/bin/python3)])

    if test "${PYTHON3}" != "no" ; then
      if $PYTHON3 -c 'import distutils.sysconfig' 2> /dev/null
      then
        AC_MSG_CHECKING([for Python3 include path])
        PYTHON3_INC=`$PYTHON3 -c "import distutils.sysconfig; \
            print(distutils.sysconfig.get_python_inc())"`
        AC_SUBST(PYTHON3_INC)
        AC_MSG_RESULT([$PYTHON3_INC])

        # check if the Python.h header is available in the include directory
        AC_MSG_CHECKING([for $PYTHON3_INC/Python.h])
        if test -f "$PYTHON3_INC/Python.h"
        then
          AC_MSG_RESULT(yes)

          AC_MSG_CHECKING([for Python3 ldflags])
          case "$host" in
            # Handle OSX Python extensions differently
            # see: http://blog.tim-smith.us/2015/09/python-extension-modules-os-x/
            #
            *-darwin* | *-macos10*)
              PYTHON3_LDFLAGS="-bundle -undefined dynamic_lookup"
              ;;
            *)
              ;;
          esac
          AC_SUBST(PYTHON3_LDFLAGS)
          AC_MSG_RESULT([$PYTHON3_LDFLAGS])

          AC_MSG_CHECKING([for Python3 extension linker])
          PYTHON3_LD=`$PYTHON3 -c "import distutils.sysconfig; \
                  print(distutils.sysconfig.get_config_var('LDSHARED'))"`
          AC_SUBST(PYTHON3_LD)
          AC_MSG_RESULT([$PYTHON3_LD])

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
          PYTHON3_CACHE_TAG=`$PYTHON3 -c "import sys; \
                  print(sys.implementation.cache_tag)"`
          AC_SUBST(PYTHON3_CACHE_TAG)
          AC_MSG_RESULT([$PYTHON3_CACHE_TAG])

          AC_MSG_CHECKING([for Python3 extension of cached and optimized bytecode])
          PYTHON3_CACHE_OPT1_EXT=`$PYTHON3 -c "import importlib.util; import sys; \
                  print('%s.pyo'%sys.implementation.cache_tag) if sys.version_info.minor<5 \
                  else print('{1}{2}'.format(*importlib.util.cache_from_source('',optimization=1).rpartition(sys.implementation.cache_tag)))"`
          AC_SUBST(PYTHON3_CACHE_OPT1_EXT)
          AC_MSG_RESULT([$PYTHON3_CACHE_OPT1_EXT])

        else
          AC_MSG_RESULT(no)
          AC_MSG_WARN([
**********************************************************************
Python.h not found!
You probably need to install python-dev, python-devel, or something similar.
**********************************************************************
])
          python3_enabled_but_failed="Python.h missing"
        fi
      else
        AC_MSG_WARN([
**********************************************************************
Failed to import distutils.sysconfig!
**********************************************************************
])
          python3_enabled_but_failed="Can't import distutils.sysconfig"
      fi
    else
      python3_enabled_but_failed="python3 executable missing"
    fi
])
