
AC_DEFUN([AX_PYTHON3_DEVEL],[

    # (AM_PATH_PYTHON) cannot be used for multiple Python version at once
    if test -z "$PYTHON3" ; then
      AC_PATH_PROGS([PYTHON3], [python3 python35 python3.5 python34 python3.4 python33 python3.3], [no])
    fi
    AC_ARG_VAR(PYTHON3, [Path to Python3 interpreter (e.g.: /usr/bin/python3)])

    if test "${PYTHON3}" != "no" ; then
      if test -z "${PYTHON3_CONFIG}" ; then
        AC_PATH_PROGS([PYTHON3_CONFIG], [python3-config python35-config python3.5-config python34-config python3.4-config python33-config python3.3-config], [no])
      fi
      AC_ARG_VAR(PYTHON3_CONFIG, [Path to Python3 config tool (e.g.: /usr/bin/python3-config)])

      if test "$PYTHON3_CONFIG" = "no"
      then
        AC_MSG_WARN([
The python3-config program was not found in the search path. Please ensure
that it is installed and its directory is included in the search path.
])
        python3_enabled_but_failed="python3-config missing"
      else

        AC_MSG_CHECKING([for Python3 include path])
        PYTHON3_INCLUDES=`$PYTHON3_CONFIG --includes 2> /dev/null`
        AC_SUBST(PYTHON3_INCLUDES)
        AC_MSG_RESULT([$PYTHON3_INCLUDES])

        AC_MSG_CHECKING([for Python3 compile flags])
        PYTHON3_CFLAGS=`$PYTHON3_CONFIG --cflags 2> /dev/null`
        AC_SUBST(PYTHON3_CFLAGS)
        AC_MSG_RESULT([$PYTHON3_CFLAGS])

        AC_MSG_CHECKING([for Python3 libs])
        PYTHON3_LIBS=`$PYTHON3_CONFIG --libs 2> /dev/null`
        AC_SUBST(PYTHON3_LIBS)
        AC_MSG_RESULT([$PYTHON3_LIBS])

        AC_MSG_CHECKING([for Python3 ld flags])
        PYTHON3_LDFLAGS=`$PYTHON3_CONFIG --ldflags 2> /dev/null`
        AC_SUBST(PYTHON3_LDFLAGS)
        AC_MSG_RESULT([$PYTHON3_LDFLAGS])

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
    fi
  else
    python3_enabled_but_failed="python3 executable missing"
  fi
])
