
AC_DEFUN([RNA_ENABLE_OSX],[

  AX_REQUIRE_DEFINED([RNA_ADD_FEATURE])

  RNA_ADD_FEATURE([universal-binary],
                  [generate universal (fat) binaries on MacOSX],
                  [no],
                  [enable_universal_binary=$enableval],
                  [enable_universal_binary=no],
                  ["-arch i386 -arch x86_64"])


  RNA_FEATURE_IF_ENABLED([universal_binary],[

    ## This option is only meant for MacOSX
    case "$host" in
        *darwin*)   ;;
        *)          AC_MSG_ERROR([
**********************************************************************
        --enable-universal-binary is intended for use on MacOSX only!
**********************************************************************
                    ])
                    ;;
    esac

    osx_arch_array=( ${enable_universal_binary} )
    if test ${#osx_arch_array[@]} -gt 1 ; then
      if test "$enable_dependency_tracking" != no ; then
        AC_MSG_WARN([disabling dependency tracking])
        AM_CONDITIONAL([AMDEP],[false]) # automake-1.9.6
        AM_CONDITIONAL([am__fastdepCC],[false])
        AMDEPBACKSLASH=
      fi
    fi

    ## set architectures, if specified by the user
    if test "x$enable_universal_binary" != "xyes"; then
      osx_arch="$enable_universal_binary"
    else
      ## otherwise, use default options
      osx_arch="-arch i386 -arch x86_64"
    fi

    CFLAGS="${CFLAGS} $osx_arch"
    CXXFLAGS="${CXXFLAGS} $osx_arch"
    LDFLAGS="${LDFLAGS} $osx_arch"

    ## trim leading/trailing empty spaces from flags
    CFLAGS="`echo ${CFLAGS} | sed -e 's/^ *//' -e 's/ *$//'`"
    CXXFLAGS="`echo ${CFLAGS} | sed -e 's/^ *//' -e 's/ *$//'`"
    LDLAGS="`echo ${CFLAGS} | sed -e 's/^ *//' -e 's/ *$//'`"

  ])


])
