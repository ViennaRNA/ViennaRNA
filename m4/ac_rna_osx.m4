AC_DEFUN([RNA_GET_MACOSX_CONFIG],[
  ## collect the scripting language interfaces we gonna build
  _osx_config=0
  AS_IF([test "x$enable_universal_binary" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [
    MacOSX Universal binary{$osx_arch}])
    _osx_config=1
  ])
  AS_IF([test "x$enable_macosx_installer" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [MacOSX Installer])
    _osx_config=1
  ])
  AS_IF([test "x$enable_macosx_sdk" = "xyes"],[
    AC_RNA_APPEND_VAR_COMMA($1, [MacOSX SDK])
    _osx_config=1
  ])
])

AC_DEFUN([RNA_ENABLE_OSX],[

  AX_REQUIRE_DEFINED([RNA_ADD_FEATURE])

  RNA_ENABLE_OSX_UNIVERSAL_BINARY
  RNA_ENABLE_OSX_SDK
  RNA_ENABLE_OSX_INSTALLER
])

AC_DEFUN([RNA_ENABLE_OSX_INSTALLER],[

  RNA_ADD_FEATURE([macosx_installer],
                  [generate MacOSX Installer Disk Image],
                  [no])

  AC_CONFIG_FILES([packaging/macosx/Makefile packaging/macosx/Distribution.xml packaging/macosx/resources/welcome.txt])

  RNA_FEATURE_IF_ENABLED([macosx_installer],[

    ## This option is only meant for MacOSX
    case "$host" in
        *darwin*)   AC_PATH_PROG([PKGBUILD], [pkgbuild],[no])
                    AC_PATH_PROG([PRODUCTBUILD], [productbuild],[no])
                    AC_PATH_PROG([HDIUTIL], [hdiutil],[no])

                    AS_IF([test "$PKGBUILD" == "no" || test "$PRODUCTBUILD" == "no" || test "$HDITOOL" == "no"],[
                        AC_MSG_ERROR([
**********************************************************************
        --enable-macosx-installer requires the programs pkgbuild,
        productbuild, and hdituil! Make sure they are in your PATH!
**********************************************************************
                        ])
                    ],[
                      AC_SUBST([MACOSX_INSTALLER], [packaging/macosx])
                    ])
                    ;;
        *)          AC_MSG_ERROR([
**********************************************************************
        --enable-macosx-installer is intended for use on MacOSX only!
**********************************************************************
                    ])
                    ;;
    esac
  ])

  AM_CONDITIONAL(WITH_MACOSX_INSTALLER, test "$enable_macosx_installer" != "no")
])

AC_DEFUN([RNA_ENABLE_OSX_SDK],[

  RNA_ADD_FEATURE([macosx_sdk],
                  [use a specific Mac OS X SDK],
                  [no],[],[],
                  ["latest"])

  RNA_ADD_FEATURE([macosx_sdk_path],
                  [specify the path to a specific Mac OS X SDK],
                  [no],
                  [enable_macosx_sdk="custom"],
                  [enable_macosx_sdk_path=auto],
                  ["auto"])

  RNA_FEATURE_IF_ENABLED([macosx_sdk],[

    ## This option is only meant for MacOSX
    case "$host" in
        *darwin*)   ;;
        *)          AC_MSG_ERROR([
**********************************************************************
        --enable-macosx-sdk is intended for use on MacOSX only!
**********************************************************************
                    ])
                    ;;
    esac

    macosx_sdk_path=
    macosx_sdk_version=
    
    # check whether a specific path was given by the user
    if test "x$enable_macosx_sdk_path" != "xauto" ; then
      if test -d "$enable_macosx_sdk_path" ; then
        macosx_sdk_path="$enable_macosx_sdk_path"
        # determine the version
        macosx_sdk_version=10.11
      else
        AC_MSG_ERROR([
**********************************************************************
  Unable to set the MacOSX SDK path!

  The path you've specified is not a directory!
**********************************************************************
        ])
      fi
    else
      # determine xcode path automatically
      if test -x /usr/bin/xcode-select ; then
        xcodepath="`xcode-select -p`"
      else
        AC_MSG_ERROR([
**********************************************************************
  /usr/bin/xcode-select is missing!
  
  Unable to automatically determine SDK path! Please specify the
  SDK path using the --macosx-sdk-path parameter, or fix your XCode
  installation
**********************************************************************
        ])
      fi

      # automatic determination of latest SDK
      if test "x$enable_macosx_sdk" = "xyes" || test "x$enable_macosx_sdk" = "xlatest" ; then
        ## check for possible SDKs in descending order
        macosx_sdk_versions="10.11 10.10 10.9 10.8 10.7 10.6 10.5"
      else
        macosx_sdk_versions="${enable_macosx_sdk}"
      fi

      for v in ${macosx_sdk_versions} ; do
        AC_MSG_CHECKING([for MacOSX SDK $v])
        if test -d "$xcodepath/Platforms/MacOSX.platform/Developer/SDKs/MacOSX$v.sdk" ; then
          macosx_sdk_path="$xcodepath/Platforms/MacOSX.platform/Developer/SDKs/MacOSX${v}.sdk"
          macosx_sdk_version="$v"
          AC_MSG_RESULT([${macosx_sdk_path}])
          break
        else
          AC_MSG_RESULT([not found])
        fi
      done
    fi

    if test -z "$macosx_sdk_path" ; then
      AC_MSG_ERROR([
**********************************************************************
  Unable to set/determine the SDK path!

  Either you specified a non-existing version, or, if you have not
  specified a version at all, your XCode installation is missing a
  version we are aware of!
**********************************************************************
                  ])
    else
      bitness=-m64

      if test "$with_macosx_version_min_required" = ""; then
          case $macosx_sdk_version in
          10.5)
              with_macosx_version_min_required="10.5";;
          *)
              with_macosx_version_min_required="10.6";;
          esac
      fi

      if test "$with_macosx_version_max_allowed" = ""; then
          with_macosx_version_max_allowed="$macosx_sdk_version"
      fi

      ## if everything seems to be correct, lets start to set configurations
      AC_MSG_WARN([setting up build for MacOSX SDK ${macosx_sdk_version}])
      AC_MSG_CHECKING([what compiler to use])
      case $macosx_sdk_version in
        10.5)
            CC="${gccprefix}gcc-4.2"
            CXX="${gccprefix}g++-4.2"
            INSTALL_NAME_TOOL=`xcrun -find install_name_tool`
            ;;
        10.6)
            # did someone copy her 10.6 sdk into xcode 4 (needed on Mountain Lion)?
            if test "$(echo $macosx_sdk_path | cut -c1-23)" = "/Applications/Xcode.app"; then
                CC="`xcrun -find gcc`"
                CXX="`xcrun -find g++`"
            else
                CC="gcc-4.2"
                CXX="g++-4.2"
            fi
            CFLAGS="${CFLAGS} $bitness"
            CXXFLAGS="${CXXFLAGS} $bitness"
            LDFLAGS="${LDFLAGS} $bitness"
            INSTALL_NAME_TOOL=`xcrun -find install_name_tool`
            #LIBTOOL=libtool
            ;;
        10.7|10.8|10.9|10.10|10.11)
            if test "$with_macosx_version_min_required" != 10.6; then
              # Use libc++ instead of libstdc++ when possible
              stdlib=-stdlib=libc++
            fi

            CC="`xcrun -find clang`"
            CXX="`xcrun -find clang++`"

            CFLAGS="${CFLAGS} $bitness"
            CXXFLAGS="${CXXFLAGS} $bitness $stdlib"
            LDFLAGS="${LDFLAGS} $bitness"

            INSTALL_NAME_TOOL=`xcrun -find install_name_tool`
            AR=`xcrun -find ar`
            NM=`xcrun -find nm`
            STRIP=`xcrun -find strip`
            RANLIB=`xcrun -find ranlib`
            ;;
      esac
      CFLAGS="${CFLAGS} -mmacosx-version-min=$with_macosx_version_min_required -isysroot $macosx_sdk_path"
      CXXFLAGS="${CXXFLAGS} -mmacosx-version-min=$with_macosx_version_min_required -isysroot $macosx_sdk_path"
      LDFLAGS="${LDFLAGS} -mmacosx-version-min=$with_macosx_version_min_required -isysroot $macosx_sdk_path"

      AC_MSG_RESULT([$CC and $CXX])
    fi

  ])
])

AC_DEFUN([RNA_ENABLE_OSX_UNIVERSAL_BINARY],[

  RNA_ADD_FEATURE([universal_binary],
                  [generate universal (fat) binaries on MacOSX],
                  [no],[],[],
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
