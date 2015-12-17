## Add compiler and linker flag for link time optimization

AC_DEFUN([RNA_ENABLE_LTO],[

  AX_REQUIRE_DEFINED([RNA_ADD_FEATURE])
  AX_REQUIRE_DEFINED([RNA_FEATURE_IF_ENABLED])
  AX_REQUIRE_DEFINED([AX_CHECK_COMPILE_FLAG])
  AX_REQUIRE_DEFINED([AX_COMPILER_VENDOR])

  RNA_ADD_FEATURE([lto],
                  [Link time optimization],
                  [yes],
                  [enable_lto=no],
                  [enable_lto=yes])

  # Allow for user-overwrite of ar/ranlib/nm command
  AC_ARG_VAR([USER_AR],[User overwrite for the ar command])
  AC_ARG_VAR([USER_RANLIB],[User overwrite for the ranlib command])
  AC_ARG_VAR([USER_NM],[User overwrite for the nm command])


  RNA_FEATURE_IF_ENABLED([lto],[

  ac_lto_supported=no

  AX_COMPILER_VENDOR
  AC_CANONICAL_HOST

  # check whether the compiler accepts LTO option
  AX_CHECK_COMPILE_FLAG([-flto], [ac_lto_supported=yes],[],[],[])

  # do we have support from the compiler?
  if test "x$ac_lto_supported" != "xno" ; then
    ## let's proceed by distinguishing the compiler we use
    ## Here we have to hack a little. Some systems do not provide the liblto plugin for
    ## ar/ranlib/nm by default. However, gcc provides some wrappers, gcc-ar, gcc-ranlib,
    ## and gcc-nm that do so. LLVM does it likewise.
    ## Therefore, we substitute the program env vars if we detected compilation with a
    ## compiler that we know uses this scheme
    AC_MSG_WARN([Trying to re-set ar/ranlib/nm to compiler specific wrappers])
    AS_CASE([$ax_cv_c_compiler_vendor],
    [gnu],[
        AC_PATH_PROG([GCC_AR], [gcc-ar], [no])
        AC_PATH_PROG([GCC_RANLIB], [gcc-ranlib], [no])
        AC_PATH_PROG([GCC_NM], [gcc-nm], [no])
        AS_IF([test "$GCC_AR" == "no" || test "$GCC_RANLIB" == "no" || test "$GCC_NM" == "no"],[
          enable_lto="no"
        ],[
          LTO_FLAGS="-flto"
          AR="$GCC_AR"
          RANLIB="$GCC_RANLIB"
          NM="$GCC_NM"
          AX_APPEND_FLAG(["$LTO_FLAGS"], [LDFLAGS])
          LTO_LDFLAGS="${LTO_FLAGS}"
        ])
      ],
    [clang],[
        ## Here we have to distinguish at least OS X, since it
        ## does not use gold plugin as Linux does
        case "$host" in
          *darwin*)   AC_MSG_WARN([Building for OS X])
                      ;;
          *linux*)    AC_MSG_WARN([Building for Linux])
                      AC_PATH_PROG([LLVM_AR], [llvm-ar], [no])
                      AC_PATH_PROG([LLVM_RANLIB], [llvm-ranlib], [no])
                      AC_PATH_PROG([LLVM_NM], [llvm-nm], [no])
                      AS_IF([test "$LLVM_AR" == "no"|| test "$LLVM_RANLIB" == "no" || test "$LLVM_NM" == "no" ],[
                        ac_lto_supported="no"
                      ],[
                        AX_CHECK_LINK_FLAG([-fuse-ld=gold], [
                          LTO_FLAGS="-flto"
                          AR="$LLVM_AR"
                          RANLIB="$LLVM_RANLIB"
                          NM="$LLVM_NM"
                          ## switch explicitely to ld-gold
                          AX_APPEND_FLAG([-fuse-ld=gold], [LDFLAGS])
                          AX_APPEND_FLAG(["$LTO_FLAGS"], [LDFLAGS])
                          ## and again for the pkg-config file
                          LTO_LDFLAGS="-fuse-ld=gold ${LTO_FLAGS}"
                         ], [
                          ac_lto_supported="no"],[], [])
                      ])
                      ;;
          *)          AC_MSG_WARN([Unknown target host, deactivating LTO support])
                      ac_lto_supported="no"
                      ;;
        esac
      ],
    [
      ac_lto_supported="no"
    ])

    ## append -flto flag if all checks went fine
    if test "x$ac_lto_supported" != "xno" ; then
      AX_APPEND_FLAG(["$LTO_FLAGS"], [CFLAGS])
      AX_APPEND_FLAG(["$LTO_FLAGS"], [CXXFLAGS])
    else
      AC_MSG_WARN([Your compiler/linker combination does not support link-time optimization (LTO)])
      enable_lto="no"
      enabled_but_failed_lto="unsupported"
    fi
  fi
  ])

  # distribute additional compiler and linker flags
  # --> set these variables instead of CFLAGS, CXXFLAGS, or LDFLAGS
  AC_SUBST([CXXFLAGS])
  AC_SUBST([CFLAGS])
  AC_SUBST([LDFLAGS])
  AC_SUBST([LTO_FLAGS])
  AC_SUBST([LTO_LDFLAGS])
  # substitute the environment variables in case they have changed above?
  AC_SUBST([AR])
  AC_SUBST([RANLIB])
  AC_SUBST([NM])
])
