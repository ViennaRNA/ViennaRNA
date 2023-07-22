#
# SVM support for Lfold -z
#

AC_DEFUN([RNA_ENABLE_SVM],[

  SVM_VERSION=3.31

  RNA_ADD_PACKAGE([svm],
                  [svm classifiers],
                  [yes],[],[],
                  [${srcdir}/src/libsvm-${SVM_VERSION}/svm.cpp ${srcdir}/src/libsvm-${SVM_VERSION}/svm.h])

  RNA_PACKAGE_IF_ENABLED([svm],[
    AC_SUBST([LIBSVM_DIR], [libsvm-${SVM_VERSION}])
    ## substitution for reference manual config
    AC_SUBST([REFDOC_PREDEF_SVM], [VRNA_WITH_SVM])
    AC_DEFINE([VRNA_WITH_SVM], [1], [Compute z-scores for RNALfold])
    SVM_LIBS="-lstdc++"
    CONFIG_SVM="#define VRNA_WITH_SVM"
    SETUPCFG_SW_SVM="True"
    SETUPCFG_SVM_MACRO="svm_macro = VRNA_WITH_SVM"
    SETUPCFG_SVM_LIBS="svm_libs  = stdc++"
  ])

  RNA_PACKAGE_IF_DISABLED([svm],[
    SETUPCFG_SW_SVM="False"
  ])

  AC_SUBST(SVM_VERSION)
  AC_SUBST(SVM_LIBS)
  AC_SUBST(CONFIG_SVM)
  AC_SUBST(SETUPCFG_SW_SVM)
  AC_SUBST(SETUPCFG_SVM_MACRO)
  AC_SUBST(SETUPCFG_SVM_LIBS)
  AM_CONDITIONAL(VRNA_AM_SWITCH_SVM, test "$with_svm" != "no")
])


#
# GSL support for RNApvmin
#

AC_DEFUN([RNA_ENABLE_GSL],[

  RNA_ADD_PACKAGE([gsl],
                  [GNU Scientific Library],
                  [yes])

  # check prerequisties for gsl support
  RNA_PACKAGE_IF_ENABLED([gsl],[
    AC_CHECK_LIB([gslcblas],[cblas_dgemm])
    AC_CHECK_LIB([gsl],[gsl_blas_dgemm])

    if test "$ac_cv_lib_gsl_gsl_blas_dgemm" != yes; then
      AC_MSG_WARN("Can't find libgsl. Falling back to default implementation.")
      enabled_but_failed_gsl="(libgsl not found)"
      with_gsl="no"
    fi
  ])

  RNA_PACKAGE_IF_ENABLED([gsl],[
    ## Check for gsl/gsl_multimin.h header first
    AC_CHECK_HEADER([gsl/gsl_multimin.h], [
      ## now, check if we can compile a program
      AC_MSG_CHECKING([whether we can compile and link programs with gsl support])
      ac_save_LIBS="$LIBS"
      LIBS="$ac_save_LIBS -lgsl -lgslcblas"

      AC_LANG_PUSH([C])

      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM(
          [[#include <stdio.h>
            #include <gsl/gsl_sf_bessel.h>
          ]],
          [[
            double x = 5.0;
            double y = gsl_sf_bessel_J0 (x);
            printf ("J0(%g) = %.18e\n", x, y);
            return 0;
          ]])
      ],[
        AC_DEFINE([VRNA_WITH_GSL], [1], [Use GNU Scientific Library])
        GSL_LIBS="-lgsl -lgslcblas"
        CONFIG_GSL="#define VRNA_WITH_GSL"
        SETUPCFG_SW_GSL="True"
        SETUPCFG_GSL_MACRO="gsl_macro = VRNA_WITH_GSL"
        SETUPCFG_GSL_LIBS="gsl_libs  = gsl gslcblas"
      ],[
        with_gsl=no
      ])
      AC_LANG_POP([C])
      LIBS="$ac_save_LIBS"
      AC_MSG_RESULT([$with_gsl])
    ], [
      AC_MSG_WARN([
==========================
Failed to find gsl/gsl_multimin.h!

You probably need to install the gsl-devel package or similar
==========================
    ])
    with_gsl=no])
  ])

  RNA_PACKAGE_IF_DISABLED([gsl],[
    SETUPCFG_SW_GSL="False"
  ])

  AC_SUBST([GSL_LIBS])
  AC_SUBST([CONFIG_GSL])
  AC_SUBST([SETUPCFG_SW_GSL])
  AC_SUBST([SETUPCFG_GSL_MACRO])
  AC_SUBST([SETUPCFG_GSL_LIBS])
  AM_CONDITIONAL(VRNA_AM_SWITCH_GSL, test "x$with_gsl" = "xyes")
])


#
# Use hash for non-redundant sampling
#

AC_DEFUN([RNA_ENABLE_NR_SAMPLE_HASH],[

  RNA_ADD_FEATURE([NRhash],
                  [Hash for non-redundant sampling datas structure],
                  [no])

  RNA_FEATURE_IF_ENABLED([NRhash],[
    AC_DEFINE([VRNA_NR_SAMPLING_HASH], [1], [Use Hash for non-redundant sampling data structure])
    CONFIG_NR_SAMPLING="#define VRNA_NR_SAMPLING_HASH"
  ])

  AC_SUBST(CONFIG_NR_SAMPLING)
])


AC_DEFUN([RNA_ENABLE_MPFR], [

  RNA_ADD_FEATURE([mpfr],
                  [Use MPFR library for aribtrary precision computations in non-redundant sampling],
                  [yes])

  RNA_FEATURE_IF_ENABLED([mpfr],[
    ## Check for mpfr.h header first
    AC_CHECK_HEADER([mpfr.h], [
      ## now, check if we can compile a program
      AC_MSG_CHECKING([whether we can compile programs with mpfr support])
      ac_save_LIBS="$LIBS"
      LIBS="$ac_save_LIBS -lmpfr -lgmp"

      AC_LANG_PUSH([C])

      AC_COMPILE_IFELSE([
        AC_LANG_PROGRAM(
          [[#include <stdio.h>
            #include <mpfr.h>
          ]],
          [[  printf ("MPFR library: %-12s\nMPFR header:  %s (based on %d.%d.%d)\n",
              mpfr_get_version (), MPFR_VERSION_STRING, MPFR_VERSION_MAJOR,
              MPFR_VERSION_MINOR, MPFR_VERSION_PATCHLEVEL);
              return 0;
          ]])
      ],[
        MPFR_LIBS="-lmpfr -lgmp"
        SETUPCFG_SW_MPFR="True"
        SETUPCFG_MPFR_MACRO="mpfr_macro  = VRNA_NR_SAMPLING_MPFR"
        SETUPCFG_MPFR_LIBS="mpfr_libs   = mpfr gmp"
        AC_DEFINE([VRNA_NR_SAMPLING_MPFR], [1], [Use MPFR for non-redundant sampling data structure operations])
      ],[
        enable_mpfr=no
      ])
      AC_LANG_POP([C])
      LIBS="$ac_save_LIBS"
      AC_MSG_RESULT([$enable_mpfr])
    ], [
      AC_MSG_WARN([
==========================
Failed to find mpfr.h!

You probably need to install the mpfr-devel package or similar
==========================
    ])
    enable_mpfr=no])
  ])

  RNA_FEATURE_IF_DISABLED([mpfr],[
    SETUPCFG_SW_MPFR="False"
  ])

  AC_SUBST(MPFR_LIBS)
  AC_SUBST(SETUPCFG_SW_MPFR)
  AC_SUBST(SETUPCFG_MPFR_MACRO)
  AC_SUBST(SETUPCFG_MPFR_LIBS)
  AM_CONDITIONAL(VRNA_AM_SWITCH_MPFR, test "x$enable_mpfr" = "xyes")
])

#
# OpenMP support
#

AC_DEFUN([RNA_ENABLE_OPENMP],[

  RNA_ADD_FEATURE([openmp],
                  [OpenMP support],
                  [yes])

  RNA_FEATURE_IF_ENABLED([openmp],[
    AC_LANG_PUSH([C])
    AX_OPENMP([],[enable_openmp="no"])
    AC_LANG_POP([C])

    AS_IF([ test "x$enable_openmp" != "xno" ],[
      OMP_CFLAGS="$OPENMP_CFLAGS"

      AC_LANG_PUSH([C++])
      AX_OPENMP([],[enable_openmp="no"])
      AC_LANG_POP([C++])

      if test "x$enable_openmp" != "xno"
      then
        RNA_CFLAGS="${RNA_CFLAGS} ${OMP_CFLAGS}"
        RNA_CXXFLAGS="${RNA_CXXFLAGS} ${OPENMP_CXXFLAGS}"
        LIBGOMPFLAG="$OPENMP_CXXFLAGS"
        CONFIG_OPENMP="#define VRNA_WITH_OPENMP"
        SETUPCFG_SW_OPENMP="True"
        SETUPCFG_OPENMP_MACRO="openmp_macro  = VRNA_WITH_OPENMP"
        SETUPCFG_OPENMP_FLAGS="openmp_cflags = $OPENMP_CXXFLAGS"
      fi
    ])
  ])

  RNA_FEATURE_IF_DISABLED([openmp],[
    SETUPCFG_SW_OPENMP="False"
  ])
  AC_SUBST(CONFIG_OPENMP)
  AC_SUBST(SETUPCFG_SW_OPENMP)
  AC_SUBST(SETUPCFG_OPENMP_MACRO)
  AC_SUBST(SETUPCFG_OPENMP_FLAGS)
  AC_SUBST(LIBGOMPFLAG)
  AC_SUBST(OPENMP_CFLAGS)
  AC_SUBST(OPENMP_CXXFLAGS)
])


AC_DEFUN([RNA_ENABLE_PTHREADS],[
  RNA_ADD_FEATURE([pthreads],
                  [Parallel input processing support],
                  [yes])

  RNA_FEATURE_IF_ENABLED([pthreads],[
    ## probe for pthreads availability
    AX_PTHREAD([
        AC_DEFINE([VRNA_WITH_PTHREADS], [1], [Use pthreads for parallel input processing])
    ], [
        enable_pthreads="no"
    ])
  ])

  AC_SUBST(PTHREAD_LIBS)
  AC_SUBST(PTHREAD_CFLAGS)
  AM_CONDITIONAL(VRNA_AM_SWITCH_PTHREADS, test "x$enable_pthreads" = "xyes")
])

#
# C11 feature support
#

AC_DEFUN([RNA_ENABLE_C11],[

  RNA_ADD_FEATURE([c11],
                  [C11 feature support (unnamed unions of unnamed structs)],
                  [yes])

  RNA_FEATURE_IF_ENABLED([c11],[
    AC_MSG_CHECKING([whether the C compiler allows unnamed unions of unnamed structs])
    # save current global flags
    AC_LANG_PUSH([C])
    AC_RUN_IFELSE([AC_LANG_SOURCE([[#include <stdlib.h>
                                    struct bla {
                                      union {
                                        struct { int a; char  b;};
                                        struct { long c; double d;};
                                      };
                                    };
                                    int main (void) { return 0;} ]])],
                                    [enable_c11=yes],
                                    [enable_c11=no],
                                    [enable_c11=no])

    AC_LANG_POP([C])
    AC_MSG_RESULT([$enable_c11])

    AS_IF([ test "x$enable_c11" != "xno" ],[
      AC_MSG_CHECKING([whether the C++ compiler allows unnamed unions of unnamed structs])
      AC_LANG_PUSH([C++])
      AC_RUN_IFELSE([AC_LANG_SOURCE([[extern "C" {
                                        #include <stdlib.h>
                                        struct bla {
                                          union {
                                            struct { int a; char  b;};
                                            struct { long c; double d;};
                                          };
                                        };
                                      };
                                      int main (void) { return 0;} ]])],
                                      [enable_c11=yes],
                                      [enable_c11=no],
                                      [enable_c11=no])

      AC_LANG_POP([C++])
      AC_MSG_RESULT([$enable_c11])
    ])
  ])

  AS_IF([ test "x$enable_c11" != "xyes" ],[
    DISABLE_C11_FEATURES=-DVRNA_DISABLE_C11_FEATURES
    CONFIG_DISABLE_C11_FEATURES="#define VRNA_DISABLE_C11_FEATURES"
    AX_APPEND_FLAG(["${DISABLE_C11_FEATURES}"], [RNA_CPPFLAGS])
  ])

  AC_SUBST(DISABLE_C11_FEATURES)
  AC_SUBST(CONFIG_DISABLE_C11_FEATURES)
])

#
#
#

AC_DEFUN([RNA_ENABLE_FLOATPF],[

  RNA_ADD_FEATURE([floatpf],
                  [Floating point precision in partition function computations],
                  [no])

  # Handle floating point precision flag
  RNA_FEATURE_IF_ENABLED([floatpf],[

    AC_DEFINE([USE_FLOAT_PF], [1], [Use floating point precision in partition function computations])
    CONFIG_FLOAT_PF="#define USE_FLOAT_PF"
    ## substitution for reference manual config
    AC_SUBST([REFDOC_PREDEF_FLOAT_PF], [USE_FLOAT_PF])
    AC_SUBST([FLOAT_PF_FLAG], [-DUSE_FLOAT_PF])
    AX_APPEND_FLAG([-DUSE_FLOAT_PF], [RNA_CPPFLAGS])
  ])

  AC_SUBST(CONFIG_FLOAT_PF)
])


#
# Warn about usage of deprecated symbols
#

AC_DEFUN([RNA_ENABLE_DEPRECATION_WARNINGS],[

  RNA_ADD_FEATURE([warn_deprecated],
                  [Warn upon usage of deprecated symbols],
                  [no])

  ## Add preprocessor define statement for deprecation warnings
  RNA_FEATURE_IF_ENABLED([warn_deprecated],[
    AC_DEFINE([VRNA_WARN_DEPRECATED], [1], [Warn upon usage of deprecated symbols])
    DEPRECATION_WARNING=-DVRNA_WARN_DEPRECATED
    AX_APPEND_FLAG([-DVRNA_WARN_DEPRECATED], [RNA_CPPFLAGS])
  ])
  AC_SUBST(DEPRECATION_WARNING)
])


#
# Colored TTY output
#

AC_DEFUN([RNA_ENABLE_COLORED_TTY],[

  RNA_ADD_FEATURE([tty_colors],
                  [Colored TTY output],
                  [yes])

  RNA_FEATURE_IF_DISABLED([tty_colors],[
    AC_DEFINE([VRNA_WITHOUT_TTY_COLORS], [1], [Do not use colors for TTY output])
    CONFIG_TTY_COLORS="#define VRNA_WITHOUT_TTY_COLORS"
  ])

  AC_SUBST(CONFIG_TTY_COLORS)
])


#
# Statically linked executables
#

AC_DEFUN([RNA_ENABLE_STATIC_BIN],[

  RNA_ADD_FEATURE([static_executables],
                  [Produce statically linked executable binaries],
                  [no])

  ## Check whether necessary static libraries are present
  RNA_FEATURE_IF_ENABLED([static_executables],[
    SAVED_LDFLAGS=$LDFLAGS
    LDFLAGS="$LDFLAGS -static"

    AC_MSG_NOTICE([Checking possiblity to build statically linked executables using C compiler])
    AC_LANG_PUSH([C])

    AC_RUN_IFELSE([AC_LANG_SOURCE([[ #include <math.h>
                                      int main (void) { return (int)log(1.);} ]])],
                                      [],
                                      [ AC_MSG_ERROR([[
############################################
Failed to statically link C program

Please make sure that static variants for
all libraries that are about to be linked
into the executables are present!

Checked compiler flag setup:
 LDFLAGS = $LDFLAGS 
 LIBS    = $LIBS
 CLFAGS  = $CFLAGS
############################################]]) ],
                                      [enable_static_executables=no])

    AC_LANG_POP([C])
    AC_MSG_NOTICE([Building statically linked C executables seems to work fine])
    LDFLAGS=$SAVED_LDFLAGS
  ])

  RNA_FEATURE_IF_ENABLED([static_executables],[
    SAVED_LDFLAGS=$LDFLAGS
    LDFLAGS="$LDFLAGS -static -lstdc++"

    AC_MSG_NOTICE([Checking possiblity to build statically linked executables using C++ compiler])
    AC_LANG_PUSH([C++])

    AC_RUN_IFELSE([AC_LANG_SOURCE([[ int main (void) { return 0;} ]])],
                                      [],
                                      [ AC_MSG_ERROR([[
############################################
Failed to statically link C++ program

Please make sure that static variants for
all libraries that are about to be linked
into the executables are present!
############################################]]) ],
                                      [enable_static_executables=no])

    AC_LANG_POP([C++])
    AC_MSG_NOTICE([Building statically linked C++ executables seems to work fine])

    LDFLAGS=$SAVED_LDFLAGS
  ])

  AM_CONDITIONAL(VRNA_AM_SWITCH_STATIC_EXECUTABLES, test "x$enable_static_executables" = "xyes")
])


#
# SIMD optimizations
#

AC_DEFUN([RNA_ENABLE_SIMD],[

  RNA_ADD_FEATURE([simd],
                  [Speed-up MFE computations using explicit SIMD instructions.],
                  [yes])

  RNA_ADD_FEATURE([sse],
                  [Deprecated switch for SIMD optimizations. Use --enable-simd/--disable-simd instead],
                  [no])

  AS_IF([test "x$enable_sse" != "xno"],[
    AC_MSG_WARN([[

############################################
Option --enable-sse is deprecated!

Please consider using the successor option --enable-simd instead.
############################################
    ]])
    enable_simd="yes"
    AC_RNA_ADD_WARNING([Deprecated option --enable-simd detected => Please use --enable-simd instead!])
  ])

  AS_IF([test "x$enable_simd" != "xno"],[
    ## Check for all supported SIMD features first
    AC_MSG_CHECKING([compiler support for AVX 512 instructions])

    ac_save_CFLAGS="$CFLAGS"
    CFLAGS="$ac_save_CFLAGS -Werror -mavx512f"
    AC_LANG_PUSH([C])

    AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM([[
                  #include <immintrin.h>
                  #include <limits.h>
                ]],
                  [[__m512i a = _mm512_set1_epi32(INT_MAX);
                    __m512i b = _mm512_set1_epi32(INT_MIN);
                    __m512i c = _mm512_set1_epi32(INT_MIN);
                    __mmask16 mask = _kand_mask16(_mm512_cmplt_epi32_mask(a, c),
                                                  _mm512_cmplt_epi32_mask(b, c));

                    b = _mm512_min_epi32(a, b);
                    int e = _mm512_mask_reduce_min_epi32(mask, b);
                ]])
    ],
    [
      AC_MSG_RESULT([yes])
      AC_DEFINE([VRNA_WITH_SIMD_AVX512], [1], [use AVX 512 implementations])
      ac_simd_capability_avx512f=yes
      SIMD_AVX512_FLAGS="-mavx512f"
    ],
    [
      AC_MSG_RESULT([no])
    ])

    AC_LANG_POP([C])
    CFLAGS="$ac_save_CFLAGS"

    AC_MSG_CHECKING([compiler support for SSE 4.1 instructions])

    ac_save_CFLAGS="$CFLAGS"
    CFLAGS="$ac_save_CFLAGS -Werror -msse4.1"
    AC_LANG_PUSH([C])

    AC_COMPILE_IFELSE(
    [
      AC_LANG_PROGRAM([[
                        #include <smmintrin.h>
                        #include <limits.h>
                      ]],
                        [[__m128i a = _mm_set1_epi32(INT_MAX);
                          __m128i b = _mm_set1_epi32(INT_MIN);
                          b = _mm_min_epi32(a, b);
                      ]])
    ],
    [
      AC_MSG_RESULT([yes])
      AC_DEFINE([VRNA_WITH_SIMD_SSE41], [1], [use SSE 4.1 implementations])
      ac_simd_capability_sse41=yes
      SIMD_SSE41_FLAGS="-msse4.1"
    ],
    [
      AC_MSG_RESULT([no])
    ])

    AC_LANG_POP([C])
    CFLAGS="$ac_save_CFLAGS"
  ])

  RNA_FEATURE_IF_ENABLED([simd],[
    SETUPCFG_SW_SIMD="True"
  ])

  RNA_FEATURE_IF_DISABLED([simd],[
    SETUPCFG_SW_SIMD="False"
  ])

  AC_SUBST(SIMD_AVX512_FLAGS)
  AC_SUBST(SIMD_SSE41_FLAGS)
  AC_SUBST(SETUPCFG_SW_SIMD)
  AM_CONDITIONAL(VRNA_AM_SWITCH_SIMD_AVX512, test "x$ac_simd_capability_avx512f" = "xyes")
  AM_CONDITIONAL(VRNA_AM_SWITCH_SIMD_SSE41, test "x$ac_simd_capability_sse41" = "xyes")
])


AC_DEFUN([RNA_ENABLE_VECTORIZE],[

  RNA_ADD_FEATURE([vectorize],
                  [Apply automatic SIMD vectorization to optimize execution speed],
                  [yes])

  RNA_FEATURE_IF_ENABLED([vectorize],[
    AC_LANG_PUSH([C])
    AX_CHECK_COMPILE_FLAG([-ftree-vectorize], [
      AX_APPEND_FLAG(["-ftree-vectorize"], [RNA_CFLAGS])
    ],[],[],[])
    AC_LANG_POP([C])
  ])
])


#
# Naview layout algorithm (Bruccoleri 1988)
#

AC_DEFUN([RNA_ENABLE_NAVIEW],[

  RNA_ADD_FEATURE([naview],
                  [Naview Layout algorithm],
                  [yes],[],[],
                  [${srcdir}/src/ViennaRNA/plotting/naview/naview.c ${srcdir}/src/ViennaRNA/plotting/naview/naview.h ])

  RNA_FEATURE_IF_ENABLED([naview],[
    AC_DEFINE([VRNA_WITH_NAVIEW_LAYOUT], [1], [Include Naview RNA structure layout algorithm])
    CONFIG_NAVIEW_LAYOUT="#define VRNA_WITH_NAVIEW_LAYOUT"
    SETUPCFG_SW_NAVIEW="True"
    SETUPCFG_NAVIEW_MACRO="naview_macro  = VRNA_WITH_NAVIEW_LAYOUT"
    NAVIEW_DIR=plotting/naview
    AX_APPEND_FLAG([-DVRNA_WITH_NAVIEW_LAYOUT], [RNA_CPPFLAGS])
    AC_CONFIG_FILES([src/ViennaRNA/plotting/naview/Makefile])
  ])

  RNA_FEATURE_IF_DISABLED([naview],[
    SETUPCFG_SW_NAVIEW="False"
  ])

  AC_SUBST(NAVIEW_DIR)
  AC_SUBST(CONFIG_NAVIEW_LAYOUT)
  AC_SUBST(SETUPCFG_NAVIEW_MACRO)
  AC_SUBST(SETUPCFG_SW_NAVIEW)
  AM_CONDITIONAL(VRNA_AM_SWITCH_NAVIEW, test "$enable_naview" != "no")
])


