#
# SVM support for Lfold -z
#

AC_DEFUN([RNA_ENABLE_SVM],[

  SVM_VERSION=3.23

  RNA_ADD_PACKAGE([svm],
                  [svm classifiers],
                  [yes],[],[],
                  [${srcdir}/src/libsvm-${SVM_VERSION}/svm.cpp ${srcdir}/src/libsvm-${SVM_VERSION}/svm.h])

  RNA_PACKAGE_IF_ENABLED([svm],[
    AC_SUBST([LIBSVM_DIR], [libsvm-${SVM_VERSION}])
    ## substitution for reference manual config
    AC_SUBST([REFDOC_PREDEF_SVM], [VRNA_WITH_SVM])
    AC_DEFINE([VRNA_WITH_SVM], [1], [Compute z-scores for RNALfold])
    CONFIG_SVM="#define VRNA_WITH_SVM"
  ])

  AC_SUBST(CONFIG_SVM)
  AM_CONDITIONAL(VRNA_AM_SWITCH_SVM, test "$with_svm" != "no")
])


#
# JSON serializer/deserializer
#

AC_DEFUN([RNA_ENABLE_JSON],[

  RNA_ADD_PACKAGE([json],
                  [json in/out support],
                  [yes],[],[],
                  [${srcdir}/src/json/json.c ${srcdir}/src/json/json.h])

  RNA_PACKAGE_IF_ENABLED([json],[
    ## substitution for reference manual config
    AC_SUBST([REFDOC_PREDEF_JSON], [VRNA_WITH_JSON_SUPPORT])
    AC_DEFINE([VRNA_WITH_JSON_SUPPORT], [1], [Add JSON support for input and output functions])
    CONFIG_JSON="#define VRNA_WITH_JSON_SUPPORT"
  ])

  AC_SUBST(CONFIG_JSON)
  AM_CONDITIONAL(VRNA_AM_SWITCH_JSON, test "$with_json" != "no")
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
    AC_DEFINE([VRNA_WITH_GSL], [1], [Use GNU Scientific Library])
    GSL_LIBS="-lgsl -lgslcblas"
    CONFIG_GSL="#define VRNA_WITH_GSL"
  ])

  AC_SUBST([GSL_LIBS])
  AC_SUBST([CONFIG_GSL])
  AM_CONDITIONAL(VRNA_AM_SWITCH_GSL, test "x$with_gsl" = "xyes")
])


#
# Boustrophedon scheme for stochastic backtracking
#

AC_DEFUN([RNA_ENABLE_BOUSTROPHEDON],[

  RNA_ADD_FEATURE([boustrophedon],
                  [Boustrophedon scheme for stochastic backtracking],
                  [yes])

  ## Add preprocessor define statement for Boustrophedon scheme in stochastic backtracking in part_func.c
  RNA_FEATURE_IF_ENABLED([boustrophedon],[
    AC_DEFINE([VRNA_WITH_BOUSTROPHEDON], [1], [Use Boustrophedon scheme for stochastic backtracking])
    CONFIG_BOUSTROPHEDON="#define VRNA_WITH_BOUSTROPHEDON"
  ])

  AC_SUBST(CONFIG_BOUSTROPHEDON)
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
      fi
    ])
  ])

  AC_SUBST(CONFIG_OPENMP)
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

  ## Add preprocessor define statement for Boustrophedon scheme in stochastic backtracking in part_func.c
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
# SSE implementations
#

AC_DEFUN([RNA_ENABLE_SSE],[

  RNA_ADD_FEATURE([sse],
                  [Speed-up MFE computations using SIMD instructions. Use one of 'sse41' and 'avx512' to choose an instruction set.],
                  [no],[],[],[sse41])

  AS_IF([test "x$enable_sse" != "xno"],[
    ## Check for all supported SIMD features first
    case $enable_sse in
        avx512)
            AC_MSG_CHECKING([compiler support for AVX 512 instructions])

            ac_save_CFLAGS="$CFLAGS"
            CFLAGS="$ac_save_CFLAGS -mavx512f"
            AC_LANG_PUSH([C])

            AC_COMPILE_IFELSE(
            [
              AC_LANG_PROGRAM([[
                          #include <immintrin.h>
                          #include <limits.h>
                        ]],
                          [[__m512i a = _mm512_set1_epi32(INT_MAX);
                            __m512i b = _mm512_set1_epi32(INT_MIN);
                            b = _mm512_min_epi32(a, b);
                        ]])
            ],
            [
              AC_MSG_RESULT([yes])
              AC_DEFINE([VRNA_WITH_SIMD_AVX512], [1], [use AVX 512 implementations])
              SIMD_FLAGS="-mavx512f"
            ],
            [
              AC_MSG_RESULT([no; using default implementation])
              enable_sse=no
              simd_failed="unable to compile for AVX512 instruction set"
            ])

            AC_LANG_POP([C])
            CFLAGS="$ac_save_CFLAGS"

            ;;
        *)
            AC_MSG_CHECKING([compiler support for SSE 4.1 instructions])

            ac_save_CFLAGS="$CFLAGS"
            CFLAGS="$ac_save_CFLAGS -msse4.1"
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
              SIMD_FLAGS="-msse4.1"
            ],
            [
              AC_MSG_RESULT([no; using default implementation])
              enable_sse=no
              simd_failed="unable to compile for SSE4.1 instruction set"
            ])

            AC_LANG_POP([C])
            CFLAGS="$ac_save_CFLAGS"
            ;;
    esac
  ])

  AS_IF([test "x$enable_sse" != "xno"],[
    AC_DEFINE([VRNA_WITH_SIMD_EXTENSIONS], [1], [use SIMD extensions])
  ])

  AC_SUBST(SIMD_FLAGS)
])


#
# Boustrophedon scheme for stochastic backtracking
#

AC_DEFUN([RNA_ENABLE_VECTORIZE],[

  RNA_ADD_FEATURE([vectorize],
                  [Apply automatic SIMD vectorization to optimize execution speed],
                  [yes])

  ## Add preprocessor define statement for Boustrophedon scheme in stochastic backtracking in part_func.c
  RNA_FEATURE_IF_ENABLED([vectorize],[
    AC_LANG_PUSH([C])
    AX_CHECK_COMPILE_FLAG([-ftree-vectorize], [
      AX_APPEND_FLAG(["-ftree-vectorize"], [RNA_CFLAGS])
    ],[],[],[])
    AC_LANG_POP([C])
  ])
])


