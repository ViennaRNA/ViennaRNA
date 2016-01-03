
#
# SVM support for Lfold -z
#

AC_DEFUN([RNA_ENABLE_SVM],[

  SVM_VERSION=3.20

  RNA_ADD_PACKAGE([svm],
                  [svm classifiers],
                  [yes],
                  [with_svm=no],
                  [with_svm=yes],
                  [${srcdir}/src/libsvm-${SVM_VERSION}/svm.cpp ${srcdir}/src/libsvm-${SVM_VERSION}/svm.h])

  RNA_PACKAGE_IF_ENABLED([svm],[
    AC_SUBST([LIBSVM_DIR], [libsvm-${SVM_VERSION}])
    AC_SUBST([WITH_SVM], [USE_SVM])
    AC_DEFINE([USE_SVM], [1], [Compute z-scores for RNALfold])
  ])

  AM_CONDITIONAL(WITH_LIBSVM, test "$with_svm" != "no")
])


#
# JSON serializer/deserializer
#

AC_DEFUN([RNA_ENABLE_JSON],[

  RNA_ADD_PACKAGE([json],
                  [json in/out support],
                  [yes],
                  [with_json=no],
                  [with_json=yes],
                  [${srcdir}/src/json/json.c ${srcdir}/src/json/json.h])

  RNA_PACKAGE_IF_ENABLED([json],[
    AC_SUBST([WITH_JSON], [WITH_JSON_SUPPORT])
    AC_DEFINE([WITH_JSON_SUPPORT], [1], [Add JSON support for input and output functions])
  ])

  AM_CONDITIONAL(WITH_JSON, test "$with_json" != "no")
])


#
# GSL support for RNApvmin
#

AC_DEFUN([RNA_ENABLE_GSL],[

  RNA_ADD_PACKAGE([gsl],
                  [GNU Scientific Library],
                  [yes],
                  [with_gsl=no],
                  [with_gsl=yes])

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
    AC_DEFINE([WITH_GSL], [1], [Use GNU Scientific Library])
  ])

  AM_CONDITIONAL(WITH_GSL, test "$with_gsl" != "no")
])


#
# Boustrophedon scheme for stochastic backtracking
#

AC_DEFUN([RNA_ENABLE_BOUSTROPHEDON],[

  RNA_ADD_FEATURE([boustrophedon],
                  [Boustrophedon scheme for stochastic backtracking],
                  [yes],
                  [enable_boustrophedon=no],
                  [enable_boustrophedon=yes])

  ## Add preprocessor define statement for Boustrophedon scheme in stochastic backtracking in part_func.c
  RNA_FEATURE_IF_ENABLED([boustrophedon],[
    AC_DEFINE([WITH_BOUSTROPHEDON], [1], [Use Boustrophedon scheme for stochastic backtracking])
  ])
])


#
# Generic Hard Constraints
#

AC_DEFUN([RNA_ENABLE_GEN_HC],[

  RNA_ADD_FEATURE([gen_hard_constraints],
                  [Generic hard constraints],
                  [no],
                  [enable_gen_hard_constraints=yes],
                  [enable_gen_hard_constraints=no])

  ## Add preprocessor define statement for generlaized hard constraints feature
  RNA_FEATURE_IF_ENABLED([gen_hard_constraints],[
    AC_DEFINE([WITH_GEN_HC], [1], [Provide generic hard constraints])
    GENERIC_HC_DEF=-DWITH_GEN_HC
  ])
  AC_SUBST(GENERIC_HC_DEF)
])


#
# OpenMP support
#

AC_DEFUN([RNA_ENABLE_OPENMP],[

  ## Add linker flag for OpenMP in pkg-config file
  RNA_FEATURE_IF_ENABLED([openmp],[
    LIBGOMPFLAG=-lgomp
  ])
  AC_SUBST(LIBGOMPFLAG)
])


#
#
#

AC_DEFUN([RNA_ENABLE_FLOATPF],[

  RNA_ADD_FEATURE([floatpf],
                  [Floating point precision in partition function computations],
                  [no],
                  [enable_floatpf=yes],
                  [enable_floatpf=no])

  # Handle floating point precision flag
  RNA_FEATURE_IF_ENABLED([floatpf],[
    AC_DEFINE([USE_FLOAT_PF], [1], [Use floating point precision in partition function computations])
  
    AC_SUBST([WITH_FLOAT_PF], [USE_FLOAT_PF])
    AC_SUBST([FLOAT_PF_FLAG], [-DUSE_FLOAT_PF])
##  AX_APPEND_FLAG([-DUSE_FLOAT_PF], [AM_CPPFLAGS])
  ])


])
