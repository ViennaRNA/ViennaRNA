dnl -*-autoconf-*-

AC_DEFUN([AC_PROG_CC_FPIC],
[ac_test_CFLAGS=${CFLAGS+set}
ac_save_CFLAGS=$CFLAGS
CFLAGS="-fpic"
AC_CACHE_CHECK(whether $CC accepts -fpic, ac_cv_prog_cc_fpic,
       [_AC_COMPILE_IFELSE([AC_LANG_PROGRAM()], [ac_cv_prog_cc_fpic=yes],
						[ac_cv_prog_cc_fpic=no])])
CFLAGS="$ac_save_CFLAGS"
if test $ac_cv_prog_cc_fpic = yes; then
  FPIC="-fpic"
fi[]dnl
]
)
