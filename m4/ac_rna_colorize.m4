# -*- autoconf -*-

AC_DEFUN([_COLORIZE_RESULT_PREPARE], [

    msg_checking= msg_result_yes= msg_result_no= msg_result_other= msg_reset=
    AS_IF([test "x${CONFIGURE_TTY}" = xyes -o -t 1], [
    msg_begin="`tput smso 2>/dev/null`"
    AS_CASE(["$msg_begin"], ['@<:@'*m],
        [msg_begin="`echo "$msg_begin" | sed ['s/[0-9]*m$//']`"
        msg_checking="${msg_begin}36m"
        AS_IF([test ${TEST_COLORS:+set}], [
        msg_result_yes=[`expr ":$TEST_COLORS:" : ".*:pass=\([^:]*\):"`]
        msg_result_no=[`expr ":$TEST_COLORS:" : ".*:fail=\([^:]*\):"`]
        msg_result_other=[`expr ":$TEST_COLORS:" : ".*:skip=\([^:]*\):"`]
        ])
        msg_result_yes="${msg_begin}${msg_result_yes:-32;1}m"
        msg_result_no="${msg_begin}${msg_result_no:-31;1}m"
        msg_result_other="${msg_begin}${msg_result_other:-34;1}m"
        msg_bold="${msg_begin}1m"
        msg_underline="${msg_begin}4m"
        msg_reset="${msg_begin}m"
        ])
    AS_UNSET(msg_begin)
    ])
    AS_REQUIRE_SHELL_FN([colorize_result],
    [AS_FUNCTION_DESCRIBE([colorize_result], [MSG], [Colorize result])],
        [AS_CASE(["$[]1"],
            [yes], [AS_ECHO(["${msg_result_yes}$[]1${msg_reset}]")],
            [no], [AS_ECHO(["${msg_result_no}$[]1${msg_reset}]")],
            [AS_ECHO(["${msg_result_other}$[]1${msg_reset}]")])])
])

AC_DEFUN([COLORIZE_RESULT], [AC_REQUIRE([_COLORIZE_RESULT_PREPARE])dnl
    AS_LITERAL_IF([$1],
    [m4_case([$1],
        [yes], [AS_ECHO(["${msg_result_yes}$1${msg_reset}"])],
        [no], [AS_ECHO(["${msg_result_no}$1${msg_reset}"])],
        [AS_ECHO(["${msg_result_other}$1${msg_reset}"])])],
    [colorize_result "$1"]) dnl
])


AC_DEFUN([AC_CHECKING],[dnl
AC_REQUIRE([_COLORIZE_RESULT_PREPARE])dnl
AS_MESSAGE([checking ${msg_checking}$1${msg_reset}...])])

AC_DEFUN([AC_MSG_RESULT], [dnl
{ _AS_ECHO_LOG([result: $1])
COLORIZE_RESULT([$1]); dnl
}])



AC_DEFUN([AC_RNA_COLOR_RESULT],[dnl
  AC_REQUIRE([_COLORIZE_RESULT_PREPARE])dnl
  AS_CASE("${$1$2}",
    [yes], [result_$2="${msg_result_yes}${$1$2}${msg_reset}"],
    [no|No*], [result_$2="${msg_result_no}${$1$2}${msg_reset}"],
    [result_$2="${msg_result_other}${$1$2}${msg_reset}"])
])


AC_DEFUN([AC_RNA_COLOR_RESULT_PACKAGE],[
  AC_RNA_COLOR_RESULT([with_], [$1])
])

AC_DEFUN([AC_RNA_COLOR_RESULT_FEATURE],[
  AC_RNA_COLOR_RESULT([enable_], [$1])
])

AC_DEFUN([AC_RNA_COLOR_RESULT_SIMPLE],[
  AC_RNA_COLOR_RESULT([], [$1])
])

AC_DEFUN([AC_RNA_STRING_APPEND_FORMAT_BOLD],[dnl
  AC_REQUIRE([_COLORIZE_RESULT_PREPARE])dnl
  [$1="${$1}${msg_bold}$2${msg_reset}"]
])

AC_DEFUN([AC_RNA_STRING_APPEND_FORMAT_UNDERLINE],[dnl
  AC_REQUIRE([_COLORIZE_RESULT_PREPARE])dnl
  [$1="${$1}${msg_underline}$2${msg_reset}"]
])

AC_DEFUN([RNA_COLOR_MAKE], [

    enable_color_make="yes"

    am__v_GEN_0="@@<:@ -t 1 @:>@ && printf \"  \\033@<:@01;34mGEN\\033@<:@00m      %s\\n\" \$$@@ || echo \"  GEN      \$$@@\";"
    am__v_CC_0="@@<:@ -t 1 @:>@ && printf \"  \\033@<:@01;34mCC\\033@<:@00m       %s\\n\" \$$@@ || echo \"  CC       \$$@@\";"
    am__v_CCLD_0="@@<:@ -t 1 @:>@ && printf \"  \\033@<:@01;34mCCLD\\033@<:@00m     %s\\n\" \$$@@ || echo \"  CCLD     \$$@@\";"
    am__v_CXX_0="@@<:@ -t 1 @:>@ && printf \"  \\033@<:@01;34mCXX\\033@<:@00m      %s\\n\" \$$@@ || echo \"  CXX      \$$@@\";"
    am__v_CXXLD_0="@@<:@ -t 1 @:>@ && printf \"  \\033@<:@01;34mCXXLD\\033@<:@00m    %s\\n\" \$$@@ || echo \"  CXXLD    \$$@@\";"

    AC_SUBST([am__v_GEN_0])
    AC_SUBST([am__v_CC_0])
    AC_SUBST([am__v_CCLD_0])
    AC_SUBST([am__v_CXX_0])
    AC_SUBST([am__v_CXXLD_0])

    AM_CONDITIONAL(WITH_COLOR_MAKE, test "x$enable_color_make" = "xyes")
])
