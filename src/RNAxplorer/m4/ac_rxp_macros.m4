# RNAxplorer 2016 Ronny Lorenz, Gregor Entzian
#

##--------------------------------##
## Several macros for AC_RXP_INIT ##
##--------------------------------##

AC_DEFUN([AC_RXP_TEST_FILE],[
AC_MSG_CHECKING([for $1])
if test -f $1 ; then
  AC_MSG_RESULT([yes])
  $2
else
  AC_MSG_RESULT([no])
  $3
fi
])

#
# RXP_PACKAGE_IF_ENABLED(package-name, action-if-enabled)
#
# Perform some action if a subpackage is enabled
# Parameters:
#        package-name:   the name of the package that is to be checked
#

AC_DEFUN([RXP_PACKAGE_IF_ENABLED],[
if test "x$with_$1" = "xyes"; then
  $2
fi
])

AC_DEFUN([RXP_PACKAGE_IF_DISABLED],[
if test "x$with_$1" = "xno"; then
  $2
fi
])

#
# RXP_FEATURE_IF_ENABLED(feature-name, action-if-enabled)
#
# Perform some action if a feature is enabled
# Parameters:
#        feature-name:   the name of the feature that is to be checked
#

AC_DEFUN([RXP_FEATURE_IF_ENABLED],[
if test "x$enable_$1" = "xyes"; then
  $2
fi
])

AC_DEFUN([RXP_FEATURE_IF_DISABLED],[
if test "x$enable_$1" = "xno"; then
  $2
fi
])


#
# RXP_ADD_PACKAGE(package-name,
#                 package-description,
#                 default-on,
#                 [action-if-not-default],
#                 [action-if-default],
#                 [files to check for])
#
# This macro handles additional package inclusion
# Parameters:
#       package-name:       a lowercase name of the optional package
#                           which is used for $with_package_name variables
#                           and --with[out]-package-name options in
#                           configure script
#                           The package_name must consist of alphanumeric
#                           characters including the dash only! Every
#                           occuring dash will be replaced by a '-' char
#                           in the --with[out]-package-name option
#
#       package-desciption: a very brief description used for the package
#                           specific help output in configure script
#
#       default-on:         package build | installed by default
#                           Values: "yes" or "no"
#
# Example: RXP_ADD_PACKAGE([foo], [the incredible Foo program], [yes], [with_foo=no], [with_foo=yes], [file1 file2])
#

AC_DEFUN([RXP_ADD_PACKAGE],[

# announce the option to include it in configure script
m4_if([$3],[yes],[
  AC_ARG_WITH(m4_translit([[$1]], [_], [-]),
            [AS_HELP_STRING([--without-m4_translit([$1], [_], [-])], [do not build | install $2])])
  AS_IF([ test "x$with_$1" != "xno"],[with_$1=yes])
  ],[
  AC_ARG_WITH(m4_translit([[$1]], [_], [-]),
            [AS_HELP_STRING([--with-m4_translit([$1], [_], [-])], [build | install $2])])
  AS_IF([ test "x$with_$1" != "xyes"],[with_$1=no])
])

# check if enabling the package makes sense at configure-time
# and deactivate it if not

RXP_PACKAGE_IF_ENABLED([$1],[
  for i in $6; do
    AC_RXP_TEST_FILE([$i],
      [with_$1=$with_$1],
      [with_$1=no])
  done
])

m4_if([$3],[yes],[
  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RXP_PACKAGE_IF_DISABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RXP_PACKAGE_IF_ENABLED([$1],[$5])])
],[
  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RXP_PACKAGE_IF_ENABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RXP_PACKAGE_IF_DISABLED([$1],[$5])])
])

])


AC_DEFUN([RXP_ADD_FEATURE],[

# announce the feature for inclusion in the configure script

m4_if([$3],[yes],[
  AC_ARG_ENABLE(m4_translit([[$1]], [_], [-]),
            [ifelse([$7], [],
                [AS_HELP_STRING([--disable-m4_translit([$1], [_], [-])], [$2])],
                [AS_HELP_STRING([--disable-m4_translit([$1], [_], [-])@<:@=ARG@:>@], [$2 @<:@default ARG=$7@:>@])]
            )])
  AS_IF([ test "x$enable_$1" = "x" ], [ enable_$1=yes ])

  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RXP_FEATURE_IF_ENABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RXP_FEATURE_IF_DISABLED([$1],[$5])])
  ],[
  AC_ARG_ENABLE(m4_translit([[$1]], [_], [-]),
            [ifelse([$7], [],
                [AS_HELP_STRING([--enable-m4_translit([$1], [_], [-])], [$2])],
                [AS_HELP_STRING([--enable-m4_translit([$1], [_], [-])@<:@=ARG@:>@], [$2 @<:@default ARG=$7@:>@])]
            )])
  AS_IF([ (test "x$enable_$1" = "x") || (test "x$enable_$1" = "xno") ], [ enable_$1=no ])
])

# check if enabling the feature makes sense at configure-time
# and deactivate it if not

RXP_FEATURE_IF_ENABLED([$1],[
  for i in $6; do
    AC_RXP_TEST_FILE([$i],
      [enable_$1=$enable_$1],
      [enable_$1=no])
  done
])

m4_if([$3],[yes],[
  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RXP_FEATURE_IF_DISABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RXP_FEATURE_IF_ENABLED([$1],[$5])])
])

])
