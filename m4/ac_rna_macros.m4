# ViennaRNA Package 2011 Ronny Lorenz
#

##--------------------------------##
## Several macros for AC_RNA_INIT ##
##--------------------------------##

AC_DEFUN([AC_RNA_APPEND_VAR_COMMA],[
  if test "x$$1" != "x" ; then
    $1="${$1}, "
  fi
  $1="${$1}$2"
])

AC_DEFUN([AC_RNA_TEST_FILE],[
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
# RNA_PACKAGE_IF_ENABLED(package-name, action-if-enabled)
#
# Perform some action if a subpackage is enabled
# Parameters:
#        package-name:   the name of the package that is to be checked
#

AC_DEFUN([RNA_PACKAGE_IF_ENABLED],[
if test "x$with_$1" = "xyes"; then
  $2
fi
])

AC_DEFUN([RNA_PACKAGE_IF_DISABLED],[
if test "x$with_$1" = "xno"; then
  $2
fi
])

#
# RNA_FEATURE_IF_ENABLED(feature-name, action-if-enabled)
#
# Perform some action if a feature is enabled
# Parameters:
#        feature-name:   the name of the feature that is to be checked
#

AC_DEFUN([RNA_FEATURE_IF_ENABLED],[
if test "x$enable_$1" = "xyes"; then
  $2
fi
])

AC_DEFUN([RNA_FEATURE_IF_DISABLED],[
if test "x$enable_$1" = "xno"; then
  $2
fi
])


#
# RNA_ADD_PACKAGE(package-name,
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
#       is-default:         package build | installed by default
#                           Values: "yes" or "no"
#
# Example: RNA_ADD_PACKAGE([foo], [the incredible Foo program], [yes], [with_bar=no], [], [file1 file2])
#

AC_DEFUN([RNA_ADD_PACKAGE],[

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

RNA_PACKAGE_IF_ENABLED([$1],[
  for i in $6; do
    AC_RNA_TEST_FILE([$i],
      [with_$1=$with_$1],
      [with_$1=no])
  done
])

m4_if([$3],[yes],[
  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RNA_PACKAGE_IF_DISABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RNA_PACKAGE_IF_ENABLED([$1],[$5])])
],[
  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RNA_PACKAGE_IF_ENABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RNA_PACKAGE_IF_DISABLED([$1],[$5])])
])

])


#
# RNA_ADD_FEATURE(feature-name,
#                 feature-description,
#                 is-default,
#                 [action-if-no-default],
#                 [action-if-default])
#
# Add a feature to the ViennaRNA Package
#

AC_DEFUN([RNA_ADD_FEATURE],[

# announce the feature for inclusion in the configure script

m4_if([$3],[yes],[
  AC_ARG_ENABLE(m4_translit([[$1]], [_], [-]),
            [ifelse([$6], [],
                [AS_HELP_STRING([--disable-m4_translit([$1], [_], [-])], [$2])],
                [AS_HELP_STRING([--disable-m4_translit([$1], [_], [-])@<:@=ARG@:>@], [$2 @<:@default ARG=$6@:>@])]
            )])
  AS_IF([ test "x$enable_$1" = "x" ], [ enable_$1=yes ])

  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RNA_FEATURE_IF_ENABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RNA_FEATURE_IF_DISABLED([$1],[$5])])
  ],[
  AC_ARG_ENABLE(m4_translit([[$1]], [_], [-]),
            [ifelse([$6], [],
                [AS_HELP_STRING([--enable-m4_translit([$1], [_], [-])], [$2])],
                [AS_HELP_STRING([--enable-m4_translit([$1], [_], [-])@<:@=ARG@:>@], [$2 @<:@default ARG=$6@:>@])]
            )])
  AS_IF([ (test "x$enable_$1" = "x") || (test "x$enable_$1" = "xno") ], [ enable_$1=no ])

  # take action depending on configure-time user settings
  m4_if(["x$4"], ["x"],[], [ RNA_FEATURE_IF_DISABLED([$1],[$4])])
  m4_if(["x$5"], ["x"],[], [ RNA_FEATURE_IF_ENABLED([$1],[$5])])
])

])

#
# RNA_FEATURE_POST
#
# Set postconditions after all features have been announced and
# initialized. This ensures, that we pass CFLAGS, CXXFLAGS,
# and LDFLAGS to the subpackages afterwards
# This macro also substitutes the automake variables
# AM_CPPFLAGS, AM_CFLAGS, AM_CXXFLAGS, and AM_LDFLAGS
#
# As a bonus, we subsitute RNALIB_CFLAGS and
# RNALIB_LIBS to pass them further to our subprojects
#

AC_DEFUN([RNA_FEATURE_POST],[

  # substitute automake variables in case we set them somewhere
  # in our autoconf mess

  AC_SUBST([AM_CPPFLAGS])
  AC_SUBST([AM_CFLAGS])
  AC_SUBST([AM_CXXFLAGS])
  AC_SUBST([AM_LDFLAGS])

  AC_SUBST([RNA_CFLAGS])
  AC_SUBST([RNA_CXXFLAGS])
  AC_SUBST([RNA_CPPFLAGS])
  AC_SUBST([RNA_LDFLAGS])

  VRNA_LIBS=" ${RNA_LDFLAGS} -L\$(top_builddir)/../../src/ViennaRNA -lRNA ${LIBGOMPFLAG}"
  VRNA_CFLAGS=" -I\$(top_srcdir)/../../src/ViennaRNA -I\$(top_srcdir)/../../src ${RNA_CPPFLAGS}"

  AC_SUBST([VRNA_LIBS])
  AC_SUBST([VRNA_CFLAGS])

  # Replace/add flags in/to ac_configure_args
  for var in CC CXX CFLAGS CXXFLAGS LDFLAGS AR RANLIB NM VRNA_CFLAGS VRNA_LIBS; do
    value=`eval echo \\${${var}}`
    if test "x$value" != "x" ; then
      AS_CASE([$ac_configure_args],
        [*${var}=*], [ac_configure_args=`AS_ECHO "$ac_configure_args" | $SED ["s|${var}=[^']*|${var}=${value}|"]`],
        [AS_VAR_APPEND([ac_configure_args],[" '${var}=${value}'"])]
      )
    fi
  done

])

