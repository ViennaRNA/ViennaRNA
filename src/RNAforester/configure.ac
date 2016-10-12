dnl Process this file with autoconf to produce a configure script.
AC_REVISION($Revision: 1.13 $)
AC_INIT(RNAforester,2.0,sschirme@techfak.uni-bielefeld.de,RNAforester)
AC_PREREQ([2.59])

AM_SILENT_RULES([yes])
AM_INIT_AUTOMAKE([1.11])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AM_CONFIG_HEADER(config.h)
		
dnl Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_LANG(C++)
AM_PROG_LEX
AC_PROG_YACC

dnl Checks for libraries.
dnl AC_CHECK_LIB(g2,g2_open_vd,,AC_MSG_WARN(Compiling without support for graphical output. Install the g2 library that comes with the RNAforester distribution or download the latest version at http://g2.sourceforge.net/))

AC_OPENMP

# Search for the RNAlib2 and set variables for Makefiles.
PKG_CHECK_MODULES([VRNA], [RNAlib2],
                  [
                    AC_DEFINE([HAVE_LIBRNA], [], [Link against RNAlib])
                    CPPFLAGS=["${VRNA_CFLAGS}  ${CPPFLAGS}"]
                    CFLAGS=["${VRNA_CFLAGS}  ${CFLAGS}"]
                    LIBS=["${VRNA_LIBS} ${LIBS}"]
                    AC_SUBST([VRNA_CFLAGS])
                    AC_SUBST([VRNA_LIBS])
                  ],
                  [])


AC_ARG_WITH([xml2-lib],
            [AS_HELP_STRING([--with-xml2-lib],[where to search for libxml2.so])],
	    [LDFLAGS="-L$withval/lib -R$withval/lib ${LDFLAGS}"
	    AC_CHECK_LIB([xml2],[main],AC_DEFINE([HAVE_LIBXML2],[1],[Define if you have libxml2]),AC_MSG_FAILURE([--with-xml2-lib was given but test for libxml2 failed]))
	    LIBS="-lxml2 ${LIBS}"
	    CPPFLAGS="-I$withval/include/libxml2/ ${CPPFLAGS}"],
	    [])

AC_ARG_WITH([xmlplusplus-lib],
            [AS_HELP_STRING([--with-xmlplusplus-lib],[where to search for libxml++-1.0.so])],
	    [LDFLAGS="-L$withval/lib -R$withval/lib ${LDFLAGS}"
	    AC_CHECK_LIB([xml++-1.0],[main],AC_DEFINE([HAVE_LIBXMLPLUSPLUS],[1],[Define if you have libxml++-1.0]),[AC_MSG_FAILURE([--with-xmlplusplus-lib was given but test for libxml++-1.0 failed])],[])
	    LIBS="-lxml++-1.0 ${LIBS}"
	    CPPFLAGS="-I$withval/include/libxml++-1.0/ ${CPPFLAGS}"],
	    [])

save_CPPFLAGS=$CPPFLAGS	
CPPFLAGS=["-Ig2-0.72/src/ ${CPPFLAGS}"]

dnl AC_CHECK_LIB(gd,main,,)
		  
dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADER(g2.h,AC_DEFINE(HAVE_LIBG2,1,Define to 1 if you have the `g2' library (-lg2).),AC_MSG_WARN(Compiling without support for graphical output. Install the g2 library that comes with the RNAforester distribution or download the latest version at http://g2.sourceforge.net/))
dnl AC_CHECK_HEADER(g2_PS.h,,AC_MSG_WARN(Compiling without support for graphical output. Install the g2 library that comes with the RNAforester distribution or download the latest version at http://g2.sourceforge.net/))

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T

CPPFLAGS="$save_CPPFLAGS -I../g2-0.72/include/ "
LDFLAGS=["-L../g2-0.72 ${LDFLAGS}"]

if ${HAVE_LIBG2}; then
  LIBS=["-lg2 ${LIBS}"]
fi	
											
AC_SUBST(CXX)
AC_SUBST(CPPFLAGS)
AC_SUBST(LIBS)
AC_SUBST(LDFLAGS)

dnl Checks for library functions.

AC_CONFIG_SUBDIRS(g2-0.72)
AC_OUTPUT(Makefile src/Makefile man/Makefile)














