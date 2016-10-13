
#
# Check for one project file
#
AC_INIT(src/g2_device.c)

PRODUCT=g2


#
# check cannonical system name,
# and set options for shared libraries (black magic)
#
AC_CANONICAL_HOST
case "$host" in
  *-linux*)
		LD_SHARED='ld'
		LD_SHARED_FLAGS='-shared -soname libg2.so.0'
		LD_SHARED_OPT=''
		LD_SHARED_EXT='so.0.$(G2_VERSION)'
		AC_DEFINE(LINUX);;
    *-aix*)
		LD_SHARED='ld'
		LD_SHARED_FLAGS='-bM:SRE -bnoentry -bexpall'
		LD_SHARED_OPT='$(LDFLAGS) -lc'
		LD_SHARED_EXT='so'
		AC_DEFINE(AIX);;
    *-osf*)
		LD_SHARED='ld'
		LD_SHARED_FLAGS='-shared'
		LD_SHARED_OPT='$(LDFLAGS) -lc'
		LD_SHARED_EXT='so'
		AC_DEFINE(OSF);;
         *)
		LD_SHARED_EXT='so'
		AC_DEFINE(UNIX);;
esac

AC_SUBST(LD_SHARED)
AC_SUBST(LD_SHARED_FLAGS)
AC_SUBST(LD_SHARED_OPT)
AC_SUBST(LD_SHARED_EXT)

#
# check C compiler, preprocesor, etc.
#
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_INSTALL
AC_CHECK_PROG(FIND, find, find, :)
AC_CHECK_PROG(MAKEDEPEND, makedepend, makedepend, :)
AC_CHECK_LIB(dl, dlopen, have_dl=yes, have_dl=no)

#
# Enable user setting 'ar' command
#
if test "x$AR" != "x" ; then
  AR="$AR"
else
  AR=ar
fi
AC_SUBST(AR)

#
# Enable user setting 'ranlib' command
#
if test "x$RANLIB" != "x" ; then
  RANLIB="$RANLIB"
else
  AC_PROG_RANLIB
fi



#
# Try to locate the X Window System include files and libraries
# and add /usr/local to include and lib path and add -lm (for testing)
#
AC_PATH_XTRA
CFLAGS="$CFLAGS $X_CFLAGS -I/usr/local/include"
LIBS="$X_LIBS -L/usr/local/lib"


DO_PS='#'
DO_FIG='#'
DO_X11='#'
DO_XDBE='#'
DO_GD='#'
DO_GIF='#'
DO_WIN32='#'


#
# check supported devices
#
AC_CHECK_LIB(m,main,DO_PS='';DO_FIG='';LIBS="$LIBS -lm";AC_DEFINE(DO_PS) AC_DEFINE(DO_FIG),break)
AC_CHECK_LIB(X11, main, DO_X11='';LIBS="$LIBS -lX11";AC_DEFINE(DO_X11))
AC_CHECK_LIB(Xext, main, DO_XDBE='';LIBS="$LIBS -lXext";AC_DEFINE(DO_XDBE))
AC_CHECK_LIB(gd, main, have_gd=yes, have_gd=no)
if test "$have_gd" = yes
then
	LIBS="$LIBS -lgd"
	AC_CHECK_FUNC(gdImagePng, DO_GD='';AC_DEFINE(DO_GD))
	AC_CHECK_FUNC(gdImageGif, DO_GIF='';AC_DEFINE(DO_GIF))
fi

#
# handling cygwin environment (and gd for cygwin == bgd)
#
AC_CYGWIN
if test "$CYGWIN" = yes
then
	LIBS="$LIBS -L/lib/w32api"
	AC_CHECK_LIB(gdi32, main, DO_WIN32='';LIBS="$LIBS -lgdi32";AC_DEFINE(DO_WIN32))
#	AC_CHECK_LIB(bgd, main, have_bgd=yes, have_bgd=no)
#	if test "$have_bgd" = yes
#	then
#		LIBS="$LIBS -lbgd"
#		AC_CHECK_FUNC(gdImagePng, DO_GD='';AC_DEFINE(DO_GD))
#	fi
fi
   
#
# Declare variables which we want substituted in the Makefile.in's
#
AC_SUBST(DO_PS)
AC_SUBST(DO_FIG)
AC_SUBST(DO_X11)
AC_SUBST(DO_XDBE)
AC_SUBST(DO_WIN32)
AC_SUBST(DO_GD)
AC_SUBST(DO_GIF)

AC_HEADER_STDC
AC_CHECK_HEADERS(limits.h)

#
# finally create makefiles using Makefile.in
#
echo
AC_OUTPUT(Makefile demo/Makefile g2_perl/Makefile.PL g2_python/Makefile)

echo
echo "Supported devices:"
if test "$DO_PS"    != "#"; then echo "    PostScript"; fi
if test "$DO_FIG"   != "#"; then echo "    FIG"; fi
if test "$DO_X11"   != "#"; then echo "    X11"; fi
if test "$DO_XDBE"  != "#"; then echo "    Xdbe"; fi
if test "$DO_WIN32" != "#"; then echo "    Win32"; fi
if test "$DO_GD"    != "#"; then echo "    GD"; fi
if test "$DO_GIF"    != "#"; then echo "    GD (GIF)"; fi

echo
echo "Unsupported devices:"
if test "$DO_PS"     = "#"; then echo "    PostScript"; fi
if test "$DO_FIG"    = "#"; then echo "    FIG"; fi
if test "$DO_X11"    = "#"; then echo "    X11"; fi
if test "$DO_XDBE"   = "#"; then echo "    Xdbe"; fi
if test "$DO_WIN32"  = "#"; then echo "    Win32"; fi
if test "$DO_GD"     = "#"; then echo "    GD"; fi
if test "$DO_GIF"    = "#"; then echo "    GD (GIF)"; fi
		   
echo
echo "Run 'make depend' to create dependencies."
echo
