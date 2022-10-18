# SYNOPSIS
#
#   AX_BLAS([ACTION-IF-FOUND[, ACTION-IF-NOT-FOUND]])
#
# DESCRIPTION
#
#   This macro looks for a library that implements the BLAS linear-algebra
#   interface (see http://www.netlib.org/blas/). On success, it sets the
#   BLAS_LIBS output variable to hold the requisite library linkages.
#
#   To link with BLAS, you should link with:
#
#       $BLAS_LIBS $LIBS
#
#   in that order. 
#
#   To call a BLAS function, use the F77_FUNC macro, defined in config.h
#
#   Many libraries are searched for, from ATLAS to CXML to ESSL. The user
#   may also use --with-blas=<lib> in order to use some specific BLAS
#   library <lib>. 
#
#   ACTION-IF-FOUND is a list of shell commands to run if a BLAS library is
#   found, and ACTION-IF-NOT-FOUND is a list of commands to run it if it is
#   not found. If ACTION-IF-FOUND is not specified, the default action will
#   define HAVE_BLAS.
#
#   This macro requires autoconf 2.50 or later.
#
#   The macro is a modified version of ACX_BLAS, from the autoconf macro
#   archive.  The original macro depends on a Fortran compiler, and this
#   version does not.  This version has not been tested extensively, but
#   it is known to work with ATLAS and vecLib.
#
# LAST MODIFICATION
#
#   2008-12-29
#
# COPYLEFT
#
#   Copyright (c) 2008 Patrick O. Perry  <patperry@stanfordalumni.org>
#   Copyright (c) 2008 Steven G. Johnson <stevenj@alum.mit.edu>
#
#   This program is free software: you can redistribute it and/or modify it
#   under the terms of the GNU General Public License as published by the
#   Free Software Foundation, either version 3 of the License, or (at your
#   option) any later version.
#
#   This program is distributed in the hope that it will be useful, but
#   WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
#   Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program. If not, see <http://www.gnu.org/licenses/>.
#
#   As a special exception, the respective Autoconf Macro's copyright owner
#   gives unlimited permission to copy, distribute and modify the configure
#   scripts that are the output of Autoconf when processing the Macro. You
#   need not follow the terms of the GNU General Public License when using
#   or distributing such scripts, even though portions of the text of the
#   Macro appear in them. The GNU General Public License (GPL) does govern
#   all other use of the material that constitutes the Autoconf Macro.
#
#   This special exception to the GPL applies to versions of the Autoconf
#   Macro released by the Autoconf Macro Archive. When you make and
#   distribute a modified version of the Autoconf Macro, you may extend this
#   special exception to the GPL to apply to your modified version as well.

AC_DEFUN([AX_BLAS], [
AC_PREREQ(2.50)
ax_blas_ok=no

AC_ARG_WITH(blas,
	[AC_HELP_STRING([--with-blas=<lib>], [use BLAS library <lib>])])
case $with_blas in
	yes | "") ;;
	no) acx_blas_ok=disable ;;
	-* | */* | *.a | *.so | *.so.* | *.o) BLAS_LIBS="$with_blas" ;;
	*) BLAS_LIBS="-l$with_blas" ;;
esac

_AX_BLAS(gemm_, [ax_blas_ok=yes; ax_blas_underscore=yes],
    [_AX_BLAS(gemm, [ax_blas_ok=yes; ax_blas_underscore=no])])

AC_SUBST(BLAS_LIBS)

if test x"$ax_blas_ok" = xyes; then
    ifelse([$1],, [
        AC_DEFINE(HAVE_BLAS,1, [Define if you have a BLAS library.])
        AH_TEMPLATE([F77_FUNC], 
                [Define to a macro mangling the given Fortan function name])
        if test x"$ax_blas_underscore" = xyes; then
            AC_DEFINE([F77_FUNC(name)], [name ## _])
        else
            AC_DEFINE([F77_FUNC(name)], [name])
        fi
        ],
        [$1])
    :
else
    ax_blas_ok=no
    $2
fi

])dnl AX_BLAS


# _AX_BLAS(GEMM, [ACTION-IF-FOUND], [ACTION-IF-NOT-FOUND])
# --------------------------------------------------------
# Look for a BLAS library in all of the standard places by checking for the 
# functions s$GEMM and d$GEMM.  On success, define BLAS_LIBS, set
# ax_blas_ok to yes, and execute ACTION-IF-FOUND.  On failure, set ax_blas_ok
# to no and execute ACTION-IF-NOT-FOUND.
AC_DEFUN([_AX_BLAS], [
AC_PREREQ(2.50)
ax_blas_ok=no

# Get fortran linker names of BLAS functions to check for.
sgemm=s$1
dgemm=d$1

ax_blas_save_LIBS="$LIBS"
LIBS="$LIBS $FLIBS"

# First, check BLAS_LIBS environment variable
if test $ax_blas_ok = no; then
if test "x$BLAS_LIBS" != x; then
	save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS"
	AC_MSG_CHECKING([for $sgemm in $BLAS_LIBS])
	AC_TRY_LINK_FUNC($sgemm, [ax_blas_ok=yes], [BLAS_LIBS=""])
	AC_MSG_RESULT($ax_blas_ok)
	LIBS="$save_LIBS"
fi
fi

# BLAS linked to by default?  (happens on some supercomputers)
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="$LIBS"
	AC_CHECK_FUNC($sgemm, [ax_blas_ok=yes])
	LIBS="$save_LIBS"
fi

# BLAS in ATLAS library? (http://math-atlas.sourceforge.net/)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(atlas, ATL_xerbla,
		[AC_CHECK_LIB(f77blas, $sgemm,
		[AC_CHECK_LIB(cblas, cblas_dgemm,
			[ax_blas_ok=yes
			 BLAS_LIBS="-lcblas -lf77blas -latlas"],
			[], [-lf77blas -latlas])],
			[], [-latlas])])
fi

# BLAS in PhiPACK libraries? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(dgemm, $dgemm,
		[AC_CHECK_LIB(sgemm, $sgemm,
			[ax_blas_ok=yes; BLAS_LIBS="-lsgemm -ldgemm -lblas"],
			[], [-lblas])],
			[], [-lblas])])
fi

# BLAS in Intel MKL library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(mkl, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lmkl"])
fi

# BLAS in Apple vecLib library?
if test $ax_blas_ok = no; then
	save_LIBS="$LIBS"; LIBS="-framework vecLib $LIBS"
	AC_CHECK_FUNC($sgemm, [ax_blas_ok=yes;BLAS_LIBS="-framework vecLib"])
	LIBS="$save_LIBS"
fi

# BLAS in Alpha CXML library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(cxml, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-lcxml"])
fi

# BLAS in Alpha DXML library? (now called CXML, see above)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(dxml, $sgemm, [ax_blas_ok=yes;BLAS_LIBS="-ldxml"])
fi

# BLAS in Sun Performance library?
if test $ax_blas_ok = no; then
	if test "x$GCC" != xyes; then # only works with Sun CC
		AC_CHECK_LIB(sunmath, acosp,
			[AC_CHECK_LIB(sunperf, $sgemm,
        			[BLAS_LIBS="-xlic_lib=sunperf -lsunmath"
                                 ax_blas_ok=yes],[],[-lsunmath])])
	fi
fi

# BLAS in SCSL library?  (SGI/Cray Scientific Library)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(scs, $sgemm, [ax_blas_ok=yes; BLAS_LIBS="-lscs"])
fi

# BLAS in SGIMATH library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(complib.sgimath, $sgemm,
		     [ax_blas_ok=yes; BLAS_LIBS="-lcomplib.sgimath"])
fi

# BLAS in IBM ESSL library? (requires generic BLAS lib, too)
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm,
		[AC_CHECK_LIB(essl, $sgemm,
			[ax_blas_ok=yes; BLAS_LIBS="-lessl -lblas"],
			[], [-lblas $FLIBS])])
fi

# Generic BLAS library?
if test $ax_blas_ok = no; then
	AC_CHECK_LIB(blas, $sgemm, [ax_blas_ok=yes; BLAS_LIBS="-lblas"])
fi

LIBS="$ax_blas_save_LIBS"

# Finally, execute ACTION-IF-FOUND/ACTION-IF-NOT-FOUND:
if test x"$ax_blas_ok" = xyes; then
    $2
    :
else
    ax_blas_ok=no
    $3
fi
])dnl _AX_BLAS
