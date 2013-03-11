pkgbindir = $(pkgdatadir)/bin
pkgbin_PROGRAMS = b2ct popt

pscript = b2mt.pl ct2b.pl dpzoom.pl mountain.pl relplot.pl rotate_ss.pl cmount.pl colorrna.pl coloraln.pl refold.pl switch.pl
pkgbin_SCRIPTS = $(pscript) Fold

EXTRA_DIST = $(pscript) Fold

LDADD = -L$(top_srcdir)/src/ViennaRNA/ -lRNA -lm

AM_CPPFLAGS = -I$(top_srcdir)/src

$(pscript) : Makefile
	$(AM_V_GEN)if test -n "@PerlCmd@"; then \
	  @PerlCmd@ -npi -e 's|^#!\S+|#!@PerlCmd@|' $@ ; \
	fi
