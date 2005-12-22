pkgbindir = $(pkgdatadir)/bin
pkgbin_PROGRAMS = b2ct popt

pscript = b2mt.pl ct2b.pl dpzoom.pl mountain.pl relplot.pl rotate_ss.pl cmount.pl colorrna.pl coloraln.pl
pkgbin_SCRIPTS = $(pscript) Fold

EXTRA_DIST = $(pscript) Fold

LDADD = -L../lib -lRNA -lm

INCLUDES = -I../H

$(pscript) : Makefile
	if test -n "@PerlCmd@"; then \
	  @PerlCmd@ -npi -e 's|^#!\S+|#!@PerlCmd@|' $@ ; \
	fi
