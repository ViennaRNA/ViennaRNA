pkgbindir = $(pkgdatadir)/bin
pkgbin_PROGRAMS = b2ct popt ct2db

pscript = b2mt.pl ct2b.pl dpzoom.pl mountain.pl relplot.pl rotate_ss.pl cmount.pl colorrna.pl coloraln.pl refold.pl switch.pl
pkgbin_SCRIPTS = $(pscript) Fold

AM_CFLAGS = $(OPENMP_CFLAGS)

GENGETOPT_CMDL =  ct2db_cmdl.c ct2db_cmdl.h

GENGETOPT_FILES =  ct2db.ggo

EXTRA_DIST = $(pscript) Fold ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

LDADD = -L$(top_builddir)/src/ViennaRNA -lRNA -lm

AM_CPPFLAGS = -I$(top_srcdir)/src

ct2db_cmdl.h ct2db_cmdl.c:  ct2db.ggo
	$(AM_V_GEN)gengetopt -i ct2db.ggo

ct2db_SOURCES=ct2db_cmdl.c ct2db.c

install-data-hook:
	$(AM_V_GEN)if test -n "@PerlCmd@"; then \
	  for i in $(pscript); \
	  do \
	    @PerlCmd@ -npi -e 's|^#!\S+|#!@PerlCmd@|' $(DESTDIR)$(pkgbindir)/$$i ; \
	  done \
	fi
