pkgbindir = $(pkgdatadir)/bin

pscript = b2mt.pl ct2b.pl dpzoom.pl mountain.pl relplot.pl rotate_ss.pl cmount.pl colorrna.pl coloraln.pl refold.pl switch.pl

bin_PROGRAMS = b2ct popt ct2db

pkgbin_SCRIPTS = $(pscript)

AM_CFLAGS = $(OPENMP_CFLAGS) @AM_CFLAGS@
AM_CXXFLAGS = $(OPENMP_CFLAGS) @AM_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@ @FLOAT_PF_FLAG@ @DEPRECATION_WARNING@
AM_LDFLAGS = @AM_LDFLAGS@

GENGETOPT_CMDL =  ct2db_cmdl.c ct2db_cmdl.h

GENGETOPT_FILES =  ct2db.ggo

EXTRA_DIST = $(pscript) ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la -lm

ct2db_cmdl.h ct2db_cmdl.c:  ct2db.ggo
	$(AM_V_GEN)gengetopt -i ct2db.ggo

ct2db_SOURCES=ct2db_cmdl.c ct2db.c

install-data-hook:
	$(AM_V_GEN)for i in $(pscript); \
	do \
	  $(PERL) -npi -e 's|^#!\S+|#!'$(PERL)'|' $(DESTDIR)$(pkgbindir)/$$i ; \
	done
