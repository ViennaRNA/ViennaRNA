pkgbindir = $(pkgdatadir)/bin

pscript = b2mt.pl ct2b.pl dpzoom.pl mountain.pl relplot.pl rotate_ss.pl cmount.pl colorrna.pl coloraln.pl refold.pl switch.pl RNAdesign.pl

bin_PROGRAMS = b2ct popt ct2db

pkgbin_SCRIPTS = $(pscript)

AM_CFLAGS = $(RNA_CFLAGS)
AM_CXXFLAGS = $(RNA_CXXFLAGS)
AM_CPPFLAGS = $(RNA_CPPFLAGS) -I$(top_srcdir)/src
AM_LDFLAGS = $(RNA_LDFLAGS)

if VRNA_AM_SWITCH_STATIC_EXECUTABLES
AM_LDFLAGS += -all-static
endif

GENGETOPT_CMDL =  ct2db_cmdl.c ct2db_cmdl.h

GENGETOPT_FILES =  ct2db.ggo

EXTRA_DIST = $(pscript) ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la -lm

SUFFIXES = .ggo _cmdl.h _cmdl.c

# this is how we make the commandline interface parser
.ggo_cmdl.c:
	$(AM_V_GEN)gengetopt -i $<

_cmdl.c_cmdl.h :

ct2db_SOURCES=ct2db_cmdl.c ct2db.c

install-data-hook:
	$(AM_V_GEN)for i in $(pscript); \
	do \
	  $(PERL) -npi -e 's|^#!\S+|#!'$(PERL)'|' $(DESTDIR)$(pkgbindir)/$$i ; \
	done
