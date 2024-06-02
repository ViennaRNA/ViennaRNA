include $(top_srcdir)/silent_rules.mk

pkgbindir = $(pkgdatadir)/bin
pkgpyexecdir = $(py3execdir)/RNA

pkgpyexec_DATA = \
    RNAconsensus.py

pscript = \
    b2mt.pl \
    dpzoom.pl \
    mountain.pl \
    relplot.pl \
    rotate_ss.pl \
    cmount.pl \
    coloraln.pl \
    refold.pl \
    switch.pl \
    RNAdesign.pl

bin_PROGRAMS = \
    b2ct \
    popt \
    ct2db

PYTHON_UTILS = RNAconsensus

bin_SCRIPTS = ${PYTHON_UTILS}
pkgpyexec_SCRIPTS = RNAconsensus.py
pkgbin_SCRIPTS = $(pscript)

BUILT_SOURCES = ${PYTHON_UTILS}

AM_CFLAGS = $(RNA_CFLAGS) $(PTHREAD_CFLAGS)
AM_CXXFLAGS = $(RNA_CXXFLAGS) $(PTHREAD_CFLAGS)
AM_CPPFLAGS = $(RNA_CPPFLAGS) -I$(top_srcdir)/src
AM_LDFLAGS = $(RNA_LDFLAGS) $(PTHREAD_LIBS)

if VRNA_AM_SWITCH_STATIC_EXECUTABLES
AM_LDFLAGS += -all-static
endif

GENGETOPT_CMDL =  ct2db_cmdl.c ct2db_cmdl.h

GENGETOPT_FILES =  ct2db.ggo


LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la -lm

# Link against stdc++ if we use SVM
if VRNA_AM_SWITCH_SVM
LDADD += $(SVM_LIBS)
endif

SUFFIXES = .ggo _cmdl.h _cmdl.c

# this is how we make the commandline interface parser
.ggo_cmdl.c:
	$(ggo_verbose)gengetopt -i $<

_cmdl.c_cmdl.h :

ct2db_SOURCES=ct2db_cmdl.c ct2db.c


${PYTHON_UTILS}: python_wrapper
	cp $< $@


EXTRA_DIST = \
    $(pscript) \
    ${GENGETOPT_FILES} \
    ${GENGETOPT_CMDL} \
    RNAconsensus.py

CLEANFILES = ${PYTHON_UTILS}


install-exec-hook:
	$(AM_V_GEN)target_path=${python3dir}; \
	for i in $(PYTHON_UTILS); \
	do \
	  $(SED) -i 's#\$${exec_prefix}#'$(DESTDIR)$(prefix)'#g' $(DESTDIR)$(bindir)/$$i ; \
	done


install-data-hook:
	$(AM_V_GEN)for i in $(pscript); \
	do \
	  $(PERL) -npi -e 's|^#!\S+|#!'$(PERL)'|' $(DESTDIR)$(pkgbindir)/$$i ; \
	done
