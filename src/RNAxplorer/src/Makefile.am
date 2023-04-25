
AM_CPPFLAGS = $(CFLAGS) $(VRNA_CFLAGS) 

# convenience library
noinst_LTLIBRARIES = libRNAxplorer_conv.la

libRNAxplorer_conv_la_SOURCES = \
  RNAwalk.c \
  RNAwalk.h \
  meshpoint.c \
  meshpoint.h \
  dist_class_sc.c \
  dist_class_sc.h \
  distorted_sampling.c \
  distorted_sampling.h \
  distorted_samplingMD.c \
  distorted_samplingMD.h \
  repellant_sampling.c \
  repellant_sampling.h \
  barrier_lower_bound.c \
  barrier_lower_bound.h \
  PathFinder.c \
  PathFinder.h \
  gradient_walker.c \
  gradient_walker.h

bin_PROGRAMS = RNAxplorer

SUFFIXES = _cmdl.c _cmdl.h .ggo

RNAxplorer_LDADD = \
    libRNAxplorer_conv.la \
    $(VRNA_LIBS)

GENGETOPT_CMDL = \
    RNAxplorer_cmdl.c RNAxplorer_cmdl.h

GENGETOPT_FILES = \
    RNAxplorer.ggo

EXTRA_DIST = \
    ${GENGETOPT_FILES} \
    ${GENGETOPT_CMDL}

# how we make the commandline interface parser
.ggo_cmdl.c:
	$(AM_V_GEN)gengetopt -i $<

_cmdl.c_cmdl.h :
## Recover from the removal of $@
	@if test -f $@; then : touch $@ ; else \
          rm -f $< ; \
          $(MAKE) $(AM_MAKEFLAGS) $< ; \
        fi


# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)


RNAxplorer_SOURCES = \
    ${GENGETOPT_CMDL} \
    RNAxplorer.c

