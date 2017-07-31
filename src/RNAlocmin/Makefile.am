AM_CPPFLAGS = $(VRNA_CFLAGS) -Wno-write-strings
AM_CXX_FLAGS = -fexceptions
AM_CFLAGS =  -fexceptions

bin_PROGRAMS = RNAlocmin

LDADD = $(VRNA_LIBS)

SUFFIXES = _cmdline.c _cmdline.h .ggo


GENGETOPT_CMDL =  RNAlocmin_cmdline.c RNAlocmin_cmdline.h

GENGETOPT_FILES =  RNAlocmin.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} \
              ${GENGETOPT_CMDL}

## Below are the rules how to create the test suite executables
# how we make the commandline interface parser
## See checkmk man page for details about *.ts files
.ggo_cmdline.c:
	$(AM_V_GEN)gengetopt -i $<

_cmdline.c_cmdline.h :

# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

RNAlocmin_SOURCES = \
  RNAlocmin_cmdline.c \
  barrier_tree.cpp barrier_tree.h \
  findpath_pk.cpp findpath_pk.h \
  flood.cpp flood.h \
  globals.cpp globals.h \
  hash_util.cpp hash_util.h \
  main.cpp \
  move_set_inside.c move_set_inside.h \
  move_set_pk.cpp move_set_pk.h \
  neighbourhood.cpp neighbourhood.h \
  pknots.cpp pknots.h \
  RNAlocmin.cpp RNAlocmin.h \
  treeplot.cpp treeplot.h

dist_man_MANS = RNAlocmin.1
