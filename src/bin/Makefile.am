AM_CFLAGS = $(OPENMP_CFLAGS) @AM_CFLAGS@
AM_CXXFLAGS = $(OPENMP_CFLAGS) @AM_CXXFLAGS@
AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@ @FLOAT_PF_FLAG@

bin_PROGRAMS = RNAfold RNAeval RNAheat RNApdist RNAdistance RNAinverse \
    RNAplot RNAsubopt RNALfold RNAcofold RNApaln RNAduplex \
    RNAalifold RNAplfold RNAup RNAaliduplex RNALalifold RNA2Dfold \
    RNAparconv RNAPKplex RNAplex RNAsnoop RNApvmin

LDADD = $(OPENMP_CFLAGS) $(top_builddir)/src/ViennaRNA/libRNA_conv.la

GENGETOPT_CMDL =  RNAfold_cmdl.c RNAfold_cmdl.h\
              RNAalifold_cmdl.c RNAalifold_cmdl.h\
              RNAcofold_cmdl.c  RNAcofold_cmdl.h\
              RNAup_cmdl.c  RNAup_cmdl.h\
              RNALfold_cmdl.c RNALfold_cmdl.h\
              RNALalifold_cmdl.c RNALalifold_cmdl.h\
              RNAplfold_cmdl.c RNAplfold_cmdl.h\
              RNAduplex_cmdl.c RNAduplex_cmdl.h\
              RNAaliduplex_cmdl.c RNAaliduplex_cmdl.h\
              RNAeval_cmdl.c RNAeval_cmdl.h\
              RNAheat_cmdl.c RNAheat_cmdl.h\
              RNAsubopt_cmdl.c RNAsubopt_cmdl.h\
              RNA2Dfold_cmdl.c RNA2Dfold_cmdl.h\
              RNAparconv_cmdl.c RNAparconv_cmdl.h\
              RNAinverse_cmdl.c RNAinverse_cmdl.h\
              RNApaln_cmdl.c RNApaln_cmdl.h\
              RNAplot_cmdl.c RNAplot_cmdl.h\
              RNAPKplex_cmdl.c RNAPKplex_cmdl.h\
              RNAplex_cmdl.c RNAplex_cmdl.h\
              RNAsnoop_cmdl.c RNAsnoop_cmdl.h\
              RNAdistance_cmdl.c RNAdistance_cmdl.h\
              RNApdist_cmdl.c RNApdist_cmdl.h \
              RNApvmin_cmdl.c RNApvmin_cmdl.h

GENGETOPT_FILES =  RNAfold.ggo RNAalifold.ggo RNAcofold.ggo\
              RNAup.ggo RNALfold.ggo RNALalifold.ggo\
              RNAplfold.ggo RNAduplex.ggo RNAaliduplex.ggo\
              RNAeval.ggo RNAheat.ggo RNAsubopt.ggo RNA2Dfold.ggo\
              RNAparconv.ggo RNAinverse.ggo RNApaln.ggo RNAplot.ggo\
              RNAPKplex.ggo RNAplex.ggo RNAsnoop.ggo RNAdistance.ggo\
              RNApdist.ggo RNApvmin.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)

# how we make the commandline interface parser
%_cmdl.h %_cmdl.c: %.ggo
	$(AM_V_GEN)gengetopt -i $<

RNAfold_SOURCES=RNAfold_cmdl.c RNAfold.c

RNAsubopt_SOURCES=RNAsubopt_cmdl.c RNAsubopt.c

RNAalifold_SOURCES=RNAalifold_cmdl.c RNAalifold.c

RNAcofold_SOURCES=RNAcofold_cmdl.c RNAcofold.c

RNAup_SOURCES=RNAup_cmdl.c RNAup.c

RNALfold_SOURCES=RNALfold_cmdl.c RNALfold.c

RNALalifold_SOURCES=RNALalifold_cmdl.c RNALalifold.c

RNAplfold_SOURCES=RNAplfold_cmdl.c RNAplfold.c

RNAduplex_SOURCES=RNAduplex_cmdl.c RNAduplex.c

RNAaliduplex_SOURCES=RNAaliduplex_cmdl.c RNAaliduplex.c

RNAeval_SOURCES=RNAeval_cmdl.c RNAeval.c

RNAheat_SOURCES=RNAheat_cmdl.c RNAheat.c

RNA2Dfold_SOURCES=RNA2Dfold_cmdl.c RNA2Dfold.c

RNAparconv_SOURCES=RNAparconv_cmdl.c RNAparconv.c

RNAinverse_SOURCES=RNAinverse_cmdl.c RNAinverse.c

RNApaln_SOURCES=RNApaln_cmdl.c RNApaln.c

RNAplot_SOURCES=RNAplot_cmdl.c RNAplot.c

RNAPKplex_SOURCES=RNAPKplex_cmdl.c RNAPKplex.c

RNAplex_SOURCES=RNAplex_cmdl.c RNAplex.c

RNAsnoop_SOURCES=RNAsnoop_cmdl.c RNAsnoop.c

RNAdistance_SOURCES=RNAdistance_cmdl.c RNAdistance.c

RNApdist_SOURCES=RNApdist_cmdl.c RNApdist.c

RNApvmin_SOURCES=RNApvmin_cmdl.c RNApvmin.c
