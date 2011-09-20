

if WITH_LIBSVM
  RNALfold_LINK = ${CXXLD} -o $@
else
  RNALfold_LINK = ${CCLD} -o $@
endif


bin_PROGRAMS = RNAfold RNAeval RNAheat RNApdist RNAdistance RNAinverse \
    RNAplot RNAsubopt RNALfold RNAcofold RNApaln RNAduplex \
    RNAalifold RNAplfold RNAup RNAaliduplex RNALalifold RNA2Dfold \
    RNAparconv RNAPKplex RNAplex RNAsnoop

LDADD = $(OPENMP_CFLAGS) ../lib/libRNA.a -lm

AM_CFLAGS = $(OPENMP_CFLAGS)

INCLUDES = -I$(srcdir)/../H -I$(srcdir)/../lib

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
              RNApdist_cmdl.c RNApdist_cmdl.h

GENGETOPT_FILES =  RNAfold.ggo RNAalifold.ggo RNAcofold.ggo\
              RNAup.ggo RNALfold.ggo RNALalifold.ggo\
              RNAplfold.ggo RNAduplex.ggo RNAaliduplex.ggo\
              RNAeval.ggo RNAheat.ggo RNAsubopt.ggo RNA2Dfold.ggo\
              RNAparconv.ggo RNAinverse.ggo RNApaln.ggo RNAplot.ggo\
              RNAPKplex.ggo RNAplex.ggo RNAsnoop.ggo RNAdistance.ggo\
              RNApdist.ggo

EXTRA_DIST =  ${GENGETOPT_FILES} ${GENGETOPT_CMDL}

# remove all generated man pages
maintainer-clean-local:
	rm -rf $(GENGETOPT_CMDL)


RNAfold.c:    RNAfold_cmdl.h

RNAfold_cmdl.h  RNAfold_cmdl.c: RNAfold.ggo
	gengetopt -i RNAfold.ggo

RNAfold_SOURCES=RNAfold_cmdl.c RNAfold.c

RNAsubopt.c:    RNAsubopt_cmdl.h

RNAsubopt_cmdl.h  RNAsubopt_cmdl.c: RNAsubopt.ggo
	gengetopt -i RNAsubopt.ggo

RNAsubopt_SOURCES=RNAsubopt_cmdl.c RNAsubopt.c


RNAalifold.c:    RNAalifold_cmdl.h

RNAalifold_cmdl.h  RNAalifold_cmdl.c: RNAalifold.ggo
	gengetopt -i RNAalifold.ggo

RNAalifold_SOURCES=RNAalifold_cmdl.c RNAalifold.c


RNAcofold.c:    RNAcofold_cmdl.h

RNAcofold_cmdl.h  RNAcofold_cmdl.c: RNAcofold.ggo
	gengetopt -i RNAcofold.ggo

RNAcofold_SOURCES=RNAcofold_cmdl.c RNAcofold.c

RNAup.c:        RNAup_cmdl.h

RNAup_cmdl.h  RNAup_cmdl.c: RNAup.ggo
	gengetopt -i RNAup.ggo

RNAup_SOURCES=RNAup_cmdl.c RNAup.c


RNALfold.c:       RNALfold_cmdl.h

RNALfold_cmdl.h  RNALfold_cmdl.c: RNALfold.ggo
	gengetopt -i RNALfold.ggo

RNALfold_SOURCES=RNALfold_cmdl.c RNALfold.c

RNALalifold.c:      RNALalifold_cmdl.h

RNALalifold_cmdl.h  RNALalifold_cmdl.c: RNALalifold.ggo
	gengetopt -i RNALalifold.ggo

RNALalifold_SOURCES=RNALalifold_cmdl.c RNALalifold.c


RNAplfold.c:      RNAplfold_cmdl.h

RNAplfold_cmdl.h  RNAplfold_cmdl.c: RNAplfold.ggo
	gengetopt -i RNAplfold.ggo

RNAplfold_SOURCES=RNAplfold_cmdl.c RNAplfold.c


RNAduplex.c:      RNAduplex_cmdl.h

RNAduplex_cmdl.h  RNAduplex_cmdl.c: RNAduplex.ggo
	gengetopt -i RNAduplex.ggo

RNAduplex_SOURCES=RNAduplex_cmdl.c RNAduplex.c


RNAaliduplex.c:   RNAaliduplex_cmdl.h

RNAaliduplex_cmdl.h  RNAaliduplex_cmdl.c: RNAaliduplex.ggo
	gengetopt -i RNAaliduplex.ggo

RNAaliduplex_SOURCES=RNAaliduplex_cmdl.c RNAaliduplex.c

RNAeval.c:   RNAeval_cmdl.h

RNAeval_cmdl.h  RNAeval_cmdl.c: RNAeval.ggo
	gengetopt -i RNAeval.ggo

RNAeval_SOURCES=RNAeval_cmdl.c RNAeval.c

RNAheat.c:   RNAheat_cmdl.h

RNAheat_cmdl.h  RNAheat_cmdl.c: RNAheat.ggo
	gengetopt -i RNAheat.ggo

RNAheat_SOURCES=RNAheat_cmdl.c RNAheat.c

RNA2Dfold.c:    RNA2Dfold_cmdl.h

RNA2Dfold_cmdl.h  RNA2Dfold_cmdl.c: RNA2Dfold.ggo
	gengetopt -i RNA2Dfold.ggo

RNA2Dfold_SOURCES=RNA2Dfold_cmdl.c RNA2Dfold.c

RNAparconv.c:    RNAparconv_cmdl.h

RNAparconv_cmdl.h  RNAparconv_cmdl.c: RNAparconv.ggo
	gengetopt -i RNAparconv.ggo

RNAparconv_SOURCES=RNAparconv_cmdl.c RNAparconv.c

RNAinverse.c:   RNAinverse_cmdl.h

RNAinverse_cmdl.h  RNAinverse_cmdl.c: RNAinverse.ggo
	gengetopt -i RNAinverse.ggo

RNAinverse_SOURCES=RNAinverse_cmdl.c RNAinverse.c

RNApaln.c:   RNApaln_cmdl.h

RNApaln_cmdl.h  RNApaln_cmdl.c: RNApaln.ggo
	gengetopt -i RNApaln.ggo

RNApaln_SOURCES=RNApaln_cmdl.c RNApaln.c

RNAplot.c:   RNAplot_cmdl.h

RNAplot_cmdl.h  RNAplot_cmdl.c: RNAplot.ggo
	gengetopt -i RNAplot.ggo

RNAplot_SOURCES=RNAplot_cmdl.c RNAplot.c

RNAPKplex.c:    RNAPKplex_cmdl.h

RNAPKplex_cmdl.h  RNAPKplex_cmdl.c: RNAPKplex.ggo
	gengetopt -i RNAPKplex.ggo

RNAPKplex_SOURCES=RNAPKplex_cmdl.c RNAPKplex.c

RNAplex.c:    RNAplex_cmdl.h

RNAplex_cmdl.h  RNAplex_cmdl.c: RNAplex.ggo
	gengetopt -i RNAplex.ggo

RNAplex_SOURCES=RNAplex_cmdl.c RNAplex.c

RNAsnoop.c:    RNAsnoop_cmdl.h

RNAsnoop_cmdl.h  RNAsnoop_cmdl.c: RNAsnoop.ggo
	gengetopt -i RNAsnoop.ggo

RNAsnoop_SOURCES=RNAsnoop_cmdl.c RNAsnoop.c

RNAdistance.c:    RNAdistance_cmdl.h

RNAdistance_cmdl.h  RNAdistance_cmdl.c: RNAdistance.ggo
	gengetopt -i RNAdistance.ggo

RNAdistance_SOURCES=RNAdistance_cmdl.c RNAdistance.c

RNApdist.c:    RNApdist_cmdl.h

RNApdist_cmdl.h  RNApdist_cmdl.c: RNApdist.ggo
	gengetopt -i RNApdist.ggo

RNApdist_SOURCES=RNApdist_cmdl.c RNApdist.c
