

if WITH_LIBSVM
  LfoldZ_prog = RNALfoldz
  LfoldZ_cmdl = RNALfoldz.ggo RNALfoldz_cmdl.c RNALfoldz_cmdl.h
  RNALfold_LINK = ${CXXLD} -o $@
  RNALfoldz_LINK = ${CXXLD} -o $@
else
  RNALfold_LINK = ${CCLD} -o $@
endif


bin_PROGRAMS = RNAfold RNAeval RNAheat RNApdist RNAdistance RNAinverse \
    RNAplot RNAsubopt RNALfold RNAcofold RNApaln RNAduplex \
    RNAalifold RNAplfold RNAup RNAaliduplex RNALalifold RNA2Dfold \
    RNAparconv RNAPKplex RNAplex RNAsnoop ${LfoldZ_prog}

LDADD = $(OPENMP_CFLAGS) ../lib/libRNA.a -lm

AM_CFLAGS = $(OPENMP_CFLAGS)

INCLUDES = -I$(srcdir)/../H -I$(srcdir)/../lib

EXTRA_DIST =  RNAfold.ggo RNAfold_cmdl.c RNAfold_cmdl.h\
              RNAalifold.ggo RNAalifold_cmdl.c RNAalifold_cmdl.h\
              RNAcofold.ggo RNAcofold_cmdl.c  RNAcofold_cmdl.h\
              RNAup.ggo RNAup_cmdl.c  RNAup_cmdl.h\
              RNALfold.ggo RNALfold_cmdl.c RNALfold_cmdl.h ${LfoldZ_cmdl}\
              RNALalifold.ggo RNALalifold_cmdl.c RNALalifold_cmdl.h\
              RNAplfold.ggo RNAplfold_cmdl.c RNAplfold_cmdl.h\
              RNAduplex.ggo RNAduplex_cmdl.c RNAduplex_cmdl.h\
              RNAaliduplex.ggo RNAaliduplex_cmdl.c RNAaliduplex_cmdl.h\
              RNAeval.ggo RNAeval_cmdl.c RNAeval_cmdl.h\
              RNAheat.ggo RNAheat_cmdl.c RNAheat_cmdl.h\
              RNAsubopt.ggo RNAsubopt_cmdl.c RNAsubopt_cmdl.h\
              RNA2Dfold.ggo RNA2Dfold_cmdl.c RNA2Dfold_cmdl.h\
              RNAparconv.ggo RNAparconv_cmdl.c RNAparconv_cmdl.h\
              RNAinverse.ggo RNAinverse_cmdl.c RNAinverse_cmdl.h\
              RNApaln.ggo RNApaln_cmdl.c RNApaln_cmdl.h\
              RNAplot.ggo RNAplot_cmdl.c RNAplot_cmdl.h\
              RNAPKplex.ggo RNAPKplex_cmdl.c RNAPKplex_cmdl.h\
              RNAplex.ggo RNAplex_cmdl.c RNAplex_cmdl.h\
              RNAsnoop.ggo RNAsnoop_cmdl.c RNAsnoop_cmdl.h

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

#if WITH_LIBSVM
RNALfoldz.c:       RNALfoldz_cmdl.h

RNALfoldz_cmdl.h  RNALfoldz_cmdl.c: RNALfoldz.ggo
	gengetopt -i RNALfoldz.ggo

RNALfoldz_SOURCES=RNALfoldz_cmdl.c RNALfoldz.c
#endif

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
