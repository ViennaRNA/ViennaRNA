
bin_PROGRAMS = RNAfold RNAeval RNAheat RNApdist RNAdistance RNAinverse \
		RNAplot RNAsubopt RNALfold RNAcofold RNApaln RNAduplex \
		RNAalifold RNAplfold RNAup RNAaliduplex RNALalifold

LDADD = ../lib/libRNA.a -lm

if USE_OPENMP
  AM_CFLAGS = -DUSE_OPENMP -fopenmp
  AM_LDFLAGS = ../lib/libRNA.a -lm -fopenmp
endif


INCLUDES = -I$(srcdir)/../H -I$(srcdir)/../lib

EXTRA_DIST =  RNAfold.ggo RNAfold_cmdl.c RNAfold_cmdl.h\
              RNAalifold.ggo RNAalifold_cmdl.c RNAalifold_cmdl.h\
              RNAcofold.ggo RNAcofold_cmdl.c  RNAcofold_cmdl.h\
              RNAup.ggo RNAup_cmdl.c  RNAup_cmdl.h\
              RNALfold.ggo RNALfold_cmdl.c RNALfold_cmdl.h\
              RNALalifold.ggo RNALalifold_cmdl.c RNALalifold_cmdl.h\
              RNAplfold.ggo RNAplfold_cmdl.c RNAplfold_cmdl.h\
              RNAduplex.ggo RNAduplex_cmdl.c RNAduplex_cmdl.h\
              RNAaliduplex.ggo RNAaliduplex_cmdl.c RNAaliduplex_cmdl.h\
              RNAeval.ggo RNAeval_cmdl.c RNAeval_cmdl.h\
              RNAheat.ggo RNAheat_cmdl.c RNAheat_cmdl.h\
              RNAsubopt.ggo RNAsubopt_cmdl.c RNAsubopt_cmdl.h

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
