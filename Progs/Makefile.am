
bin_PROGRAMS = RNAfold RNAeval RNAheat RNApdist RNAdistance RNAinverse \
		RNAplot RNAsubopt RNALfold RNAcofold RNApaln RNAduplex \
		RNAalifold RNAplfold RNAup RNAaliduplex

LDADD = ../lib/libRNA.a -lm

INCLUDES = -I$(srcdir)/../H -I$(srcdir)/../lib

EXTRA_DIST =  RNAfold.ggo RNAfold_cmdl.c RNAfold_cmdl.h\
              RNAalifold.ggo RNAalifold_cmdl.c RNAalifold_cmdl.h\
              RNAcofold.ggo RNAcofold_cmdl.c  RNAcofold_cmdl.h\
              RNAup.ggo RNAup_cmdl.c  RNAup_cmdl.h\
              RNALfold.ggo RNALfold_cmdl.c RNALfold_cmdl.h

RNAfold.c:    RNAfold_cmdl.h

RNAfold_cmdl.h  RNAfold_cmdl.c: RNAfold.ggo
	gengetopt -i RNAfold.ggo

RNAfold_SOURCES=RNAfold_cmdl.c RNAfold.c

RNAalifold.c:    RNAalifold_cmdl.h

RNAalifold_cmdl.h  RNAalifold_cmdl.c: RNAalifold.ggo
	gengetopt -i RNAalifold.ggo

RNAalifold_SOURCES=RNAalifold_cmdl.c RNAalifold.c

RNAcofold.c:    RNAcofold_cmdl.h

RNAcofold_cmdl.h  RNAcofold_cmdl.c: RNAcofold.ggo
	gengetopt -i RNAcofold.ggo

RNAcofold_SOURCES=RNAcofold_cmdl.c RNAcofold.c

RNAup_cmdl.h  RNAup_cmdl.c: RNAup.ggo
	gengetopt -i RNAup.ggo

RNAup_SOURCES=RNAup_cmdl.c RNAup.c

RNALfold_cmdl.h  RNALfold_cmdl.c: RNALfold.ggo
	gengetopt -i RNALfold.ggo

RNALfold_SOURCES=RNALfold_cmdl.c RNALfold.c
