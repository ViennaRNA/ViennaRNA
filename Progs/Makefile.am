
bin_PROGRAMS = RNAfold RNAeval RNAheat RNApdist RNAdistance RNAinverse \
		RNAplot RNAsubopt RNALfold RNAcofold RNApaln RNAduplex \
		RNAalifold RNAplfold RNAup RNAaliduplex

LDADD = ../lib/libRNA.a -lm

INCLUDES = -I$(srcdir)/../H -I$(srcdir)/../lib

EXTRA_DIST=RNAfold.ggo RNAfold_cmdl.c RNAfold_cmdl.h

RNAfold.c:    RNAfold_cmdl.h

RNAfold_cmdl.h  RNAfold_cmdl.c: RNAfold.ggo
	gengetopt -i RNAfold.ggo

RNAfold_SOURCES=RNAfold_cmdl.c RNAfold.c
