bin_PROGRAMS = RNAsubopt popt

RNAsubopt_SOURCES = RNAsubopt.c list.c subopt.c

RNAsubopt_LDADD = fold.o -L../lib -lRNA -lm

noinst_HEADERS = list.h

LDADD = -L../lib -lRNA -lm

INCLUDES = -I../H

fold.o:	../lib/fold.c
	$(COMPILE) -I../lib -DSUBOPT -c ../lib/fold.c