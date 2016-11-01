AM_CPPFLAGS = -I$(top_srcdir)/src

if WITH_LIBRNA_API3
AM_CFLAGS = @VRNA_CFLAGS@
LDADD = @VRNA_LIBS@
else
AM_CFLAGS = @VRNA2_CFLAGS@
LDADD = @VRNA2_LIBS@
endif

bin_PROGRAMS = Kinfold
SUBDIRS = Example

Kinfold_SOURCES = baum.c cache.c globals.c main.c nachbar.c \
		  baum.h cache_util.h globals.h   nachbar.h \
		  cmdline.c cmdline.h


EXTRA_DIST = Laplace/extract_data.pl \
             Laplace/laplace.sh \
             Laplace/README \
             Laplace/seq1 \
             Laplace/seq1_ref.ps \
             Laplace/to_boxplot.R \
             kinfold.ggo

globals.o: cmdline.h

cmdline.h cmdline.c: kinfold.ggo
	gengetopt -i kinfold.ggo

dist_man_MANS = Kinfold.1
