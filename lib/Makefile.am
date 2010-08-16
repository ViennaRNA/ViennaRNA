lib_LIBRARIES = libRNA.a

libRNA_a_SOURCES  = dist_vars.c fold_vars.c part_func.c read_epars.c\
                    treedist.c energy_par.c inverse.c ProfileDist.c \
                    RNAstruct.c utils.c fold.c naview.c PS_dot.c \
                    stringdist.c params.c subopt.c list.c Lfold.c \
                    cofold.c part_func_co.c ProfileAln.c duplex.c \
                    alifold.c alipfold.c aln_util.c LPfold.c \
                    part_func_up.c ribo.c findpath.c loop_energies.c \
                    MEA.c aliLfold.c mm.c 2Dfold.c 2Dpfold.c

libRNA_a_CFLAGS = $(AM_CFLAGS) @FPIC@

if USE_OPENMP
  AM_CFLAGS = -DUSE_OPENMP -fopenmp
  AM_LDFLAGS = ../lib/libRNA.a -lm -fopenmp
endif

noinst_HEADERS = intloops.h list.h

EXTRA_DIST = circfold.inc alicircfold.inc

INCLUDES = -I$(srcdir)/../H
