lib_LIBRARIES = libRNA.a


if WITH_LIBSVM
  SVM_SOURCE = ../@libsvm_dir@/svm.cpp
endif

libRNA_a_SOURCES  = dist_vars.c fold_vars.c part_func.c read_epars.c\
                    treedist.c energy_par.c inverse.c ProfileDist.c \
                    RNAstruct.c utils.c fold.c naview.c PS_dot.c \
                    stringdist.c params.c subopt.c list.c Lfold.c \
                    cofold.c part_func_co.c ProfileAln.c duplex.c \
                    alifold.c alipfold.c aln_util.c LPfold.c \
                    part_func_up.c ribo.c findpath.c convert_epars.c\
                    MEA.c aliLfold.c mm.c 2Dfold.c 2Dpfold.c ${SVM_SOURCE}

libRNA_a_CFLAGS = $(AM_CFLAGS) @FPIC@

if USE_OPENMP
  AM_CFLAGS = -DUSE_OPENMP -fopenmp
  AM_LDFLAGS = ../lib/libRNA.a -lm -fopenmp
endif

noinst_HEADERS  = intl11.h intl11dH.h intl21.h intl21dH.h \
                  intl22.h intl22dH.h list.h\
                  1.8.4_epars.h 1.8.4_intloops.h

EXTRA_DIST = circfold.inc alicircfold.inc model_avg.inc model_sd.inc

INCLUDES = -I$(srcdir)/../H

if WITH_LIBSVM
  libRNA_INCLUDES = -I@libsvm_include@
endif
