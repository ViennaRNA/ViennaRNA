lib_LIBRARIES = libRNA.a


if WITH_LIBSVM
  SVM_SRC     = ../@LIBSVM_DIR@/svm.cpp
  SVM_INC     = -I${srcdir}/../@LIBSVM_DIR@
  SVM_H       = ../@LIBSVM_DIR@/svm.h
  SVM_UTILS   = svm_utils.c
endif

libRNA_a_SOURCES  = dist_vars.c fold_vars.c part_func.c read_epars.c\
                    treedist.c energy_par.c inverse.c ProfileDist.c \
                    RNAstruct.c utils.c fold.c naview.c PS_dot.c \
                    stringdist.c params.c subopt.c list.c Lfold.c \
                    cofold.c part_func_co.c ProfileAln.c duplex.c \
                    alifold.c alipfold.c aln_util.c LPfold.c \
                    part_func_up.c ribo.c findpath.c convert_epars.c\
                    MEA.c aliLfold.c mm.c 2Dfold.c 2Dpfold.c plex_functions.c\
                    plot_layouts.c ali_plex.c c_plex.c plex.c snofold.c\
                    snoop.c move_set.c gquad.c ${SVM_SRC} ${SVM_UTILS}

libRNA_a_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS) @FPIC@


noinst_HEADERS  = intl11.h intl11dH.h intl21.h intl21dH.h \
                  intl22.h intl22dH.h list.h\
                  1.8.4_epars.h 1.8.4_intloops.h

EXTRA_DIST = circfold.inc alicircfold.inc model_avg.inc model_sd.inc ${SVM_H}

AM_CPPFLAGS = -I$(srcdir)/../H ${SVM_INC} @FPIC@

