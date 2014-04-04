lib_LIBRARIES = libRNA.a


if WITH_LIBSVM
  SVM_SRC     = ../@LIBSVM_DIR@/svm.cpp
  SVM_INC     = -I$(top_srcdir)/src/@LIBSVM_DIR@
  SVM_H       = ../@LIBSVM_DIR@/svm.h
  SVM_UTILS   = svm_utils.c
  SVM_UTILS_H = svm_utils.h
endif

pkginclude_HEADERS =  energy_par.h pair_mat.h RNAstruct.h dist_vars.h \
                      fold.h part_func.h stringdist.h edit_cost.h fold_vars.h \
                      profiledist.h treedist.h energy_const.h inverse.h PS_dot.h \
                      utils.h params.h subopt.h cofold.h duplex.h alifold.h aln_util.h \
                      part_func_co.h ProfileAln.h LPfold.h Lfold.h part_func_up.h ribo.h \
                      findpath.h loop_energies.h data_structures.h MEA.h mm.h \
                      naview.h 2Dfold.h 2Dpfold.h read_epars.h convert_epars.h PKplex.h\
                      plot_layouts.h ali_plex.h plex.h snofold.h snoop.h move_set.h\
                      gquad.h constraints.h eval.h perturbation_fold.h centroid.h\
                      structure_utils.h model.h ${SVM_UTILS_H}

libRNA_a_SOURCES  = dist_vars.c fold_vars.c part_func.c read_epars.c\
                    treedist.c energy_par.c inverse.c ProfileDist.c \
                    RNAstruct.c utils.c fold.c naview.c PS_dot.c \
                    stringdist.c params.c subopt.c list.c Lfold.c \
                    cofold.c part_func_co.c ProfileAln.c duplex.c \
                    alifold.c alipfold.c aln_util.c LPfold.c \
                    part_func_up.c ribo.c findpath.c convert_epars.c\
                    MEA.c aliLfold.c mm.c 2Dfold.c 2Dpfold.c plex_functions.c\
                    plot_layouts.c ali_plex.c c_plex.c plex.c snofold.c\
                    snoop.c move_set.c gquad.c eval.c constraints.c data_structures.c\
                    perturbation_fold.c centroid.c structure_utils.c model.c \
                    ${SVM_SRC} ${SVM_UTILS}

libRNA_a_CFLAGS = $(AM_CFLAGS) $(OPENMP_CFLAGS) @FPIC@


noinst_HEADERS  = intl11.h intl11dH.h intl21.h intl21dH.h \
                  intl22.h intl22dH.h list.h\
                  1.8.4_epars.h 1.8.4_intloops.h

EXTRA_DIST = $(pkginclude_HEADERS) circfold.inc alicircfold.inc model_avg.inc model_sd.inc ${SVM_H}

AM_CPPFLAGS = -I$(top_srcdir)/src ${SVM_INC} @FPIC@

