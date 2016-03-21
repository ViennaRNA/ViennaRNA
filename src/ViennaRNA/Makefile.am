AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = $(OPENMP_CFLAGS)
AM_CXXFLAGS = $(OPENMP_CFLAGS)
AM_CPPFLAGS =  -I$(top_srcdir)/src ${SVM_INC} ${JSON_INC} @FLOAT_PF_FLAG@ @DEPRECATION_WARNING@

# the library names to build
lib_LTLIBRARIES = libRNA.la

# convenience library
noinst_LTLIBRARIES = libRNA_conv.la

# where to install the headers on the system
libRNA_ladir = $(includedir)/ViennaRNA
libRNA_la_SOURCES = 

# ldflags to avoid version information
libRNA_la_LIBADD = libRNA_conv.la
libRNA_la_LDFLAGS = -avoid-version -static

# Dummy C++ source to cause C++ linking.
if WITH_LIBSVM
nodist_EXTRA_libRNA_la_SOURCES = dummy.cxx
else
nodist_EXTRA_libRNA_la_SOURCES = 
endif

if WITH_LIBSVM
  SVM_SRC     = ../@LIBSVM_DIR@/svm.cpp
  SVM_INC     = -I$(top_srcdir)/src/@LIBSVM_DIR@
  SVM_H       = ../@LIBSVM_DIR@/svm.h
  SVM_UTILS   = svm_utils.c
  SVM_UTILS_H = svm_utils.h
endif

if WITH_JSON
  JSON_SRC    = ../json/json.c
  JSON_H      = ../json/json.h
  JSON_INC    = -I$(top_srcdir)/src/json
endif

pkginclude_HEADERS =  energy_par.h pair_mat.h RNAstruct.h dist_vars.h mfe.h \
                      fold.h part_func.h stringdist.h edit_cost.h fold_vars.h \
                      profiledist.h treedist.h energy_const.h inverse.h PS_dot.h \
                      utils.h params.h subopt.h cofold.h duplex.h alifold.h aln_util.h \
                      part_func_co.h ProfileAln.h LPfold.h Lfold.h part_func_up.h ribo.h \
                      findpath.h data_structures.h MEA.h mm.h loop_energies.h \
                      exterior_loops.h hairpin_loops.h interior_loops.h multibranch_loops.h \
                      naview.h 2Dfold.h 2Dpfold.h read_epars.h convert_epars.h PKplex.h\
                      plot_layouts.h ali_plex.h plex.h snofold.h snoop.h move_set.h\
                      gquad.h constraints.h eval.h perturbation_fold.h centroid.h\
                      structure_utils.h string_utils.h model.h file_formats.h dp_matrices.h \
                      plot_structure.h plot_aln.h ligand.h boltzmann_sampling.h \
                      constraints_SHAPE.h constraints_hard.h constraints_soft.h \
                      ${SVM_UTILS_H} ${SVM_H} ${JSON_H}

libRNA_conv_la_SOURCES  = dist_vars.c data_structures.c part_func.c read_epars.c\
                    treedist.c energy_par.c inverse.c ProfileDist.c \
                    RNAstruct.c utils.c mfe.c fold.c naview.c PS_dot.c \
                    stringdist.c params.c subopt.c list.c Lfold.c \
                    cofold.c part_func_co.c ProfileAln.c duplex.c \
                    alifold.c alipfold.c aln_util.c LPfold.c string_utils.c \
                    part_func_up.c ribo.c findpath.c convert_epars.c\
                    MEA.c aliLfold.c mm.c 2Dfold.c 2Dpfold.c plex_functions.c \
                    plot_layouts.c ali_plex.c c_plex.c plex.c snofold.c \
                    snoop.c move_set.c gquad.c eval.c constraints.c ligand.c \
                    perturbation_fold.c centroid.c structure_utils.c model.c \
                    file_formats.c dp_matrices.c plot_structure.c plot_aln.c \
                    exterior_loops.c hairpin_loops.c interior_loops.c multibranch_loops.c \
                    boltzmann_sampling.c ${SVM_SRC} ${SVM_UTILS} ${JSON_SRC} \
                    constraints_SHAPE.c constraints_hard.c constraints_soft.c

# static library for subpackages
BUILT_SOURCES = libRNA.a

libRNA.a: libRNA.la
	$(LIBTOOL) --config > libtoolconfig.tmp
## ksh requires a path on the sourced file.
	. ./libtoolconfig.tmp; cp "$$objdir/libRNA.a" ./
	rm -f libtoolconfig.tmp

clean-local:
	-rm -f libRNA.a

noinst_HEADERS  = intl11.h intl11dH.h intl21.h intl21dH.h \
                  intl22.h intl22dH.h list.h\
                  1.8.4_epars.h 1.8.4_intloops.h ${SVM_H} ${JSON_H}

EXTRA_DIST = $(pkginclude_HEADERS) circfold.inc alicircfold.inc model_avg.inc model_sd.inc ${SVM_H} ${JSON_H}
