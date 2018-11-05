AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = $(RNA_CFLAGS) $(SIMD_FLAGS) $(PTHREAD_CFLAGS)
AM_CXXFLAGS = $(RNA_CXXFLAGS) $(SIMD_FLAGS) $(PTHREAD_CFLAGS)
AM_CPPFLAGS = $(RNA_CPPFLAGS) ${SVM_INC} -I$(top_srcdir)/src ${JSON_INC}
AM_LDFLAGS = $(RNA_LDFLAGS) $(PTHREAD_LIBS)

SUBDIRS = static

vrna_constraintsdir     = $(pkgincludedir)/constraints
vrna_paramsdir          = $(pkgincludedir)/params
vrna_utilsdir           = $(pkgincludedir)/utils
vrna_iodir              = $(pkgincludedir)/io
vrna_searchdir          = $(pkgincludedir)/search
vrna_plottingdir        = $(pkgincludedir)/plotting
vrna_loopsdir           = $(pkgincludedir)/loops
vrna_datastructuresdir  = $(pkgincludedir)/datastructures

# the library names to build
lib_LTLIBRARIES = libRNA.la

# convenience library
noinst_LTLIBRARIES =  libRNA_conv.la \
                      libRNA_constraints.la \
                      libRNA_datastructures.la \
                      libRNA_utils.la \
                      libRNA_plotting.la \
                      libRNA_params.la \
                      libRNA_loops.la \
                      libRNA_eval.la \
                      libRNA_special_const.la

# where to install the headers on the system
libRNA_ladir = $(includedir)/ViennaRNA
libRNA_la_SOURCES = 

# ldflags to avoid version information
libRNA_la_LIBADD = libRNA_conv.la
libRNA_la_LDFLAGS = -avoid-version -static $(LTO_LDFLAGS) $(PTHREAD_LIBS)

libRNA_conv_la_LIBADD = \
    libRNA_constraints.la \
    libRNA_datastructures.la \
    libRNA_utils.la \
    libRNA_plotting.la \
    libRNA_params.la \
    libRNA_eval.la \
    libRNA_loops.la

# Dummy C++ source to cause C++ linking.
if VRNA_AM_SWITCH_SVM
nodist_EXTRA_libRNA_la_SOURCES = dummy.cxx
else
nodist_EXTRA_libRNA_la_SOURCES = 
endif

if VRNA_AM_SWITCH_SVM
  SVM_SRC         = ../@LIBSVM_DIR@/svm.cpp
  SVM_INC         = -I$(top_srcdir)/src/@LIBSVM_DIR@
  SVM_H           = ../@LIBSVM_DIR@/svm.h
  SVM_UTILS       = utils/svm_utils.c
  SVM_UTILS_H     = utils/svm.h
  SVM_UTILS_H_OLD = svm_utils.h
endif

if VRNA_AM_SWITCH_JSON
  JSON_SRC    = ../json/json.c
  JSON_H      = ../json/json.h
  JSON_INC    = -I$(top_srcdir)/src/json
endif

pkginclude_HEADERS =  \
    pair_mat.h \
    RNAstruct.h \
    dist_vars.h \
    mfe.h \
    mfe_window.h \
    fold.h \
    part_func.h \
    part_func_window.h \
    stringdist.h \
    edit_cost.h \
    fold_vars.h \
    profiledist.h \
    treedist.h \
    inverse.h \
    subopt.h \
    cofold.h \
    duplex.h \
    alifold.h \
    part_func_co.h \
    concentrations.h \
    ProfileAln.h \
    LPfold.h \
    Lfold.h \
    part_func_up.h \
    ribo.h \
    findpath.h \
    data_structures.h \
    char_stream.h \
    stream_output.h \
    fold_compound.h \
    MEA.h \
    mm.h \
    loop_energies.h \
    exterior_loops.h \
    hairpin_loops.h \
    interior_loops.h \
    multibranch_loops.h \
    2Dfold.h \
    2Dpfold.h \
    PKplex.h\
    ali_plex.h \
    plex.h \
    snofold.h \
    snoop.h \
    move_set.h\
    gquad.h \
    constraints.h \
    eval.h \
    perturbation_fold.h \
    centroid.h\
    utils.h \
    structure_utils.h \
    string_utils.h \
    aln_util.h \
    model.h \
    dp_matrices.h \
    PS_dot.h \
    plot_structure.h \
    plot_aln.h \
    plot_utils.h \
    plot_layouts.h \
    naview.h \
    boltzmann_sampling.h \
    equilibrium_probs.h \
    constraints_ligand.h \
    constraints_SHAPE.h \
    constraints_hard.h \
    constraints_soft.h \
    alphabet.h \
    sequence.h \
    grammar.h \
    structured_domains.h \
    unstructured_domains.h \
    file_utils.h \
    file_formats.h \
    file_formats_msa.h \
    energy_par.h \
    energy_const.h \
    params.h \
    read_epars.h \
    convert_epars.h \
    commands.h \
    units.h \
    combinatorics.h \
    neighbor.h \
    walk.h \
    ${SVM_UTILS_H_OLD} \
    ${SVM_H} \
    ${JSON_H}

vrna_constraints_HEADERS = \
    constraints/basic.h \
    constraints/ligand.h \
    constraints/SHAPE.h \
    constraints/hard.h \
    constraints/soft.h


vrna_loops_HEADERS = \
    loops/all.h \
    loops/external.h \
    loops/hairpin.h \
    loops/internal.h \
    loops/multibranch.h


vrna_utils_HEADERS = \
    utils/basic.h \
    utils/strings.h \
    utils/structures.h \
    utils/alignments.h \
    utils/higher_order_functions.h \
    utils/cpu.h \
    ${SVM_UTILS_H}


vrna_search_HEADERS = \
    search/BoyerMoore.h


vrna_plotting_HEADERS = \
    plotting/alignments.h \
    plotting/layouts.h \
    plotting/probabilities.h \
    plotting/structures.h \
    plotting/utils.h \
    plotting/naview.h


vrna_io_HEADERS = \
    io/utils.h \
    io/file_formats.h \
    io/file_formats_msa.h


vrna_params_HEADERS = \
    params/default.h \
    params/constants.h \
    params/basic.h \
    params/io.h \
    params/convert.h


vrna_datastructures_HEADERS = \
    datastructures/basic.h \
    datastructures/lists.h \
    datastructures/char_stream.h \
    datastructures/stream_output.h \
    datastructures/hash_tables.h


libRNA_conv_la_SOURCES = \
    fold_compound.c \
    dist_vars.c \
    part_func.c \
    part_func_wrappers.c \
    pf_fold.c \
    treedist.c \
    inverse.c \
    ProfileDist.c \
    RNAstruct.c \
    mfe.c \
    mfe_window.c \
    mfe_wrappers.c \
    mfe_window_wrappers.c \
    fold.c \
    stringdist.c \
    subopt.c \
    Lfold.c \
    cofold.c \
    part_func_co.c \
    concentrations.c \
    ProfileAln.c \
    duplex.c \
    alifold.c \
    alipfold.c \
    LPfold.c \
    part_func_up.c \
    ribo.c \
    findpath.c \
    MEA.c \
    mm.c \
    2Dfold.c \
    2Dpfold.c \
    plex_functions.c \
    ali_plex.c \
    c_plex.c \
    plex.c \
    snofold.c \
    snoop.c \
    move_set.c \
    gquad.c \
    perturbation_fold.c \
    centroid.c \
    model.c \
    dp_matrices.c \
    boltzmann_sampling.c \
    equilibrium_probs.c \
    neighbor.c \
    walk.c \
    ${SVM_SRC} \
    ${JSON_SRC} \
    alphabet.c \
    sequence.c \
    unstructured_domains.c \
    grammar.c

libRNA_eval_la_SOURCES = \
    eval.c \
    eval_wrappers.c

libRNA_constraints_la_SOURCES = \
    constraints/constraints.c \
    constraints/SHAPE.c \
    constraints/hard.c \
    constraints/soft.c \
    constraints/ligand.c

libRNA_utils_la_SOURCES = \
    utils/utils.c \
    utils/string_utils.c \
    utils/structure_utils.c \
    utils/structure_tree.c \
    utils/msa_utils.c \
    utils/higher_order_functions.c \
    utils/cpu.c \
    io/io_utils.c \
    io/file_formats.c \
    io/file_formats_msa.c \
    search/BoyerMoore.c \
    commands.c \
    units.c \
    combinatorics.c \
    ${SVM_UTILS}

libRNA_plotting_la_SOURCES = \
    plotting/alignments.c \
    plotting/layouts.c \
    plotting/probabilities.c \
    plotting/structures.c \
    plotting/plot_utils.c \
    plotting/naview.c

libRNA_loops_la_SOURCES = \
    loops/external.c \
    loops/external_bt.c \
    loops/external_pf.c \
    loops/hairpin.c \
    loops/hairpin_bt.c \
    loops/hairpin_pf.c \
    loops/internal.c \
    loops/internal_bt.c \
    loops/internal_pf.c \
    loops/multibranch.c \
    loops/multibranch_bt.c \
    loops/multibranch_pf.c

libRNA_params_la_SOURCES = \
    params/io.c \
    params/default.c \
    params/params.c \
    params/convert.c

libRNA_datastructures_la_SOURCES = \
    datastructures/basic_datastructures.c \
    datastructures/lists.c \
    datastructures/char_stream.c \
    datastructures/stream_output.c \
    datastructures/hash_tables.c

libRNA_special_const_la_SOURCES = \
    special_const.c

# static library for subpackages
all-local:      libRNA.a

libRNA.a: libRNA.la
	$(AM_V_GEN)( \
    $(LIBTOOL) --config > libtoolconfig.tmp && \
    ( \. ./libtoolconfig.tmp; cp "$$objdir/libRNA.a" ./ ) && \
    rm -f libtoolconfig.tmp \
  )

clean-local:
	-rm -f libRNA.a

noinst_HEADERS  = params/intl11.h \
                  params/intl11dH.h \
                  params/intl21.h \
                  params/intl21dH.h \
                  params/intl22.h \
                  params/intl22dH.h \
                  params/1.8.4_epars.h \
                  params/1.8.4_intloops.h \
                  list.h\
                  ${SVM_H} \
                  ${JSON_H}\
                  special_const.h

nodist_pkginclude_HEADERS = vrna_config.h

EXTRA_DIST =  $(pkginclude_HEADERS) \
              loops/external_hc.inc \
              loops/external_sc.inc \
              loops/external_sc_pf.inc \
              loops/hairpin_hc.inc \
              loops/hairpin_sc.inc \
              loops/hairpin_sc_pf.inc \
              loops/internal_hc.inc \
              loops/internal_sc.inc \
              loops/internal_sc_pf.inc \
              loops/multibranch_hc.inc \
              loops/multibranch_sc.inc \
              loops/multibranch_sc_pf.inc \
              params/svm_model_avg.inc \
              params/svm_model_sd.inc \
              ${SVM_H} \
              ${JSON_H} \
              color_output.inc \
              special_const.h
