AUTOMAKE_OPTIONS = subdir-objects

include $(top_srcdir)/silent_rules.mk

AM_CFLAGS = $(RNA_CFLAGS) $(PTHREAD_CFLAGS)
AM_CXXFLAGS = $(RNA_CXXFLAGS) $(PTHREAD_CFLAGS)
AM_CPPFLAGS = $(RNA_CPPFLAGS) -I$(top_srcdir)/src ${JSON_INC}
AM_LDFLAGS = $(RNA_LDFLAGS) $(PTHREAD_LIBS)

SUBDIRS = static @NAVIEW_DIR@
DIST_SUBDIRS = static @NAVIEW_DIR@

vrna_backtrackdir       = $(pkgincludedir)/backtrack
vrna_combinatoricsdir   = $(pkgincludedir)/combinatorics
vrna_constraintsdir     = $(pkgincludedir)/constraints
vrna_datastructuresdir  = $(pkgincludedir)/datastructures
vrna_evaldir            = $(pkgincludedir)/eval
vrna_grammardir         = $(pkgincludedir)/grammar
vrna_inversedir         = $(pkgincludedir)/inverse
vrna_iodir              = $(pkgincludedir)/io
vrna_landscapedir       = $(pkgincludedir)/landscape
vrna_loopsdir           = $(pkgincludedir)/loops
vrna_mfedir             = $(pkgincludedir)/mfe
vrna_paramsdir          = $(pkgincludedir)/params
vrna_partfuncdir        = $(pkgincludedir)/partfunc
vrna_plottingdir        = $(pkgincludedir)/plotting
vrna_probabilitiesdir   = $(pkgincludedir)/probabilities
vrna_probingdir         = $(pkgincludedir)/probing
vrna_puzzlerdir         = $(pkgincludedir)/plotting/RNApuzzler
vrna_samplingdir        = $(pkgincludedir)/sampling
vrna_searchdir          = $(pkgincludedir)/search
vrna_sequencesdir       = $(pkgincludedir)/sequences
vrna_structuresdir      = $(pkgincludedir)/structures
vrna_suboptdir          = $(pkgincludedir)/subopt
vrna_utilsdir           = $(pkgincludedir)/utils
vrna_zscoredir          = $(pkgincludedir)/zscore

# the library names to build
lib_LTLIBRARIES = libRNA.la

# convenience library
noinst_LTLIBRARIES =  \
    libRNA_backtrack.la \
    libRNA_combinatorics.la \
    libRNA_concentrations.la \
    libRNA_constraints.la \
    libRNA_conv.la \
    libRNA_datastructures.la \
    libRNA_eval.la \
    libRNA_grammar.la \
    libRNA_inverse.la \
    libRNA_landscape.la \
    libRNA_mfe.la \
    libRNA_params.la \
    libRNA_partfunc.la \
    libRNA_plotting.la \
    libRNA_probabilities.la \
    libRNA_probing.la \
    libRNA_sampling.la \
    libRNA_sequences.la \
    libRNA_special_const.la \
    libRNA_structures.la \
    libRNA_subopt.la \
    libRNA_utils.la

if VRNA_AM_SWITCH_SVM
noinst_LTLIBRARIES += libRNA_svm.la
endif

# where to install the headers on the system
libRNA_ladir = $(includedir)/ViennaRNA
libRNA_la_SOURCES = 

# ldflags to avoid version information
libRNA_la_LIBADD = libRNA_conv.la
libRNA_la_LDFLAGS = -avoid-version -static $(LTO_LDFLAGS) $(PTHREAD_LIBS)

libRNA_concentrations_la_CPPFLAGS = $(AM_CPPFLAGS) $(DLIB_CPPFLAGS)

libRNA_conv_la_LIBADD = \
    libRNA_backtrack.la \
    libRNA_combinatorics.la \
    libRNA_concentrations.la \
    libRNA_constraints.la \
    libRNA_datastructures.la \
    libRNA_eval.la \
    libRNA_grammar.la \
    libRNA_inverse.la \
    libRNA_landscape.la \
    libRNA_mfe.la \
    libRNA_params.la \
    libRNA_partfunc.la \
    libRNA_plotting.la \
    libRNA_probabilities.la \
    libRNA_probing.la \
    libRNA_sampling.la \
    libRNA_sequences.la \
    libRNA_structures.la \
    libRNA_subopt.la \
    libRNA_utils.la

if VRNA_AM_SWITCH_NAVIEW
  libRNA_plotting_la_LIBADD = @NAVIEW_DIR@/libRNA_naview.la
endif


if VRNA_AM_SWITCH_SIMD_SSE41
noinst_LTLIBRARIES += libRNA_utils_sse41.la
libRNA_conv_la_LIBADD += libRNA_utils_sse41.la
libRNA_utils_sse41_la_CFLAGS = $(SIMD_SSE41_FLAGS)
endif

if VRNA_AM_SWITCH_SIMD_AVX512
noinst_LTLIBRARIES += libRNA_utils_avx512.la
libRNA_conv_la_LIBADD += libRNA_utils_avx512.la
libRNA_utils_avx512_la_CFLAGS = $(SIMD_AVX512_FLAGS)
endif

# Dummy C++ source to cause C++ linking.
if VRNA_AM_SWITCH_SVM
nodist_EXTRA_libRNA_la_SOURCES = dummy.cxx
else
nodist_EXTRA_libRNA_la_SOURCES = 
endif

if VRNA_AM_SWITCH_SVM
  libRNA_conv_la_LIBADD += libRNA_svm.la

  SVM_SRC         = ../@LIBSVM_DIR@/svm.cpp
  SVM_INC         = -I$(top_srcdir)/src/@LIBSVM_DIR@
  SVM_H           = ../@LIBSVM_DIR@/svm.h
  SVM_UTILS       = utils/svm_utils.c
  SVM_UTILS_H     = utils/svm.h
  SVM_UTILS_H_OLD = svm_utils.h

  libRNA_svm_la_CPPFLAGS = $(CPPFLAGS)  -I$(top_srcdir)/src ${SVM_INC} -U_OPENMP ## disable OpenMP support for libSVM
  libRNA_utils_la_CPPFLAGS = $(AM_CPPFLAGS) ${SVM_INC}
endif

JSON_SRC    = ../json/json.c
JSON_H      = ../json/json.h
JSON_INC    = -I$(top_srcdir)/src/json


CEPHES_SRC  = ../cephes/kn.c \
							../cephes/const.c \
							../cephes/mtherr.c \
							../cephes/expn.c
CEPHES_H    = ../cephes/mconf.h


DEPRECATED_INCLUDE_FILES = \
    aln_util.h \
    alphabet.h \
    boltzmann_sampling.h \
    centroid.h\
    char_stream.h \
    commands.h \
    combinatorics.h \
    constraints.h \
    constraints_hard.h \
    constraints_ligand.h \
    constraints_SHAPE.h \
    constraints_soft.h \
    constraints/SHAPE.h \
    convert_epars.h \
    data_structures.h \
    energy_const.h \
    energy_par.h \
    equilibrium_probs.h \
    eval.h \
    exterior_loops.h \
    file_formats.h \
    file_formats_msa.h \
    file_utils.h \
    findpath.h \
    grammar.h \
    gquad.h \
    hairpin_loops.h \
    interior_loops.h \
    inverse.h \
    Lfold.h \
    loop_energies.h \
    LPfold.h \
    MEA.h \
    mfe.h \
    mfe_window.h \
    multibranch_loops.h \
    neighbor.h \
    pair_mat.h \
    params.h \
    part_func.h \
    part_func_window.h \
    pf_multifold.h \
    PKplex.h\
    plot_aln.h \
    plot_layouts.h \
    plot_structure.h \
    plot_utils.h \
    PS_dot.h \
    read_epars.h \
    sequence.h \
    stream_output.h \
    string_utils.h \
    structure_utils.h \
    subopt.h \
    subopt_zuker.h \
    units.h \
    utils.h \
    walk.h


if VRNA_AM_SWITCH_NAVIEW
DEPRECATED_INCLUDE_FILES += naview.h 
endif

pkginclude_HEADERS =  \
    2Dfold.h \
    2Dpfold.h \
    alifold.h \
    ali_plex.h \
    cofold.h \
    concentrations.h \
    dist_vars.h \
    dp_matrices.h \
    duplex.h \
    edit_cost.h \
    fold_compound.h \
    fold.h \
    fold_vars.h \
    heat_capacity.h \
    mm.h \
    model.h \
    move_set.h\
    part_func_co.h \
    part_func_up.h \
    perturbation_fold.h \
    pk_plex.h\
    plex.h \
    ProfileAln.h \
    profiledist.h \
    ribo.h \
    RNAstruct.h \
    snofold.h \
    snoop.h \
    stringdist.h \
    structured_domains.h \
    treedist.h \
    unstructured_domains.h \
    wrap_dlib.h \
    ${JSON_H} \
    ${CEPHES_H} \
    ${DEPRECATED_INCLUDE_FILES}


if VRNA_AM_SWITCH_SVM
pkginclude_HEADERS += \
    ${SVM_UTILS_H_OLD} \
    ${SVM_H}
endif

vrna_constraints_HEADERS = \
    constraints/basic.h \
    constraints/ligand.h \
    constraints/hard.h \
    constraints/SHAPE.h \
    constraints/soft.h \
    constraints/soft_special.h


vrna_probing_HEADERS = \
    probing/basic.h \
    probing/SHAPE.h


vrna_loops_HEADERS = \
    loops/all.h \
    loops/external.h \
    loops/hairpin.h \
    loops/internal.h \
    loops/multibranch.h


vrna_eval_HEADERS = \
    eval/basic.h \
    eval/structures.h \
    eval/exterior.h \
    eval/hairpin.h \
    eval/internal.h \
    eval/multibranch.h \
    eval/gquad.h


vrna_mfe_HEADERS = \
    mfe/global.h \
    mfe/local.h \
    mfe/exterior.h \
    mfe/internal.h \
    mfe/multibranch.h \
    mfe/gquad.h


vrna_backtrack_HEADERS = \
    backtrack/global.h \
    backtrack/exterior.h \
    backtrack/hairpin.h \
    backtrack/internal.h \
    backtrack/multibranch.h \
    backtrack/gquad.h


vrna_partfunc_HEADERS =   \
    partfunc/global.h     \
    partfunc/local.h      \
    partfunc/multifold.h  \
    partfunc/exterior.h   \
    partfunc/internal.h   \
    partfunc/multibranch.h \
    partfunc/gquad.h


vrna_probabilities_HEADERS =   \
    probabilities/basepairs.h     \
    probabilities/structures.h


vrna_subopt_HEADERS = \
    subopt/basic.h    \
    subopt/wuchty.h   \
    subopt/zuker.h \
    subopt/gquad.h


vrna_sampling_HEADERS = \
    sampling/basic.h

vrna_combinatorics_HEADERS = \
    combinatorics/basic.h


vrna_inverse_HEADERS = \
    inverse/basic.h


vrna_grammar_HEADERS = \
    grammar/basic.h \
    grammar/mfe.h \
    grammar/partfunc.h


vrna_utils_HEADERS = \
    utils/basic.h \
    utils/strings.h \
    utils/structures.h \
    utils/alignments.h \
    utils/higher_order_functions.h \
    utils/cpu.h \
    utils/units.h \
    utils/log.h


if VRNA_AM_SWITCH_SVM
vrna_utils_HEADERS += \
    ${SVM_UTILS_H}


vrna_zscore_HEADERS = \
    zscore/basic.h
endif


vrna_search_HEADERS = \
    search/BoyerMoore.h


vrna_plotting_HEADERS = \
    plotting/alignments.h \
    plotting/layouts.h \
    plotting/probabilities.h \
    plotting/structures.h \
    plotting/utils.h


vrna_puzzler_HEADERS = \
    plotting/RNApuzzler/RNApuzzler.h \
    plotting/RNApuzzler/RNAturtle.h


vrna_io_HEADERS = \
    io/utils.h \
    io/file_formats.h \
    io/file_formats_msa.h \
    io/commands.h


vrna_params_HEADERS = \
    params/default.h \
    params/constants.h \
    params/basic.h \
    params/io.h \
    params/convert.h \
    params/salt.h \
    params/ribosum.h


vrna_structures_HEADERS = \
    structures/benchmark.h \
    structures/centroid.h \
    structures/dotbracket.h \
    structures/helix.h \
    structures/mea.h \
    structures/metrics.h \
    structures/pairtable.h \
    structures/problist.h \
    structures/shapes.h \
    structures/tree.h \
    structures/utils.h

vrna_sequences_HEADERS = \
    sequences/alphabet.h \
    sequences/sequence.h \
    sequences/alignments.h \
    sequences/utils.h


vrna_datastructures_HEADERS = \
    datastructures/array.h \
    datastructures/basic.h \
    datastructures/lists.h \
    datastructures/char_stream.h \
    datastructures/stream_output.h \
    datastructures/string.h \
    datastructures/hash_tables.h \
    datastructures/heap.h \
    datastructures/sparse_mx.h \
    datastructures/dp_matrices.h


vrna_landscape_HEADERS = \
    landscape/paths.h \
    landscape/findpath.h \
    landscape/neighbor.h \
    landscape/walk.h \
    landscape/move.h


libRNA_conv_la_SOURCES = \
    fold_compound.c \
    dist_vars.c \
    treedist.c \
    ProfileDist.c \
    RNAstruct.c \
    stringdist.c \
    ProfileAln.c \
    duplex.c \
    mm.c \
    2Dfold.c \
    2Dpfold.c \
    plex_functions.c \
    ali_plex.c \
    c_plex.c \
    plex.c \
    snofold.c \
    snoop.c \
    perturbation_fold.c \
    model.c \
    ${JSON_SRC} \
    unstructured_domains.c \
    heat_capacity.c \
    ${CEPHES_SRC}


if VRNA_AM_SWITCH_SVM
libRNA_svm_la_SOURCES = \
    ${SVM_SRC}
endif


libRNA_eval_la_SOURCES = \
    eval/eval_exterior.c \
    eval/eval_hairpin.c \
    eval/eval_internal.c \
    eval/eval_multibranch.c \
    eval/eval_gquad.c \
    eval/exp_eval_exterior.c \
    eval/exp_eval_hairpin.c \
    eval/exp_eval_internal.c \
    eval/exp_eval_multibranch.c \
    eval/eval_structures.c \
    eval/eval_covariance.c \
    eval/eval_wrappers.c


libRNA_concentrations_la_SOURCES = \
    concentrations.c \
    ../@DLIB_DIR@/dlib/all/source.cpp \
    wrap_dlib.cpp


libRNA_constraints_la_SOURCES = \
    constraints/constraints.c \
    constraints/hard.c \
    constraints/soft.c \
    constraints/soft_cb_multi.c \
    constraints/sc_cb_mod.c \
    constraints/sc_cb_mod_parser.c \
    constraints/sc_cb_mod_wrappers.c \
    constraints/ligand.c


libRNA_probing_la_SOURCES = \
    probing/probing.c \
    probing/SHAPE.c


libRNA_structures_la_SOURCES = \
    structures/structure_benchmark.c \
    structures/centroid.c \
    structures/structure_dotbracket.c \
    structures/structure_helix.c \
    structures/mea.c \
    structures/structure_metrics.c \
    structures/structure_pairtable.c \
    structures/structure_problist.c \
    structures/structure_shapes.c \
    structures/structure_tree.c \
    structures/structure_utils.c


libRNA_sequences_la_SOURCES = \
    sequences/alphabet.c \
    sequences/sequence.c \
    sequences/msa.c \
    sequences/seq_utils.c


libRNA_utils_la_SOURCES = \
    utils/utils.c \
    utils/string_utils.c \
    utils/higher_order_functions.c \
    utils/cpu.c \
    utils/units.c \
    utils/log.c \
    io/io_utils.c \
    io/file_formats.c \
    io/file_formats_msa.c \
    search/BoyerMoore.c \
    io/commands.c


if VRNA_AM_SWITCH_SVM
libRNA_utils_la_SOURCES += \
    zscore/zscore.c \
    ${SVM_UTILS}
endif


if VRNA_AM_SWITCH_SIMD_SSE41
libRNA_utils_sse41_la_SOURCES = \
    utils/higher_order_functions_sse41.c
endif

if VRNA_AM_SWITCH_SIMD_AVX512
libRNA_utils_avx512_la_SOURCES = \
    utils/higher_order_functions_avx512.c
endif

libRNA_plotting_la_SOURCES = \
    plotting/alignments.c \
    plotting/layouts.c \
    plotting/probabilities.c \
    plotting/structures.c \
    plotting/plot_utils.c \
    plotting/svg.c \
    plotting/eps.c \
    plotting/gml.c \
    plotting/xrna.c \
    plotting/ssv.c \
    plotting/RNApuzzler/RNApuzzler.c \
    plotting/RNApuzzler/RNAturtle.c


libRNA_mfe_la_SOURCES = \
    mfe/mfe.c \
    mfe/mfe_window.c \
    mfe/mfe_wrappers.c \
    mfe/mfe_window_wrappers.c \
    mfe/mfe_exterior.c \
    mfe/mfe_exterior_window.c \
    mfe/mfe_internal.c \
    mfe/mfe_multibranch.c \
    mfe/mfe_gquad.c \
    mfe/alifold.c \
    mfe/cofold.c \
    mfe/Lfold.c \
    mfe/fold.c


libRNA_backtrack_la_SOURCES = \
    backtrack/bt_exterior.c \
    backtrack/bt_exterior_f5.c \
    backtrack/bt_exterior_f3.c \
    backtrack/bt_hairpin.c \
    backtrack/bt_internal.c \
    backtrack/bt_multibranch.c \
    backtrack/bt_gquad.c


libRNA_partfunc_la_SOURCES = \
    partfunc/partfunc.c \
    partfunc/pf_window.c \
    partfunc/pf_wrappers.c \
    partfunc/pf_multifold.c \
    partfunc/pf_exterior.c \
    partfunc/pf_internal.c \
    partfunc/pf_multibranch.c \
    partfunc/pf_gquad.c \
    partfunc/alipfold.c \
    partfunc/part_func_up.c \
    partfunc/part_func_co.c \
    partfunc/pf_fold.c


libRNA_probabilities_la_SOURCES = \
    probabilities/equilibrium_probs.c \
    probabilities/probs_structures.c


libRNA_subopt_la_SOURCES = \
    subopt/subopt.c \
    subopt/subopt_zuker.c \
    subopt/subopt_gquad.c


libRNA_sampling_la_SOURCES = \
    sampling/boltzmann_sampling.c \
    sampling/bs_wrappers.c


libRNA_params_la_SOURCES = \
    params/io.c \
    params/default.c \
    params/params.c \
    params/convert.c \
    params/salt.c \
    params/ribosum.c


libRNA_inverse_la_SOURCES = \
    inverse/inverse.c


libRNA_combinatorics_la_SOURCES = \
    combinatorics/combinatorics.c


libRNA_grammar_la_SOURCES = \
    grammar/grammar.c \
    grammar/gr_extension_mfe.c \
    grammar/gr_extension_pf.c


libRNA_datastructures_la_SOURCES = \
    datastructures/array.c \
    datastructures/basic_datastructures.c \
    datastructures/lists.c \
    datastructures/char_stream.c \
    datastructures/stream_output.c \
    datastructures/string.c \
    datastructures/hash_tables.c \
    datastructures/heap.c \
    datastructures/sparse_mx.c \
    datastructures/dp_matrices.c


libRNA_landscape_la_SOURCES = \
    move_set.c \
    landscape/move.c \
    landscape/findpath.c \
    landscape/neighbor.c \
    landscape/local_neighbors.c \
    landscape/walk.c

libRNA_special_const_la_SOURCES = \
    params/special_const.c \
    io/sanitize.c

RNAPUZZLER_H = \
    plotting/RNApuzzler/headers/AABB_struct.h \
    plotting/RNApuzzler/headers/boundingBoxes_struct.h \
    plotting/RNApuzzler/headers/configArc_struct.h \
    plotting/RNApuzzler/headers/config_struct.h \
    plotting/RNApuzzler/headers/configtree_struct.h \
    plotting/RNApuzzler/headers/tBaseInformation_struct.h

RNAPUZZLER_INC = \
    plotting/RNApuzzler/includes/coordinates.inc \
    plotting/RNApuzzler/includes/definitions.inc \
    plotting/RNApuzzler/includes/vector_math.inc \
    plotting/RNApuzzler/includes/boundingBoxes.inc \
    plotting/RNApuzzler/includes/configtree.inc \
    plotting/RNApuzzler/includes/drawingconfig.inc \
    plotting/RNApuzzler/includes/intersectLevelBoundingBoxes.inc \
    plotting/RNApuzzler/includes/intersectLevelLines.inc \
    plotting/RNApuzzler/includes/intersectLevelTreeNodes.inc \
    plotting/RNApuzzler/includes/postscriptArcs.inc \
    plotting/RNApuzzler/includes/boundingWedge.inc \
    plotting/RNApuzzler/includes/calcDeltas.inc \
    plotting/RNApuzzler/includes/handleAncestorIntersections.inc \
    plotting/RNApuzzler/includes/handleConfigChanges.inc \
    plotting/RNApuzzler/includes/handleSiblingIntersections.inc \
    plotting/RNApuzzler/includes/intersectionType.inc \
    plotting/RNApuzzler/includes/optimize.inc \
    plotting/RNApuzzler/includes/resolveExteriorChildIntersections.inc \
    plotting/RNApuzzler/includes/resolveIntersections.inc \
    plotting/RNApuzzler/includes/rotationAngle.inc \
    plotting/RNApuzzler/includes/svgArcs.inc

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
                  constraints/sc_cb_intern.h \
                  ${RNAPUZZLER_H} \
                  list.h \
                  ${SVM_H} \
                  ${JSON_H} \
                  params/special_const.h \
                  ${CEPHES_H} \
                  io/sanitize.h

nodist_pkginclude_HEADERS = vrna_config.h

EXTRA_DIST =  $(pkginclude_HEADERS) \
              constraints/hc_depot.inc \
              constraints/exterior_hc.inc \
              constraints/exterior_sc.inc \
              constraints/exterior_sc_pf.inc \
              constraints/hairpin_hc.inc \
              constraints/hairpin_sc.inc \
              constraints/hairpin_sc_pf.inc \
              constraints/internal_hc.inc \
              constraints/internal_sc.inc \
              constraints/internal_sc_pf.inc \
              constraints/multibranch_hc.inc \
              constraints/multibranch_sc.inc \
              constraints/multibranch_sc_pf.inc \
              params/svm_model_avg.inc \
              params/svm_model_sd.inc \
              sampling/data_structures_nonred.inc \
              plotting/ps_helpers.inc \
              ${RNAPUZZLER_INC} \
              landscape/local_neighbors.inc \
              ${SVM_H} \
              ${JSON_H} \
              intern/color_output.h \
              intern/gquad_helpers.h \
              intern/grammar_dat.h \
              intern/unistd_win.h \
              params/special_const.h \
              io/sanitize.h

if VRNA_AM_SWITCH_SVM
EXTRA_DIST += intern/zscore_dat.h
endif
