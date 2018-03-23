AUTOMAKE_OPTIONS = subdir-objects

AM_CFLAGS = $(RNA_CFLAGS) $(SIMD_CFLAGS)
AM_CXXFLAGS = $(RNA_CXXFLAGS)
AM_CPPFLAGS = $(RNA_CPPFLAGS) -I$(top_srcdir)/src ${SVM_INC} ${JSON_INC}
AM_LDFLAGS = $(RNA_LDFLAGS)

# the library names to build
lib_LTLIBRARIES = libRNA.la

# convenience library
noinst_LTLIBRARIES =  libRNA_conv.la \
                      libRNA_constraints.la \
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
libRNA_la_LDFLAGS = -avoid-version -static $(LTO_LDFLAGS)

libRNA_conv_la_LIBADD = \
    libRNA_constraints.la \
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
  SVM_SRC     = ../@LIBSVM_DIR@/svm.cpp
  SVM_INC     = -I$(top_srcdir)/src/@LIBSVM_DIR@
  SVM_H       = ../@LIBSVM_DIR@/svm.h
  SVM_UTILS   = svm_utils.c
  SVM_UTILS_H = svm_utils.h
endif

if VRNA_AM_SWITCH_JSON
  JSON_SRC    = ../json/json.c
  JSON_H      = ../json/json.h
  JSON_INC    = -I$(top_srcdir)/src/json
endif

pkginclude_HEADERS =  \
    energy_par.h \
    pair_mat.h \
    RNAstruct.h \
    dist_vars.h \
    mfe.h \
    fold.h \
    part_func.h \
    stringdist.h \
    edit_cost.h \
    fold_vars.h \
    profiledist.h \
    treedist.h \
    energy_const.h \
    inverse.h \
    PS_dot.h \
    utils.h \
    params.h \
    subopt.h \
    cofold.h \
    duplex.h \
    alifold.h \
    aln_util.h \
    part_func_co.h \
    concentrations.h \
    ProfileAln.h \
    LPfold.h \
    Lfold.h \
    part_func_up.h \
    ribo.h \
    findpath.h \
    data_structures.h \
    fold_compound.h \
    MEA.h \
    mm.h \
    loop_energies.h \
    exterior_loops.h \
    hairpin_loops.h \
    interior_loops.h \
    multibranch_loops.h \
    naview.h \
    2Dfold.h \
    2Dpfold.h \
    read_epars.h \
    convert_epars.h \
    PKplex.h\
    plot_layouts.h \
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
    structure_utils.h \
    string_utils.h \
    model.h \
    file_formats.h \
    dp_matrices.h \
    plot_structure.h \
    plot_aln.h \
    plot_utils.h \
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
    file_formats_msa.h \
    file_utils.h \
    char_stream.h \
    commands.h \
    units.h \
    combinatorics.h \
    neighbor.h \
    walk.h \
    ${SVM_UTILS_H} \
    ${SVM_H} \
    ${JSON_H}

libRNA_conv_la_SOURCES = \
    fold_compound.c \
    dist_vars.c \
    data_structures.c \
    part_func.c \
    pf_fold.c \
    treedist.c \
    inverse.c \
    ProfileDist.c \
    RNAstruct.c \
    mfe.c \
    fold.c \
    stringdist.c \
    subopt.c \
    list.c \
    mfe_window.c \
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
    ${SVM_SRC} ${SVM_UTILS} ${JSON_SRC} \
    alphabet.c \
    sequence.c \
    unstructured_domains.c

libRNA_eval_la_SOURCES = \
    eval.c \
    eval_wrappers.c

libRNA_constraints_la_SOURCES = \
    constraints.c \
    constraints_SHAPE.c \
    constraints_hard.c \
    constraints_soft.c \
    constraints_ligand.c

libRNA_utils_la_SOURCES = \
    utils.c \
    aln_util.c \
    string_utils.c \
    structure_utils.c \
    file_utils.c \
    file_formats.c \
    file_formats_msa.c \
    char_stream.c \
    commands.c \
    units.c \
    combinatorics.c

libRNA_plotting_la_SOURCES = \
    plot_utils.c \
    plot_aln.c \
    plot_layouts.c \
    plot_structure.c \
    naview.c \
    PS_dot.c

libRNA_loops_la_SOURCES = \
    exterior_loops.c \
    exterior_loops_bt.c \
    exterior_loops_pf.c \
    hairpin_loops.c \
    hairpin_loops_bt.c \
    hairpin_loops_pf.c \
    interior_loops.c \
    interior_loops_bt.c \
    interior_loops_pf.c \
    multibranch_loops.c \
    multibranch_loops_bt.c \
    multibranch_loops_pf.c

libRNA_params_la_SOURCES = \
    read_epars.c \
    energy_par.c \
    params.c \
    convert_epars.c

libRNA_special_const_la_SOURCES = \
    special_const.c

# static library for subpackages
BUILT_SOURCES = libRNA.a

libRNA.a: libRNA.la
	$(AM_V_GEN)( \
    $(LIBTOOL) --config > libtoolconfig.tmp && \
    ( \. ./libtoolconfig.tmp; cp "$$objdir/libRNA.a" ./ ) && \
    rm -f libtoolconfig.tmp \
  )

clean-local:
	-rm -f libRNA.a

noinst_HEADERS  = intl11.h intl11dH.h intl21.h intl21dH.h \
                  intl22.h intl22dH.h list.h\
                  1.8.4_epars.h 1.8.4_intloops.h ${SVM_H} ${JSON_H}\
                  special_const.h

nodist_pkginclude_HEADERS = vrna_config.h

EXTRA_DIST =  $(pkginclude_HEADERS) \
              exterior_loops_hc.inc \
              exterior_loops_sc.inc \
              exterior_loops_sc_pf.inc \
              hairpin_loops.inc \
              interior_loops.inc \
              multibranch_loops.inc \
              model_avg.inc \
              model_sd.inc \
              ${SVM_H} \
              ${JSON_H} \
              color_output.inc \
              special_const.h
