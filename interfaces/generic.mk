swig_verbose = $(swig_verbose_@AM_V@)
swig_verbose_ = $(swig_verbose_@AM_DEFAULT_V@)
swig_verbose_0 = @echo "  SWIG     $@";

SWIG_main_src = $(srcdir)/../RNA.i

SWIG_tmaps = \
  $(srcdir)/../ptr2array.i

SWIG_misc_src = \
  $(srcdir)/../aln_utils.i \
  $(srcdir)/../basic_algorithms.i \
  $(srcdir)/../boltzmann_sampling.i \
  $(srcdir)/../combinatorics.i \
  $(srcdir)/../compare.i \
  $(srcdir)/../commands.i \
  $(srcdir)/../constraints.i \
  $(srcdir)/../constraints_hard.i \
  $(srcdir)/../constraints_soft.i \
  $(srcdir)/../constraints_SHAPE.i \
  $(srcdir)/../constraints_ligand.i \
  $(srcdir)/../data_structures.i \
  $(srcdir)/../duplex.i \
  $(srcdir)/../eval.i \
  $(srcdir)/../file_formats.i \
  $(srcdir)/../fold_compound.i \
  $(srcdir)/../grammar.i \
  $(srcdir)/../inverse.i \
  $(srcdir)/../loops.i \
  $(srcdir)/../mfe.i \
  $(srcdir)/../mfe_window.i \
  $(srcdir)/../model_details.i \
  $(srcdir)/../move.i \
  $(srcdir)/../neighbor.i \
  $(srcdir)/../params.i \
  $(srcdir)/../part_func.i \
  $(srcdir)/../paths.i \
  $(srcdir)/../pf_window.i \
  $(srcdir)/../plotting.i \
  $(srcdir)/../sequence.i \
  $(srcdir)/../subopt.i \
  $(srcdir)/../structure_utils.i \
  $(srcdir)/../utils.i \
  $(srcdir)/../walk.i

SWIG_module_name = RNA

SWIG_wrapper = RNA_wrap.cpp

SWIG_src = $(SWIG_main_src) $(SWIG_misc_src) $(SWIG_tmaps)
