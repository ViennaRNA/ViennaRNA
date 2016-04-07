
SWIG_main_src = $(srcdir)/../RNA.i

SWIG_tmaps = \
  $(srcdir)/../ptr2array.i

SWIG_misc_src = \
  $(srcdir)/../compare.i \
  $(srcdir)/../constraints.i \
  $(srcdir)/../constraints_hard.i \
  $(srcdir)/../constraints_soft.i \
  $(srcdir)/../constraints_SHAPE.i \
  $(srcdir)/../constraints_ligand.i \
  $(srcdir)/../eval.i \
  $(srcdir)/../fold_compound.i \
  $(srcdir)/../inverse.i \
  $(srcdir)/../mfe.i \
  $(srcdir)/../model_details.i \
  $(srcdir)/../params.i \
  $(srcdir)/../part_func.i \
  $(srcdir)/../plotting.i \
  $(srcdir)/../subopt.i \
  $(srcdir)/../utils.i \
  $(srcdir)/../file_formats.i

SWIG_module_name = RNA

SWIG_wrapper = RNA_wrap.cpp

SWIG_src = $(SWIG_main_src) $(SWIG_misc_src) $(SWIG_tmaps)
