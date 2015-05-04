
SWIG_main_src = $(srcdir)/../RNA.i

SWIG_tmaps = \
  $(srcdir)/../tmaps.i \
  $(srcdir)/../ptr2array.i

SWIG_misc_src = \
  $(srcdir)/../essentials.i \
  $(srcdir)/../plotting.i \
  $(srcdir)/../utils.i

SWIG_module_name = RNA

SWIG_wrapper = RNA_wrap.c
