
SWIG_main_src = \
    $(srcdir)/../RNAxplorer.i

SWIG_misc_src = \
    $(srcdir)/../barrier_lower_bound.i \
    $(srcdir)/../distorted_sampling.i \
    $(srcdir)/../distorted_samplingMD.i \
    $(srcdir)/../meshpoint.i \
    $(srcdir)/../RNAwalk.i \
    $(srcdir)/../paths.i


SWIG_module_name = RNAxplorer

SWIG_wrapper = RNAxplorer_wrap.cpp

SWIG_src = \
    $(SWIG_main_src) \
    $(SWIG_misc_src)
