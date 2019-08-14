include ../generic.mk

LANG_SPECIFIC_INTERFACE_FILES = \
  $(srcdir)/tmaps.i \
  $(srcdir)/globals-md.i \
  $(srcdir)/callbacks-fc.i \
  $(srcdir)/callbacks-sc.i \
  $(srcdir)/callbacks-ud.i \
  $(srcdir)/callbacks-subopt.i \
  $(srcdir)/callbacks-boltzmann-sampling.i \
  $(srcdir)/callbacks-mfe-window.i \
  $(srcdir)/callbacks-pf-window.i \
  $(builddir)/version.i

INTERFACE_FILES = $(SWIG_src) \
                  $(LANG_SPECIFIC_INTERFACE_FILES)

pkgpy2execdir = $(py2execdir)/RNA

pkgpy2exec_LTLIBRARIES = _RNA.la
pkgpy2exec_DATA = RNA/__init__.py

_RNA_la_SOURCES = $(INTERFACE_FILES) \
                  $(SWIG_wrapper)

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py \
  $(LANG_SPECIFIC_INTERFACE_FILES)

_RNA_la_CFLAGS = $(RNA_CFLAGS) $(PTHREAD_CFLAGS)
_RNA_la_CXXFLAGS = $(RNA_CXXFLAGS) $(PTHREAD_CFLAGS) -Wno-unused-label -Wno-unused-variable
_RNA_la_CPPFLAGS = -I$(PYTHON2_INC) $(RNA_CPPFLAGS) -I$(top_srcdir)/src

_RNA_la_LIBADD = \
    $(LDADD) \
    $(top_builddir)/src/ViennaRNA/libRNA_conv.la

if VRNA_AM_SWITCH_SVM
_RNA_la_LIBADD += -lstdc++
_RNA_la_CPPFLAGS += -DVRNA_WITH_SVM
RNA_CPPFLAGS += -DVRNA_WITH_SVM
endif

if VRNA_AM_SWITCH_MPFR
_RNA_la_LIBADD += $(MPFR_LIBS)
endif

_RNA_la_LDFLAGS = -Xcompiler $(PTHREAD_CFLAGS) \
                  $(PTHREAD_LIBS) \
                  -avoid-version \
                  -module \
                  -shared \
                  -export-dynamic \
                  -shrext \
                  $(PYTHON2_SO) \
                  $(PYTHON2_LDFLAGS) 

##
## SWIG interface wrapper code generation
##

if HAS_SWIG

$(SWIG_wrapper): $(SWIG_src) $(INTERFACE_FILES)
	$(swig_verbose)$(SWIG) -I$(top_srcdir)/src \
                  -I$(srcdir) -I$(builddir) \
                  $(RNA_CPPFLAGS) \
                  -python -c++ \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

endif

# We rename RNA.py to __init__.py so that it can be installed
# in an "RNA" directory and be treated as a package.
RNA/__init__.py: $(SWIG_module_name).py
	$(AM_V_GEN)( test -d RNA || $(MKDIR_P) RNA ) && \
  cp `test -f RNA.py || echo '$(srcdir)/'`RNA.py RNA/__init__.py

RNA/_RNA$(PYTHON2_SO): _RNA.la
	$(AM_V_GEN)( test -d RNA || $(MKDIR_P) RNA ) && \
  ( $(LIBTOOL) --config > libtoolconfig.tmp ) && \
  ( . ./libtoolconfig.tmp; cp $$objdir/_RNA$(PYTHON2_SO) RNA ) && \
  rm -f libtoolconfig.tmp

CLEANFILES = \
    RNA/_RNA$(PYTHON2_SO) \
    RNA/__init__.py \
    RNA/__init__.pyc

all-local:  RNA/_RNA$(PYTHON2_SO)

clean-local:
	-rm -rf RNA

# Remove the .la file - RNA.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pkgpy2execdir)/_RNA.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pkgpy2exec_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pkgpy2execdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNA.py

distclean-local:
	-rm -f $(SWIG_wrapper) RNA.py

