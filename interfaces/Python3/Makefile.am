include ../generic.mk

LANG_SPECIFIC_INTERFACE_FILES = \
  $(srcdir)/tmaps.i \
  $(srcdir)/file_py3.i \
  $(srcdir)/globals-md.i \
  $(srcdir)/callbacks-fc.i \
  $(srcdir)/callbacks-sc.i \
  $(srcdir)/callbacks-ud.i \
  $(srcdir)/callbacks-subopt.i \
  $(srcdir)/callbacks-mfe-window.i \
  $(builddir)/version.i

INTERFACE_FILES = $(SWIG_src) \
                  $(LANG_SPECIFIC_INTERFACE_FILES)

pkgpyexecdir = $(py3execdir)/RNA
pkgpycmoddir = $(py3execdir)
pkgpycachedir = $(pkgpyexecdir)/__pycache__

pkgpycmod_LTLIBRARIES = _RNA.la
pkgpyexec_DATA =  RNA/__init__.py
pkgpycache_DATA = RNA/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
                  RNA/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@

_RNA_la_SOURCES = $(INTERFACE_FILES) \
                  $(SWIG_wrapper)

BUILT_SOURCES = _RNA$(PYTHON3_SO)

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py \
  $(LANG_SPECIFIC_INTERFACE_FILES)

_RNA_la_CFLAGS = $(RNA_CFLAGS) $(PYTHON3_CFLAGS)
_RNA_la_CXXFLAGS = $(RNA_CXXFLAGS)
_RNA_la_CPPFLAGS = $(RNA_CPPFLAGS) -I$(top_srcdir)/src $(PYTHON3_INCLUDES)

_RNA_la_LIBADD = $(LDADD) $(top_builddir)/src/ViennaRNA/libRNA_conv.la

if VRNA_AM_SWITCH_SVM
_RNA_la_LIBADD += -lstdc++
_RNA_la_CPPFLAGS += -DVRNA_WITH_SVM
RNA_CPPFLAGS += -DVRNA_WITH_SVM
endif

_RNA_la_LDFLAGS = -avoid-version -module -shared -export-dynamic -shrext $(PYTHON3_SO) $(PYTHON3_LDFLAGS)

if VRNA_AM_SWITCH_LTO
_RNA_la_CFLAGS += $(LTO_CFLAGS)
_RNA_la_CXXFLAGS += $(LTO_CXXFLAGS)
_RNA_la_LDFLAGS += $(LTO_LDFLAGS)
endif


if HAS_SWIG

$(SWIG_wrapper): $(SWIG_src)
	$(swig_verbose)$(SWIG) -I$(top_srcdir)/src \
                  -I$(srcdir) -I$(builddir) \
                  -DPY3 \
                  $(_RNA_la_CPPFLAGS) \
                  -python -modern -py3 -c++ \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

endif

# We rename RNA.py to __init__.py so that it can be installed
# in an "RNA" directory and be treated as a package.
RNA/__init__.py: $(SWIG_module_name).py
	$(AM_V_GEN)( test -d RNA || $(MKDIR_P) RNA ) && \
  cp `test -f RNA.py || echo '$(srcdir)/'`RNA.py RNA/__init__.py

# We "import _RNA" first so that if we fail to import the glue library
# we don't generate a broken .pyc or .pyo.
RNA/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc: RNA/__init__.py _RNA$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNA:$$PYTHONPATH" $(PYTHON3) -c "import _RNA" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -c "import RNA" )

RNA/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@: RNA/__init__.py _RNA$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNA:$$PYTHONPATH" $(PYTHON3) -O -c "import _RNA" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -O -c "import RNA" )

_RNA$(PYTHON3_SO): _RNA.la
	$(AM_V_GEN)( $(LIBTOOL) --config > libtoolconfig.tmp ) && \
  ( . ./libtoolconfig.tmp; cp $$objdir/_RNA$(PYTHON3_SO) . ) && \
  rm -f libtoolconfig.tmp

CLEANFILES = \
    _RNA$(PYTHON3_SO) \
    RNA/__init__.py \
    RNA/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
    RNA/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@\
    $(srcdir)/*.pyc

clean-local:
	-rm -rf RNA _RNA$(PYTHON3_SO)


# Remove the .la file - RNA.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pkgpycmoddir)/_RNA.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pkgpycmod_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pkgpycmoddir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNA.py

distclean-local:
	-rm -f $(SWIG_wrapper) RNA.py

