include ../generic.mk

INTERFACE_FILES = $(top_srcdir)/interfaces/RNAxplorer.i

pkgpyexecdir = $(py3execdir)/RNAxplorer
pkgpycachedir = $(pkgpyexecdir)/__pycache__

pkgpyexec_LTLIBRARIES = _RNAxplorer.la
pkgpyexec_DATA =  RNAxplorer/__init__.py
pkgpycache_DATA = RNAxplorer/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
                  RNAxplorer/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@

_RNAxplorer_la_SOURCES = $(INTERFACE_FILES) \
                         $(SWIG_wrapper)

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py

_RNAxplorer_la_CFLAGS = $(RNAxplorer_CFLAGS)
_RNAxplorer_la_CXXFLAGS = $(RNAxplorer_CXXFLAGS) -Wno-unused-label -Wno-unused-variable
_RNAxplorer_la_CPPFLAGS = -I$(PYTHON3_INC) $(RNAxplorer_CPPFLAGS) -I$(top_srcdir)/src

_RNAxplorer_la_LIBADD = \
    $(LDADD) \
    $(top_builddir)/src/libRNAxplorer_conv.la

_RNAxplorer_la_LDFLAGS  = $(RNALIB_LIBS) \
                          -avoid-version \
                          -module \
                          -shared \
                          -export-dynamic \
                          -shrext \
                          $(PYTHON3_SO) \
                          $(PYTHON3_LDFLAGS)


$(SWIG_wrapper): $(SWIG_src) $(INTERFACE_FILES)
	$(swig_verbose)$(SWIG) \
                  -I$(top_srcdir)/src \
                  -I$(srcdir) \
                  -I$(builddir) \
                  -DPY3 \
                  $(RNAxplorer_CPPFLAGS) \
                  -python \
                  -modern \
                  -py3 \
                  -c++ \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

# We rename RNAxplorer.py to __init__.py so that it can be installed
# in an "RNAxplorer" directory and be treated as a package.
RNAxplorer/__init__.py: $(SWIG_module_name).py
	$(AM_V_GEN)( test -d RNAxplorer || $(MKDIR_P) RNAxplorer ) && \
  cp `test -f RNAxplorer.py || echo '$(srcdir)/'`RNAxplorer.py RNAxplorer/__init__.py

# We "import _RNAxplorer" first so that if we fail to import the glue library
# we don't generate a broken .pyc or .pyo.
RNAxplorer/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc: RNAxplorer/__init__.py RNAxplorer/_RNAxplorer$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNAxplorer:$$PYTHONPATH" $(PYTHON3) -c "import _RNAxplorer" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -c "import RNAxplorer" )

RNAxplorer/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@: RNAxplorer/__init__.py RNAxplorer/_RNAxplorer$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNAxplorer:$$PYTHONPATH" $(PYTHON3) -O -c "import _RNAxplorer" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -O -c "import RNAxplorer" )

RNAxplorer/_RNAxplorer$(PYTHON3_SO): _RNAxplorer.la
	$(AM_V_GEN)( test -d RNAxplorer || $(MKDIR_P) RNAxplorer ) && \
  ( $(LIBTOOL) --config > libtoolconfig.tmp ) && \
  ( . ./libtoolconfig.tmp; cp $$objdir/_RNAxplorer$(PYTHON3_SO) RNAxplorer ) && \
  rm -f libtoolconfig.tmp

CLEANFILES = \
    RNAxplorer/_RNAxplorer$(PYTHON3_SO) \
    RNAxplorer/__init__.py \
    RNAxplorer/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
    RNAxplorer/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@\
    $(srcdir)/*.pyc

all-local:  RNAxplorer/_RNAxplorer$(PYTHON3_SO)

clean-local:
	-rm -rf RNAxplorer

# Remove the .la file - RNAxplorer.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pkgpyexecdir)/_RNAxplorer.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pkgpyexec_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pkgpyexecdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNAxplorer.py

distclean-local:
	-rm -f $(SWIG_wrapper) RNAxplorer.py

