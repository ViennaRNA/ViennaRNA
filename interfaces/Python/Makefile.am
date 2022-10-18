include ../generic.mk

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py

INTERFACE_FILES = $(top_srcdir)/interfaces/RNAxplorer.i

pkgpy3execdir = $(py3execdir)/RNAxplorer
pkgpy3cmoddir = $(py3execdir)

pkgpy3cmod_LTLIBRARIES = _RNAxplorer.la
pkgpy3exec_DATA = RNAxplorer/__init__.py

_RNAxplorer_la_SOURCES = $(INTERFACE_FILES) \
                         $(SWIG_wrapper)

BUILT_SOURCES = _RNAxplorer$(PYTHON3_SO)


_RNAxplorer_la_CFLAGS = $(RNAxplorer_CFLAGS) $(PYTHON3_CFLAGS)
_RNAxplorer_la_CXXFLAGS = $(RNAxplorer_CXXFLAGS)
_RNAxplorer_la_CPPFLAGS = -I$(PYTHON3_INC) $(PYTHON3_CPPFLAGS) -I$(top_srcdir)/src

_RNAxplorer_la_LIBADD = $(LDADD) $(top_builddir)/src/libRNAxplorer_conv.la $(PYTHON3_LDFLAGS)

_RNAxplorer_la_LDFLAGS = $(RNALIB_LIBS) -avoid-version -module -shared -export-dynamic -shrext $(PYTHON3_SO)


$(SWIG_wrapper): $(SWIG_src)
	$(swig_verbose)$(SWIG) -I$(top_srcdir)/src \
                  -I$(srcdir) -I$(builddir) \
                  -python -modern -py3 -c++ \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

# We rename RNAxplorer.py to __init__.py so that it can be installed
# in an "RNAxplorer" directory and be treated as a package.
RNAxplorer/__init__.py: $(SWIG_module_name).py
	-test -d RNAxplorer || $(MKDIR_P) RNAxplorer
	cp `test -f RNAxplorer.py || echo '$(srcdir)/'`RNAxplorer.py RNAxplorer/__init__.py

_RNAxplorer$(PYTHON3_SO): _RNAxplorer.la
	$(AM_V_GEN)( $(LIBTOOL) --config > libtoolconfig.tmp ) && \
  ( . ./libtoolconfig.tmp; cp $$objdir/_RNAxplorer$(PYTHON3_SO) . ) && \
  rm -f libtoolconfig.tmp

clean-local:
	-rm -rf RNAxplorer _RNAxplorer$(PYTHON3_SO)

# Remove the .la file - RNAxplorer.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pkgpy3cmoddir)/_RNAxplorer.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pkgpy3cmod_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pkgpy3cmoddir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNAxplorer.py

distclean-local:
	-rm -f $(SWIG_wrapper) RNAxplorer.py

