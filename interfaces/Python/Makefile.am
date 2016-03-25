include ../generic.mk

AM_CFLAGS += $(PYTHON2_CFLAGS)

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py


pkgpy2execdir = $(py2execdir)/RNA

pkgpy2exec_LTLIBRARIES = _RNA.la
pkgpy2exec_DATA = RNA/__init__.py

_RNA_la_SOURCES = $(top_srcdir)/interfaces/RNA.i $(SWIG_wrapper)

_RNA_la_CPPFLAGS = -I$(top_srcdir)/src $(PYTHON2_INCLUDES) $(SWIG_PYTHON_CPPFLAGS) @FLOAT_PF_FLAG@ @DEPRECATION_WARNING@
_RNA_la_LIBADD = $(LDADD) $(top_builddir)/src/ViennaRNA/libRNA_conv.la $(PYTHON2_LIBS)
if WITH_LIBSVM
_RNA_la_LIBADD += -lstdc++
endif

_RNA_la_LDFLAGS = -avoid-version -module -shared -export-dynamic -shrext "$(PYTHON_SO)"
if WITH_LTO
_RNA_la_LDFLAGS += -Wc,$(LTO_LDFLAGS)
endif

$(SWIG_wrapper): $(SWIG_src)
	$(AM_V_GEN)$(SWIG) -I$(top_srcdir)/src @FLOAT_PF_FLAG@ \
                  -python -modern  -c++ \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

# We rename RNA.py to __init__.py so that it can be installed
# in an "RNA" directory and be treated as a package.
RNA/__init__.py: $(SWIG_module_name).py
	-test -d RNA || $(MKDIR_P) RNA
	cp `test -f RNA.py || echo '$(srcdir)/'`RNA.py RNA/__init__.py

RNA/_RNA$(PYTHON2_SO): _RNA.la
	$(MKDIR_P) RNA
	$(LIBTOOL) --config > libtoolconfig.tmp
## ksh requires a path on the sourced file.
	. ./libtoolconfig.tmp; cp $$objdir/_RNA$(PYTHON2_SO) RNA
	rm -f libtoolconfig.tmp

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

