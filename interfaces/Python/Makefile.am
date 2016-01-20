include ../generic.mk

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py


pyexec_LTLIBRARIES = _RNA.la
pkgpylibdir = $(pyexecdir)/RNA
pkgpylib_DATA = RNA/__init__.py RNA/__init__.pyc

_RNA_la_SOURCES = $(top_srcdir)/interfaces/RNA.i $(SWIG_wrapper)

_RNA_la_CPPFLAGS = -I$(top_srcdir)/src $(PYTHON_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) @FLOAT_PF_FLAG@ @DEPRECATION_WARNING@
_RNA_la_LIBADD = $(LDADD) $(top_builddir)/src/ViennaRNA/libRNA_conv.la $(PYTHON_LIBS) $(PYTHON_EXTRA_LIBS) -lstdc++
_RNA_la_LDFLAGS = -avoid-version -module -shared -export-dynamic $(PYTHON_EXTRA_LDFLAGS) -Wc,$(LTO_LDFLAGS)

$(SWIG_wrapper): $(SWIG_src)
	$(AM_V_GEN)$(SWIG) -I$(top_srcdir)/src @FLOAT_PF_FLAG@ \
                  $(AX_SWIG_PYTHON_OPT) \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

$(SWIG_module_name).py: $(SWIG_wrapper)

# We rename RNA.py to __init__.py so that it can be installed
# in an "RNA" directory and be treated as a package.
RNA/__init__.py: $(SWIG_module_name).py
	-test -d RNA || mkdir RNA
##	$(MKDIR_P) RNA
	cp `test -f RNA.py || echo '$(srcdir)/'`RNA.py RNA/__init__.py

# We "import _RNA" first so that if we fail to import the glue library
# we don't generate a broken .pyc or .pyo.
RNA/__init__.pyc: RNA/__init__.py _RNA.la
	PYTHONPATH=".libs:$$PYTHONPATH" $(PYTHON) -c "import _RNA"
	PYTHONPATH="RNA:.libs:$$PYTHONPATH" $(PYTHON) -c "import RNA"

clean-local:
	-rm -rf RNA

# Remove the .la file - RNA.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pyexecdir)/_RNA.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pyexec_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pyexecdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNA.py

distclean-local:
	-rm -f $(SWIG_wrapper) RNA.py

