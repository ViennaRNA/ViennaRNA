sphinx_verbose = $(sphinx_verbose_@AM_V@)
sphinx_verbose_ = $(sphinx_verbose_@AM_DEFAULT_V@)
sphinx_verbose_0 = @echo "  SPHINX    $@";

doxy2swig_verbose = $(doxy2swig_verbose_@AM_V@)
doxy2swig_verbose_ = $(doxy2swig_verbose_@AM_DEFAULT_V@)
doxy2swig_verbose_0 = @echo "  DOX2SWIG $@";

include ../generic.mk

if WITH_REFERENCE_MANUAL_XML
doxy2swig_pydoc_file = documentation.i
doxy2swig_pydoc_define = -DDOXY2SWIG
else
doxy2swig_pydoc_file = 
doxy2swig_pydoc_define = 
endif


LANG_SPECIFIC_INTERFACE_FILES = \
  $(srcdir)/tmaps.i \
  $(srcdir)/file_py3.i \
  $(srcdir)/globals-md.i \
  $(srcdir)/callbacks-fc.i \
  $(srcdir)/callbacks-sc.i \
  $(srcdir)/callbacks-ud.i \
  $(srcdir)/callbacks-subopt.i \
  $(srcdir)/callbacks-boltzmann-sampling.i \
  $(srcdir)/callbacks-melting.i \
  $(srcdir)/callbacks-landscape.i \
  $(srcdir)/callbacks-mfe-window.i \
  $(srcdir)/callbacks-pf-window.i \
  $(builddir)/version.i \
  $(doxy2swig_pydoc_file)


INTERFACE_FILES = $(SWIG_src) \
                  $(LANG_SPECIFIC_INTERFACE_FILES)

pkgpyexecdir = $(py3execdir)/RNA
pkgpycachedir = $(pkgpyexecdir)/__pycache__

pkgpyexec_LTLIBRARIES = _RNA.la
pkgpyexec_DATA =  RNA/__init__.py
pkgpycache_DATA = RNA/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
                  RNA/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@

_RNA_la_SOURCES = $(INTERFACE_FILES) \
                  $(SWIG_wrapper)

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(SWIG_module_name).py \
  $(LANG_SPECIFIC_INTERFACE_FILES)

_RNA_la_CFLAGS = $(RNA_CFLAGS) $(PTHREAD_CFLAGS)
_RNA_la_CXXFLAGS = $(RNA_CXXFLAGS) $(PTHREAD_CFLAGS) -Wno-unused-label -Wno-unused-variable
_RNA_la_CPPFLAGS = -I$(PYTHON3_INC) $(RNA_CPPFLAGS) -I$(top_srcdir)/src

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
                  $(PYTHON3_SO) \
                  $(PYTHON3_LDFLAGS)

if HAS_SWIG

if WITH_REFERENCE_MANUAL_XML

documentation.i: $(srcdir)/../classify.csv $(srcdir)/../typedefs.csv
	$(doxy2swig_verbose)test -d $(top_builddir)/doc/xml || ( test -d $(top_srcdir)/doc/xml && $(LN_S) $(top_srcdir)/doc/xml $(top_builddir)/doc/xml ) && \
	$(PYTHON3) $(top_srcdir)/doc/doxy2swig/doxy2swig.py \
      -b $(srcdir)/../classify.csv \
      -m $(srcdir)/../typedefs.csv \
      -p "vrna_" \
      -p "VRNA_" \
      --prefix-name="RNA" \
      -s "_t" \
      -s "_s" \
      -s "_e" \
      -s "_f" \
      --constructor-suffix="_t" \
      --constructor-suffix="_s" \
      -t \
      -a \
      -c \
      -w 100 \
      -q \
      -d $(PACKAGE_VERSION) \
      $(top_builddir)/doc/xml/index.xml \
      documentation.i

endif WITH_REFERENCE_MANUAL_XML

$(SWIG_wrapper): $(SWIG_src) $(INTERFACE_FILES)
	$(swig_verbose)$(SWIG) \
                  -I$(top_srcdir)/src \
                  -I$(srcdir) \
                  -I$(builddir) \
                  -DPY3 \
                  $(doxy2swig_pydoc_define) \
                  $(RNA_CPPFLAGS) \
                  -python \
                  -py3 \
                  -c++ \
                  -keyword \
                  -w511,509 \
                  -o $(SWIG_wrapper) \
                  $(SWIG_main_src)


$(SWIG_module_name).py: $(SWIG_wrapper)

endif HAS_SWIG

if WITH_PYTHON_DOC

##--------------------------------------------------##
## prepare Python documention                       ##
##--------------------------------------------------##

PYDOC_FILES_HTML = doc/build/* 

##--------------------------------------------------##
## prepare all rules in case reference manuals are  ##
## going to be (re)generated                        ##
##--------------------------------------------------##

doc/build/index.html: $(wildcard $(srcdir)/doc/source/*) RNA/__init__.py RNA/_RNA$(PYTHON3_SO)
	$(sphinx_verbose)$(SPHINXBUILD) -b html doc/source doc/build >sphinx.log 2>&1


$(PYDOC_FILES_HTML): doc/build/index.html

endif WITH_PYTHON_DOC

# We rename RNA.py to __init__.py so that it can be installed
# in an "RNA" directory and be treated as a package.
RNA/__init__.py: $(SWIG_module_name).py
	$(AM_V_GEN)( test -d RNA || $(MKDIR_P) RNA ) && \
  cp `test -f RNA.py || echo '$(srcdir)/'`RNA.py RNA/__init__.py

# We "import _RNA" first so that if we fail to import the glue library
# we don't generate a broken .pyc or .pyo.
RNA/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc: RNA/__init__.py RNA/_RNA$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNA:$$PYTHONPATH" $(PYTHON3) -c "import _RNA" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -c "import RNA" )

RNA/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@: RNA/__init__.py RNA/_RNA$(PYTHON3_SO)
	$(AM_V_GEN)( PYTHONPATH="RNA:$$PYTHONPATH" $(PYTHON3) -O -c "import _RNA" ; \
  PYTHONPATH=".:$$PYTHONPATH" $(PYTHON3) -O -c "import RNA" )

RNA/_RNA$(PYTHON3_SO): _RNA.la
	$(AM_V_GEN)( test -d RNA || $(MKDIR_P) RNA ) && \
  ( $(LIBTOOL) --config > libtoolconfig.tmp ) && \
  ( . ./libtoolconfig.tmp; cp $$objdir/_RNA$(PYTHON3_SO) RNA ) && \
  rm -f libtoolconfig.tmp


##--------------------------------------------------##
## Tell autoconf/automake to include the necessary  ##
## files in the distribution archive as well as in  ##
## the installation routine                         ##
##--------------------------------------------------##
html_DATA = $(PYDOC_FILES_HTML)

if WITH_PYTHON_DOC

vrna_htmlpydoc_inst: installdirs doc/build
	@$(NORMAL_INSTALL)
	test -d $(builddir)/doc/build || $(LN_S) $(srcdir)/doc/build $(builddir)/doc/build
	{ for i in $(builddir)/doc/build/*; do cp -r "$$i" "$(DESTDIR)$(pyhtmldir)/"; done }
	{ for i in "$(DESTDIR)$(pyhtmldir)"/*; do test -f "$$i" && chmod 644 "$$i"; done }
	chmod 755 "$(DESTDIR)$(pyhtmldir)"
	{ test ! -h $(builddir)/doc/build || rm -f $(builddir)/doc/build; }

vrna_htmlpydoc_uninst:
	@$(NORMAL_UNINSTALL)
	rm -rf "$(DESTDIR)$(pyhtmldir)"

else

vrna_htmlpydoc_inst:
vrna_htmlpydoc_uninst:

endif


install-htmlDATA: vrna_htmlpydoc_inst
uninstall-htmlDATA: vrna_htmlpydoc_uninst


CLEANFILES = \
    RNA/_RNA$(PYTHON3_SO) \
    RNA/__init__.py \
    RNA/__pycache__/__init__.@PYTHON3_CACHE_TAG@.pyc \
    RNA/__pycache__/__init__.@PYTHON3_CACHE_OPT1_EXT@\
    $(srcdir)/*.pyc

all-local:  RNA/_RNA$(PYTHON3_SO)

clean-local:
	-rm -rf RNA

# Remove the .la file - RNA.la is never linked against (it's a module)
# and Python doesn't use libltdl.  Note that the library gets installed by
# install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(pkgpyexecdir)/_RNA.la

# Because we don't install the .la file, "make uninstall" doesn't work and
# we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(pkgpyexec_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(pkgpyexecdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) RNA.py documentation.i version.i
	-rm -rf doc/build

distclean-local:
	-rm -f $(SWIG_wrapper) RNA.py documentation.i version.i
