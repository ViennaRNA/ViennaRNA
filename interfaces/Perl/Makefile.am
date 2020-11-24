include ../generic.mk


PERL_MODULE = $(SWIG_module_name).pm

perlextdir = $(prefix)$(PERL_ARCH_RELATIVE_INSTALL_DIR)/auto/RNA
perlext_LTLIBRARIES = RNA.la

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

RNA_la_SOURCES =  $(INTERFACE_FILES) \
                  $(SWIG_wrapper)

perlmoddir = $(prefix)$(PERL_LIB_RELATIVE_INSTALL_DIR)
nobase_dist_perlmod_DATA = \
  $(PERL_MODULE) \
  RNA/Benchmark.pm \
  RNA/Design.pm \
  RNA/Files.pm \
  RNA/Reliability.pm \
  RNA/Utils.pm \
  RNA/Params.pm

RNA_la_CFLAGS = $(RNA_CFLAGS) $(PTHREAD_CFLAGS)
RNA_la_CXXFLAGS = $(RNA_CXXFLAGS) $(PTHREAD_CFLAGS) -Wno-unused-label -Wno-unused-variable
RNA_la_CPPFLAGS = -I$(PERL_EXT_INC) $(RNA_CPPFLAGS) -I$(top_srcdir)/src

RNA_la_LIBADD = \
    $(LDADD) \
    $(top_builddir)/src/ViennaRNA/libRNA_conv.la

if VRNA_AM_SWITCH_SVM
RNA_la_LIBADD += -lstdc++
RNA_CPPFLAGS += -DVRNA_WITH_SVM
RNA_la_CPPFLAGS += -DVRNA_WITH_SVM
endif

if VRNA_AM_SWITCH_MPFR
RNA_la_LIBADD += $(MPFR_LIBS)
endif

RNA_la_LDFLAGS =

RNA_la_LDFLAGS += -Xcompiler $(PTHREAD_CFLAGS) \
                  $(PTHREAD_LIBS) \
                  -avoid-version \
                  -module \
                  -shared \
                  -export-dynamic \
                  -shrext \
                  $(PERL_SO)

SOURCES = $(srcdir)/RNA.pod

EXTRA_DIST = \
  $(SWIG_wrapper) \
  RNA.pod \
  $(PERL_MODULE) \
  RNA/t/Design.t \
  $(LANG_SPECIFIC_INTERFACE_FILES)

##
## SWIG interface wrapper code generation
##

if HAS_SWIG

$(SWIG_wrapper) $(PERL_MODULE): $(INTERFACE_FILES)
	$(swig_verbose)$(SWIG) -I$(top_srcdir)/src \
                  -I$(srcdir) -I$(builddir) \
                  $(RNA_CPPFLAGS) \
                  -perl5 -proxy -const -c++ \
                  -module RNA -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

endif

##
## copy shared library for perl-tests
##

auto/RNA/RNA$(PERL_SO): RNA.la
	$(AM_V_GEN)( $(MKDIR_P) auto/RNA ) && \
  ( \
    $(LIBTOOL) --config > libtoolconfig.tmp && \
    ( . ./libtoolconfig.tmp; cp "$$objdir/RNA$(PERL_SO)" auto/RNA ) && \
    rm -f libtoolconfig.tmp \
  )

all-local:  auto/RNA/RNA$(PERL_SO)

clean-local:
	-rm -rf auto

## Remove the .la file - RNA.la is never linked against (it's a module)
## and Perl doesn't use libltdl.  Note that the library gets installed by
## install-data, so that's where we need to hook.
install-data-hook:
	rm -f $(DESTDIR)$(perlextdir)/RNA.la

## Because we don't install the .la file, "make uninstall" doesn't work and
## we need to remove the file ourselves.
uninstall-local:
	eval `grep '^dlname=' $(perlext_LTLIBRARIES)` ; \
	rm -f $(DESTDIR)$(perlextdir)/"$$dlname"

## remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) $(PERL_MODULE)

distclean-local:
	-rm -f $(SWIG_wrapper) $(PERL_MODULE)
