include ../generic.mk

PERL_MODULE = $(SWIG_module_name).pm

perlextdir = $(PERL_EXT_LIB)/auto/RNA
perlext_LTLIBRARIES = RNA.la

RNA_la_SOURCES = $(SWIG_src) RNA_wrap.c

perlmoddir = $(PERL_EXT_LIB)
nobase_dist_perlmod_DATA = \
  $(PERL_MODULE) \
  RNA/Benchmark.pm \
  RNA/Design.pm \
  RNA/Files.pm \
  RNA/Reliability.pm \
  RNA/Utils.pm \
  test.pl \
  RNA/t/Design.t

RNA_la_CPPFLAGS = -I$(top_srcdir)/src $(PERL_EXT_CPPFLAGS) $(SWIG_PERL_CPPFLAGS) $(PERLXS_CFLAGS) @FLOAT_PF_FLAG@ @DEPRECATION_WARNING@
RNA_la_LIBADD = $(LDADD) $(top_builddir)/src/ViennaRNA/libRNA_conv.la $(PERLXS_LIBS) -lstdc++
RNA_la_LDFLAGS = -avoid-version -module -shared -export-dynamic $(PERLXS_LDFLAGS) -Wc,$(LTO_LDFLAGS)

SOURCES = RNA.pod

BUILT_SOURCES = auto/RNA/RNA.so

EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(PERL_MODULE) \
  RNAfold.pl \
  RNAfold.cgi

##
## SWIG interface wrapper code generation
##

$(SWIG_wrapper) $(PERL_MODULE): $(SWIG_src)
	$(AM_V_GEN)$(SWIG) -I$(top_srcdir)/src @FLOAT_PF_FLAG@ \
                  -perl5 -proxy -const \
                  -module RNA -o $(SWIG_wrapper) \
                  $(SWIG_main_src)

##
## copy shared library for perl-tests
##

auto/RNA/RNA.so: RNA.la
	$(MKDIR_P) auto/RNA
	$(LIBTOOL) --config > libtoolconfig.tmp
## ksh requires a path on the sourced file.
	. ./libtoolconfig.tmp; cp "$$objdir/RNA.so" auto/RNA
	rm -f libtoolconfig.tmp

##
## Perl Test Suite
##

AM_TESTS_ENVIRONMENT = PERL5LIB='$(srcdir)/.:${PERL5LIB}'; export PERL5LIB;
LOG_COMPILER = $(PERL)
AM_LOG_FLAGS = -Mstrict -w
TESTS = test.pl \
        RNA/t/Design.t

clean-local:
	-rm -rf auto test.coord test_dp.ps test_ss.ps test_ss_a.ps test.par

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
