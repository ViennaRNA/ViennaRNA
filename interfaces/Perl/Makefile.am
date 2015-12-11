include ../generic.mk

PERL_MODULE = $(SWIG_module_name).pm

SOURCES = Makefile.PL RNA.pod
EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(PERL_MODULE) \
  RNAfold.pl \
  test.pl \
  RNAfold.cgi \
  Makefile.PL \
  RNA/Benchmark.pm \
  RNA/Design.pm \
  RNA/Files.pm \
  RNA/Reliability.pm \
  RNA/Utils.pm \
  RNA/t/Design.t

AM_CPPFLAGS = -I$(top_srcdir)/src

if WITH_LARGE_PF
  USE_LARGE_PF = -DLARGE_PF
else
  USE_LARGE_PF =
endif


##
## SWIG interface wrapper code generation
##

$(SWIG_wrapper) $(PERL_MODULE): $(SWIG_src)
	$(AM_V_GEN)swig $(INCLUDES) $(USE_LARGE_PF) \
                  -perl5 -proxy -const \
                  -module RNA -o $(SWIG_wrapper) \
                  $(SWIG_main_src)


##
## Perl interface building
##

Perl_interface: $(SWIG_wrapper) Makefile.perl
	$(AM_V_GEN)$(MAKE) -f Makefile.perl



#
# MAD HACK #1
#
# Create Makefile.perl via Perls MakeMaker and patch its uninstall
# targets such that files in prefixed install targets are reported
# by the uninstall target too
#
Makefile.perl: Makefile.PL Makefile
	$(AM_V_GEN)@PerlCmd@ Makefile.PL PREFIX=$(prefix) INSTALL_BASE= CC=$(CC)
	$(SED) -i \
  's#$$(NOECHO) $$(UNINSTALL) $$(PERL_ARCHLIB)/auto/$$(FULLEXT)/.packlist#$$(NOECHO) $$(UNINSTALL) $$(DESTINSTALLARCHLIB)/auto/$$(FULLEXT)/.packlist#g' \
  Makefile.perl
	$(SED) -i \
  's#$$(NOECHO) $$(UNINSTALL) $$(SITEARCHEXP)/auto/$$(FULLEXT)/.packlist#$$(NOECHO) $$(UNINSTALL) $$(DESTINSTALLSITEARCH)/auto/$$(FULLEXT)/.packlist#g' \
  Makefile.perl
	$(SED) -i \
  's#$$(NOECHO) $$(UNINSTALL) $$(VENDORARCHEXP)/auto/$$(FULLEXT)/.packlist#$$(NOECHO) $$(UNINSTALL) $$(DESTINSTALLVENDORARCH)/auto/$$(FULLEXT)/.packlist#g' \
  Makefile.perl

##
## Misc make rules
##

all-local: Perl_interface
	

install-exec-local: all-local
	$(MAKE) -f Makefile.perl pure_install

#
# MAD HACK #2
# Perform actual uninstall for files that won't get uninstalled but
# only reported by MakeMakers original uninstall target
#
uninstall-local:
	`$(MAKE) -f Makefile.perl uninstall | $(GREP) unlink | $(SED) -e 's#unlink#rm -f#'`

check-local: all-local
	$(MAKE) -f Makefile.perl test

mostlyclean-local: Makefile.perl
	$(MAKE) -f Makefile.perl clean

# remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) Makefile.perl.old $(PERL_MODULE)

distclean-local:
	-rm -f $(SWIG_wrapper) Makefile.perl.old $(PERL_MODULE)
