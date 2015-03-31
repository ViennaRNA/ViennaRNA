include ../generic.mk

PERL_MODULE = $(SWIG_module_name).pm

SOURCES = Makefile.PL RNA.pod
EXTRA_DIST = \
  $(SWIG_wrapper) \
  $(PERL_MODULE) \
  RNAfold.pl \
  test.pl \
  RNAfold.cgi \
  Makefile.PL

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

Makefile.perl: Makefile.PL Makefile
	$(AM_V_GEN)@PerlCmd@ Makefile.PL PREFIX=$(prefix) INSTALL_BASE= 


##
## Misc make rules
##

all-local: Perl_interface
	

install-exec-local: all-local
	$(MAKE) -f Makefile.perl install

##uninstall-local: all-local
##	$(MAKE) -f Makefile.perl uninstall

check-local: all-local
	$(MAKE) -f Makefile.perl test

mostlyclean-local: Makefile.perl
	$(MAKE) -f Makefile.perl clean

# remove all generated content
maintainer-clean-local:
	-rm -f $(SWIG_wrapper) Makefile.perl.old $(PERL_MODULE)

