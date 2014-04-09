if WITH_LARGE_PF
  USE_LARGE_PF = -DLARGE_PF
else
  USE_LARGE_PF =
endif

SOURCES = Makefile.PL RNA.pod

EXTRA_DIST = RNA.pm RNA_wrap.c RNAfold.pl test.pl RNAfold.cgi Makefile.PL
AM_CPPFLAGS = -I$(top_srcdir)/src


Makefile.perl: Makefile.PL Makefile
	$(AM_V_GEN)@PerlCmd@ Makefile.PL PREFIX=$(prefix) INSTALL_BASE= 

# the perl command repairs some errors in the SWIG CVS version
RNA_wrap.c: ../RNA.i ../tmaps.i
	$(AM_V_GEN)swig -perl5 -shadow $(USE_LARGE_PF) -I../ RNA.i

Perl_interface: RNA_wrap.c Makefile.perl
	$(AM_V_GEN)$(MAKE) -f Makefile.perl

all-local: Perl_interface
	

install-exec-local: all-local
	$(MAKE) -f Makefile.perl install

check-local: all-local
	$(MAKE) -f Makefile.perl test

mostlyclean-local: Makefile.perl
	$(MAKE) -f Makefile.perl clean

# remove all generated content
maintainer-clean-local:
	-rm -f RNA_wrap.c Makefile.perl.old RNA.pm

