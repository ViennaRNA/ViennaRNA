SOURCES = RNA.i ptr2array.i Makefile.PL 

EXTRA_DIST = RNA.pm RNA_wrap.c RNAfold.pl test.pl RNAfold.cgi RNA.i Makefile.PL

Makefile.perl: Makefile.PL Makefile
	@PerlCmd@ Makefile.PL

# the perl command repairs some errors in the SWIG CVS version
RNA_wrap.c: RNA.i
	swig -perl5 -shadow RNA.i
#	@PerlCmd@ -i -npe 's/\A\s+\;\s+\Z//; s/SWIGTYPE_p_unsigned_short/SWIGTYPE_p_short/' RNA_wrap.c

all-local: RNA_wrap.c Makefile.perl
	$(MAKE) -f Makefile.perl

install-exec-local: all-local
	$(MAKE) -f Makefile.perl install

check-local: all-local
	$(MAKE) -f Makefile.perl test

mostlyclean-local: Makefile.perl
	$(MAKE) -f Makefile.perl clean
