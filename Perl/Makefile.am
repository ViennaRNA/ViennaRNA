SOURCES = RNA.i ptr2array.i Makefile.PL 

EXTRA_DIST = RNA.pm RNA_wrap.doc RNA_wrap.c RNAfold.pl test.pl RNAfold.cgi

Makefile.perl: Makefile.PL Makefile
	@PerlCmd@ Makefile.PL

all-local: Makefile.perl
	$(MAKE) -f Makefile.perl

install-exec-local: all-local
	$(MAKE) -f Makefile.perl install

check-local: all-local
	$(MAKE) -f Makefile.perl test

mostlyclean-local: Makefile.perl
	$(MAKE) -f Makefile.perl clean
