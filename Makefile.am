
if MAKE_PERL_EXT
   PLDIR = Perl
endif

if MAKE_CLUSTER
  CLDIR = Cluster
endif

if MAKE_FORESTER
   FDIR = RNAforester
endif

if MAKE_KINFOLD
   KINDIR = Kinfold
endif

SUBDIRS = lib Progs man Utils H $(PLDIR) $(CLDIR) $(KINDIR) $(FDIR)

pkgdata_DATA = vienna1.8.4.par default.par

EXTRA_DIST = $(pkgdata_DATA) INSTALL.configure Readseq @SVM_SOURCE_ARCHIVE@
