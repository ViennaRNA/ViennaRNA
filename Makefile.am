
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

SUBDIRS = lib Progs man Utils H $(PLDIR) $(CLDIR) $(FDIR) $(KINDIR)

#DIST_SUBDIRS = lib Progs man Utils AliFold H Perl Cluster

pkgdata_DATA = vienna13.par default.par

EXTRA_DIST = $(pkgdata_DATA) INSTALL.configure Readseq
