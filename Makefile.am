PLDIR = 
CLDIR = 
if MAKE_PERL_EXT
   PLDIR = Perl
endif

if MAKE_CLUSTER
  CLDIR = Cluster
endif

SUBDIRS = lib Progs man Subopt Utils AliFold H $(PLDIR) $(CLDIR)
DIST_SUBDIRS = lib Progs man Subopt Utils AliFold H Perl Cluster

pkgdata_DATA = vienna13.par default.par

EXTRA_DIST = $(pkgdata_DATA)
