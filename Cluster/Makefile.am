bin_PROGRAMS = AnalyseDists AnalyseSeqs
man_MANS = AnalyseDists.1 AnalyseSeqs.1

AnalyseDists_SOURCES = split.c \
          distance_matrix.c \
          cluster.c \
          treeplot.c \
          AD_main.c

AnalyseSeqs_SOURCES = statgeom.c \
          distance_matrix.c \
          cluster.c \
          AS_main.c \
          split.c \
          treeplot.c \
          PS3D.c
noinst_HEADERS = cluster.h split.h treeplot.h distance_matrix.h statgeom.h \
	utils.h PS3D.h StrEdit_CostMatrix.h  

LDADD = utils.o -lm

utils.o:
	$(COMPILE) -c ../lib/utils.c

EXTRA_DIST = $(man_MANS)
