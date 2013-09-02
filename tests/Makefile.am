
TESTS = check_fold

check_PROGRAMS = check_fold
check_fold_SOURCES = check_fold.c $(top_builddir)/src/ViennaRNA/fold.h
check_fold_CFLAGS = @CHECK_CFLAGS@ -I $(top_builddir)/src/  $(OPENMP_CFLAGS)
check_fold_LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@

