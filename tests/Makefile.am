
TESTS = check_fold

check_PROGRAMS = check_fold
check_fold_SOURCES = check_constraints.c check_fold.c check_loop_energies.c check_utils.c
check_fold_CFLAGS = @CHECK_CFLAGS@ -I $(top_builddir)/src/  $(OPENMP_CFLAGS)
check_fold_LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@

