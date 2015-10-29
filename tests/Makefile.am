AM_CPPFLAGS = -I$(top_srcdir)/src

TESTS = check_check

check_PROGRAMS = check_check
check_check_SOURCES = check_check_main.c check_constraints.c check_fold.c check_loop_energies.c check_utils.c
check_check_CFLAGS = @CHECK_CFLAGS@ -I $(top_builddir)/src/  $(OPENMP_CFLAGS)
check_check_LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@

