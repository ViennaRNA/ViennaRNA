AM_CPPFLAGS = -I$(top_srcdir)/src

EXTRA_DIST =  energy_evaluation.ts \
              constraints.ts \
              fold.ts \
              utils.ts

TESTS = energy_evaluation \
        constraints \
        fold \
        utils

check_PROGRAMS = ${TESTS}

## Below are the rules how to create the test suite executables
## See checkmk man page for details about *.ts files

# Check structure prediction interface

fold.c:  fold.ts
	checkmk fold.ts > fold.c

fold_SOURCES = fold.c
fold_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS)
fold_LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@

# Check energy evaluation interface

energy_evaluation.c:  energy_evaluation.ts
	checkmk energy_evaluation.ts > energy_evaluation.c

energy_evaluation_SOURCES = energy_evaluation.c
energy_evaluation_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS)
energy_evaluation_LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@

# Check constraints interface

constraints.c:  constraints.ts
	checkmk constraints.ts > constraints.c

constraints_SOURCES = constraints.c
constraints_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS)
constraints_LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@

# Check utilities

utils.c:  utils.ts
	checkmk utils.ts > utils.c

utils_SOURCES = utils.c
utils_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS)
utils_LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@
