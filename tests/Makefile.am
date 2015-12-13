AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS) @FLOAT_PF_FLAG@
LDADD = $(top_builddir)/src/ViennaRNA/libRNA.a @CHECK_LIBS@

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
%.c: %.ts
	checkmk $< > $@
