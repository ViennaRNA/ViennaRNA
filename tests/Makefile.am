AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS) @FLOAT_PF_FLAG@
LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la @CHECK_LIBS@

SUFFIXES = .c .ts

CHECKMK_FILES = \
              energy_evaluation.ts \
              constraints.ts \
              fold.ts \
              utils.ts

CHECK_CFILES = \
              energy_evaluation.c \
              constraints.c \
              fold.c \
              utils.c

TESTS = energy_evaluation \
        constraints \
        fold \
        utils

check_PROGRAMS = ${TESTS}

## Below are the rules how to create the test suite executables
## See checkmk man page for details about *.ts files
.ts.c:
	$(AM_V_GEN)checkmk $< > $@

EXTRA_DIST =  ${CHECKMK_FILES} ${CHECK_CFILES}
