AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = $(RNA_CPPFLAGS) -I$(top_srcdir)/src
AM_CFLAGS = $(RNA_CFLAGS) @CHECK_CFLAGS@ -I$(top_builddir)/src/
AM_CXXFLAGS = $(RNA_CXXFLAGS)
AM_LDFLAGS = $(RNA_LDFLAGS)

LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la @CHECK_LIBS@

SUFFIXES = .c .ts
TEST_EXTENSIONS = .pl .t .py .py3

if WITH_UNIT_TESTS

AM_TESTS_ENVIRONMENT =  \
## Some environment initializations are kept in a separate shell
## file 'tests-env.sh', which can make it easier to also run tests
## from the command line.
                        . $(builddir)/test-env.sh; \
## On Solaris, prefer more POSIX-compliant versions of the standard
## tools by default.
                        if test -d /usr/xpg4/bin; then \
                          PATH=/usr/xpg4/bin:$$PATH; export PATH; \
                        fi;

## With this, the test scripts will be able to print diagnostic
## messages to the original standard error stream, even if the test
## driver redirects the stderr of the test scripts to a log file
## before executing them.
AM_TESTS_FD_REDIRECT = 9>&2

######################################
## unit tests for Perl5 interface   ##
######################################
if WITH_PERL_TESTS

BUILT_SOURCES = perl5/test-RNA-Design.t

PL_LOG_COMPILER = $(PERL)
T_LOG_COMPILER = $(PERL)
AM_PL_LOG_FLAGS = -Mstrict -w \
                  -I$(top_builddir)/interfaces/Perl \
                  -I$(top_srcdir)/interfaces/Perl
AM_T_LOG_FLAGS = -Mstrict -w \
                  -I$(top_builddir)/interfaces/Perl \
                  -I$(top_srcdir)/interfaces/Perl

# add test scripts here
PERL_TESTS =  perl5/test-RNA.pl \
              perl5/test-RNA-callbacks.pl \
              perl5/test-RNA-combinatorics.pl \
              perl5/test-RNA-Design.t \
              perl5/test-RNA-mfe_eval.pl \
              perl5/test-RNA-mfe_window.pl \
              perl5/test-RNA-constraints.pl


# add test script output files here
PERL_TEST_OUTPUT =  test.coord \
                    test_dp.ps \
                    test_ss.ps \
                    test_ss_a.ps \
                    test.par \
                    test-RNA-mfe_eval.pl.out \
                    outputFile_test.txt

perl5/test-RNA-Design.t: $(top_srcdir)/interfaces/Perl/RNA/t/Design.t
	cp $(top_srcdir)/interfaces/Perl/RNA/t/Design.t $(builddir)/perl5/test-RNA-Design.t

endif

######################################  
## unit tests for Python2 interface ##
######################################
if WITH_PYTHON_TESTS

PY_LOG_COMPILER = $(PYTHON2)
AM_PY_LOG_FLAGS =

# add test scripts here
PYTHON2_TESTS = python/test-RNA-combinatorics.py \
                python/test-RNA-constraints.py \
                python/test-RNA-mfe_eval.py \
                python/test-RNA-mfe_window.py \
                python/test-RNA-sc-callbacks.py \
                python/test-RNA-subopt.py \
                python/test-RNA-utils.py \
                python/test-RNA.py

# add test script output files here
PYTHON2_TEST_OUTPUT = test-RNA-mfe_eval.py.out

endif

######################################
## unit tests for Python3 interface ##
######################################
if WITH_PYTHON3_TESTS

PY3_LOG_COMPILER = $(PYTHON3)
AM_PY3_LOG_FLAGS =

# add test scripts here
PYTHON3_TESTS = python3/test-RNA-combinatorics.py3 \
                python3/test-RNA-constraints.py3 \
                python3/test-RNA-mfe_eval.py3 \
                python3/test-RNA-mfe_window.py3 \
                python3/test-RNA-sc-callbacks.py3 \
                python3/test-RNA-subopt.py3 \
                python3/test-RNA-utils.py3 \
                python3/test-RNA.py3

# add test script output files here
PYTHON3_TEST_OUTPUT = test-RNA-mfe_eval.py3.out

endif

######################################
## unit tests for C-library         ##
######################################
if WITH_CHECK

CHECKMK_FILES = \
              energy_evaluation.ts \
              constraints.ts \
              constraints_soft.ts \
              fold.ts \
              utils.ts \
              eval_structure.ts \
              walk.ts \
              neighbor.ts

CHECK_CFILES = \
              energy_evaluation.c \
              constraints.c \
              constraints_soft.c \
              fold.c \
              utils.c \
              eval_structure.c \
              walk.c \
              neighbor.c

LIBRARY_TESTS = energy_evaluation \
                constraints \
                constraints_soft \
                fold \
                utils \
                eval_structure \
                walk \
                neighbor

check_PROGRAMS = ${LIBRARY_TESTS}

endif

########################################
## unit tests for executable programs ##
########################################
if WITH_EXECUTABLE_TESTS

EXECUTABLE_TESTS  = \
                  RNAfold/general.sh \
                  RNAfold/partfunc.sh \
                  RNAfold/special.sh \
                  RNAfold/long.sh

endif

TESTS = \
      ${EXECUTABLE_TESTS} \
      ${LIBRARY_TESTS} \
      ${PERL_TESTS} \
      ${PYTHON2_TESTS} \
      ${PYTHON3_TESTS}

## Below are the rules how to create the test suite executables
## See checkmk man page for details about *.ts files
.ts.c:
	$(AM_V_GEN)checkmk $< > $@

endif

EXTRA_DIST =  data \
              RNAfold/results \
              ${CHECKMK_FILES} ${CHECK_CFILES} \
              ${PERL_TESTS} \
              ${PYTHON2_TESTS} \
              ${PYTHON3_TESTS} \
              ${EXECUTABLE_TESTS}

clean-local:
	-rm -rf ${PERL_TEST_OUTPUT} \
                $(PYTHON2_TEST_OUTPUT) \
                $(PYTHON3_TEST_OUTPUT) \
                *.pyc \
                __pycache__
