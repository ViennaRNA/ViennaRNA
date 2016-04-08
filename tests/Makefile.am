AUTOMAKE_OPTIONS = subdir-objects

AM_CPPFLAGS = -I$(top_srcdir)/src
AM_CFLAGS = @CHECK_CFLAGS@ -I$(top_builddir)/src/  $(OPENMP_CFLAGS) @FLOAT_PF_FLAG@
LDADD = $(top_builddir)/src/ViennaRNA/libRNA_conv.la @CHECK_LIBS@

SUFFIXES = .c .ts
TEST_EXTENSIONS = .pl .t .py .py3

if WITH_UNIT_TESTS

AM_TESTS_ENVIRONMENT =  PYTHONPATH='.:${PYTHONPATH}'; \
                        export PYTHONPATH;

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
              perl5/test-RNA-Design.t \
              perl5/test-RNA-mfe_eval.pl

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
PYTHON2_TESTS = python/test-RNA-callbacks.py \
                python/test-RNA-mfe_eval.py \
                python/test-RNA-constraints.py \
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
PYTHON3_TESTS = python3/test-RNA-callbacks.py3 \
                python3/test-RNA-mfe_eval.py3 \
                python3/test-RNA-constraints.py3 \
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
              fold.ts \
              utils.ts \
              eval_structure.ts 

CHECK_CFILES = \
              energy_evaluation.c \
              constraints.c \
              fold.c \
              utils.c \
              eval_structure.c

LIBRARY_TESTS = energy_evaluation \
                constraints \
                fold \
                utils \
                eval_structure 

check_PROGRAMS = ${LIBRARY_TESTS}

endif

TESTS = ${LIBRARY_TESTS} ${PERL_TESTS} ${PYTHON2_TESTS} ${PYTHON3_TESTS}

## Below are the rules how to create the test suite executables
## See checkmk man page for details about *.ts files
.ts.c:
	$(AM_V_GEN)checkmk $< > $@

endif

EXTRA_DIST =  data \
              ${CHECKMK_FILES} ${CHECK_CFILES} \
              ${PERL_TESTS} \
              ${PYTHON2_TESTS} \
              ${PYTHON3_TESTS}

clean-local:
	-rm -rf ${PERL_TEST_OUTPUT} \
                $(PYTHON2_TEST_OUTPUT) \
                $(PYTHON3_TEST_OUTPUT) \
                *.pyc \
                __pycache__
