pkgexampledir = $(pkgdatadir)/examples
pkgexampledatadir = $(pkgexampledir)/data
pkgperlexampledir = $(pkgexampledir)/perl5
pkgpythonexampledir = $(pkgexampledir)/python

examples_c = \
    callback_subopt.c \
    example1.c \
    example_old.c \
    fold_compound_md.c \
    fold_compound_mfe.c \
    helloworld_mfe.c \
    helloworld_mfe_comparative.c \
    helloworld_probabilities.c \
    soft_constraints_up.c

examples_perl = \
    Perl5/helloworld_mfe.pl \
    Perl5/helloworld_mfe_comparative.pl \
    Perl5/helloworld_nondefault.pl \
    Perl5/maximum_matching.pl \
    Perl5/oo_example1.pl

examples_python = \
    Python/boltzmann_sampling.py \
    Python/helloworld_mfe_comparative.py \
    Python/helloworld_mfe.py \
    Python/helloworld_nondefault.py \
    Python/maximum_matching.py \
    Python/oo_example1.py \
    Python/subopt.py \
    Python/RNAfold_MEA.py

data_files = \
    files/alignment_clustal.aln \
    files/alignment_fasta.fa \
    files/alignment_stockholm.stk \
    files/alignment_maf.maf

pkgexample_DATA = $(examples_c)
pkgperlexample_DATA = $(examples_perl)
pkgpythonexample_DATA = $(examples_python)
pkgexampledata_DATA = $(data_files)


EXTRA_DIST = \
    $(examples_c) \
    $(examples_perl) \
    $(examples_python) \
    $(data_files)
