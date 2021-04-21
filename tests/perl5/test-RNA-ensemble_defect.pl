#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 4;
use Data::Dumper;
use FileHandle;

use RNA;
use warnings;


######################### End of black magic.

##########################################################################################################################
#starting with ensemble_defect test

my $seq     = "AGGAAACCUUAAUUGGUUA";
my $structpk= ".((...))(([[..))]].";

####################################################
##test_ensemble_defect
print "test_ensemble_defect\n";
my $fc = new RNA::fold_compound($seq);
(my $ss, my $gfe) = $fc->pf();

my $ed = $fc->ensemble_defect($structpk);
is($ed, 0.6140797258673892);

my $pt = RNA::ptable($structpk);
$ed = $fc->ensemble_defect($pt);
is($ed, 0.6140797258673892);

$ed = $fc->ensemble_defect($structpk, RNA::BRACKETS_ANY);
is($ed, 0.7279171755397522);

$pt = RNA::ptable($structpk, RNA::BRACKETS_ANY);
$ed = $fc->ensemble_defect($pt);
is($ed, 0.7279171755397522);

undef $pt;
undef $fc;
