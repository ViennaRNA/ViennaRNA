#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 6;
use Data::Dumper;
use FileHandle;

use RNA;
use warnings;


######################### End of black magic.

##########################################################################################################################
#starting with ensemble_defect test

my $seq       = "AGGAAACCUUAAUUGGUUA";
my $structpk  = ".((...))(([[..))]].";
my $seq2      = "AAAAAAAA";
my $struct2   = "(......)";
my $struct3   = "..(..)..";

####################################################
##test_ensemble_defect
print "test_ensemble_defect\n";
my $fc = new RNA::fold_compound($seq);
(my $ss, my $gfe) = $fc->pf();

my $ed = $fc->ensemble_defect($structpk);
is($ed, 0.614080983833787, "Ensemble Defect");

my $pt = RNA::ptable($structpk);
$ed = $fc->ensemble_defect($pt);
is($ed, 0.614080983833787, "Ensemble Defect (pair table)");

$ed = $fc->ensemble_defect($structpk, RNA::BRACKETS_ANY);
is($ed, 0.7279184335061499, "Ensemble Defect (pk)");

$pt = RNA::ptable($structpk, RNA::BRACKETS_ANY);
$ed = $fc->ensemble_defect($pt);
is($ed, 0.7279184335061499, "Ensemble Defect (pk, pair table)");


$fc = new RNA::fold_compound($seq2);
($ss, $gfe) = $fc->pf();

$ed = $fc->ensemble_defect($struct2) * length($seq2);
is($ed, 2.0, "Ensemble Defect (Sanity check 1)");

$ed = $fc->ensemble_defect($struct3) * length($seq2);
is($ed, 2.0, "Ensemble Defect (Sanity check 2)");

undef $pt;
undef $fc;
