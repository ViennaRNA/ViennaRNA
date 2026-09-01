#!/usr/bin/perl
#

use RNA;

my $seq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGA";

# create fold_compound data structure
# (required for all subsequently applied  algorithms)
my $fc = new RNA::fold_compound($seq);

# compute MFE and MFE structure
my ($mfe_struct, $mfe) = $fc->mfe();

# rescale Boltzmann factors for partition function computation
$fc->exp_params_rescale($mfe);

# compute partition function
my ($pp, $pf) = $fc->pf();

# compute centroid structure
my ($centroid_struct, $dist) = $fc->centroid();

# compute free energy of centroid structure
my $centroid_en = $fc->eval_structure($centroid_struct);

# print only the centroid structure, its free energy
# in kcal/mol and the distance to the ensemble
printf("%s\n%s (%6.2f d = %.2f)\n", $seq, $centroid_struct, $centroid_en, $dist);
