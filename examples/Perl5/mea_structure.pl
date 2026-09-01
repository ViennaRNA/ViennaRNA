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

# compute MEA structure
my ($mea_struct, $mea) = $fc->MEA(1.0);

# compute free energy of MEA structure
my $mea_en = $fc->eval_structure($mea_struct);

# print only the MEA structure, its free energy
# in kcal/mol and the maximum expected accuracy
printf("%s\n%s (%6.2f MEA = %.2f)\n", $seq, $mea_struct, $mea_en, $mea);
