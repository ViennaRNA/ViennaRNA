#!/usr/bin/perl

use warnings;
use strict;

use RNA;

my $seq1 = "CGCAGGGAUACCCGCG";

# compute minimum free energy (mfe) and corresponding structure
my ($ss, $mfe) = RNA::fold($seq1);

# print output
printf "%s [ %6.2f ]\n", $ss, $mfe;
