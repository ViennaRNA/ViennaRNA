#!/usr/bin/perl
#

use strict;
use warnings;
use Data::Dumper;
use Test::More tests => 2;

use RNA;

my @content   = (1, 1, 1, 1);
my @content2  = (3, 1, 2);
my @perm_real   = ( [0, 3, 2, 1],
                  [0, 3, 1, 2],
                  [0, 2, 3, 1],
                  [0, 2, 1, 3],
                  [0, 1, 3, 2],
                  [0, 1, 2, 3]);
my @perm_real2  = (
                  [1, 0, 0, 0, 2, 2],
                  [1, 0, 0, 2, 0, 2],
                  [1, 0, 0, 2, 2, 0],
                  [1, 0, 2, 0, 0, 2],
                  [1, 0, 2, 0, 2, 0],
                  [1, 0, 2, 2, 0, 0],
                  [1, 2, 0, 0, 0, 2],
                  [1, 2, 0, 0, 2, 0],
                  [1, 2, 0, 2, 0, 0],
                  [1, 2, 2, 0, 0, 0]);
my %reshash;
my %goldhash;

map { $goldhash{join(",", @{$_})} = 1 } @perm_real;

print "test_enumerate_necklaces (4 sequences, 1 strand per sequences)\n";
my $permutations = RNA::enumerate_necklaces(\@content);
map { $reshash{join(",", @{$_})} = 1 } @{$permutations};

is_deeply(\%goldhash, \%reshash);

%reshash  = ();
%goldhash = ();

map { $goldhash{join(",", @{$_})} = 1 } @perm_real2;

print "test_enumerate_necklaces (3 sequences, different number of strands per sequence)\n";
$permutations = RNA::enumerate_necklaces(\@content2);
map { $reshash{join(",", @{$_})} = 1 } @{$permutations};

is_deeply(\%goldhash, \%reshash);
