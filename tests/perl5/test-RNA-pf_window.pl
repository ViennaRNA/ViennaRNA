#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)

use strict;
use Test::More tests => 18;
use Data::Dumper;
use FileHandle;

use RNA;
use warnings;


my $longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU";
my $md;
my $fc;
my %data;
my $bpp;
my $up;
my $max_pos;
my $min_v;
my $max_v;
my $min_entries;
my $max_entries;


sub pf_window_callback {

  my ($v, $v_size, $i, $maxsize, $what, $data) = @_;

  if ($what & RNA::PROBS_WINDOW_UP) {
    my %d_up = ('i' => $i, 'up' => $v);
    push @{$data->{'unpaired_probs'}}, \%d_up;
  } else {
    for (my $j = 0; $j <= $#$v; $j++) {
      next if !defined($v->[$j]);
      next if $v->[$j] < 0.01;
      my %d_bpp = ('i' => $i, 'j' => $j, 'p' => $v->[$j]);
      push @{$data->{'pair_probs'}}, \%d_bpp;
    }
  }
}


print "test_pfl_fold\n";
$bpp = RNA::pfl_fold($longseq, 200, 150, 0.01);

# sanity check for base pair probabilities
$max_pos = 0;
$min_v = 2.;
$max_v = -1;
for my $bp (@{$bpp}) {
  $min_v = $bp->{p} if $min_v > $bp->{p};
  $max_v = $bp->{p} if $max_v < $bp->{p};
  $max_pos = $bp->{i} if $max_pos < $bp->{i};
  $max_pos = $bp->{j} if $max_pos < $bp->{j};
}
ok(scalar(@{$bpp}) == 640);
ok($min_v >= 0.01);
ok($max_v <= 1.0);
ok($max_pos == 300);


print "test_pfl_fold_up\n";
$up = RNA::pfl_fold_up($longseq, 10, 200, 150);

# sanity check for unpaired probabilities
$min_entries = 0;
$max_entries = 0;
$min_v = 2.;
$max_v = -1;
for my $u (@{$up}) {
  next if ! defined($u);
  $min_entries = scalar(@{$u}) if $min_entries < scalar(@{$u});
  $max_entries = scalar(@{$u}) if $max_entries < scalar(@{$u});
  for (@{$u}) {
    next if !defined($_);
    $min_v = $_ if $_ < $min_v;
    $max_v = $_ if $_ > $max_v;
  };
}
ok(scalar(@{$up}) == 301); # one entry for each of the 300nt of longseq plus the 0th element
ok($min_entries == 11); # 10 + 0th (undef) entry
ok($max_entries == 11); # 10 + 0th (undef) entry
ok($min_v >= 0);
ok($max_v <= 1.0);


print "test_probs_window\n";
$md = new RNA::md();
$md->{max_bp_span} = 150;
$md->{window_size} = 200;
$fc = new RNA::fold_compound($longseq, $md, RNA::OPTION_WINDOW);
%data = ( 'pair_probs' => [], 'unpaired_probs' => [] );
$fc->probs_window(10, RNA::PROBS_WINDOW_BPP | RNA::PROBS_WINDOW_UP, \&pf_window_callback, \%data);

# sanity check for base pair probabilities
$max_pos = 0;
$min_v = 2.;
$max_v = -1;
for my $bp (@{$data{'pair_probs'}}) {
  $min_v = $bp->{p} if $min_v > $bp->{p};
  $max_v = $bp->{p} if $max_v < $bp->{p};
  $max_pos = $bp->{i} if $max_pos < $bp->{i};
  $max_pos = $bp->{j} if $max_pos < $bp->{j};
}
ok(scalar(@{$data{'pair_probs'}}) == 640);
ok($min_v >= 0.01);
ok($max_v <= 1.0);
ok($max_pos == 300);

# sanity check for unpaired probabilities
$min_entries = 0;
$max_entries = 0;
$min_v = 2.;
$max_v = -1;
for my $u (@{$data{'unpaired_probs'}}) {
  next if ! defined($u);
  $min_entries = scalar(@{$u->{'up'}}) if $min_entries < scalar(@{$u->{'up'}});
  $max_entries = scalar(@{$u->{'up'}}) if $max_entries < scalar(@{$u->{'up'}});
  for (@{$u->{'up'}}) {
    next if !defined($_);
    $min_v = $_ if $_ < $min_v;
    $max_v = $_ if $_ > $max_v;
  };
}
ok(scalar(@{$data{'unpaired_probs'}}) == 300); # one entry for each of the 300nt of longseq
ok($min_entries == 11); # 10 + 0th (undef) entry
ok($max_entries == 11); # 10 + 0th (undef) entry
ok($min_v >= 0);
ok($max_v <= 1.0);
