#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)

use strict;
use warnings;
use Test::More tests => 23;
use Data::Dumper;
use FileHandle;

use RNA;
use RNAHelpers qw(:Paths);


my $datadir = getDataDirPath();
my $shortseq = "UGGGAAUAGUCUCUUCCGAGUCUCGCGGGCGACGGGCGAUCUUCGAAAGUGGAAUCCGUA";
my $longseq = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUCGCUACUUAGGCGAUCCCUGCCAAUGAGGGUCAAGGAGUUGAAUUAUCGGGCCACAUCGACGUGGCCUUUACGGCCAGGUAAUUCAAAGGCCUCAAGUCCU";
my $md;
my $fc;
my $fc2;
my %data;
my @data;
my $bpp;
my $up;
my $max_pos;
my $min_v;
my $max_v;
my $min_entries;
my $max_entries;


sub getShapeDataFromFile {
    my $filepath = shift(@_);
    my @retVec;
    push(@retVec,-999.0);# data list is 1-based, so we add smth. at pos 0
    my $count=1;
    open(my $fh, "<", $filepath) || die "Couldn't open '".$filepath."' for reading because: ".$!;

    foreach my $line (<$fh>)
    {

        my @sp= split(/ /,$line);
        my $pos = $sp[0];
        my $value = $sp[1];

        if($pos == $count)
        {
            push(@retVec,$value +0.00); # make float from string
        } else {
            for(my $i=$pos;$i<=$count;$i++)
            {
                push(@retVec,-999.0);
            }
            push(@retVec,$value);
            $count=$pos;
        }
        $count+=1;

    }
    return @retVec;
}


sub getShapeSequenceFromFile {
    my $f = shift(@_);
    my $retSeq ="";
    open(my $fh, "<", $f) || die "Couldn't open '".$f."' for reading because: ".$!;
    foreach my $line (<$fh>)
    {
        chomp $line;
        return $line; #return the first line
    }
}


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


sub bpp_callback {
  my ($v, $v_size, $i, $maxsize, $what, $data) = @_;

  if ($what & RNA::PROBS_WINDOW_BPP) {
    for (my $j = 0; $j <= $#$v; $j++) {
      next if !defined($v->[$j]);
      next if $v->[$j] < 0.01;
      my %d_bpp = ('i' => $i, 'j' => $j, 'p' => $v->[$j]);
      push @{$data}, \%d_bpp;
    }
  }
}


sub up_callback {
  my ($v, $v_size, $i, $maxsize, $what, $data) = @_;

  if ($what & RNA::PROBS_WINDOW_UP) {
    my %d_up = ('i' => $i, 'up' => $v);
    push @{$data}, \%d_up;
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


print "test_pfl_SHAPE\n";
my @benchmark_set = ( "Lysine_riboswitch_T._martima", "TPP_riboswitch_E.coli" );

for my $b (@benchmark_set) {
  my $seq           = getShapeSequenceFromFile($datadir . $b . ".db");
  my @reactivities  = getShapeDataFromFile($datadir . $b . ".shape_2rows");
  %data = ( 'pair_probs' => [], 'unpaired_probs' => [] );

  $md = new RNA::md();
  $md->{'max_bp_span'} = length($seq);
  $md->{'window_size'} = length($seq);

  # compute pairing probabilities using local partition function implementation
  # with L = W = len(seq)
  $fc = new RNA::fold_compound($seq, $md, RNA::OPTION_WINDOW);
  $fc->sc_add_SHAPE_deigan(\@reactivities, 1.9, -0.7, RNA::OPTION_WINDOW);
  $fc->probs_window(0, RNA::PROBS_WINDOW_BPP, \&pf_window_callback, \%data);

  # compute pairing probabilities using global parition function
  my $fc2 = new RNA::fold_compound($seq);
  $fc2->sc_add_SHAPE_deigan(\@reactivities, 1.9, -0.7, RNA::OPTION_DEFAULT);
  $fc2->pf();
  my $bpp = $fc2->bpp();

  # check for differences in pairing probabilities between local and global implementation
  # Hint: There must not be any!
  my $equal_probabilities = 1;
  foreach my $prob (@{$data{'pair_probs'}}) {
      my $p = $prob->{'p'};
      my $i = $prob->{'i'};
      my $j = $prob->{'j'};
      my $p2 = $bpp->[$i][$j];
      ($equal_probabilities = 0, last) if ! (sprintf("%g", $p) eq sprintf("%g", $p2));
  }
  ok($equal_probabilities == 1);
}


# Compute partition function and base pair probabilities both, using the implementation
# for local structure prediction and global structure prediction.
# When comparing both results, equilibrium probabilities must not have changed!
print "test_probs_window_full\n";
@data = ();
$md = new RNA::md();
$md->{'max_bp_span'} = length($shortseq);
$md->{'window_size'} = length($shortseq);

# compute pairing probabilities using local partition function implementation
# with L = W = length(seq)
$fc = new RNA::fold_compound($shortseq, $md, RNA::OPTION_WINDOW);
$fc->probs_window(0, RNA::PROBS_WINDOW_BPP, \&bpp_callback, \@data);

# compute pairing probabilities using global parition function
$fc2 = new RNA::fold_compound($shortseq);
$fc2->pf();
$bpp = $fc2->bpp();

# check for differences in pairing probabilities between local and global implementation
# Hint: There must not be any!
my $equal_probabilities = 1;
foreach my $prob (@data) {
  my $p = $prob->{'p'};
  my $i = $prob->{'i'};
  my $j = $prob->{'j'};
  my $p2 = $bpp->[$i][$j];
  ($equal_probabilities = 0, last) if ! (sprintf("%g", $p) eq sprintf("%g", $p2));
}
ok($equal_probabilities == 1);


# Compute partition function and base pair probabilities both, constrained
# and unconstrained, where the constraint simply shifts the free energy base
# line by -1 kcal/mol per nucleotide.
# When comparing both results, equilibrium probabilities must not have changed,
# except for free energy of the ensemble!
print "test_pfl_sc\n";
$fc = new RNA::fold_compound($longseq, undef, RNA::OPTION_WINDOW);

@data = ();
# unconstrained partition function
$fc->probs_window(0, RNA::PROBS_WINDOW_BPP, \&bpp_callback, \@data);

# add constraints
foreach my $i (1..length($longseq)) {
  $fc->sc_add_up($i, -1.0, RNA::OPTION_WINDOW);
}

foreach my $i (1..length($longseq)) {
  foreach my $j (($i + 1)..length($longseq)) {
    $fc->sc_add_bp($i, $j, -2, RNA::OPTION_WINDOW);
  }
}

my @data2 = ();
# constrained partition function
$fc->probs_window(0, RNA::PROBS_WINDOW_BPP, \&bpp_callback, \@data2);

$equal_probabilities = 1;
foreach my $d (@data) {
  my ($item) = grep {
    ($_->{'i'} == $d->{'i'}) && ($_->{'j'} == $d->{'j'})
  } @data2;
  ($equal_probabilities = 0, last) if ! defined($item);
  ($equal_probabilities = 0, last) if ! (sprintf("%g", $d->{'p'}) eq sprintf("%g", $item->{'p'}));
}
ok($equal_probabilities == 1);


# Compute unpaired probabilities both, constrained and unconstrained, where the
# constraint simply shifts the free energy base line by -1 kcal/mol per nucleotide.
# When comparing both results, equilibrium probabilities must not have changed!
print "test_probs_window_up\n";

my $ulength = 45;

$md = new RNA::md();
$md->{'max_bp_span'} = 150;
$md->{'window_size'} = 200;

$fc = new RNA::fold_compound($longseq, $md, RNA::OPTION_WINDOW);

@data = ();

# unconstrained unpaired probabilities
$fc->probs_window($ulength, RNA::PROBS_WINDOW_UP, \&up_callback, \@data);

# add constraints
foreach my $i (1..length($longseq)) {
  $fc->sc_add_up($i, -1.0, RNA::OPTION_WINDOW);
}

foreach my $i (1..length($longseq)) {
  foreach my $j (($i + 1)..length($longseq)) {
    $fc->sc_add_bp($i, $j, -2, RNA::OPTION_WINDOW);
  }
}

@data2 = ();

# constrained unpaired probabilities
$fc->probs_window($ulength, RNA::PROBS_WINDOW_UP, \&up_callback, \@data2);

$equal_probabilities = 1;
foreach my $i (1..length($longseq)) {
  foreach my $u (1..$ulength) {
    last if $i - $u + 1 <= 0;
    ($equal_probabilities = 0, last) if ! (sprintf("%g", $data[$i - 1]->{'up'}[$u]) eq sprintf("%g", $data2[$i - 1]->{'up'}[$u]));
  }
}
ok($equal_probabilities == 1);
