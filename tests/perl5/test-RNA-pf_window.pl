#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)

use strict;
use warnings;
use Test::More tests => 37;
use Data::Dumper;
use FileHandle;

use RNA;
use RNAHelpers qw(:Paths :Messages);


my $datadir = getDataDirPath();

my $kT = 0.61632077549999997;
# maximum allowed difference beteen compared probabilties
my $allowed_diff = 1e-7;
# some sequences to work with
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


sub up_split_callback {
  my ($v, $v_size, $i, $maxsize, $what, $data) = @_;

  if ($what & RNA::PROBS_WINDOW_UP) {
    $what = $what & ~RNA::PROBS_WINDOW_UP;
    # Non-split case:
    $data->[$i] = $v if $what == RNA::ANY_LOOP;
    # all the cases where probability is split into different loop contexts
    $data->{'ext'}->[$i] = $v if $what == RNA::EXT_LOOP;
    $data->{'hp'}->[$i] = $v if $what == RNA::HP_LOOP;
    $data->{'int'}->[$i] = $v if $what == RNA::INT_LOOP;
    $data->{'mb'}->[$i] = $v if $what == RNA::MB_LOOP;
  }
}


MsgChecking("whether pfl_fold() returns reasonable results");
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


MsgChecking("whether pfl_fold_up() returns reasonable results");
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


MsgChecking("whether we can do the same with a single call to probs_window()");
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


MsgChecking("whether probs_window() with SHAPE reactivity data and full-length window does the same as regular pf()");
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
MsgChecking("whether probs_window() with full-length window does the same as regular pf()");
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
MsgChecking("whether lowering the entire ensemble by -1.0 kcal/mol per nucleotide using soft constraints doesn't change equilibrium pairing probabilities");
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
MsgChecking("whether lowering the entire ensemble by -1.0 kcal/mol per nucleotide using soft constraints doesn't change equilibrium unpaired probabilities");
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


# Compute unpaired probabilities both, split into different loop contexts and full probability
# for all loop contexts. This check verifies that the sum of the individual loop types actually
# adds up to that obtained for any loop
MsgChecking("whether individual loop context unpaired probabilities sum up properly");
%data = ( 'ext' => [], 'hp' => [], 'int' => [], 'mb' => [] );
@data = ();

$md = new RNA::md();
$md->{'max_bp_span'} = 150;
$md->{'window_size'} = 200;

$fc = new RNA::fold_compound($longseq, $md, RNA::OPTION_WINDOW);
$fc->probs_window($ulength, RNA::PROBS_WINDOW_UP | RNA::PROBS_WINDOW_UP_SPLIT, \&up_split_callback, \%data);
$fc->probs_window($ulength, RNA::PROBS_WINDOW_UP, \&up_split_callback, \@data);

# sum up unpaired probabilities of individual loop contexts
my $max_diff = 0;
foreach my $i (1..length($longseq)) {
  foreach my $u (1..$ulength) {
    next if !defined($data{'ext'}->[$i][$u]);
    my $diff = abs($data{'ext'}->[$i][$u] +
                   $data{'hp'}->[$i][$u] +
                   $data{'int'}->[$i][$u] +
                   $data{'mb'}->[$i][$u] -
                   $data[$i][$u]);

    $max_diff = $diff if $diff > $max_diff;
  }
}
ok($max_diff < $allowed_diff);


MsgChecking("whether unpaired probabilities for u = 1 are identical to what we get using regular partition function method");
@data = ();
$ulength = 1;
$md = new RNA::md();
$md->{'max_bp_span'} = 150;
$md->{'window_size'} = 200;

$fc = new RNA::fold_compound($shortseq, $md, RNA::OPTION_WINDOW);
$fc->probs_window($ulength, RNA::PROBS_WINDOW_UP, \&up_split_callback, \@data);

# compute pairing probability and from that
# unpaired probability for individual nucleotides
$fc2 = new RNA::fold_compound($shortseq);
$fc2->pf();
$bpp = $fc2->bpp();

my @p_single = (0) x (length($shortseq) + 1);

foreach my $i (1..length($shortseq)) {
  foreach my $j ($i..length($shortseq)) {
    $p_single[$i] += $bpp->[$i][$j];
    $p_single[$j] += $bpp->[$i][$j];
  }
}

$max_diff = 0;
foreach my $i (1..length($shortseq)) {
  my $diff = abs((1. - $p_single[$i]) - $data[$i][1]);
  $max_diff = $diff if $diff > $max_diff;
}
ok($max_diff < $allowed_diff);


MsgChecking("whether unpaired probabilties for short RNAs are identical to what we get from exhaustive enumeration");
RNA::init_rand();
my $randseq_length = 20;
$ulength = 5;
$md = new RNA::md();
$md->{'max_bp_span'} = $randseq_length;
$md->{'window_size'} = $randseq_length;
# turn off dangles due to scaling issues that prevent us
# from comparing against exhaustive results
$md->{'dangles'} = 0;

# repreat the test 3 times
foreach my $trial (0..2) {
  %data = ( 'ext' => [], 'hp' => [], 'int' => [], 'mb' => [] );
  my $randseq = RNA::random_string($randseq_length, "ACGU");

  $fc = new RNA::fold_compound($randseq, $md, RNA::OPTION_WINDOW);
  $fc->probs_window($ulength, RNA::PROBS_WINDOW_UP | RNA::PROBS_WINDOW_UP_SPLIT, \&up_split_callback, \%data);

  # to not run into problems with undef values, fill up with 0 everything thats missing
  foreach my $i (1..$randseq_length) {
    foreach my $u (1..$ulength) {
      $data{'ext'}->[$i]->[$u] = 0 if ! defined $data{'ext'}->[$i]->[$u];
      $data{'hp'}->[$i]->[$u] = 0 if ! defined $data{'hp'}->[$i]->[$u];
      $data{'int'}->[$i]->[$u] = 0 if ! defined $data{'int'}->[$i]->[$u];
      $data{'mb'}->[$i]->[$u] = 0 if ! defined $data{'mb'}->[$i]->[$u];
    }
  }

  # compute unpaired probabilities from exhaustive enumeration
  my $pf = 0.0;
  my @pu_ext = map { [ (0) x ($ulength + 1) ] } (0..$randseq_length);
  my @pu_hp = map { [ (0) x ($ulength + 1) ] } (0..$randseq_length);
  my @pu_int = map { [ (0) x ($ulength + 1) ] } (0..$randseq_length);
  my @pu_mb = map { [ (0) x ($ulength + 1) ] } (0..$randseq_length);

  # 1st, compute partition functions for loop contexts using subopt()
  $fc = new RNA::fold_compound($randseq, $md);
  foreach my $s (@{$fc->subopt(5000)}) {
    $pf += exp(-$s->{energy} / $kT);
    my $ss = RNA::db_to_element_string($s->{structure});
    foreach my $u (1..$ulength) {
      foreach my $j (1..$randseq_length) {
        next if $j - $u < 0;
        my $sub = substr($ss, $j - $u, $u), "\n";
        $pu_ext[$j][$u] += exp(-$s->{energy}/$kT) if $sub eq ("e" x $u);
        $pu_hp[$j][$u] += exp(-$s->{energy}/$kT) if $sub eq ("h" x $u);
        $pu_int[$j][$u] += exp(-$s->{energy}/$kT) if $sub eq ("i" x $u);
        $pu_mb[$j][$u] += exp(-$s->{energy}/$kT) if $sub eq ("m" x $u);
      }
    }
  }

  # 2nd, get unpaired probabilities
  foreach my $i (1..$randseq_length) {
    foreach my $u (1..$ulength) {
      $pu_ext[$i][$u] = (defined $pu_ext[$i][$u]) ? $pu_ext[$i][$u] / $pf : 0.;
      $pu_hp[$i][$u]  = (defined $pu_hp[$i][$u])  ? $pu_hp[$i][$u] / $pf  : 0.;
      $pu_int[$i][$u] = (defined $pu_int[$i][$u]) ? $pu_int[$i][$u] / $pf : 0.;
      $pu_mb[$i][$u]  = (defined $pu_mb[$i][$u])  ? $pu_mb[$i][$u] / $pf  : 0.;
    }
  }

  # 3rd, compute maximum difference for each position
  my @diff_ext = (0) x ($randseq_length + 1);
  my @diff_hp = (0) x ($randseq_length + 1);
  my @diff_int = (0) x ($randseq_length + 1);
  my @diff_mb = (0) x ($randseq_length + 1);

  foreach my $i (1..$randseq_length) {
    foreach my $u (1..$ulength) {
      my $diff;
      $diff = abs($pu_ext[$i]->[$u] - $data{'ext'}->[$i]->[$u]);
      $diff_ext[$i] = $diff if $diff > $diff_ext[$i];
      $diff = abs($pu_hp[$i]->[$u] - $data{'hp'}->[$i]->[$u]);
      $diff_hp[$i] = $diff if $diff > $diff_hp[$i];
      $diff = abs($pu_int[$i]->[$u] - $data{'int'}->[$i]->[$u]);
      $diff_int[$i] = $diff if $diff > $diff_int[$i];
      $diff = abs($pu_mb[$i]->[$u] - $data{'mb'}->[$i]->[$u]);
      $diff_mb[$i] = $diff if $diff > $diff_mb[$i];
    }
  }

  # Finally, compute average difference
  my $diff_ext_avg = 0;
  my $diff_hp_avg = 0;
  my $diff_int_avg = 0;
  my $diff_mb_avg = 0;
  $diff_ext_avg += $_ foreach @diff_ext;
  $diff_hp_avg += $_ foreach @diff_hp;
  $diff_int_avg += $_ foreach @diff_int;
  $diff_mb_avg += $_ foreach @diff_mb;

  $diff_ext_avg /= $randseq_length;
  $diff_hp_avg /= $randseq_length;
  $diff_int_avg /= $randseq_length;
  $diff_mb_avg /= $randseq_length;

  # compute variance
  my $diff_ext_var = 0;
  my $diff_hp_var = 0;
  my $diff_int_var = 0;
  my $diff_mb_var = 0;
  $diff_ext_var += ($_ - $diff_ext_avg)**2 foreach @diff_ext;
  $diff_hp_var += ($_ - $diff_hp_avg)**2 foreach @diff_hp;
  $diff_int_var += ($_ - $diff_int_avg)**2 foreach @diff_int;
  $diff_mb_var += ($_ - $diff_mb_avg)**2 foreach @diff_mb;

  printf("\tTrial %d - Difference between pflfold unpaired probs and exhaustive enumeration\n", $trial);
  printf("\tExterior loop    (avg, var)\t%g\t%g\n", $diff_ext_avg, $diff_ext_var);
  printf("\tHairpin loop     (avg, var)\t%g\t%g\n", $diff_hp_avg, $diff_hp_var);
  printf("\tInterior loop    (avg, var)\t%g\t%g\n", $diff_int_avg, $diff_int_var);
  printf("\tMultibranch loop (avg, var)\t%g\t%g\n", $diff_mb_avg, $diff_mb_var);

  ok($diff_ext_avg < $allowed_diff);
  ok($diff_hp_avg < $allowed_diff);
  ok($diff_int_avg < $allowed_diff);
  ok($diff_mb_avg < $allowed_diff);
}
