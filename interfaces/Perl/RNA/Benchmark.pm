package RNA::Benchmark;

use strict;
use Exporter;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(CompareStructures
                  TruePositiveRate
                  PositivePredictiveValue
                  MathewsCorrelationCoefficient
                  FMeasure);

%EXPORT_TAGS = (  ALL => [qw(&CompareStructures &TruePositiveRate &PositivePredictiveValue &MathewsCorrelationCoefficient &FMeasure)],
                  MEASURES => [qw(&TruePositiveRate &PositivePredictiveValue &MathewsCorrelationCoefficient &FMeasure)],
                  MATCHERS => [qw(&CompareStructures)]);

=head1 AUTHOR

Ronny Lorenz (ronny@tbi.univie.ac.at)

=head1 NAME

RNA::Benchmark - A set of subroutines to benchmark RNA secondary structure
prediction performance

=head1 DESCRIPTION

This package provides various subroutines for benchmarking the prediction
performance of RNA secondary structure prediction tools.
The basic data structure for comparing secondary structures is a pair table,
which can be generated from dot-bracket annotation via the make_pair_table()
subroutine of the RNA::Utils package.

=head1 METHODS

=cut


=head2 CompareStructures($pt_gold, $pt_other, $fuzzy = 0, $verbose = 0)

Compare a secondary structure to a 'gold standard' to assess the number of
false positives (FP) and false negatives (FN)

Takes a trusted secondary structure ($pt_gold), and a predicted secondary
structure ($pt_other) provided as references to pair_tables to compute the
number of base pairs in the trusted structure (true positives, or TP), the
number of base pairs in the predicted structure that are not in the trusted
structure (false positives, or FP), and the number of base pairs in the
trusted structure but not in the predicted one (false negatives, or FN).
The maximal number of base pairs that are not part of the trusted structure,
the true negatives (TN), are estimates as n*(n-1)/2 for an RNA sequence of
length n.

The optional parameter $fuzzy allows for base pair slippage. Hence, for any
base pair (i,j) in the gold standard set of pairs ($pt_gold), a base pair
(p, q) in the second structure ($pt_other) is considered a true positive, if
i - $fuzzy <= p <= i + $fuzzy, and j - $fuzzy <= q <= j + $fuzzy.

The optional parameter $verbose is mainly for debugging purposes, but may be
used to produce progress output of the comparison.

The subroutine returns a set of four values ($TP, $FP, $TN, $FN), representing
the true positives $TP, the true negatives $TN, the false negatives $FN, and
the false positives $FP.

=cut

sub CompareStructures{
  my $pt_gold   = shift;
  my $pt_other  = shift;
  my $fuzzy     = shift;
  my $verbose   = shift;

  $fuzzy    = 0 if !defined($fuzzy);
  $verbose  = 0 if !defined($verbose);

  my ($TP, $TN, $FN, $FP) = (0, 0, 0, 0);

  # number of false positive, but compatible pairs
  my $compatible = 0;

  # count the number of base pairs in $pt_gold, and $pt_other
  my $bps_gold  = 0;
  my $bps_other = 0;

  for my $i (1 .. $pt_gold->[0] ){
    next if $pt_gold->[$i] < $i;
    $bps_gold++;
  }
  for my $i (1 .. $pt_other->[0] ){
    next if $pt_other->[$i] < $i;
    $bps_other++;
  }

  # got through $pt_other and compare its pairs to $pt_gold
  for my $i ( 1 .. $pt_other->[0] ) {
    next if ( $pt_other->[$i] < $i );

    my $j                = $pt_other->[$i];
    my $is_true_positive = 0;
    my $is_inconsistent  = 0;
    my $is_contradicting = 0;

    #############################################################
    # TRUE POSITIVES
    #############################################################

    # let's see if position i matches in gold position j +/- fuzzy
    for my $add ( 0 .. $fuzzy ) {
      if ( $pt_gold->[$i] == $j + $add and $j + $add <= $pt_other->[0] ) {
        $is_true_positive = 1;
        $bps_gold--;
        print STDERR "BP ($i, $j) in Reference matches BP ($i, ", ( $j + $add ), ") in Prediction.\n" if ($verbose);
        last;
      }
      if ( $pt_gold->[$i] == $j - $add and $j - $add >= 1 ) {
        $is_true_positive = 1;
        $bps_gold--;
        print STDERR "BP ($i, $j) in Reference matches BP ($i, ", ( $j - $add ), ") in Prediction.\n" if ($verbose);
        last;
      }
    }

    # let's see if position j matches in gold position i +/- fuzzy
    if ( $is_true_positive == 0 ) {
      for my $add ( 0 .. $fuzzy ) {
        if ( $pt_gold->[$j] == $i + $add and $i + $add <= $pt_other->[0] ) {
          $is_true_positive = 1;
          $bps_gold--;
          print STDERR "BP ($i, $j) in Reference matches BP ($j, ", ( $i + $add ), ") in Prediction.\n" if ($verbose);
          last;
        }
        if ( $pt_gold->[$j] == $i - $add and $i - $add >= 1 ) {
          $is_true_positive = 1;
          $bps_gold--;
          print STDERR "BP ($i, $j) in Reference matches BP ($j, ", ( $i - $add ), ") in Prediction.\n" if ($verbose);
          last;
        }
      }
    }

    #############################################################
    # FALSE POSITIVES
    #############################################################

    # let's check if base-pair i * j is inconsistent with reference structure
    # base i or base j in reference is paired to something else
    my $tmp_string1 = '';
    if ( $is_true_positive == 0 ) {
      if ( $pt_gold->[$i] != 0 ) {
        $is_inconsistent = 1;
        $tmp_string1     = "Reference: ($i, $pt_gold->[$i]) ";
      }
      if ( $pt_gold->[$j] != 0 ) {
        $is_inconsistent = 1;
        $tmp_string1 .= "Reference: ($j, $pt_gold->[$j]) ";
      }
    }

    # let's check if a base-pair is contradicting
    # there is a pair k * l in ther reference that k < i < l < j
    my $tmp_string2 = '';
    for my $k ( 1 .. $pt_gold->[0] ) {
      next if ( $pt_gold->[$k] < $k);
      my $l = $pt_gold->[$k];
      if ( $k < $i and $i < $l and $l < $j ) {
        $tmp_string2      = "crossing pairs: $k < $i < $l < $j";
        $is_contradicting = 1;
        last;
      }
    }

    $TP++ if ( $is_true_positive == 1 );

    if ($is_true_positive == 0){
      if ( $is_inconsistent == 1 or $is_contradicting == 1 ) {
        print STDERR "BP ($i, $j) is a false positive." if ( $verbose );
        print STDERR "\t$tmp_string1" if ( $verbose and $tmp_string1 ne '');
        print STDERR "\t$tmp_string2"  if ( $verbose and $tmp_string2 ne '');
        print STDERR "\n" if ( $verbose );
      } else {
        print STDERR "BP ($i, $j) is a false positive, but compatible with Reference.\n" if ( $verbose );
        $compatible++;
      }
      $FP++;
    }
  }

  print "Reference contains ", $bps_gold, " BPs that were not matched to Prediction\n" if ( $verbose );

  $FN = $bps_gold;

  # make an educated guess for the number of true negatives
  $TN = 0.5 * ($pt_gold->[0] * ($pt_gold->[0] - 1));

  return ($TP, $FP, $TN, $FN, $compatible);
}

=head2 TruePositiveRate($TP, $FP, $FN)

Compute the true positive rate (TPR), a.k.a. sensitivity

Takes the number of true positives ($TP), false positives ($FP),
and false negatives ($FN) to return the True Positive Rate.

=cut

sub TruePositiveRate{
  my ($TP, $FP, $FN) = @_;
  return $TP / ($TP + $FN) if ($TP + $FN > 0);
  return 0;
}

=head2 PositivePredictiveValue($TP, $FP, $FN)

Compute the positive predictive value (PPV), a.k.a. selectivity

Takes the number of true positives ($TP), false positives ($FP),
and false negatives ($FN) to return the Positve Predictive Value.

=cut

sub PositivePredictiveValue{
  my ($TP, $FP, $FN) = @_;
  return $TP / ($TP + $FP) if ($TP + $FP > 0);
  return 0;
}

=head2 MathewsCorrelationCoefficient($TP, $FP, $FN, $TN)

Compute the Mathews Correlation Coefficient (MCC)

Takes the number of true positives ($TP), false positives ($FP),
false negatives ($FN), and true negatives ($TN) to return the
Positve Predictive Value.

=cut

sub MathewsCorrelationCoefficient{
  my ($TP, $FP, $FN, $TN) = @_;
  return ($TP * $TN - $FP * $FN) / sqrt(($TP + $FP) * ($TP + $FN) * ($TN + $FP) * ($TN + $FN));
}

=head2 FMeasure($PPV, $TPR)

Compute the F1-Measure

Takes the positive predictive value ($PPV), and the true positive
rate ($TPR) to compute the F1-Measure.

=cut

sub FMeasure{
  my ($ppv, $tpr) = @_;
  return 2 * ($ppv * $tpr) / ($ppv + $tpr) if ($ppv + $tpr > 0);
  return 0;
}

=head2  Bootstrap($data, $intervals = [95], $n = 1000)

I<Compute confidence intervals via bootstrapping>

This subroutine computes a set of confidence intervals of mean
values using the bootstrapping method.

The input data B<$data> has to be a reference to an array of data points,
which themself may either be a I<SCALAR>, or a I<HASH> of scalars.
The additional, optional arguments B<$interval>, and B<$n> may be used to
specify the confidence intervals to compute, and the number of iterations
of the bootstrapping method, respectively. Here, B<$intervals> has to be
a reference to an array of scalar integer values in the range of [1 .. 99],
and B<$n> a positive integer number.

The subroutine computes the arithmetic mean of the input data and each
of the requested confidence intervals. The confidence interval data is
then stored in a I<HASH> with the following keys:

=over

=item B<INTERVAL> ... The interval (as given in the input arguments)

=item B<FROM> ... Lower bound of the confidence interval

=item B<TO> ... Upper bound of the confidence interval

=back

The function returns a pair of two values

B<($mean, $interval_data)>

whose data type entirely depends on the input data types passed to it.
While the data type of the entries in B<$data> affects both output values,
the number of intervals passed via the B<$intervals> argument influences
only the second variable of the returned pair:

=over

=item * I<$intervals is a single number>

B<$interval_data> is a I<HASH> reference pointing to the computed data.

=item * I<$intervals is an array of numbers>:

B<$interval_data> is an I<ARRAY> reference pointing to a list of
I<HASH>es that contain the individual confidence interval data.

=item * I<$data is an array of single numbers>:

B<$mean> is a single I<SCALAR> value representing the arithmetic mean
of the input data.

=item * I<$data is an array of hashes of numbers>:

B<$mean> is a I<HASH> reference having the same keys as the input data.
Each associated value represents the arithmetic mean of the particular
input data subset.

Furthermore, the keys B<FROM>, and B<TO> within the I<HASH> reference(s)
of B<$interval_data> do not represent I<SCALAR> values, but are I<HASH>
references themself, where the keys are, again, the same as provided by
the input data.

=back

=cut

sub Bootstrap{
  my $default_key = "some_very_unlikely_hash_key";
  my $data        = shift;
  my $intervals   = shift;
  my $bootstrap_n = shift;

  # first prepare the input data
  my @input_data = ();
  if(ref($data->[0]) eq 'HASH'){
    @input_data = @{$data};
  } else {
    foreach my $d (@{$data}){
      my %h;
      $h{$default_key} = $d;
      push @input_data, \%h;
    }
  }

  # second prepare the confidence intervals
  my @confidence_intervals;

  if(defined($intervals)){
    if(ref($intervals) eq 'ARRAY'){
      @confidence_intervals = @{$intervals};
    } elsif(int($intervals) == $intervals){
      push @confidence_intervals, $intervals;
    }
  }

  push @confidence_intervals, 95 if not @confidence_intervals;

  # last we assign a default number of repetitions if not specified directly
  $bootstrap_n  = 1000 if not defined($bootstrap_n);

  # now, lets start the fun!

  my $sample_n  = scalar(@input_data);
  my @keys      = keys(%{$input_data[0]});

  my %bootstrap_means;
  my %input_mean;

  # compute simple mean value
  foreach my $k (@keys) {
    for(my $j = 0; $j < $sample_n; $j++){
      $input_mean{$k} += $input_data[$j]->{$k};
    }
    $input_mean{$k} /= $sample_n;
  }

  # generate bootstrap distribution of sample means
  for(my $n = 0; $n < $bootstrap_n; $n++){
    for(my $j = 0; $j < $sample_n; $j++){
      my $pos = int(rand($sample_n - 1));
      foreach my $k (@keys) {
        if(!exists($bootstrap_means{$k})){
          my @a = ();
          $bootstrap_means{$k} = \@a;
        }
        if(defined($bootstrap_means{$k}->[$n])){
          $bootstrap_means{$k}->[$n] += $input_data[$pos]->{$k};
        } else {
          $bootstrap_means{$k}->[$n] = $input_data[$pos]->{$k};
        }
      }
    }
    foreach my $k (@keys) {
      $bootstrap_means{$k}->[$n] /= $sample_n;
    }
  }

  # sort the bootstrapped mean values
  my %sorted_bootstrap_means;
  foreach my $k (@keys) {
    my @sorted_a = sort(@{$bootstrap_means{$k}});
    $sorted_bootstrap_means{$k} = \@sorted_a;
  }

  # compose the output
  my @output = ();
  for(my $i = 0; $i < @confidence_intervals; $i++){
    my %dat;
    my $percentile = (100 - $confidence_intervals[$i]) / 2.;
    my $pos1 = int(($percentile / 100.) * $bootstrap_n - 0.5);
    my $pos2 = int(((100. - $percentile) / 100.) * $bootstrap_n - 0.5);

    $dat{INTERVAL}  = $confidence_intervals[$i];

    if((scalar(@keys) == 1) and ($keys[0] eq $default_key)){
      $dat{FROM}      = $sorted_bootstrap_means{$default_key}->[$pos1];
      $dat{TO}        = $sorted_bootstrap_means{$default_key}->[$pos2];
    } else {
      my %from;
      my %to;

      $dat{FROM} = \%from;
      $dat{TO} = \%to;

      foreach my $k (@keys) {
        $dat{FROM}->{$k} = $sorted_bootstrap_means{$k}->[$pos1];
        $dat{TO}->{$k}   = $sorted_bootstrap_means{$k}->[$pos2];
      }
    }
    push @output, \%dat;
  }

  my ($out_mean, $out);
  if((scalar(@keys) == 1) and ($keys[0] eq $default_key)){
    $out_mean = $input_mean{$default_key};
  } else {
    $out_mean = \%input_mean;
  }
  if(scalar(@confidence_intervals) == 1){
    $out = $output[0];
  } else {
    $out = \@output;
  }

  return ($out_mean, $out);
}

1;

=head1  SEE ALSO

=over

=item * RNA::Utils ... For converting secondary structures to pair_table format

=item * RNA::Files ... For parsing secondary structures from different input file formats

=back

=cut
