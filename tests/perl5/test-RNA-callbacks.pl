#!/usr/bin/perl
#

use strict;
use warnings;
use Data::Dumper;
use Test::More tests => 1;
#BEGIN { plan tests => 0; }

use RNA;

is(1,1); # at least one test should be here
my $a = new RNA::fold_compound("GGGGAAAACCCC");

my %b = ('test' => "something");
my %c = ('what' => "theheck");

sub bla {
  my ($d,$data) = @_;
  print "about to start MFE recursions\n" if $d == RNA::STATUS_MFE_PRE;
  print "finished MFE recursions\n" if $d == RNA::STATUS_MFE_POST;
  print Dumper($data),"\n";
}

sub blubb {
  my ($i,$j,$k,$l,$d,$data) = @_;
  return -1000 if $d == RNA::DECOMP_PAIR_HP;
  return 0;
}

sub bt {
  my ($i,$j,$k,$l,$d,$data) = @_;
  my @bp;

  # The backtracking callback must return an array of base pairs
  # Here, the base pairs may be given in one of the three ways:
  # either as a hash, with keys 'i' and 'j', or an array @a with
  # $a[0] = i, $a[1] = j, or an object of type RNA::basepair
  #
  # below are examples for all the possibilities

  # 1. We create an array of hash references, where
  # each hash must have an 'i' and a 'j' key that specify the
  # coordinates of the base pair (i,j)
  my %bp;
  $bp{'i'} = $i+1;
  $bp{'j'} = $j-1;
  #push @bp, \%bp if $d == RNA::DECOMP_PAIR_HP;

  # 2. We create an array of array references, where
  # the base pair arrays have the i and j positions of
  # base pair (i,j) at position 0 and 1, respectively
  #push @bp, [($i+1), ($j-1)] if $d == RNA::DECOMP_PAIR_HP;

  # 3. We create an array of RNA::basepair objects
  my $pair = new RNA::basepair();
  $pair->{"i"} = $i+1;
  $pair->{"j"} = $j-1;
  push @bp, $pair if $d == RNA::DECOMP_PAIR_HP;

  return @bp;
}

$a->add_auxdata(\%b, undef);
$a->add_callback(\&bla);
$a->sc_add_data(\%c, undef);
$a->sc_add_f(\&blubb);
$a->sc_add_bt(\&bt);

my ($s,$mfe) = $a->mfe();
printf("%s %6.2f\n", $s, $mfe);

$a->DESTROY();

#
# test subopt callback
#
sub print_result{
  my $structure = shift;
  my $energy    = shift;
  printf("%s [%6.2f]\n", $structure, $energy) if defined($structure);
}

$RNA::uniq_ML = 1;  # for subopt calls
$a = new RNA::fold_compound("GGGGAAAACCCC");
$a->subopt_cb(500, \&print_result);


$a->DESTROY();
