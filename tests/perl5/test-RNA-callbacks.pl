#!/usr/bin/perl
#

use strict;
use warnings;
use Data::Dumper;
use Test;
BEGIN { plan tests => 0; }

use RNA;

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

$a->add_auxdata(\%b, undef);
$a->add_callback(\&bla);
$a->sc_add_data(\%c, undef);
$a->sc_add_f(\&blubb);

my ($s,$mfe) = $a->mfe();
printf("%s %6.2f\n", $s, $mfe);

$a->DESTROY();
