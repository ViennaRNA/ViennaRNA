#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 4;

use RNA;
use warnings;
use RNAHelpers qw(:Paths :Messages);


my $datadir = getDataDirPath();

my $seq_con     = "CCCAAAAGGGCCCAAAAGGG";
my $str_con     = "..........(((....)))";
my $str_con_def = "(((....)))(((....)))";
my $seq_long    = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUC";
my $fc; #fold compound reference
my $ss=""; #return string
my $mfe=0;
my $energy=0;
my @data = ();
my $d = 0;
my $failed = 0;
my $e   = 0;
my $s = "";
my $c = 0;

sub mfe_window_callback {

  my ($start, $end, $structure, $energy, $data) = @_;
  my %Lfold_hit = ();
  $Lfold_hit{'structure'} = $structure;
  $Lfold_hit{'start'}     = $start;
  $Lfold_hit{'end'}       = $end;
  $Lfold_hit{'energy'}    = $energy;

  push @{$data}, \%Lfold_hit;
}



##################################
##test_sc_set_up

MsgChecking("whether setting unpaired soft constraints to all nucleotides at once works");
my $seq_sc  =      "CCCAAAAGGG";
$fc         = new RNA::fold_compound($seq_sc);
($ss, $mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,"(((....)))");

$fc->sc_init();

my @m= (0.0,-5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);  #E 1 0 1 -5 ,  position 1 gets -5 if unpaired , vector starts with 0 and not 1

$fc->sc_set_up(\@m);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n", $ss, $mfe);
is(sprintf ("%6.2f", $mfe), sprintf("%6.2f", -5.70));
undef @m;
##################################   We have no argin typemap for multidimensional lists yet, so this doesn't work
#test_sc_set_bp
#MsgChecking("whether setting base pair soft constraints to all nucleotides at once works");
#
#@m = ();
#push @m, [(0) x (length($seq_sc) + 1)] for 0 .. length($seq_sc);
#
# add energy of -5 to basepair 1-9 if formed, prefered structure should now be ((.....))., with a energy of -4.90
#$m[1][9] = -5.0; # base 1-9 should get -5.0 if basepair
#$m[9][1] = -5.0; # base 1-9 should get -5.0 if basepair
#
#$fc = new RNA::fold_compound($seq_sc);
#$fc->sc_set_bp(\@m);
#($ss, $mfe) = $fc->mfe();
#print("%s [%6.2f]\n", $ss, $mfe);
#is(sprintf ("%6.2f", $mfe), sprintf("%6.2f", -4.90));


##################################
##test_sc_mfe_window_add_bp

MsgChecking("whether adding base pair soft constraints works in conjunction with sliding-window MFE prediction");

$fc = new RNA::fold_compound($seq_long, undef, RNA::OPTION_WINDOW);

# add twice -10.0 kcal/mol on base pair (55,60)
$fc->sc_add_bp(55, 60, -10.0, RNA::OPTION_WINDOW);
$fc->sc_add_bp(55, 60, -10.0, RNA::OPTION_WINDOW);

# allow pair (55,60) to actually form (might be non-canonical)
$fc->hc_add_bp(55, 60, RNA::CONSTRAINT_CONTEXT_NO_REMOVE | RNA::CONSTRAINT_CONTEXT_ALL_LOOPS);

@data = ();
$mfe = $fc->mfe_window_cb(\&mfe_window_callback, \@data);

$failed = 0;
foreach my $hit (@data) {
    if ($hit->{'start'} <= 55 && $hit->{'end'} >= 60) {
        # compose actual dot bracket string (including potential 5' dangle nucleotide
        if ($hit->{'start'} > 1) {
          $d = 1;
          $ss = ".".$hit->{'structure'};
        } else {
          $d = 0;
          $ss = $hit->{'structure'};
        }

        # get corresponding subsequence
        $s = substr($seq_long, $hit->{'start'} - 1 - $d, $hit->{'end'} - $hit->{'start'} + 2);

        # re-evaluate free energy of subsequence/hit
        $e = RNA::energy_of_struct($s, $ss);

        # energy difference between both must be -20.0 kcal/mol (if the constrained base pair is present)
        $failed = 1 if ! (sprintf("%6.2f", $hit->{'energy'}) eq sprintf("%6.2f", $e - 20.0));
    }
}
ok($failed == 0);

##################################
##test_sc_mfe_window_add_bp

MsgChecking("whether adding unpaired soft constraints works in conjunction with sliding-window MFE prediction");

$fc = new RNA::fold_compound($seq_long, undef, RNA::OPTION_WINDOW);

# add twice -5.0 kcal/mol per unpaired nucleotide in segment [55,60]
for my $i (50..60) {
    $fc->sc_add_up($i, -5.0, RNA::OPTION_WINDOW);
}

# force segment [50,60] to stay unpaired
for my $i (50..60) {
    $fc->hc_add_up($i, RNA::CONSTRAINT_CONTEXT_ALL_LOOPS);
}

@data = ();
$mfe = $fc->mfe_window_cb(\&mfe_window_callback, \@data);

$failed = 0;

foreach my $hit (@data) {
    if ($hit->{'start'} <= 55 && $hit->{'end'} >= 60) {
        # count number of unpiared nucleotides with bonus
        $c = 0;
        foreach my $i (50..60) {
            $c++ if substr($hit->{'structure'}, $i - $hit->{'start'}, 1) eq ".";
        }

        # compose actual dot bracket string (including potential 5' dangle nucleotide
        if ($hit->{'start'} > 1) {
          $d = 1;
          $ss = ".".$hit->{'structure'};
        } else {
          $d = 0;
          $ss = $hit->{'structure'};
        }

        # get corresponding subsequence
        $s = substr($seq_long, $hit->{'start'} - 1 - $d, $hit->{'end'} - $hit->{'start'} + 1 + $d);

        # re-evaluate free energy of subsequence/hit
        $e = RNA::energy_of_struct($s, $ss);

        # energy difference between both must be c * -5.0 kcal/mol
        $failed = 1 if ! (sprintf("%6.2f", $hit->{'energy'}) eq sprintf("%6.2f", $e + ($c * -5.0)));
    }
}
ok($failed == 0);

undef $fc;
