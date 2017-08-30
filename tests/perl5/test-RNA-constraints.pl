#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 8;

use RNA;
use warnings;
use RNAHelpers qw(:Paths);


my $datadir = getDataDirPath();

my $seq_con     = "CCCAAAAGGGCCCAAAAGGG";
my $str_con     = "..........(((....)))";
my $str_con_def = "(((....)))(((....)))";
my $seq_long    = "AUUUCCACUAGAGAAGGUCUAGAGUGUUUGUCGUUUGUCAGAAGUCCCUAUUCCAGGUACGAACACGGUGGAUAUGUUCGACGACAGGAUCGGCGCACUACGUUGGUAUCAUGUCCUCCGUCCUAACAAUUAUACAUCGAGAGGCAAAAUUUCUAAUCCGGGGUCAGUGAGCAUUGCCAUUUUAUAACUCGUGAUCUCUC";
my $fc; #fold compound reference
my $ss=""; #return string
my $mfe=0;
my $energy=0;


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
## test_constraints_add
##################################

my $hc_file = $datadir . "hc.txt";
print "test_constraints_add";
$fc = new RNA::fold_compound($seq_con);
$fc->constraints_add($hc_file);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,$str_con);

$fc->hc_init();
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,$str_con_def);

##################################
#sc.txt = E 3 8 1 -5
##################################
my $sc_file = $datadir . "sc.txt";
$fc->sc_init();
$fc->constraints_add($sc_file);
($ss,my $mfeNew) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfeNew);
is(sprintf ("%6.2f",$mfe), sprintf ("%6.2f",$mfeNew+5));

##################################
##test_hc_add_up:
##################################
print "test_hc_add_up\n";

$fc = new RNA::fold_compound($seq_con);
$fc->hc_add_up(1,RNA::CONSTRAINT_CONTEXT_ALL_LOOPS);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,".((....)).(((....)))");

##################################
## test_hc_add_bp_nonspecific
##################################
print "test_hc_add_bp_nonspecific";

$fc= new RNA::fold_compound("GGGCCCCCCCCCCCCCCCCC");
$fc->hc_add_bp_nonspecific(20,-1); # force the last base to pair with some bases upstream
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,"(((..............)))");

##################################
## test_hc_add_bp
##################################
print "test_hc_add_bp";

$fc= new RNA::fold_compound($seq_con);
$fc->hc_add_bp(1,20,RNA::CONSTRAINT_CONTEXT_ENFORCE | RNA::CONSTRAINT_CONTEXT_ALL_LOOPS);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,"(((..............)))");

##################################
##  test_hc_add_from_db
##################################
print "test_hc_add_from_db";

$fc = new RNA::fold_compound($seq_con);
$fc->hc_add_from_db("xxx.................");
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,$str_con);

##################################
##  test_hc_mfe_window (base pairs)
##################################
print "test hc_mfe_window_bp\n";

$fc = new RNA::fold_compound($seq_long, undef, RNA::OPTION_WINDOW);
$fc->hc_add_bp(1, 10, RNA::CONSTRAINT_CONTEXT_ALL_LOOPS | RNA::CONSTRAINT_CONTEXT_ENFORCE);
$fc->hc_add_bp(101, 110, RNA::CONSTRAINT_CONTEXT_ALL_LOOPS | RNA::CONSTRAINT_CONTEXT_ENFORCE);
my @data = ();
$mfe = $fc->mfe_window_cb(\&mfe_window_callback, \@data);

my $everythingFine = 1;
foreach my $hit (@data) {
    if (($hit->{'start'} <= 101) && ($hit->{'end'} >= 110)) {
        # must contain base pair (101,110)
        my $pt = RNA::ptable($hit->{'structure'});
        $everythingFine = 0 if $pt->[101 - $hit->{'start'} + 1] != (110 - $hit->{'start'} + 1);
    }
    if (($hit->{'start'} == 1) && ($hit->{'end'} >= 10)) {
        # must contain base pair (101,110)
        # must contain base pair (101,110)
        my $pt = RNA::ptable($hit->{'structure'});
        $everythingFine = 0 if $pt->[1] != 10;
    }
}

ok($everythingFine == 1);

undef $fc;
