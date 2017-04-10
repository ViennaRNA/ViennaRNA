#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 2;

use RNA;
use warnings;
use RNApath;


my $datadir = RNApath::getDataDirPath();

my $seq_con     = "CCCAAAAGGGCCCAAAAGGG";
my $str_con     = "..........(((....)))";
my $str_con_def = "(((....)))(((....)))";
my $fc; #fold compound reference
my $ss=""; #return string
my $mfe=0;
my $energy=0;



##################################
##test_sc_set_up

print "test_sc_set_up";
my $seq_sc  =      "CCCAAAAGGG";
$fc         = new RNA::fold_compound($seq_sc);
($ss, $mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,"(((....)))");

$fc->sc_init();

my @m= (0.0,-5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);  #E 1 0 1 -5 ,  position 1 gets -5 if unpaired , vector starts with 0 and not 1

$fc->sc_set_up(\@m);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is(sprintf ("%6.2f",$mfe), sprintf("%6.2f",-5.70));
undef @m;
##################################   !!!geht noch nicht
#test_sc_set_bp
print "test_sc_set_bp";

undef $fc;
