#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 1;
use Data::Dumper;
use FileHandle;

use RNA;
use warnings;
use RNAHelpers qw(:Paths);


my $datadir = getDataDirPath();

my $seq_con     = "CCCAAAAGGGCCCAAAAGGG";
my $str_con     = "..........(((....)))";
my $str_con_def = "(((....)))(((....)))";
my $fc; #fold compound reference
my $ss=""; #return string
my $mfe=0;
my $energy=0;


##test_sc_add_hi_motif
print "test_sc_add_hi_motif";
my $fc= new RNA::fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG");
#struct =          ".(((((..((((((((((((((((((....(((((((............)))))))........)))))))))))))...)))))))))).............."
#structWithMotif=  "((((...((((((((......)))))...)))...))))...((.((((((((((.((((((.((.((((....)))).))..)))))).)))))))))))).."
my $r=$fc->sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC","(...((((&)...)))...)",-9.22);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($r,1);

##################################
# test theophylline ligand binding interface
$fc = new RNA::fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG");
($ss, $mfe) = $fc->mfe();
printf "%s [ %6.2f ]\n", $ss, $mfe;

$fc->sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC", "(...((((&)...)))...)", -9.22);
($ss, $mfe) = $fc->mfe();
printf "%s [ %6.2f ]\n", $ss, $mfe;

$fc->sc_remove();

$fc->sc_add_hi_motif("GAAAAAU", "(.....)", -19);
($ss, $mfe) = $fc->mfe();
printf "%s [ %6.2f ]\n", $ss, $mfe;

undef $fc;
