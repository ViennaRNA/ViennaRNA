#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 11;
use Data::Dumper;
use FileHandle;

use RNA;
use warnings;


my $seq1  =                  "CGCAGGGAUACCCGCG";
my $struct1=                "(((.(((...))))))";
my $seq1Dimer =         "CGCAGGGA&ACCCGCG";
my $struct1Dimer=        "(((.(((..))))))";
my $struct11 =                 "(((.((.....)))))";
my $struct1_pt =         [length($struct1),16,15,14,0,13,12,11,0,0,0,7,6,5,3,2,1];
my $fc; #fold compound reference
my $ss=""; #return string
my $mfe=0;
my $energy=0;

my $outFile="test-RNA-mfe_eval.pl.out";
my $fh;

##################################
print "test_mfe\n";
$fc= new RNA::fold_compound($seq1);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,$struct1);
undef $fc;
##################################
print "test_mfe_Dimer\n";
$fc= new RNA::fold_compound($seq1Dimer);
($ss,$mfe) = $fc->mfe_dimer();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,$struct1Dimer);
undef $fc;
##################################
print "test_eval_structure\n";
$fc = new RNA::fold_compound($seq1);
$energy= $fc->eval_structure($struct1);
printf("%s [%6.2f] \n",$struct1,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",-5.60));
undef $fc;
##################################

print "test_eval_structure_pt\n";
$fc=new RNA::fold_compound($seq1);
$energy= $fc->eval_structure_pt($struct1_pt) / 100.; #/100 for dcal
printf("%s [%6.2f] \n",$struct1,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",-5.60));
undef $fc;
##################################
print "test_eval_structure_verbose";
$fc=new RNA::fold_compound($seq1);
open($fh, ">", $outFile) || die "Couldn't open '".$outFile."' for reading because: ".$!;
$energy = $fc->eval_structure_verbose($struct1,$fh);
$energy = $fc->eval_structure_verbose($struct1,undef); # print to stdout
printf("%s [%6.2f] \n",$struct1,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",-5.60));
undef $fc;
undef $fh;
##################################
print "test_eval_structure_pt_verbose\n";
open($fh, ">", $outFile) || die "Couldn't open '".$outFile."' for reading because: ".$!;
$fc=new RNA::fold_compound($seq1);

$energy = $fc->eval_structure_pt_verbose($struct1_pt,$fh) / 100.;  # / 100 for dcal
$energy = $fc->eval_structure_pt_verbose($struct1_pt,undef) / 100.;  # / 100 for dcal, print to stdout
printf("%s [%6.2f] \n",$struct1,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",-5.60));
undef $fc;
undef $fh;
##################################
print "test_eval_covar_structure\n";
my $s1="CCCCAAAACGGG";
my $s2="CCCGAAAAGGGG";
my $s3="CCCCAAAAGGGG";
my @ali = ($s1,$s2,$s3);
my $covarStructure = "((((....))))";

$fc = new RNA::fold_compound(\@ali);
my $ps_energy=$fc->eval_covar_structure($covarStructure);
printf("%s [%6.2f] \n",$covarStructure,$ps_energy);
ok($ps_energy != 0);
undef $fc;
##################################
print "test_eval_loop_pt\n";
$fc= new RNA::fold_compound($seq1);
$energy= $fc->eval_loop_pt(6,$struct1_pt) / 100.; #/100 for dcal
printf("%s [%6.2f] \n",$struct1,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",-3.3));
undef $fc;
# ##################################
#
print "test_eval_move_del\n";
$fc = new RNA::fold_compound($seq1);
$energy = $fc->eval_move($struct1,-7,-11);  # remove basepair (6,11) ,  $energy change should be 2.10
printf("%s moveset (-7,-11) --> [%6.2f] \n",$struct1,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",2.10));
# ##################################
#
print "test_eval_move_ins\n";
$fc = new RNA::fold_compound($seq1);
$energy = $fc->eval_move($struct11,7,11);  # add basepair (6,11) ,  $energy change should be -2.10
printf("%s moveset (7,11) --> [%6.2f] \n",$struct11,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",-2.10));
# ##################################
#
print "test_eval_move_pt_del\n";
$fc = new RNA::fold_compound($seq1);
$energy = $fc->eval_move_pt($struct1_pt,-7,-11) /100;  # remove basepair (6,11) ,  $energy change should be 2.10
printf("%s moveset (-7,-11) --> [%6.2f] \n",$struct1,$energy);
is(sprintf("%6.2f",$energy), sprintf("%6.2f",2.10));
##################################


