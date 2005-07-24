#!/usr/bin/perl -Iblib/arch -Iblib/lib

# Last changed Time-stamp: <2005-07-23 16:04:56 ivo>

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test;
use lib qw|blib/arch blib/lib|;

BEGIN {  plan tests => 19; }

use RNA;
use warnings;

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

my $seq1  ="CGCAGGGAUACCCGCG";
my $struc1="(((.((....)).)))";
my $seq2  ="GCGCCCAUAGGGACGC";
my $struc2="((((((...))).)))";
# calculate a hamming distance (from util.c)
ok(RNA::hamming($seq1, $seq2), 16);

# check a global variable
ok($RNA::temperature, 37);

# fold a sequence

# old obsolete way of calling fold()
my $struct = $seq1;  # wierd way of allocating space
my $mfe=RNA::fold($seq1, $struct);
ok($struct, $struc1);

# new better interface
($struct, $mfe) = RNA::fold($seq1);
ok($struct eq $struc1);
# check energy
ok(RNA::energy_of_struct($seq1,$struc1), $mfe);

# check constrained folding
$RNA::fold_constrained = 1;
my($struct3, $cmfe) = RNA::fold($seq1, '....xx....xx...');
ok($struct3, '(((..........)))');
ok(RNA::energy_of_struct($seq1,$struct3), $cmfe);
$RNA::fold_constrained = 0;

# test cofold
$RNA::cut_point = length($seq1)+1;
my($costruct, $comfe) = RNA::cofold($seq1 . $seq2);
ok($costruct, '(((.((....)).)))((((((...))).)))');
$cmfe = RNA::energy_of_struct($seq1 . $seq2, $costruct);
ok(abs($comfe-$cmfe)<1e-5);
$RNA::cut_point=-1;

# pf_fold
my $f = RNA::pf_fold($seq1, $struct);
ok(($f<$mfe)&&($mfe-$f<0.8));

# tree distance
my $xstruc = RNA::expand_Full($struc1);
my $T1 = RNA::make_tree($xstruc);
$xstruc = RNA::expand_Full($struc2);
my $T2 = RNA::make_tree($xstruc);
$RNA::edit_backtrack = 1;
my $tree_dist = RNA::tree_edit_distance($T1, $T2); 
# print RNA::get_aligned_line(0), RNA::get_aligned_line(1),"\n";
ok($tree_dist,4);

# check access to a C array
#ok(RNA::ptrvalue($RNA::iindx,3),108);
ok(RNA::intP_getitem($RNA::iindx,3),108);

#RNA::shortP_setitem($RNA::xsubi, 0, 171);
#RNA::shortP_setitem($RNA::xsubi, 1, 42);
#RNA::shortP_setitem($RNA::xsubi, 2, 93);
RNA::memmove($RNA::xsubi, pack('S3', 171,42,93));
ok(pack('S3', 171,42,93), RNA::cdata($RNA::xsubi, 6));

# get a bp prob in two different ways
my $p1 = RNA::get_pr(2,15);
my $ii = RNA::intP_getitem($RNA::iindx, 2);
my $p2 = RNA::doubleP_getitem($RNA::pr, $ii-15);
ok(($p1<0.999) && ($p1>0.99) && (abs($p1-$p2)<1.2e-7));

my $bpf = RNA::Make_bp_profile(length($seq1));
my @bpf = unpack("f*",RNA::cdata($bpf, length($seq1)*4*3));
ok (($bpf[2*3]+$bpf[2*3+1]>.99999)&&$bpf[2*3+2]==0 &&
    ($bpf[2*3+1]>=$p1));

my $pack = RNA::pack_structure($struc1);
ok (RNA::unpack_structure($pack), $struc1);


RNA::parse_structure($struc1);
ok(($RNA::loops==2) && ($RNA::pairs==5)&&($RNA::unpaired==6) &&
  (RNA::intP_getitem($RNA::loop_degree,1)==2));


RNA::PS_rna_plot($seq1, $struc1, "test_ss.ps");
my $anote = "2 15 1 gmark\n" . "3 cmark\n";
RNA::PS_rna_plot_a($seq1, $struc1, "test_ss_a.ps", undef, $anote);
RNA::PS_dot_plot($seq1, "test_dp.ps");
RNA::ssv_rna_plot($seq1, $struct, "test.coord");
# print "$seq1, $struct, $mfe, $f\n";
print "please check the two postscript files test_ss.ps and test_dp.ps\n";
RNA::write_parameter_file("test.par");

$RNA::symbolset = "GC";
my $start = RNA::random_string(length $struc1, $RNA::symbolset);
my ($sinv, $cost) = RNA::inverse_fold($start, $struc1);
my ($ss, $en) = RNA::fold($sinv);
ok($ss, $struc1);

RNA::free_pf_arrays();
RNA::free_arrays();

$RNA::subopt_sorted = 1;
$RNA::noLonelyPairs = 1;
my $solution = RNA::subopt($seq1, undef, 500, undef);

printf "%d suboptimals\n", $solution->size();
for (0..$solution->size()) {
  # the last access should produce a "value out of range" warning
  printf "%s %6.2f\n",  $solution->get($_)->{structure},
  			$solution->get($_)->{energy}
	if defined  $solution->get($_);
}
$RNA::cut_point = 3;
my $e =  RNA::energy_of_struct("GCGC", "(())");
ok(int($e*100+0.5), 70);

