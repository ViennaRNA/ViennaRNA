# Before `make install' is performed this script should be runnable with
# `make test'. After `make install' it should work as `perl test.pl'
# Last changed Time-stamp: <2001-07-20 12:57:10 ivo>

######################### We start with some black magic to print on failure.
# Change 1..1 below to 1..last_test_to_print .
# (It may become useful if the test is moved to ./t subdirectory.)

BEGIN { $| = 1; print "1..13\n"; }
END {print "not ok 1\n" unless $loaded;}
use RNA;
$loaded = 1;
print "ok 1\n";

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

$seq1  ="CGCAGGGAUACCCGCG"; 
$struc1="(((.((....)).)))";
$seq2  ="GCGCCCAUAGGGACGC";
$struc2="((((((...))).)))";
# calculate a hamming distance (from util.c)
if (RNA::hamming($seq1, $seq2) == 16) 
    {print "ok 2\n"; } else { print "not ok 2\n"; }
# check a global variable
if ($RNA::temperature == 37)
    {print "ok 3\n"; } else { print "not ok 3\n"; }
# fold a sequence
RNA::initialize_fold(length($seq1));
$struct = $seq1;  # wierd way of allocating space
$mfe=RNA::fold($seq1, $struct);
if ($struct eq $struc1)
    {print "ok 4\n"; } else { print "not ok 4\n"; }
# check energy
if (RNA::energy_of_struct($seq1,$struc1) == $mfe)
    {print "ok 5\n"; } else { print "not ok 5\n"; }
# pf_fold
$f = RNA::pf_fold($seq1, $struct);
if (($f<$mfe)&&($mfe-$f<0.8)) 
    {print "ok 6\n"; } else { print "not ok 6\n"; }
# tree distance
$xstruc = RNA::expand_Full($struc1);
$T1 = RNA::make_tree($xstruc);
$xstruc = RNA::expand_Full($struc2);
$T2 = RNA::make_tree($xstruc);
$tree_dist = RNA::tree_edit_distance($T1, $T2);  
if ($tree_dist==4)
    {print "ok 7\n"; } else { print "not ok 7\n"; }
# check access to a C array
if (RNA::ptrvalue($RNA::iindx,3)==108)
    {print "ok 8\n"; } else { print "not ok 8\n";}
# get a bp prob in two different ways
$p1 = RNA::get_pr(2,15);
$ii = RNA::ptrvalue($RNA::iindx, 2);
$p2 = RNA::ptrvalue($RNA::pr, $ii-15);
if (($p1<0.999) && ($p1>0.99) && (abs($p1-$p2)<1.2e-7))
    {print "ok 9\n"; } else { print "not ok 9 $p1 $p2\n" ;}

$bpf = RNA::Make_bp_profile(length($seq1));
if (1) {
print "omitting test 10; deref_any is currently not functional\n";
} else { 
$bpfi = RNA::ptrcast(RNA::deref_any($bpf, 2), 'float *');
if ((RNA::ptrvalue($bpfi, 0)+RNA::ptrvalue($bpfi,1)>.99999)&&
    (RNA::ptrvalue($bpfi, 1)>=$p1)) {
   print "ok 10\n"; }
else { print "not ok 10 $p1 $p2 ",RNA::ptrvalue($bpfi, 1, "float"),"\n" ;}
}
$pack = RNA::pack_structure($struc1);
if (RNA::unpack_structure($pack) eq $struc1) {
   print "ok 11\n";
} else {
   print "not ok 11\n";
}

RNA::parse_structure($struc1);
if (($RNA::loops==2) && ($RNA::pairs==5)&&($RNA::unpaired==6)&&
    (RNA::ptrvalue($RNA::loop_degree,1)==2)) {
   print "ok 12\n";
} else {
   print "not ok 12\n";
}

RNA::PS_rna_plot($seq1, $struc1, "test_ss.ps");
$anote = "2 15 1 gmark\n" . "3 cmark\n";
RNA::PS_rna_plot_a($seq1, $struc1, "test_ss_a.ps", undef, $anote);
RNA::PS_dot_plot($seq1, "test_dp.ps");
RNA::ssv_rna_plot($seq1, $struct, "test.coord");
# print "$seq1, $struct, $mfe, $f\n";
print "please check the two postscript files test_ss.ps and test_dp.ps\n";
RNA::write_parameter_file("test.par");

$RNA::symbolset = "GC";
RNA::inverse_fold($seq1, $struc1);
$ss = $seq1;
RNA::fold($seq1, $ss);
print 'not ' if ($ss ne $struc1);
print "ok 13\n";
RNA::free_arrays();

$RNA::sorted = 1;
$RNA::noLonelyPairs = 1;
$solution = RNA::subopt($seq1, undef, 500, undef);

printf "%d suboptimals\n", $solution->size();
for (0..$solution->size()-1) {
  printf "%s %6.2f\n",  $solution->get($_)->{structure},
  			$solution->get($_)->{energy};
}

$RNA::cut_point = 3;
$e =  RNA::energy_of_struct("GCGC", "(())");
print  ((int($e*100)==-340)?'':'not ');
print "ok 14\n";
