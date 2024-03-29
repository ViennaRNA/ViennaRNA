#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 56;
use Data::Dumper;
use FileHandle;

use RNA;
use warnings;
use RNAHelpers qw(:VERSION);

######################### End of black magic.

# Insert your test code below (better if it prints "is 13"
# (correspondingly "not is 13") depending on the success of chunk 13
# of the test code):

my $seq1  ="CGCAGGGAUACCCGCG";
my $struc1="(((.(((...))))))";
my $seq2  ="GCGCCCAUAGGGACGC";
my $struc2="((((((...))).)))";


# Version number
is($RNA::VERSION, $RNAHelpers::VERSION, "Version number (\$RNA::VERSION)");

# calculate a hamming distance (from util.c)
is(RNA::hamming($seq1, $seq2), 16, "Hamming distance (RNA::hamming)");

is(RNA::bp_distance($struc1, $struc2), 6, "Base pair distance (RNA::bp_distance)");

# check a global variable
is($RNA::temperature, 37, "Global default temperature (\$RNA::temperature)");

# fold a sequence

# old obsolete way of calling fold()
my $struct = reverse $seq1;  # wierd way of allocating space
my $mfe=RNA::fold($seq1, $struct);
is($struct, $struc1, "MFE folding - obsolete way (RNA::fold)");

# new better interface
($struct, $mfe) = RNA::fold($seq1);
is($struct , $struc1, "MFE folding - new way (RNA::fold)");
# check energy
is(RNA::energy_of_struct($seq1,$struc1), $mfe, "Compare MFE prediction against energy evaluation (RNA::energy_of_struct)");

# check constrained folding
$RNA::fold_constrained = 1;
my($struct3, $cmfe) = RNA::fold($seq1, '....xx....xx...');
is($struct3, '(((..........)))', "Constrained MFE prediction - structure (RNA::fold)");
is(RNA::energy_of_struct($seq1,$struct3), $cmfe, "Constrained MFE prediction - energy (RNA::energy_of_struct)");
$RNA::fold_constrained = 0;

# test cofold
$RNA::cut_point = length($seq1)+1;
my($costruct, $comfe) = RNA::cofold($seq1 . $seq2);
is($costruct, '(((.(((...))))))((((((...))).)))', "MFE cofold prediction - structure (RNA::cofold)");
$cmfe = RNA::energy_of_struct($seq1 . $seq2, $costruct);
ok(abs($comfe-$cmfe) < 1e-5, "MFE cofold prediction - energy (RNA::energy_of_struct)");
my ($x,$ac,$bc,$fcab,$cf) = RNA::co_pf_fold($seq1. $seq2, $struct);

ok(($cf<$comfe)&&($comfe-$cf<1.3), "Partition function cofold - Ensemble energy (RNA::co_pf_fold)");
#test concentration computation
my $fcaa;
my $fcbb;
my ($usel1,$usel2, $usel3);
($x,$usel1, $usel2, $fcaa, $usel3)= RNA::co_pf_fold($seq1. $seq1, $struct);
$RNA::cut_point = length($seq2)+1;
($x,$usel1, $usel2, $fcbb, $usel3)= RNA::co_pf_fold($seq2. $seq2, $struct);
my ($AB,$AA,$BB,$A,$B)=RNA::get_concentrations($fcab, $fcaa, $fcbb,$ac, $bc, 1e-5, 1e-5);

$AB/=2e-5;
$AA/=2e-5;
$BB/=2e-5;
$A/=2e-5;
$B/=2e-5;

ok((abs($AB-0.0)+abs($AA-0.00578)+abs($BB-0.01100)+abs($A-0.48843)+abs($B-0.47801))<0.0001, "Partition function cofold - concentrations (RNA::get_concentrations)");
$RNA::cut_point=-1;

# pf_fo ld
my ($str,$f) = RNA::pf_fold($seq1, $struct);
ok(($f<$mfe)&&($mfe-$f<0.8), "Partition function - MFE vs. ensemble energy comparison");

# tree distance
my $xstruc = RNA::expand_Full($struc1);
my $T1 = RNA::make_tree($xstruc);
$xstruc = RNA::expand_Full($struc2);
my $T2 = RNA::make_tree($xstruc);
$RNA::edit_backtrack = 1;
my $tree_dist = RNA::tree_edit_distance($T1, $T2);
# print RNA::get_aligned_line(0), RNA::get_aligned_line(1),"\n";
is($tree_dist,2, "Tree edit distance (RNA::tree_edit_distance)");

# check access to a C array
#is(RNA::ptrvalue($RNA::iindx,3),108);
is(RNA::intP_getitem($RNA::iindx,3),108, "Access to global iindx array (\$RNA::iindx)");
# memmove does not work in current swig versions
# RNA::memmove($RNA::xsubi, pack('S3', 171,42,93));
# use shortP_setitem instead
RNA::ushortP_setitem($RNA::xsubi, 0, 171);
RNA::ushortP_setitem($RNA::xsubi, 1, 42);
RNA::ushortP_setitem($RNA::xsubi, 2, 93);
is(RNA::cdata($RNA::xsubi, 6),pack('S3', 171,42,93), "Access to global xsubi array (\$RNA::xsubi)");
#my $foo = pack("S2",1,2);
#print "$foo\n",


# get a bp prob in two different ways
my $p1 = RNA::get_pr(2,15);
my $ii = RNA::intP_getitem($RNA::iindx, 2);
my $p2 = (RNA::pf_float_precision() != 0) ? RNA::floatP_getitem($RNA::pr, $ii-15) : RNA::doubleP_getitem($RNA::pr, $ii-15);
ok(($p1<0.999) && ($p1>0.99) && (abs($p1-$p2)<1.2e-7), "Access to global pr array for base pair probabilities (RNA::get_pr and via RNA::intP_getitem)");


my $bpf = RNA::Make_bp_profile(length($seq1));
my @bpf = unpack("f*",RNA::cdata($bpf, length($seq1)*4*3));
ok (($bpf[2*3]+$bpf[2*3+1]>.99999)&&$bpf[2*3+2]==0 &&
    ($bpf[2*3+1]>=$p1), "Base pair profile (RNA::Make_bp_profile)");

my $pack = RNA::pack_structure($struc1);
is (RNA::unpack_structure($pack), $struc1, "Structure compression (RNA::pack_structure and RNA::unpack_structure)");


RNA::parse_structure($struc1);
is($RNA::loops,2, "Structure parsing - (RNA::parse_structure and \$RNA::loops)");
is($RNA::pairs,6, "Structure parsing - (RNA::parse_structure and \$RNA::pairs)");
is($RNA::unpaired,4, "Structure parsing - (RNA::parse_structure and \$RNA::unpaired)");
is($RNA::loop_degree->[1],2, "Structure parsing - (RNA::parse_structure and \$RNA::loop_degree)");


RNA::PS_rna_plot($seq1, $struc1, "test_ss.ps");
my $anote = "2 15 1 gmark\n" . "3 cmark\n";
RNA::PS_rna_plot_a($seq1, $struc1, "test_ss_a.ps", undef, $anote);
RNA::PS_dot_plot($seq1, "test_dp.ps");
RNA::ssv_rna_plot($seq1, $struct, "test.coord");
# print "$seq1, $struct, $mfe, $f\n";
#print "please check the two postscript files test_ss.ps and test_dp.ps\n";
RNA::write_parameter_file("test.par");

$RNA::symbolset = "GC";
my $start = RNA::random_string(length $struc1, "GC");
my ($sinv, $cost) = RNA::inverse_fold($start, $struc1);
my ($ss, $en) = RNA::fold($sinv);
is($ss, $struc1, "Inverse folding (RNA::inverse)");

RNA::free_pf_arrays();
RNA::free_arrays();
RNA::free_co_arrays();

$RNA::subopt_sorted = 1;
$RNA::noLonelyPairs = 1;
my $solution = RNA::subopt($seq1, undef, 500, undef);

printf "%d suboptimals\n", $solution->size();
for (0..$solution->size()) {
  # the last access should produce a "value out of range" warning
  printf "%s %6.2f\n",  $solution->get($_)->{structure},
  			$solution->get($_)->{energy}
	if defined  $solution->get($_)->{structure};
}

# test native array output of subopt()
$solution = RNA::subopt($seq1, 500);

printf "%d suboptimals\n", scalar(@{$solution});
foreach my $s (@{$solution}){
  printf("%s %6.2f\n", $s->{structure}, $s->{energy});
}

$solution = "";

$RNA::cut_point = 3;
my $e =  RNA::energy_of_struct("GCGC", "(())");
is(int($e*100+0.5), 70, "Energy evaluation - dimer (RNA::energy_of_struct)");

my $duplex = RNA::duplexfold($seq1, $seq2);

is($duplex->{structure}, ".(((.....(((.&)))))).", "Duplex MFE prediction (RNA::duplexfold)");
undef $duplex;

my @align = ("GCCAUCCGAGGGAAAGGUU", "GAUCGACAGCGUCU-AUCG", "CCGUCUUUAUGAGUCCGGC");
my ($css, $cen) = RNA::alifold(\@align);
is($css,"(((.(((...)))..))).", "Comparative MFE prediction - structure (RNA::alifold)");
is(RNA::aln_consensus_mis(\@align), "SMBHBHYDRBGDVWmVKBB", "Comparative MFE prediction - most informative sequence (RNA::aln_consens_mis)");
RNA::free_alifold_arrays();

# check the move_set.h functions
$RNA::cut_point=-1;
my $struc1_move = "(..............)";
# move_standard( sequence, start structure, move_type(GRADIENT, FIRST, ADAPTIVE), verbosity, shifts, noLP)
my ($s,$energy) = RNA::move_standard($seq1, $struc1_move, 0, 1, 0, 0);
print("energy = $energy, s = $s\n");
is($s, "................", "Move set - gradient descent (RNA::move_standard)");

$struc1_move = "(..............)";
($s,$energy) = RNA::move_standard($seq1, $struc1_move, 1, 1, 0, 0);
print("energy = $energy, s = $s\n");
is($s, "(((.((....)).)))", "Move set - First (RNA::move_standard)");

# test simple_xy_coordinates
my $coords = RNA::simple_xy_coordinates($struc1);

foreach my $c (@{$coords}){
  print $c->{X}, ",", $c->{Y}, "\n";
}

#TODO use xsubi random generator, set seed and check here
#my $struc1_move = "(..............)";
#RNA::move_standard($seq1, $struc1_move, 2, 0, 0, 0);
#is("(((.(((...))))))", $struc1_move);
#print STDERR join(',', unpack('S3', RNA::cdata($RNA::xsubi, 6))), "\n";

#
# Check new scripting language interface
#

# check model details structure
my $md = RNA::md->new(); # default values
is(int($md->{dangles}), 2, "Model details - default dangles (RNA::md->new()->{dangles})");
is($md->{temperature}, 37.0, "Model details - default temperature (RNA::md->new()->{temperature})");

$RNA::dangles     = 0;
is(int($RNA::dangles), 0, "Global model settings - set dangles (\$RNA::dangles)");

$RNA::temperature = 40.1;
is($RNA::temperature, 40.1, "Global model settings - set temperature (\$RNA::temperature)");

$md = RNA::md->new(); # global values
is(int($md->{dangles}), 0, "Model details - globally set dangles (\$RNA::dangles and RNA::md->new()->{dangles})");
is($md->{temperature}, 40.1, "Model details - globally set temperature (\$RNA::temperature and RNA::md->new()->{temperature})");

# reset globals to default
$RNA::dangles = 2;
$RNA::temperature = 37.0;

# check parameter structures
my $params = RNA::param->new();
is($params->{temperature}, 37.0, "Energy parameters - default temperature (RNA::param->new()->{temperature})");

$params = RNA::param->new($md);
is($params->{temperature}, 40.1, "Energy parameters - initialized with non-defaut temperature (RNA::param->new(\$md)->{temperature})");

my $pf_params = RNA::exp_param->new();
is($pf_params->{temperature}, 37.0, "Boltzmann factors - default temperature (RNA::exp_param->new()->{temperature})");

$pf_params = RNA::exp_param->new($md);
is($pf_params->{temperature}, 40.1, "Boltzmann factors - initialized with non-defaut temperature (RNA::exp_param->new(\$md)->{temperature})");

undef $md;

##########################################################################################################################
#starting with fold_compound test (Mario Koestl)

$seq1  =        "CGCAGGGAUACCCGCG";
my $struct1=    "(((.(((...))))))";
my $struct11=   "(((.((.....)))))";
$seq2  =        "GCGCCCAUAGGGACGC";
my $struct2=    "((((((...))).)))";
my $seq3  =     "GCGCACAUAGUGACGC";
$struct3=       "(..(((...)))...)";
@align = ($seq1,$seq2,$seq3);

####################################################
##test_create_fold_compound_Single
my $fc = new RNA::fold_compound($seq1);
is($fc->{type}, RNA::FC_TYPE_SINGLE, "Fold compound - default type (new RNA::fold_compound(\$seq)->{type})");
undef $fc;
####################################################
##test_create_fold_compound_Align:
$fc = new RNA::fold_compound(\@align);
is($fc->{type},RNA::FC_TYPE_COMPARATIVE, "Fold compound - alignment type (new RNA::fold_compound(\$alignment)->{type})");
undef $fc;
####################################################
##test_centroid
print "test_centroid\n";
$fc = new RNA::fold_compound(\@align);
$fc->pf();
my ($sc,$dist) = $fc->centroid();
printf "\n%s [ %6.2f ]\n", $sc, $dist;
ok($sc, "Fold compound method - centroid structure (\$fc->centroid())");
ok($dist, "Fold compound method - centroid distance (\$fc->centroid())");
undef $fc;
####################################################
##test_pf:
print "test_pf";
$fc = new RNA::fold_compound($seq1);
($ss, my $gfe) = $fc->pf();
printf "\n%s [ %6.2f ]\n", $ss, $gfe;
ok($ss, "Fold compound method - partition function (\$fc->pf())");
my $bp_dis = $fc->mean_bp_distance();
printf "\n%s [ %6.2f ]\n", $seq1, $bp_dis;
ok($bp_dis, "Fold compound method - mean base pair distance (\$fc->centroid())");

####################################################
## test_pf_dimer:
print "test_pf_dimer\n";
$fc = new RNA::fold_compound($seq1 ."&". $seq2);
($costruct, $comfe) = $fc->mfe_dimer();
is($costruct, "(((.(((...))))))((((((...))).)))");
$cmfe = $fc->eval_structure($costruct);
ok(abs($comfe-$cmfe) < 1e-5);

($x,$ac,$bc,$fcab,$cf) = $fc->pf_dimer();
ok(($cf < $comfe) and ($comfe - $cf < 1.3));


####################################################
## test_filename_sanitize:

print "test_filename_sanitize_simple\n";

my $fn = "bla/bla??_foo\\bar\"r<u>m:ble";
my $fs = RNA::filename_sanitize($fn);
is($fs, "blabla_foobarrumble");

$fs = RNA::filename_sanitize($fn, '-');
is($fs, "bla-bla--_foo-bar-r-u-m-ble");

print "test_filename_sanitize_special_names\n";
$fn = "??";
$fs = RNA::filename_sanitize($fn);
is($fs, "");

$fs = RNA::filename_sanitize($fn, '.');
is($fs, "");

print "test_filename_sanitize_long_names\n";
$fn = ("A" x 120).("B" x 120).("C" x 10)."DEFGHIJ.svg";
$fs = RNA::filename_sanitize($fn);
is($fs, ("A" x 120).("B" x 120).("C" x 10)."D.svg");

$fn = "A.".("A" x 120).("B" x 120).("C" x 10)."DEFGHIsvg";
$fs = RNA::filename_sanitize($fn);
is($fs, "A.".("A" x 120).("B" x 120).("C" x 10)."DEF");

# $md = new RNA::md();
# $md->{noGU} = 1;
# my $fc = new RNA::fold_compound("GGGGGGGGGGGGAAAUUUUUUCCCCCC", $md);
# print "SEPP ",RNA::vrna_mfe($fc, undef), "\n";
# undef $fc;
# undef $md;
#
#$RNA::noGU = 1;
#$md = new RNA::md("global");
#$fc = new RNA::fold_compound("GGGGGGGGGGGGAAAUUUUUUCCCCCC", $md);
#print RNA::vrna_mfe($fc, undef), "\n";
#
#undef $fc;
#undef $md;
#$md = new RNA::md();
#$fc = new RNA::fold_compound("GGGGGGGGGGGGAAAUUUUUUCCCCCC", $md);
#print RNA::vrna_mfe($fc, undef), "\n";
#
#$RNA::noGU = 0;
#undef $fc;
#undef $md;
#$md = new RNA::md();
#$fc = new RNA::fold_compound("GGGGGGGGGGGGAAAUUUUUUCCCCCC", $md);
#($ss, $mfe) = $fc->mfe();
#printf "%s [ %6.2f ]\n", $ss, $mfe;
#undef $fc;
#undef $md;
