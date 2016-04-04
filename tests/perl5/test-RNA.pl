#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test;
use Data::Dumper;
use FileHandle;
BEGIN { plan tests => 50; }

use RNA;
use warnings;

######################### End of black magic.

# Insert your test code below (better if it prints "ok 13"
# (correspondingly "not ok 13") depending on the success of chunk 13
# of the test code):

my $seq1  ="CGCAGGGAUACCCGCG";
my $struc1="(((.(((...))))))";
my $seq2  ="GCGCCCAUAGGGACGC";
my $struc2="((((((...))).)))";


# calculate a hamming distance (from util.c)
ok(RNA::hamming($seq1, $seq2), 16);

ok(RNA::bp_distance($struc1, $struc2), 6);

# check a global variable
ok($RNA::temperature, 37);

# fold a sequence

# old obsolete way of calling fold()
my $struct = reverse $seq1;  # wierd way of allocating space
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
ok($costruct, '(((.(((...))))))((((((...))).)))');
$cmfe = RNA::energy_of_struct($seq1 . $seq2, $costruct);
ok(abs($comfe-$cmfe)<1e-5);
my ($x,$ac,$bc,$fcab,$cf) = RNA::co_pf_fold($seq1. $seq2, $struct);
ok(($cf<$comfe)&&($comfe-$cf<1.3));
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
ok((abs($AB-0.0)+abs($AA-0.00578)+abs($BB-0.01100)+abs($A-0.48843)+abs($B-0.47801))<0.0001);
$RNA::cut_point=-1;

# pf_fo ld
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
ok($tree_dist,2);

# check access to a C array
#ok(RNA::ptrvalue($RNA::iindx,3),108);
ok(RNA::intP_getitem($RNA::iindx,3),108);

# memmove does not work in current swig versions
# RNA::memmove($RNA::xsubi, pack('S3', 171,42,93));
# use shortP_setitem instead
RNA::ushortP_setitem($RNA::xsubi, 0, 171);
RNA::ushortP_setitem($RNA::xsubi, 1, 42);
RNA::ushortP_setitem($RNA::xsubi, 2, 93);
ok(RNA::cdata($RNA::xsubi, 6),pack('S3', 171,42,93));

# get a bp prob in two different ways
my $p1 = RNA::get_pr(2,15);
my $ii = RNA::intP_getitem($RNA::iindx, 2);
my $p2 = (RNA::pf_float_precision() != 0) ? RNA::floatP_getitem($RNA::pr, $ii-15) : RNA::doubleP_getitem($RNA::pr, $ii-15);
ok(($p1<0.999) && ($p1>0.99) && (abs($p1-$p2)<1.2e-7));

my $bpf = RNA::Make_bp_profile(length($seq1));
my @bpf = unpack("f*",RNA::cdata($bpf, length($seq1)*4*3));
ok (($bpf[2*3]+$bpf[2*3+1]>.99999)&&$bpf[2*3+2]==0 &&
    ($bpf[2*3+1]>=$p1));

my $pack = RNA::pack_structure($struc1);
ok (RNA::unpack_structure($pack), $struc1);


RNA::parse_structure($struc1);
ok(($RNA::loops==2) && ($RNA::pairs==6)&&($RNA::unpaired==4) &&
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
my $start = RNA::random_string(length $struc1, "GC");
my ($sinv, $cost) = RNA::inverse_fold($start, $struc1);
my ($ss, $en) = RNA::fold($sinv);
ok($ss, $struc1);

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
my $solution = RNA::subopt($seq1, 500);

printf "%d suboptimals\n", scalar(@{$solution});
foreach my $s (@{$solution}){
  printf("%s %6.2f\n", $s->{structure}, $s->{energy});
}

$solution = "";

$RNA::cut_point = 3;
my $e =  RNA::energy_of_struct("GCGC", "(())");
ok(int($e*100+0.5), 70);

my $duplex = RNA::duplexfold($seq1, $seq2);

ok($duplex->{structure}, ".(((.....(((.&)))))).");
undef $duplex;

my @align = ("GCCAUCCGAGGGAAAGGUU", "GAUCGACAGCGUCU-AUCG", "CCGUCUUUAUGAGUCCGGC");
my ($css, $cen) = RNA::alifold(\@align);
ok($css,"(((.(((...)))..))).");
ok(RNA::consens_mis(\@align), "SMBHBHYDRBGDVWmVKBB");
RNA::free_alifold_arrays();

# check the move_set.h functions
$RNA::cut_point=-1;
my $struc1_move = "(..............)";
# move_standar( sequence, start structure, move_type(GRADIENT, FIRST, ADAPTIVE), verbosity, shifts, noLP)
RNA::move_standard($seq1, $struc1_move, 0, 0, 0, 0);
ok($struc1_move, "................");

$struc1_move = "(..............)";
RNA::move_standard($seq1, $struc1_move, 1, 0, 0, 0);
ok($struc1_move, "(((.((....)).)))");

# test simple_xy_coordinates
my $coords = RNA::simple_xy_coordinates($struc1);

foreach my $c (@{$coords}){
  print $c->{X}, ",", $c->{Y}, "\n";
}

#TODO use xsubi random generator, set seed and check here
#my $struc1_move = "(..............)";
#RNA::move_standard($seq1, $struc1_move, 2, 0, 0, 0);
#ok("(((.(((...))))))", $struc1_move);
#print STDERR join(',', unpack('S3', RNA::cdata($RNA::xsubi, 6))), "\n";

#
# Check new scripting language interface
#

# check model details structure
my $md = RNA::md->new(); # default values
ok(int($md->{dangles}), 2);
ok($md->{temperature}, 37.0);

$RNA::dangles     = 0;
$RNA::temperature = 40.1;
$md = RNA::md->new("global"); # global values
ok(int($md->{dangles}), 0);
ok(int($RNA::dangles), 0);
ok($md->{temperature}, 40.1);

# reset globals to default
$RNA::dangles = 2;
$RNA::temperature = 37.0;

# check parameter structures
my $params = RNA::param->new();
ok($params->get_temperature(), 37.0);

$params = RNA::param->new($md);
ok($params->get_temperature(), 40.1);

my $pf_params = RNA::exp_param->new();
ok($pf_params->get_temperature(), 37.0);

$pf_params = RNA::exp_param->new($md);
ok($pf_params->get_temperature(), 40.1);

undef $md;

my $fc = new RNA::fold_compound($seq1);
($ss, $mfe) = $fc->mfe();
printf "%s [ %6.2f ]\n", $ss, $mfe;
ok($ss eq $struc1);

undef $fc;

# test theophylline ligand binding interface
$RNA::noLonelyPairs = 0;
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

##########################################################################################################################
#starting with fold_compound test (Mario Koestl)


my $sq1  =  "CGCAGGGAUACCCGCG";
my $struct1="(((.(((...))))))";
my $sq2  ="  GCGCCCAUAGGGACGC";
my $struct2="((((((...))).)))";
my $sq4  ="  GCGCACAUAGUGACGC";
my $struct4="(..(((...)))...)";
my @alignment = ("CGCAGGGAUACCCGCG","GCGCCCAUAGGGACGC","GCGCACAUAGUGACGC");

####################################################
$fc = new RNA::fold_compound($sq1);
$mfe = sprintf "%6.2f", $fc->eval_structure($struct1);
ok($mfe,sprintf "%6.2f",-5.60 ); #as computed with RNAeval
printf "\n%s [ %6.2f ]\n", $struct1, $mfe;
undef $fc;

###################################################
#		12345678
# seqquence 	CCAAAAGG:
# structure	((....))
   
#$fc = new RNA::fold_compound("CCAAAAGG");
#my @pairtable = (length("((....))"),7,6,0,0,0,0,2,1);   #pairtable[0] = length of structure, 0 if no basepair, 


#$mfe = sprintf "%6.2f", $fc->eval_structure_pt(RNA::vrna_ptable("((....))"));
#ok($mfe,sprintf "%6.2f",0.80*100 ); #as computed with RNAeval, * 100 because eval_structure_pt results 10kal/mol, and not kcal/mol
#printf "\n%s [ %6.2f ]\n", "((....))", $mfe;
#undef $fc; 


##################Test for alignments####################################
#my @align2 = ("GCCAUCCGAGGGAAAGGUU", "GAUCGACAGCGUCU-AUCG", "CCGUCUUUAUGAGUCCGGC");
#my ($css2, $cen2) = RNA::alifold(\@align);

#my $fc2 = new RNA::fold_compound($sq1);
#$fc = new RNA::fold_compound(@align2);

#print "",$fc2->testFunction(\@align2);
#print "",$fc2->testFunction2(\@align2);
#(my $cs,$mfe) = $fc->eval_covar_structure();
#$mfe = sprintf("%6.2f",$mfe);
#ok($mfe,sprintf("%6.2f",-1.93) ); #as computed with RNAAlifold
#printf "Fold_compound Test\ns\n%s [ %6.2f ]\n", $ss, $mfe;
#undef $fc;

######################################################3
#my $sq1  =  "CGCAGGGAUACCCGCG";
#my $struct1="(((.(((...))))))";
$fc = new RNA::fold_compound($sq1);
my $fh= FileHandle->new("outputFile_test.txt","w");
$mfe = sprintf "%6.2f", $fc->eval_structure_verbose($struct1,$fh);
ok($mfe,sprintf "%6.2f",-5.60 ); #as computed with RNAeval
printf "\n%s [ %6.2f ]\n", $struct1, $mfe;
undef $fc;
###################################################
#my $sq1  =  "CGCAGGGAUACCCGCG";
#my $struct1="(((.(((...))))))";
#	      0123456789012345
#  $str_del= "(((.((.....)))))" -> m1 = -6, m2 = -10  , energyStr = -3.5 , energy_move = 2.1
#  $str_ins= "(((.(((...))))))" -> m1 = 6, m2 = 10  , energyStr =   -5.6, energy_move = -2.1

#	     01234567890123456 
my $str=     "(((.(((...))))))";
my $str_del= "(((.((.....)))))";
my $m1=-7;
my $m2 = -11;
$fc = new RNA::fold_compound($sq1);
my $e_move =  $fc->eval_move($str,$m1,$m2); #delete the basePair
ok(sprintf "%6.2f", $e_move, 2.1 ); #as computed with RNAeval
printf "\n%s  with Moveset: (%i,%i) = [ %6.2f ] \n", $str,$m1,$m2, $e_move;

$m1=7;
$m2=11;
$e_move = sprintf "%6.2f", $fc->eval_move($str_del,$m1,$m2); # add the basepair
ok(sprintf "%6.2f", $e_move,-2.1); 
printf "\n%s  with Moveset: (%i,%i) = [ %6.2f ] \n", $str_del,$m1,$m2, $e_move;

undef $fc;
#############################################



####################################




#$md = new RNA::md();
#$md->{noGU} = 1;
#$fc = new RNA::fold_compound("GGGGGGGGGGGGAAAUUUUUUCCCCCC", $md);
#print RNA::vrna_mfe($fc, undef), "\n";
#undef $fc;
#undef $md;
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
