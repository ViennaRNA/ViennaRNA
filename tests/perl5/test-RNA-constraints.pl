#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 14;
use Data::Dumper;
use FileHandle;

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



sub getShapeDataFromFile
{
    my $filepath = shift(@_);
    my @retVec;
    push(@retVec,-999.0);# data list is 1-based, so we add smth. at pos 0
    my $count=1;
    open(my $fh, "<", $filepath) || die "Couldn't open '".$filepath."' for reading because: ".$!;

    foreach my $line (<$fh>)
    {

        my @sp= split(/ /,$line);
        my $pos = $sp[0];
        my $value = $sp[1];

        if($pos == $count)
        {
            push(@retVec,$value +0.00); # make float from string
        }else
        {
            for(my $i=$pos;$i<=$count;$i++)
            {
                push(@retVec,-999.0);
            }
            push(@retVec,$value);
            $count=$pos;
        }
        $count+=1;

    }
    return @retVec;
}

sub getShapeSequenceFromFile
{
    my $f = shift(@_);
    my $retSeq ="";
    open(my $fh, "<", $f) || die "Couldn't open '".$f."' for reading because: ".$!;
    foreach my $line (<$fh>)
    {
        return $line; #return the first line
    }
}
##################################
##test_constraints_add

#seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
#str_con_def=    "(((....)))(((....)))"
#hc.txt=         "P 1 0 2"
#str_con=        "..........(((....)))"

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

#sc.txt = E 3 8 1 -5
my $sc_file = $datadir . "sc.txt";
$fc->sc_init();
$fc->constraints_add($sc_file);
($ss,my $mfeNew) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfeNew);
is(sprintf ("%6.2f",$mfe), sprintf ("%6.2f",$mfeNew+5));
##################################
##test_hc_add_up:
#seq_con  =      "CCCAAAAGGGCCCAAAAGGG"
#str_con_def=    "(((....)))(((....)))"
#str_con=    "..........(((....)))"
print "test_hc_add_up\n";
$fc = new RNA::fold_compound($seq_con);
$fc->hc_add_up(1,RNA::CONSTRAINT_CONTEXT_ALL_LOOPS);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,".((....)).(((....)))");
##################################
##test_hc_add_bp_nonspecific
print "test_hc_add_bp_nonspecific";
#GGGCCCCCCCCCCCCCCCCC
#(((......)))........
$fc= new RNA::fold_compound("GGGCCCCCCCCCCCCCCCCC");
$fc->hc_add_bp_nonspecific(20,-1); # force the last base to pair with some bases upstream
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,"(((..............)))");
##################################
##def test_hc_add_bp
print "test_hc_add_bp";
#$seq_con  =      "CCCAAAAGGGCCCAAAAGGG";
# $str_con_def=    "(((....)))(((....)))";
$fc= new RNA::fold_compound($seq_con);
$fc->hc_add_bp(1,20,RNA::CONSTRAINT_CONTEXT_ENFORCE | RNA::CONSTRAINT_CONTEXT_ALL_LOOPS);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,"(((..............)))");
##################################
##test_hc_add_from_db
print "test_hc_add_from_db";
#seq_con  =      "CCCAAAAGGGCCCAAAAGGG";
#str_con_def=    "(((....)))(((....)))";
#hc.txt=    "xxx.................";
#str_con=    "..........(((....)))";
$fc = new RNA::fold_compound($seq_con);
$fc->hc_add_from_db("xxx.................");
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,$str_con);
##################################
##test_sc_set_up
print "test_sc_set_up";
#        "1234567890
my $seq_sc  =      "CCCAAAAGGG";
$fc = new RNA::fold_compound($seq_sc);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ss,"(((....)))");

$fc->sc_init();

my @m= (0.0,-5.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);  #E 1 0 1 -5 ,  position 1 gets -5 if unpaired , vector starts with 0 and not 1

$fc->sc_set_up(\@m);
($ss,$mfeNew) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfeNew);
is(sprintf ("%6.2f",$mfeNew), sprintf("%6.2f",-5.70));
undef @m;
##################################   !!!geht noch nicht
#test_sc_set_bp
print "test_sc_set_bp";



#my @mm= (
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0],
#[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]);
#
##print Dumper(\@mm);
#$mm[1][9] = -5.0;
#$mm[9][1] = -5.0;
#
#$seq_sc  =      "CCCAAAAGGG";
#$fc = new RNA::fold_compound($seq_sc);
#$fc->sc_set_bp(\@mm);
#($ss,$mfeNew) = $fc->mfe();
#printf("%s [%6.2f] \n",$ss,$mfeNew);
## mfe unconstrained is -2.5
#is(sprintf ("%6.2f",$mfeNew), sprintf("%6.2f",-4.90));


##################################
##test_sc_add_deigan
print "test_sc_add_deigan\n";
my $seq  =  getShapeSequenceFromFile($datadir . "Lysine_riboswitch_T._martima.db");
my @reactivities = getShapeDataFromFile($datadir . "Lysine_riboswitch_T._martima.shape_2rows");

$fc= new RNA::fold_compound($seq);
print "@reactivities \n";

$fc->sc_add_SHAPE_deigan(\@reactivities,1.9,-0.7,RNA::OPTION_DEFAULT);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is(sprintf ("%6.2f",$mfe), sprintf("%6.2f",-121.55));
##################################
##test_sc_add_SHAPE_deigan2
print "test_sc_add_SHAPE_deigan2";
my $seq_ribo  =  getShapeSequenceFromFile($datadir . "TPP_riboswitch_E.coli.db");
my @reactivities_ribo = getShapeDataFromFile($datadir . "TPP_riboswitch_E.coli.shape_2rows");

$fc= new RNA::fold_compound($seq_ribo);
print "@reactivities_ribo \n";

$fc->sc_add_SHAPE_deigan(\@reactivities_ribo,1.9,-0.7,RNA::OPTION_DEFAULT);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is(sprintf ("%6.2f",$mfe), sprintf("%6.2f",-52.61));
##################################
# I just added completely randomly choosen sequences and shape data, this test only checks if the function can be called, not if it results correct results.
##test_sc_add_SHAPE_deigan_ali
print "test_sc_add_SHAPE_deigan_ali";
# shape data from
my $shapeSeq1 = "CCCAAAAGGG";
my $shapeSeq2 = "AAUAAAAAUU";

#my @shapeData1 = (-999.0,0.04,1.12,1.1,-999.0,0.05,0.5,0.3,-999.0,1.4);
#my @shapeData2 = (-999.0,-999.0,-999.0,1.23,1.4,0.05,0.5,0.3,-999.0,1.4);
my @shapeAli = ($shapeSeq1,$shapeSeq2);
$fc= new RNA::fold_compound(\@shapeAli);

my @assoc = (-1,1,2);
my $ret = $fc->sc_add_SHAPE_deigan_ali(\@shapeAli, \@assoc,1.8,-0.6);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($ret,1);
##################################
##test_sc_add_SHAPE_zarringhalam
print("test_sc_add_SHAPE_zarringhalam");
$seq_ribo  =  getShapeSequenceFromFile($datadir . "TPP_riboswitch_E.coli.db");
@reactivities_ribo = getShapeDataFromFile($datadir . "TPP_riboswitch_E.coli.shape_2rows");

$fc= new RNA::fold_compound($seq_ribo);
print "@reactivities_ribo \n";

$fc->sc_add_SHAPE_zarringhalam(\@reactivities_ribo,0.5,0.5,"O");
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is(sprintf ("%6.2f",$mfe), sprintf("%6.2f",-5.28));

##################################

##test_sc_add_hi_motif
print "test_sc_add_hi_motif";
$fc= new RNA::fold_compound("GGUGAUACCAGAUUUCGCGAAAAAUCCCUUGGCAGCACCUCGCACAUCUUGUUGUCUGAUUAUUGAUUUUUCGCGAAACCAUUUGAUCAUAUGACAAGAUUGAG");
#struct =          ".(((((..((((((((((((((((((....(((((((............)))))))........)))))))))))))...)))))))))).............."
#structWithMotif=  "((((...((((((((......)))))...)))...))))...((.((((((((((.((((((.((.((((....)))).))..)))))).)))))))))))).."
my $r=$fc->sc_add_hi_motif("GAUACCAG&CCCUUGGCAGC","(...((((&)...)))...)",-9.22);
($ss,$mfe) = $fc->mfe();
printf("%s [%6.2f] \n",$ss,$mfe);
is($r,1);
##################################
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
