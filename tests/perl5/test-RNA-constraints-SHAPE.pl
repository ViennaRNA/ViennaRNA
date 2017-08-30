#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)
use strict;
use Test::More tests => 4;
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
        chomp $line;
        return $line; #return the first line
    }
}

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

undef $fc;
