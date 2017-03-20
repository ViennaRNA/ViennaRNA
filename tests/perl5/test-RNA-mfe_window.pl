#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)

use strict;
use Test::More tests => 21;
use Data::Dumper;
use FileHandle;

use RNA;
use warnings;


my $seq1 = "CGCAGGGAUACCCGCG";
my $s1="CCCCAAAACGGG";
my $s2="CCCGAAAAGGGG";
my $s3="CCCCAAAAGGGG";
my @ali = ($s1,$s2,$s3);
my $fc;
my $mfe;
my @data = ();


sub mfe_window_callback {

  my ($start, $end, $structure, $energy, $data) = @_;
  my %Lfold_hit = ();
  $Lfold_hit{'structure'} = $structure;
  $Lfold_hit{'start'}     = $start;
  $Lfold_hit{'end'}       = $end;
  $Lfold_hit{'energy'}    = $energy;

  push @{$data}, \%Lfold_hit;
}


print "test_mfe_window\n";
$fc= new RNA::fold_compound($seq1, undef, RNA::OPTION_MFE | RNA::OPTION_WINDOW);
$mfe = $fc->mfe_window();
is(sprintf("%6.2f",$mfe), sprintf("%6.2f",-5.60));
undef $fc;

print "test_Lfold_cb\n";
@data = ();
$mfe = RNA::Lfold_cb($seq1, 150, \&mfe_window_callback, \@data);
ok(scalar(@data) == 2);
is($data[0]->{'structure'}, "(((...))).");
is(sprintf("%6.2f", $data[0]->{'energy'}), sprintf("%6.2f", -2.60));
is($data[1]->{'structure'}, "(((.(((...))))))");
is(sprintf("%6.2f", $data[1]->{'energy'}), sprintf("%6.2f", -5.60));


print "test_mfe_window_cb\n";
$fc= new RNA::fold_compound($seq1, undef, RNA::OPTION_MFE | RNA::OPTION_WINDOW);
@data = ();
$mfe = $fc->mfe_window_cb(\&mfe_window_callback, \@data);
ok(scalar(@data) == 2);
is($data[0]->{'structure'}, "(((...))).");
is(sprintf("%6.2f", $data[0]->{'energy'}), sprintf("%6.2f", -2.60));
is($data[1]->{'structure'}, "(((.(((...))))))");
is(sprintf("%6.2f", $data[1]->{'energy'}), sprintf("%6.2f", -5.60));
undef $fc;


print "test_aliLfold_cb\n";
@data = ();
$mfe = RNA::aliLfold_cb(\@ali, 150, \&mfe_window_callback, \@data);
ok(scalar(@data) == 2);
is($data[0]->{'structure'}, "(((.....)))");
is(sprintf("%6.2f", $data[0]->{'energy'}), sprintf("%6.2f", -1.30));
is($data[1]->{'structure'}, "(((......)))");
is(sprintf("%6.2f", $data[1]->{'energy'}), sprintf("%6.2f", -2.70));


print "test_mfe_window_cb (comparative)\n";
$fc= new RNA::fold_compound(\@ali, undef, RNA::OPTION_MFE | RNA::OPTION_WINDOW);
@data = ();
$mfe = $fc->mfe_window_cb(\&mfe_window_callback, \@data);
ok(scalar(@data) == 2);
is($data[0]->{'structure'}, "(((.....)))");
is(sprintf("%6.2f", $data[0]->{'energy'}), sprintf("%6.2f", -1.30));
is($data[1]->{'structure'}, "(((......)))");
is(sprintf("%6.2f", $data[1]->{'energy'}), sprintf("%6.2f", -2.70));
undef $fc;
