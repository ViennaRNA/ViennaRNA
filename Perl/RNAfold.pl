#!/usr/local/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2000-11-22 17:42:02 ivo>

use RNA;
use Getopt::Long;
use strict;

 Getopt::Long::config("no_ignore_case");

use vars qw/$opt_debug $opt_v $ParamFile $pf $ns_bases/;
my $sfact=1.07;
&usage() unless GetOptions("p|p1" => \$pf,
			   "p0" => sub {$pf=1; $RNA::do_backtrack=0},
			   "C"   => \$RNA::fold_constrained,
			   "T=f" => \$RNA::temperature,
			   "4" => sub {$RNA::tetra_loop = 0},
			   "d|d0" => sub {$RNA::dangles=0},
			   "d2" => sub {$RNA::dangles=2},,
			   "noGU" => \$RNA::noGU,
			   "noCloseGU" => \$RNA::no_closingGU,
			   "noLP" => \$RNA::noLonelyPairs,
			   "e=i" => \$RNA::energy_set,
			   "P=s" => \$ParamFile,
			   "nsp=s" => \$ns_bases,
			   "S=f" => \$sfact);

 RNA::read_parameter_file($ParamFile) if ($ParamFile);

if ($ns_bases) {
   $RNA::nonstandards = "";
   foreach my $p ( split(/,/, $ns_bases) ) {
	if ($p =~ s/^-//) {
	    $RNA::nonstandards .= reverse($p)
	    }
	$RNA::nonstandards .= $p;
    }
    print "$RNA::nonstandards\n";
}

my $istty = (-t STDIN) && (-t STDOUT);
if (($RNA::fold_constrained)&&($istty)) {
   print "Input constraints using the following notation:\n";
   print "| : paired with another base\n";
   print ". : no constraint at all\n";
   print "x : base must not pair\n";
   print "< : base i is paired with a base j<i\n";
   print "> : base i is paired with a base j>i\n";
   print "matching brackets ( ): base i pairs base j\n";
} 
	
if ($istty) {
   print "\nInput string (upper or lower case); @ to quit\n";
   for (1..8) { print "....,....$_";}
   print "\n";
}
my $fname;
while (<>) {	# main loop: continue until end of file
   my ($string, $structure, $cstruc);
   # skip comment lines and get filenames
   if (/^>\s*(\S*)/) {
      $fname = $1;
      next;
   }
   last if (/^@/);
   
   if (/(\S+)/) {
      $string = $1;
   } else {
      next;
   }
   
   $string = uc($string);
   my $length = length($string);
   printf("length = %d\n", $length) if ($istty);
   
   if ($RNA::fold_constrained) {
      $_ = <>;
      $cstruc = $1 if (/(\S+)/);
      die("constraint string has wrong length")
	  if (length($cstruc)!=$length);
      $structure = $cstruc;
   } else {
      $structure = $string; # wierd way to allocate space
   }
   my $min_en = RNA::fold($string, $structure);
   print "$string\n$structure";
   if ($istty) {
      printf("\n minimum free energy = %6.2f kcal/mol\n", $min_en);
   } else {
      printf(" (%6.2f)\n", $min_en);
   }
   my $ffname = ($fname) ? ($fname . '_ss.ps') : 'rna.ps';
   &RNA::PS_rna_plot($string, $structure, $ffname);
   
   if ($pf) {

      # recompute with dangles as in pf_fold()
      $RNA::dangles=2 if ($RNA::dangles);
      $min_en = RNA::energy_of_struct($string, $structure); 
      
      my $kT = ($RNA::temperature+273.15)*1.98717/1000.; # in Kcal 
      $RNA::pf_scale = exp(-($sfact*$min_en)/$kT/$length);
      print STDERR "scaling factor $RNA::pf_scale\n" if ($length>2000);
      
      $structure = $cstruc if ($RNA::fold_constrained);
      my $energy = RNA::pf_fold($string, $structure);
      
      if ($RNA::do_backtrack) {
	 print $structure;
	 printf(" [%6.2f]\n", $energy) if (!$istty);
	 print "\n";
      }
      if (($istty)||(!$RNA::do_backtrack)) {
	 printf(" free energy of ensemble = %6.2f kcal/mol\n", $energy);
	 printf(" frequency of mfe structure in ensemble %g\n",
		exp(($energy-$min_en)/$kT));
      }
	
      if ($RNA::do_backtrack) {
	 $ffname = ($fname)?($fname . "_dp.ps"):"dot.ps";
	 &RNA::PS_dot_plot($string, $ffname);
      }
   }
   undef $fname;
} 

 RNA::free_pf_arrays() if ($pf);
 RNA::free_arrays();

sub usage()
{
   die("usage: " . 
       "RNAfold [-p[0]] [-C] [-T temp] [-4] [-d] [-noGU] [-noCloseGU]\n" .
       "               [-e e_set] [-P paramfile] [-nsp pairs] [-S scale]");
}


# End of file
