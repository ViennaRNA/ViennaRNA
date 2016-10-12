#!env perl
package RNA::Params;


use strict;
use Exporter;
use warnings;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw();
%EXPORT_TAGS = ();

#
# below are some global variables
#

my $C_header = <<EOF
/*
    Automatically generated from parameter files shipped with
    the RNAstructure program package.

    See also http://rna.urmc.rochester.edu/NNDB/index.html for
    a detailed description of the parameters we use here.

    (c) Ronny Lorenz 2016
*/

/*
     Current free energy parameters are summarized in:

     D.H. Turner and D.H Mathews
     "NNDB: the nearest neighbor parameter database for predicting stability of
      nucleic acid secondary structure"
     Nucleic Acids Res. 2010 Jan; 38(Database issue): D280-D282.

     D.H. Mathews, M.D. Disney, J.L. Childs, S.J. Schroeder, M. Zuker, D.H. Turner
     "Incorporating chemical modification constraints into a dynamic programming
      algorithm for prediction of RNA secondary structure"
     Proc. Natl Acad. Sci. USA. 2004; 101:7287-7292

     D.H. Mathews, J. Sabina, M. Zuker, D.H. Turner
     "Expanded sequence dependence of thermodynamic parameters improves
     prediction of RNA secondary structure"
     JMB, 288, pp 911-940, 1999

     Enthalpies taken from:

     A. Walter, D Turner, J Kim, M Lyttle, P M"uller, D Mathews, M Zuker
     "Coaxial stacking of helices enhances binding of oligoribonucleotides.."
     PNAS, 91, pp 9218-9222, 1994

     D.H. Turner, N. Sugimoto, and S.M. Freier.
     "RNA Structure Prediction",
     Ann. Rev. Biophys. Biophys. Chem. 17, 167-192, 1988.

     John A.Jaeger, Douglas H.Turner, and Michael Zuker.
     "Improved predictions of secondary structures for RNA",
     PNAS, 86, 7706-7710, October 1989.

     L. He, R. Kierzek, J. SantaLucia, A.E. Walter, D.H. Turner
     "Nearest-Neighbor Parameters for GU Mismatches...."
     Biochemistry 1991, 30 11124-11132

     A.E. Peritz, R. Kierzek, N, Sugimoto, D.H. Turner
     "Thermodynamic Study of Internal Loops in Oligoribonucleotides..."
     Biochemistry 1991, 30, 6428--6435
*/
EOF
;

my $C_defaults = <<EOF
#include "ViennaRNA/energy_const.h"

#define NST 0     /* Energy for nonstandard stacked pairs */
#define DEF -50   /* Default terminal mismatch, used for I */
                  /* and any non_pairing bases */
#define NSM 0     /* terminal mismatch for non standard pairs */

#define PUBLIC

PUBLIC double Tmeasure = 37+K0;  /* temperature of param measurements */

EOF
;

my $C_footer = <<EOF
#include "intl11.h"
#include "intl11dH.h"
#include "intl21.h"
#include "intl21dH.h"
#include "intl22.h"
#include "intl22dH.h"
EOF
;

my $C_footer_D = <<EOF
#include "intl11_D.h"
#include "intl11dH_D.h"
#include "intl21_D.h"
#include "intl21dH_D.h"
#include "intl22_D.h"
#include "intl22dH_D.h"
EOF
;

my $C_footer_RD = <<EOF
#include "intl11_RD.h"
#include "intl11dH_RD.h"
#include "intl21_RD.h"
#include "intl21dH_RD.h"
#include "intl22_RD.h"
#include "intl22dH_RD.h"
EOF
;

my @identifiers = ( 
                    { 'tag' => "stack", 'write' => \&write_stack, 'write_c' => \&write_stack_C, 'cord' => 2, 'hybrid_gen' => \&hybrid_stack},
                    { 'tag' => "mismatch_hairpin", 'write' => \&write_mm, 'write_c' => \&write_mm_C, 'cord' => 7, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "mismatch_interior", 'write' => \&write_mm, 'write_c' => \&write_mm_C, 'cord' => 6, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "mismatch_interior_1n", 'write' => \&write_mm, 'write_c' => \&write_mm_C, 'cord' => 9, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "mismatch_interior_23", 'write' => \&write_mm, 'write_c' => \&write_mm_C, 'cord' => 10, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "mismatch_multi", 'write' => \&write_mm, 'write_c' => \&write_mm_C, 'cord' => 8, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "mismatch_exterior", 'write' => \&write_mm, 'write_c' => \&write_mm_C, 'cord' => 11, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "dangle5", 'write' => \&write_dangle, 'write_c' => \&write_dangle_C, 'cord' => 12, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "dangle3", 'write' => \&write_dangle, 'write_c' => \&write_dangle_C, 'cord' => 13, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "int11", 'write' => \&write_int11, 'write_c' => \&write_int11_C, 'cord' => 17, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "int21", 'write' => \&write_int21, 'write_c' => \&write_int21_C, 'cord' => 18, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "int22", 'write' => \&write_int22, 'write_c' => \&write_int22_C, 'cord' => 19, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "hairpin", 'write' => \&write_loop, 'write_c' => \&write_loop_C, 'cord' => 3, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "bulge", 'write' => \&write_loop, 'write_c' => \&write_loop_C, 'cord' => 4, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "interior", 'write' => \&write_loop, 'write_c' => \&write_loop_C, 'cord' => 5, 'hybrid_gen' => \&average_arrays},
                    { 'tag' => "miscloop", 'write' => \&write_misc, 'write_c' => \&write_misc_C, 'cord' => 1, 'hybrid_gen' => \&average_hashes},
                    { 'tag' => "Hexaloops", 'write' => \&write_tloop, 'write_c' => \&write_tloop_C, 'cord' => 14, 'hybrid_gen' => \&hybrid_tloop},
                    { 'tag' => "Tetraloops", 'write' => \&write_tloop, 'write_c' => \&write_tloop_C, 'cord' => 15, 'hybrid_gen' => \&hybrid_tloop},
                    { 'tag' => "Triloops", 'write' => \&write_tloop, 'write_c' => \&write_tloop_C, 'cord' => 16, 'hybrid_gen' => \&hybrid_tloop}
                  );

my %C_varnames  = (
                    'stack' => "stack",
                    'mismatch_hairpin' => "mismatchH",
                    'mismatch_interior' => "mismatchI",
                    'mismatch_interior_1n' => "mismatch1nI",
                    'mismatch_interior_23' => "mismatch23I",
                    'mismatch_multi' => "mismatchM",
                    'mismatch_exterior' => "mismatchExt",
                    'dangle5' => "dangle5_",
                    'dangle3' => "dangle3_",
                    'int11' => "int11_",
                    'int21' => "int21_",
                    'int22' => "int22_",
                    'hairpin' => "hairpin",
                    'bulge' => "bulge",
                    'interior' => "internal_loop",
                    'Triloops' => "Triloop",
                    'Tetraloops' => "Tetraloop",
                    'Hexaloops' => "Hexaloop"
                  );

# recognizes input file suffixes
my @suffixes = (".dh", ".dat");

my %file2parser = (
                    'tstack' => \&rd_mismatch_rot,
                    'tstackm' => \&rd_mismatch_rot,
                    'tstackh' => \&rd_mismatch,
                    'tstacki1n' => \&rd_mismatch,
                    'tstacki23' => \&rd_mismatch,
                    'tstacki' => \&rd_mismatch,
                    'stack' => \&rd_stack,
                    'int11' => \&rd_int11,
                    'int21' => \&rd_int21,
                    'int22' => \&rd_int22,
                    'dangle' => \&rd_dangle,
                    'miscloop' => \&rd_miscloop,
                    'tloop' => \&rd_tloop,
                    'triloop' => \&rd_tloop,
                    'hexaloop' => \&rd_tloop,
                    'loop' => \&rd_loop
                  );

# this hash stores the different parameter types a file usually contains
# Note: the parameter types must be in the same order as they are returned
# by the corresponding retrieval function, specified above!
my %file2params = (
                    'tstack' => ["mismatch_exterior"],
                    'tstackm' => ["mismatch_multi"],
                    'tstackh' => ["mismatch_hairpin"],
                    'tstacki1n' => ["mismatch_interior_1n"],
                    'tstacki23' => ["mismatch_interior_23"],
                    'tstacki' => ["mismatch_interior"],
                    'stack' => ["stack"],
                    'int11' => ["int11"],
                    'int21' => ["int21"],
                    'int22' => ["int22"],
                    'dangle' => ["dangle5", "dangle3"],
                    'miscloop' => ["miscloop"],
                    'tloop' => ["Tetraloops"],
                    'triloop' => ["Triloops"],
                    'hexaloop' => ["Hexaloops"],
                    'loop' => ["hairpin", "bulge", "interior"]
                  );

__PACKAGE__->main( @ARGV ) unless caller();


sub main {
  my $class = shift;

  my $sourcecode;
  my $parfile;
  my $outfile;
  my @input;
  my $dna;
  my $rna;
  my $hybrid;
  my $verbose;
  my $suffix = "";
  my $help = 0;

  GetOptions("c"          => \$sourcecode,  # generate C source code
             "p"          => \$parfile,     # generate Parameter file
             "o=s"        => \$outfile,     # put output in a file
             "i=s"        => \@input,       # Input files to be processed
             "d|dna"      => \$dna,         # print dna parameters
             "r|rna"      => \$rna,         # print rna parameters
             "h|hybrid"   => \$hybrid,      # print hybrid parameters
             "suffix=s"   => \$suffix,      # suffix for descriptor/variable names
             "v|verbose"  => \$verbose,     # be verbose
             "help|?"     => \$help
  ) or pod2usage(2);

  pod2usage(1) if $help;
  pod2usage(1) if not (defined($parfile) or defined($sourcecode));
  pod2usage(1) if @input < 1;
  pod2usage(1) if not (defined($rna) or defined($dna) or defined($hybrid));


  my $test = $class->new(('verbose' => defined($verbose), 'suffix' => $suffix));

  foreach my $f (@input) {  
    $test->addFiles($f);
  }

  my $fh = \*STDOUT;
  if(defined($outfile)){
    open($fh, ">", $outfile) || die("Can't open output file $outfile!");
  }

  if(defined($parfile)){
    $test->print_parameters("RNA", $fh) if defined($rna);
    $test->print_parameters("DNA", $fh) if defined($dna);
    $test->print_parameters("HYBRID", $fh) if defined($hybrid);
  }

  if(defined($sourcecode)){
    $test->print_parameters_C("RNA") if defined($rna);
    $test->print_parameters_C("DNA") if defined($dna);
    $test->print_parameters_C("HYBRID") if defined($hybrid);
  }

  close($fh) if defined($outfile);

  exit(0);

}



sub new{
  my $class   = shift;
  my %options = @_;

  my $self  = {
    'files'             =>  undef,
    'verbose'           => 0,
    'parameters_rna'    => undef,
    'parameters_dna'    => undef,
    'parameters_hybrid' => undef,
    'suffix'            => "",
  };

  bless $self, $class;

  foreach my $k (keys(%{$self})){
    $self->{$k} = $options{$k} if exists($options{$k});
  }

  # init all remaining attributes
  $self->_init();

  return $self;
}


sub addFiles {
  my $self = shift;

  foreach my $f (@_) {
    my @files;
    my @more_files;
    # split comma separated lists
    @files = split(/,/,$f);

    # check whether any of the input files is actually an entire directory
    # and if found, add all the files within this directory
    foreach my $ff (@files){
      if( -d $ff ){
        print STDERR "adding files in directory $f\n" if $self->{'verbose'} != 0;
        opendir(my $dh, $ff) || die "Can't open directory $ff\n";
        while(readdir $dh){
          next if $_ eq ".";
          next if $_ eq "..";
          next if not -f $ff."/".$_;
          # only add files that match our expected pattern
          push @more_files, "$ff/$_" if $_ =~ /.dat$/;
          push @more_files, "$ff/$_" if $_ =~ /.dh$/;
        }
        closedir $dh;
      } else {
        push @more_files, $ff;
      }
    }

    # process the files we've just added
    foreach my $ff (@more_files){
      $self->process_file($ff);
    }

    push @{$self->{'files'}}, @more_files;
  }

}

sub process_file {
  my $self    = shift;
  my $infile  = shift;

  my($file_id, $file_dir, $file_suffix) = fileparse($infile, @suffixes);

  my $en = ($file_suffix eq ".dh") ? "enthalpies" : "energies";

  # remove possible dna prefix
  my $is_dna = $file_id =~ s/^dna// || 0;

  if(exists($file2parser{$file_id}) and exists($file2params{$file_id})){

    print STDERR "processing \"$infile\"...\n" if $self->{'verbose'} != 0;

    # parse file (note: max. 3 return values here)
    my @params = $file2parser{$file_id}($self, $infile);

    print STDERR "parameters parsed have unexpected count", next if scalar(@params) != scalar(@{$file2params{$file_id}});

    # insert the parameters into our parameter storage containers
    for(my $i=0; $i < @{$file2params{$file_id}}; $i++){
      $self->_insert_parameters($file2params{$file_id}[$i], $en, $params[$i], $is_dna);
    }

  } else {
    print STDERR "ignoring $infile...\n" if $self->{'verbose'};
  }

}

sub print_parameters {
  my $self = shift;
  my $na = shift;
  my $fh = shift;
  my $par;

  $fh = $fh || \*STDOUT;

  # print ViennaRNA parameter file header
  print $fh "## RNAfold parameter file v2.0\n";

  if($na eq "RNA"){
    print $fh "## RNA parameters\n\n";
    $par = \%{$self->{'parameters_rna'}};
    @{$self->{_pnames}} = @{$self->{_pnames_rna}};
    @{$self->{_bnames}} = @{$self->{_bnames_rna}};
    $self->{'suffix'} = "";
  } elsif($na eq "DNA"){
    print $fh "## DNA parameters\n\n";
    $par = \%{$self->{'parameters_dna'}};
    @{$self->{_pnames}} = @{$self->{_pnames_dna}};
    @{$self->{_bnames}} = @{$self->{_bnames_dna}};
    $self->{'suffix'} = "_D";
  } elsif($na eq "HYBRID"){
    print $fh "## RNA/DNA hybrid parameters\n\n";
    $self->_generate_hybrid_parameters();
    $par = \%{$self->{'parameters_hybrid'}};
    @{$self->{_pnames}} = @{$self->{_pnames_hybrid}};
    @{$self->{_bnames}} = @{$self->{_bnames_hybrid}};
    $self->{'suffix'} = "_RD";
  } else {
    print $fh "## Unrecognized parameter set\n\n";
  }

  # merge energies and enthalpies for special loop types
  $self->_merge_e_dh();

  if(defined($par)){
    foreach my $k (@identifiers){
      if(exists($par->{$k->{'tag'}})){
        foreach my $type (sort(keys(%{$par->{$k->{'tag'}}}))){
          my ($id, $data) = ($k->{'tag'}, $par->{$k->{'tag'}}{$type});
          my $string = $k->{'write'}($self, $id, $type, $data);
          print $fh $string, "\n" if defined($string);
        }
      }
    }
  }
  print $fh "# END\n";
}


sub print_parameters_C {
  my $self = shift;
  my $na = shift;

  my $par;
  my $header = "";
  my $footer = "";
  my $epar_fh;
  my $file_suffix = "";

  $header = $C_header."\n";
  $header .= "/* RNAfold parameters v2.x */\n";

  if($na eq "RNA"){
    $par = \%{$self->{'parameters_rna'}};
    $header .= "/* RNA parameters */\n\n";
    $footer = $C_footer."\n";
    @{$self->{_pnames}} = @{$self->{_pnames_rna}};
    @{$self->{_bnames}} = @{$self->{_bnames_rna}};
    $self->{'suffix'} = "";
  } elsif($na eq "DNA"){
    $file_suffix = "_D";
    $par = \%{$self->{'parameters_dna'}};
    $header .= "/* DNA parameters */\n\n";
    $footer = $C_footer_D."\n";
    @{$self->{_pnames}} = @{$self->{_pnames_dna}};
    @{$self->{_bnames}} = @{$self->{_bnames_dna}};
    $self->{'suffix'} = "_D";
  } elsif($na eq "HYBRID"){
    $file_suffix = "_RD";
    $self->_generate_hybrid_parameters();
    $par = \%{$self->{'parameters_hybrid'}};
    $header .= "/* RNA/DNA hybrid parameters */\n\n";
    $footer = $C_footer_RD."\n";
    @{$self->{_pnames}} = @{$self->{_pnames_hybrid}};
    @{$self->{_bnames}} = @{$self->{_bnames_hybrid}};
    $self->{'suffix'} = "_RD";
  } else {
    print STDERR "Unrecognized parameter set in RNA::Params->print_parameters_C()!\n";
  }

  # merge energies and enthalpies for special loop types
  $self->_merge_e_dh();

  if(defined($par)){
    open($epar_fh, ">", "energy_par".$file_suffix.".c") || die "Can't open file energy_par".$file_suffix.".c for writing!";

    # print file header
    print $epar_fh $header, "\n";

    # print some default values
    print $epar_fh $C_defaults, "\n";

    foreach my $k (sort { $a->{'cord'} <=> $b->{'cord'} } @identifiers){

      if(exists($par->{$k->{'tag'}})){
        foreach my $type (sort(keys(%{$par->{$k->{'tag'}}}))){
          my ($id, $data) = ($k->{'tag'}, $par->{$k->{'tag'}}{$type});

          # create separate files for special interior loop parameters
          if(($k->{'tag'} eq "int11") || ($k->{'tag'} eq "int21") || ($k->{'tag'} eq "int22")){
            my $intloop_fh;
            my $fname = $k->{'tag'};
            $fname =~ s/int/intl/;
            $fname .= "dH" if $type eq "enthalpies";
            $fname .= $file_suffix;
            $fname .= ".h";
            open($intloop_fh, ">", $fname) || die("Can't open file $fname for writing!");
            print $intloop_fh $k->{'write_c'}($self,$id, $type, $data), "\n";
            close $intloop_fh;
          } else {
            my $string = $k->{'write_c'}($self, $id, $type, $data);
            print $epar_fh $string if defined($string);
          }
        }
      }
    }

    # print file footer
    print $epar_fh $footer, "\n";

    close($epar_fh);
  }
}


sub _init{
  my $self = shift;

  $self->{_rev} = ();
  @{$self->{_pair}} = ( [ 0, 0, 0, 0, 0, 0, 0, 0],
                        [ 0, 0, 0, 0, 5, 0, 0, 5],
                        [ 0, 0, 0, 1, 0, 0, 0, 0],
                        [ 0, 0, 2, 0, 3, 0, 0, 0],
                        [ 0, 6, 0, 4, 0, 0, 0, 6],
                        [ 0, 0, 0, 0, 0, 0, 2, 0],
                        [ 0, 0, 0, 0, 0, 1, 0, 0],
                        [ 0, 6, 0, 0, 5, 0, 0, 0]);
  @{$self->{_pnames_rna}} = ('NP', 'CG', 'GC', 'GU', 'UG', 'AU', 'UA', 'NN');
  @{$self->{_bnames_rna}} = ('E','A','C','G','U');
  @{$self->{_pnames_dna}} = ('NP', 'CG', 'GC', 'GT', 'TG', 'AT', 'TA', 'NN');
  @{$self->{_bnames_dna}} = ('E','A','C','G','T');
  @{$self->{_pnames_hybrid}} = ('NP', 'CG', 'GC', 'GT', 'UG', 'AT', 'UA', 'NN');
  @{$self->{_bnames_hybrid}} = ('E','A','C','G','U/T');
  $self->{_pnames} = $self->{_pnames_rna};
  $self->{_bnames} = $self->{_bnames_rna};
  %{$self->{_num}} = ('A' => 1, 'C' => 2, 'G' => 3, 'U' => 4, 'T' => 4);

  # In mfold order pairs are AU CG GC UA GU UG
  # Vtype translate pair numbers from mfold to Vienna
  @{$self->{_Vtype}} = (undef,5,1,2,6,3,4);

  # prepare @rev array
  for my $a (1..6) {
    for my $b (1..6) {
      $self->{_rev}[$self->{_pair}[$a][$b]] = $self->{_pair}[$b][$a];
    }
  }

  %{$self->{'parameters_rna'}} = ();
  %{$self->{'parameters_dna'}} = ();
  %{$self->{'parameters_hybrid'}} = ();

}

sub max {
  if (!defined($_[0])) {
    warn 'undefined values in max';
    return 666;
  }
  my $max = shift;
  foreach (@_) {
    $max = $_ if $max < $_;
  }
  return $max;
}

# add entries for unknown base and nonstandard pairs
# i.e. entries for pairtype 7 and base type 0
# add_missing($aref, dim1, dim2, ...)
sub add_missing {
  my $aref = shift;
  return unless ref($aref);
  my $d = shift;
  my @dims = @_;

  for (1..$d) {
    add_missing($$aref[$_], @dims);
  }

  if (@dims == 0) {
    $$aref[7] = max(@{$aref}[1..$d]) if $d==6;
    $$aref[0] = max(@{$aref}[1..$d]) if $d==4;
    return;
  }
  foreach (@dims) {$_=7 if $_==6};
  my @i = map 0, 1..scalar(@dims);
  $i[0]=1 if $dims[0]>4;

  while ($i[-1]<=$dims[-1]) {
    my $indices = '[' . join('][', @i) . ']';
    my $code;
    $code = '$$aref[7]'. $indices . '= max(map($$aref[$_]'.$indices.', 1..6))'
      if $d>4;
    $code = '$$aref[0]'. $indices . '= max(map($$aref[$_]'.$indices.', 1..4))'
      if $d==4;
#    print "$code\n";
    eval $code;

    for my $d (0..$#i) {
      $i[$d]++;
      last if ($i[$d]<=$dims[$d]);
      $i[$d]= $dims[$d]==4 ? 0 : 1 if $d<$#i;
    }
  }
}

sub _insert_parameters {
  my $self    = shift;
  my $id      = shift;
  my $type    = shift;
  my $par     = shift;
  my $is_dna  = shift;

  if($is_dna){
    print STDERR "found DNA \"", $id, " ", $type, "\"\n" if $self->{'verbose'} != 0;
    %{$self->{'parameters_dna'}{$id}} = () if not exists($self->{'parameters_dna'}{$id});
    $self->{'parameters_dna'}{$id}{$type} = $par;
  } else {
    print STDERR "found RNA \"", $id, " ", $type, "\"\n" if $self->{'verbose'} != 0;
    %{$self->{'parameters_rna'}{$id}} = () if not exists($self->{'parameters_rna'}{$id});
    $self->{'parameters_rna'}{$id}{$type} = $par;
  }
}

# merge energies and enthalpies for some loop types
sub _merge_e_dh {
  my $self = shift;

  my @special_hp = ('miscloop', 'Hexaloops', 'Tetraloops', 'Triloops');
  my @parameter_sets = ('parameters_hybrid', 'parameters_dna', 'parameters_rna');

  for my $t (@special_hp){
    for my $set (@parameter_sets){
      if(exists($self->{$set}{$t})){
        if(exists($self->{$set}{$t}{'energies'}) and exists($self->{$set}{$t}{'enthalpies'})){
          my $energies = $self->{$set}{$t}{'energies'};
          my $enthalpies = $self->{$set}{$t}{'enthalpies'};
          #delete $self->{$set}{$t}{'energies'};
          #delete $self->{$set}{$t}{'enthalpies'};
          delete $self->{$set}{$t}{'both'} if exists($self->{$set}{$t}{'both'});
          $self->{$set}{$t}{'both'} = {'energies' => $energies, 'enthalpies' => $enthalpies};
        }
      }
    }
  }
}



sub _generate_hybrid_parameters {
  my $self = shift;
  print STDERR "generating hybrid parameters\n" if $self->{'verbose'};

  my $par_rna = \%{$self->{'parameters_rna'}};
  my $par_dna = \%{$self->{'parameters_dna'}};
  my $par_hybrid = \%{$self->{'parameters_hybrid'}};
  
  foreach my $k (@identifiers){
    if(exists($par_rna->{$k->{'tag'}}) and exists($par_dna->{$k->{'tag'}})){
      foreach my $type (sort(keys(%{$par_rna->{$k->{'tag'}}}))){
        if(exists($par_rna->{$k->{'tag'}}{$type}) and exists($par_dna->{$k->{'tag'}}{$type})){  # make sure that data is available for both, RNA and DNA
          my ($id, $data_rna, $data_dna) = ($k->{'tag'}, $par_rna->{$k->{'tag'}}{$type}, $par_dna->{$k->{'tag'}}{$type});
          my $bla = $k->{'hybrid_gen'}( $data_rna, $data_dna, $type );
          if(not exists($par_hybrid->{$k->{'tag'}})){
            my %h;
            $par_hybrid->{$k->{'tag'}} = \%h;
          }
          $par_hybrid->{$k->{'tag'}}{$type} = $bla;
        } else {
          die "$type not available in both parameter sets\n";
        }
      }
    }
  }
}

sub hybrid_stack {

  my $p_rna = shift;
  my $p_dna = shift;
  my $type  = shift;

  my @stacking_energies = (
    #  NP   CG    GC    GU    UG    AU    UA    NS
    undef,
    [  undef, -1.7, -2.1,    0,    0, -0.9, -0.9,    0],
    [  undef, -2.9, -2.7,    0,    0, -1.1, -1.3,    0],
    [  undef,    0,    0,    0,    0,    0,    0,    0],
    [  undef,    0,    0,    0,    0,    0,    0,    0],
    [  undef, -1.8, -2.1,    0,    0, -0.9, -1.0,    0],
    [  undef, -1.6, -1.5,    0,    0, -0.2, -0.6,    0],
    [  undef,    0,    0,    0,    0,    0,    0,    0]
  );

  my @stacking_enthalpies = (
    #  NP    CG     GC     GU     UG     AU     UA     NS
    undef,
    [  undef, -16.3,  -9.3,     0,     0,  -7.0,  -9.0,     0],
    [  undef, -12.8,  -8.0,     0,     0,  -7.8,  -5.5,     0],
    [  undef,     0,     0,     0,     0,     0,     0,     0],
    [  undef,     0,     0,     0,     0,     0,     0,     0],
    [  undef,  -9.1,  -5.9,     0,     0,  -8.3,  -7.8,     0],
    [  undef, -10.4,  -8.6,     0,     0, -11.5,  -7.8,     0],
    [  undef,     0,     0,     0,     0,     0,     0,     0]
  );

  return \@stacking_energies   if $type eq "energies";
  return \@stacking_enthalpies if $type eq "enthalpies";
}

sub hybrid_tloop {
  my %tloop = ();
  return \%tloop;
}

sub average_arrays {
  my $a = shift;
  my $b = shift;

  my @foobar = ();

  # check whether both $a, and $b are array references
  if((ref $a eq 'ARRAY') and (ref $b eq 'ARRAY')){
    # loop over both arrays (take bigger array as loopcount
    my $stop = (@{$a} > @{$b}) ? $#$a : $#$b;
    for(my $i = 0; $i <= $stop; $i++){
      # skip entry if one of the array values is undef
      if(defined($a->[$i]) and defined($b->[$i])){
        # iterate further if array values are array references again
        if(ref $a->[$i] eq 'ARRAY'){
          if(ref $b->[$i] eq 'ARRAY'){
            $foobar[$i] = average_arrays($a->[$i], $b->[$i]);
          } else {
            die "arrays did not match!\n";
          }
        } else { # actually average values
          if(($a->[$i] eq ".") or ($b->[$i] eq ".")){
            $foobar[$i] = ".";
          } else {
            $foobar[$i] = ($a->[$i] + $b->[$i]) / 2.;
          }
        }
      }
    }
  } elsif((ref $a eq 'HASH') and (ref $b eq 'HASH')) {
    print "passed variables are hash references\n"; 
  } else {
    die "passed variables are no array references";
  }

  return \@foobar;
}

sub average_hashes {
  my $a = shift;
  my $b = shift;

  my %foobar = ();

  # check whether both $a, and $b are array references
  if((ref $a eq 'HASH') and (ref $b eq 'HASH')) {
    # loop over all hash keys
    my @keys = keys(%{$a});
    foreach my $k (@keys){

      die "key $k not in both hashes\n" if not exists $b->{$k};
      my $val = undef;
      # loop over arrays if values are array reference
      if((ref $a->{$k} eq 'ARRAY') and (ref $b->{$k} eq 'ARRAY')){
        my @dat = ();
        my $stop = (@{$a->{$k}} > @{$b->{$k}}) ? $#{$a->{$k}} : $#{$b->{$k}};
        for(my $i = 0; $i <= $stop; $i++){
          if(defined($a->{$k}->[$i]) and defined($b->{$k}->[$i])){
            $dat[$i] = ($a->{$k}->[$i] + $b->{$k}->[$i]) / 2.;
          }
        }
        $val = \@dat;
      } elsif((ref $a->{$k} eq '') and (ref $b->{$k} eq '')){
        $val = ($a->{$k} + $b->{$k}) / 2.;
      } else {
        print $a->{$k}, " = ", (ref $a->{$k}), " ", $b->{$k}, " = ", (ref $b->{$k}), "\n";
        print STDERR "do not know what to do with input\n";
      }
      $foobar{$k} = $val;
    }
  } else {
    die "passed variables are no array references";
  }

  return \%foobar;
}



#############################
# Parsing subroutines below #
#############################

sub skip_comments {
  my $fh = shift;
  while (<$fh>) { last if /^\s+3. \<---* 5./; }
}

#############################
# rd_stack($filename)
#
# Read stacking energies/enthalpies from $filename
#
#############################
sub rd_stack {
  my $self     = shift;
  my $filename = shift;
  my ($a,$c) = (1,1);
  my @stack;
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  skip_comments($fh);
  while (<$fh>) {
    chomp;
    my @F = split;
    if (@F != 16) { # must have 16 elements on line
      skip_comments($fh); # skip to next block of data
      $a++ if @stack; $c=1;
      next;
    }

    my ($b, $d) = (1,1);
    foreach my $e (@F) {
      warn "unexpected token $e in line $. file $ARGV" unless $e =~ /\d*\.\d*/;
      my $t1 = $self->{_pair}[$a][$b];
      my $t2 = $self->{_pair}[$d][$c];
      $stack[$t1][$t2] = $e if $t1 && $t2;
      $d++; if ($d>4) {$b++; $d=1;}
    }
    $c++;
  }

  close($fh);
  warn "incorrect number of entries $#stack $a $c" unless ($#stack == 6);

  add_missing(\@stack, 6,6);
  return \@stack;
}

#############################
# rd_mismatch($filename, [$reverse])
#
# Read terminal mismatch energies/enthalpies from $filename
# if $reverse is passed as well, orientation of read parameters
# will be reversed
#
#############################
sub rd_mismatch_rot {
  return rd_mismatch(@_, 1);
}

sub rd_mismatch {
  my $self      = shift;
  my $filename  = shift;
  my $reverse   = shift || 0;
  my ($a,$c) = (1,1);
  my @mm = ();
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  skip_comments($fh);

  $reverse = 0 if not defined($reverse);
  while (<$fh>) {
    chomp;
    my @F = split;
    if (@F != 16) { # must have 16 elements on line
      skip_comments($fh); # skip to next block of data
      $a++ if @mm; $c=1;
      next;
    }

    my ($b, $d) = (1,1);
    foreach my $e (@F) {
      warn "unexpected token $e in line $. file $filename\n" unless $e =~ /\d*\.\d*/;
      if($reverse){
        my $type = $self->{_pair}[$b][$a];
        $mm[$type][$d][$c] = $e if $type;
      } else {
        my $type = $self->{_pair}[$a][$b];
        $mm[$type][$c][$d] = $e if $type;
      }
      # print STDERR "$e $type $c $d\n";
      $d++; if ($d>4) {$b++; $d=1;}
    }
    $c++;
  }

  close($fh);
  warn "incorrect number of entries $#mm" unless ($#mm == 6);

  add_missing(\@mm, 6,4,4);
  return \@mm;
}

#############################
# rd_dangle($filename)
#
# Read dangling end energies/enthalpies from $filename
#
#############################
sub rd_dangle {
  my $self      = shift;
  my $filename  = shift;
  my $a = 1;
  my @d3 = ();
  my @d5 = ();
  my @dang = ();
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  skip_comments($fh);
  while (<$fh>) {
    chomp;
    my @F = split;
    if (@F != 16) { # must have 16 elements on line
      skip_comments($fh); # skip to next block of data
      $a++ if @dang;
      if ($a== 5) {
        $a=1;
        if (@d3 ==0) {@d3 = @dang; @dang = ();}
      }
      next;
    }

    my ($b, $d) = (1,1);
    foreach my $e (@F) {
      warn "unexpected token $e" unless $e =~ /\d*\.\d*/;
      my $type = $self->{_rev}[$self->{_pair}[$a][$b]];
      $dang[$type][$d] = $e if $type;
      #      print STDERR "$e $type $c $d\n";
      $d++; if ($d>4) {$b++; $d=1;}
    }
  }

  @d5 = @dang;
  close($fh);

  warn "incorrect number of entries $#d3 $#d5" unless $#d3 ==6 && $#d5==6;
  add_missing(\@d3, 6,4);
  add_missing(\@d5, 6,4);
  return \@d5, \@d3;
}

#############################
# rd_int11($filename)
#
# Read 1-1 interior loop energies/enthalpies from $filename
#
#############################
sub rd_int11 {
  my $self = shift;
  my $filename  = shift;
  my ($p1, $x) = (1,1);
  my @int11 = ();
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  skip_comments($fh);
  while (<$fh>) {
    chomp;
    my @F = split;
    if (@F != 24) { # must have 16 elements on line
      skip_comments($fh); # skip to next block of data
      $p1++ if @int11; $x=1;
      next;
    }

    my $typ1 = $self->{_Vtype}[$p1];
    my ($p2, $y) = (1,1);
    foreach my $e (@F) {
      warn "unexpected token $e file $ARGV\n" unless $e =~ /\d*\.\d*/;
      my $typ2 = $self->{_rev}[$self->{_Vtype}[$p2]];
      $int11[$typ1][$typ2][$x][$y] = $e;
      $y++; if ($y>4) {$p2++; $y=1;}
    }
    $x++;
  }

  close($fh);

  warn "incorrect number of entries $#int11" unless ($#int11 == 6);
  add_missing(\@int11,6,6,4,4);
  return \@int11;
}

#############################
# rd_int21($filename)
#
# Read 2-1 interior loop energies/enthalpies from $filename
#
#############################
sub rd_int21 {
  my $self = shift;
  my $filename  = shift;
  my ($p1, $x, $y) = (1,1,1);
  my @int21 = ();
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  skip_comments($fh);
  while (<$fh>) {
    chomp;
    my @F = split;
    if (@F != 24) { # must have 16 elements on line
      skip_comments($fh); # skip to next block of data
      $y++ if @int21;
      if ($y>4) {$y=1; $p1++}
      $x=1;
      next;
    }

    my $typ1 = $self->{_Vtype}[$p1];
    my ($p2, $z) = (1,1);
    foreach my $e (@F) {
      warn "unexpected token $e in line $. file $ARGV\n" unless $e =~ /\d*\.\d*/;
      my $typ2 = $self->{_rev}[$self->{_Vtype}[$p2]];
      $int21[$typ1][$typ2][$x][$y][$z] = $e;
      $z++; if ($z>4) {$p2++; $z=1;}
    }
    $x++;
  }

  close($fh);

  warn "incorrect number of entries $#int21" unless ($#int21 == 6);
  add_missing(\@int21,6,6,4,4,4);
  return \@int21;
}

#############################
# rd_int22($filename)
#
# Read 2-2 interior loop energies/enthalpies from $filename
#
#############################
sub rd_int22 {
  my $self = shift;
  my $filename  = shift;
  my ($p1, $p2, $x1, $x2) = (1,1,1,1);
  my @int22 = ();
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  skip_comments($fh);
  my ($typ1, $typ2);
  while (<$fh>) {
    chomp;
    my @F = split;
    if (/^\s+5. ---*\> 3./) {
      $_=<$fh>;
      die "$.: $_ can't recognize pair" unless /([AUGCT]) .* ([AUGCT])/;
      my ($si, $sp) = ($1,$2);
      $_=<$fh>;
      die $_, "can't recognize pair" unless /([AUGCT]) .* ([AUGCT])/;
      my ($sj, $sq) = ($1,$2);
      $typ1 = $self->{_pair}[$self->{_num}{$si}][$self->{_num}{$sj}];
      $typ2 = $self->{_pair}[$self->{_num}{$sq}][$self->{_num}{$sp}];;
      $x1 = $x2 = 1;
    }
    next unless (@F==16)&&/\d*\.\d*/;
    my ($y1, $y2) = (1,1);
    foreach my $e (@F) {
      warn "unexpected token $e in line $. file $ARGV\n" unless $e =~ /\d*\.\d*/;
      $int22[$typ1][$typ2][$x1][$y1][$y2][$x2] = $e;
      $y2++; if ($y2>4) {$y1++; $y2=1;}
    }
    $x2++; if ($x2>4) {$x1++; $x2=1;}
  }

  close($fh);

  warn "incorrect number of entries $#int22\n" unless ($#int22 == 6);
  add_missing(\@int22,6,6,4,4,4,4);
  return \@int22;
}

#############################
# rd_miscloop($filename)
#
# Read Misc loop energies/enthalpies from $filename
#
#############################
sub rd_miscloop {
  my $self = shift;
  my $filename  = shift;
  my $rflag=0;
  my %miscloop = ();
  my @keys = ('extrapolation', 'max_ninio', 'ninio', 'multi', 'multi2',
              'termAU', 'GGG', 'cslope', 'cinter', 'ch3', 'init', 'GAIL');
  my @regs = ('hairpin loops', 'maximum', 'array', 'offset.*penalty', 'offset.*penalty',
              'terminal', 'GGG', 'slope', 'intercept', 'hairpin',
              'intermolec', 'GAIL', '^unkown$');
  my ($prev, $pprev) = ('', '', 0);
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  while (<$fh>) {
    if ($rflag) {
      my $r = shift @regs;
      unless ($pprev =~ /$r/i) {
        warn "didn't find expected comment, skipping: ", $pprev;
        $rflag=0; unshift @regs, $r;
        next;
      }
      my $k = shift @keys;
      my @v = split;
      $miscloop{$k} = (@v==1) ? $v[0] : [@v];
      $rflag=0;
    }
    $rflag = 1 if /^-->/;
    $pprev = $prev;
    $prev = $_;
  }

  close($fh);

  return \%miscloop;
}

#############################
# rd_loop($filename)
#
# Read other loop energies/enthalpies from $filename
#
#############################
sub rd_loop {
  my $self = shift;
  my $filename  = shift;
  my @hl = ('.');
  my @bl = ('.');
  my @il = ('.');
  my $s=1;
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  while (<$fh>) {
    my @F = split;
    next unless @F==4 && $F[0] =~ /^\d+$/;
    warn "error parsing loop data" if $F[0] != $s; $s++;
    push  @il, $F[1];
    push  @bl, $F[2];
    push  @hl, $F[3];
  }

  close($fh);

  return \@hl, \@bl, \@il;
}

#############################
# rd_tloop($filename)
#
# Read tloop energies/enthalpies from $filename
#
#############################
sub rd_tloop {
  my $self = shift;
  my $filename  = shift;
  #my @tl = ();
  my %tl = ();
  my $fh;
  open($fh, "<".$filename) or die("Could not open file $filename\n");
  while (<$fh>) {
    my @F = split;
    next unless @F==2 && $F[0] =~ /^[AUGCT]+$/ && $F[1] =~ /^-?\d+\.?\d*$/;
    $tl{$F[0]} = $F[1];
    #push  @tl, \@F;
  }

  close($fh);
  #return \%tl;
  return \%tl;
}

#############################
# Output subroutines below  #
#############################


sub write_stack {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @stack = @{$_[0]};

  my $string = "";
  $string .= "# ".$id."_".$type if $type eq "enthalpies";
  $string .= "# ".$id if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "\n";
  {
    local $" = '     ';
    $string .= "/*   @{$self->{_pnames}}[1..7]          */\n";
  }
  for my $type1 (1..7) {
    for my $type2 (1..7) {
      if (defined $stack[$type1][$type2]) {
        $string .= sprintf(" %6.0f", $stack[$type1][$type2]*100);
        if ($stack[$type1][$type2] != $stack[$type2][$type1]) {
          warn "entry $type1 $type2 is unsymetric!";
        }
      }
    }
    $string .= "    /* $self->{_pnames}[$type1] */\n";
  }
  return $string;
}

sub write_stack_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @stack = @{$_[0]};

  my $string = "PUBLIC int ".$C_varnames{$id};
  $string .= "dH" if $type eq "enthalpies";
  $string .= "37" if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "[NBPAIRS+1][NBPAIRS+1] = \n";

  {
    local $" = '      ';
    $string .= "/*     @{$self->{_pnames}}        */\n";
  }
  my @rows;
  for my $type1 (0..7) {
    my @cols;
    for my $type2 (0..7) {
      if (defined $stack[$type1][$type2]) {
        push @cols, sprintf(" %6.0f", $stack[$type1][$type2]*100);
        if ($stack[$type1][$type2] != $stack[$type2][$type1]) {
          warn "entry $type1 $type2 is unsymetric!";
        }
      } else { # missing data
        push @cols, "    INF";
      }
    }
    push @rows, "{" . join(",", @cols) ."} /* $self->{_pnames}[$type1] */\n";
  }

  $string .= "{" . join(",", @rows) . "};\n";
  return $string;
}

sub write_mm {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @mm = @{$_[0]};

  my $string = "";
  $string .= "# ".$id."_".$type if $type eq "enthalpies";
  $string .= "# ".$id if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "\n";
  for my $p (1 .. 7) {
    for my $i (0 .. 4) {
      for my $j (0..4) {
        if (defined($mm[$p][$i][$j])) {
          $string .= sprintf(" %6.0f", $mm[$p][$i][$j]*100);
        } else { $string .= "    . ";}
      }
      $string .= "    /* $self->{_pnames}[$p],$self->{_bnames}[$i] */\n";
    }
  }
  return $string;
}

sub write_mm_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @mm = @{$_[0]};

  my $string = "PUBLIC int ".$C_varnames{$id};
  $string .= "dH" if $type eq "enthalpies";
  $string .= "37" if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "[NBPAIRS+1][5][5] = \n";

  my @rows;
  for my $p (0 .. 7) {
    my @level1;
    for my $i (0 .. 4) {
      my @level2;
      for my $j (0 .. 4) {
        if (defined($mm[$p][$i][$j])) {
          push @level2, sprintf(" %6.0f", $mm[$p][$i][$j]*100);
        } else {
          push @level2, "    INF";
        }
      }
      push @level1, "{" . join (",", @level2) ."} /* $self->{_bnames}[$i] */\n";
    }
    push @rows, "/* $self->{_pnames}[$p] */\n {" . join(" ,", @level1) ." }\n";
  }

  $string .= "{" . join(",", @rows) . "};\n";
  return $string;
}

sub write_dangle {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @stack = @{$_[0]};

  my $string = "";
  $string .= "# ".$id."_".$type if $type eq "enthalpies";
  $string .= "# ".$id if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "\n";

  my @d = @{$_[0]};
  {
    local $" = '     ';
    $string .= "/*   @{$self->{_bnames}}          */\n";
  }
  # print "   INF   INF   INF   INF   INF   /* NP */\n";
  for my $p (1 .. 7) {
    for my $i (0 .. 4) {
      if (defined($d[$p][$i])) {
        $string .= sprintf(" %6.0f", $d[$p][$i]*100);
      } else {$string .= "    . ";}
    }
    $string .= "    /* $self->{_pnames}[$p] */\n";
  }

  return $string;
}

sub write_dangle_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @stack = @{$_[0]};

  my $string = "PUBLIC int ".$C_varnames{$id};
  $string .= "dH" if $type eq "enthalpies";
  $string .= "37" if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "[NBPAIRS+1][5] = \n";

  my @d = @{$_[0]};
  {
    local $" = '       ';
    $string .= "/*      @{$self->{_bnames}}        */\n";
  }

  my @rows;
  for my $p (0 .. 7) {
    my @cols;
    for my $i (0 .. 4) {
      if (defined($d[$p][$i])) {
        push @cols, sprintf(" %6.0f", $d[$p][$i]*100);
      } else {
        push @cols, "    INF";
      }
    }
    push @rows, "{" . join(",", @cols) . "} /* $self->{_pnames}[$p] */\n";
  }
  $string .= "{" . join(",", @rows) . "};\n";
  return $string;
}

sub write_int11 {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @int11 = @{$_[0]};

  my $string = "";
  $string .= "# ".$id."_".$type if $type eq "enthalpies";
  $string .= "# ".$id if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "\n";

  for my $p1 (1 .. 7) {
    for my $p2 (1 .. 7) {
      for my $i (0 .. 4) {
        for my $j (0..4) {
          if (defined($int11[$p1][$p2][$i][$j])) {
            $string .= sprintf(" %6.0f", $int11[$p1][$p2][$i][$j]*100);
            warn "Unsymmetric entry [$p1][$p2][$i][$j]"
              if $int11[$p1][$p2][$i][$j] != $int11[$p2][$p1][$j][$i];
          } else {$string .= "    . ";}
        }
        $string .= "    /* $self->{_pnames}[$p1],$self->{_pnames}[$p2],$self->{_bnames}[$i] */\n";
      }
    }
  }
  return $string;
}

sub write_int11_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @int11 = @{$_[0]};

  my $string = "PUBLIC int ".$C_varnames{$id};
  $string .= "dH" if $type eq "enthalpies";
  $string .= "37" if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "[NBPAIRS+1][NBPAIRS+1][5][5] = \n";

  my @rows;
  for my $p1 (0 .. 7) {
    my @level1;
    for my $p2 (0 .. 7) {
      my @level2;
      for my $i (0 .. 4) {
        my @level3;
        for my $j (0..4) {
          if (defined($int11[$p1][$p2][$i][$j])) {
            push @level3, sprintf(" %6.0f", $int11[$p1][$p2][$i][$j] * 100);
            warn "Unsymmetric entry [$p1][$p2][$i][$j]"
              if $int11[$p1][$p2][$i][$j] != $int11[$p2][$p1][$j][$i];
          } else {
            push @level3, "    INF";
          }
        }
        push @level2, "{" . join(",", @level3) . "} /* $self->{_pnames}[$p1],$self->{_pnames}[$p2],$self->{_bnames}[$i] */";
      }
      push @level1, "{" . join("\n  ,", @level2) . "\n  }";
    }
    push @rows, "{" . join("\n ,", @level1) . "\n }";
  }
  $string .= "{" . join("\n,", @rows) . "};\n";
  return $string;
}

sub write_int21 {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @int21 = @{$_[0]};

  my $string = "";
  $string .= "# ".$id."_".$type if $type eq "enthalpies";
  $string .= "# ".$id if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "\n";

  for my $p1 (1 .. 7) {
    for my $p2 (1 .. 7) {
      for my $i (0 .. 4) {
        for my $j (0..4) {
          for my $k (0..4) {
            if (defined($int21[$p1][$p2][$i][$j][$k])) {
              $string .= sprintf(" %6.0f", $int21[$p1][$p2][$i][$j][$k]*100);
            } else {$string .= "    . ";}
          }
          $string .= "    /* $self->{_pnames}[$p1],$self->{_pnames}[$p2],$self->{_bnames}[$i],$self->{_bnames}[$j] */\n";
        }
      }
    }
  }
  return $string;
}

sub write_int21_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @int21 = @{$_[0]};

  my $string = "PUBLIC int ".$C_varnames{$id};
  $string .= "dH" if $type eq "enthalpies";
  $string .= "37" if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "[NBPAIRS+1][NBPAIRS+1][5][5][5] = \n";

  my @rows;
  for my $p1 (0 .. 7) {
    my @level1;
    for my $p2 (0 .. 7) {
      my @level2;
      for my $i (0 .. 4) {
        my @level3;
        for my $j (0..4) {
          my @level4;
          for my $k (0..4) {
            if (defined($int21[$p1][$p2][$i][$j][$k])) {
              push @level4, sprintf(" %6.0f", $int21[$p1][$p2][$i][$j][$k]*100);
            } else {
              push @level4, "    INF";
            }
          }
          push @level3, "{" . join(",", @level4) . "} /* $self->{_pnames}[$p1],$self->{_pnames}[$p2],$self->{_bnames}[$i],$self->{_bnames}[$j] */";
        }
        push @level2, "{" . join("\n   ,", @level3) . "\n   }";
      }
      push @level1, "{" . join("\n  ,", @level2) . "\n  }";
    }
    push @rows, "{" . join("\n ,", @level1) . "\n }";
  }
  $string .= "{" . join("\n,", @rows) . "};\n";
  return $string;
}

sub write_int22 {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @int22 = @{$_[0]};

  my $string = "";
  $string .= "# ".$id."_".$type if $type eq "enthalpies";
  $string .= "# ".$id if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "\n";

  for my $p1 (1 .. 6) {
    for my $p2 (1 .. 6) {
      for my $i (1 .. 4) {
        for my $j (1 .. 4) {
          for my $k (1..4) {
            for my $l (1..4) {
              if (defined($int22[$p1][$p2][$i][$j][$k][$l])) {
                $string .= sprintf(" %6.0f", $int22[$p1][$p2][$i][$j][$k][$l]*100);
                warn "Unsymmetric entry [$p1][$p2][$i][$j][$k][$l]"
                  if $int22[$p1][$p2][$i][$j][$k][$l] !=
                     $int22[$p2][$p1][$k][$l][$i][$j];
              } else {$string .= "    . ";}
            }
            $string .= "    /* $self->{_pnames}[$p1],$self->{_pnames}[$p2],$self->{_bnames}[$i],$self->{_bnames}[$j],$self->{_bnames}[$k] */\n";
          }
        }
      }
    }
  }
  return $string;
}

sub write_int22_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @int22 = @{$_[0]};

  my $string = "PUBLIC int ".$C_varnames{$id};
  $string .= "dH" if $type eq "enthalpies";
  $string .= "37" if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "[NBPAIRS+1][NBPAIRS+1][5][5][5][5] = \n";

  my @rows;
  for my $p1 (0 .. 7) {
    my @level1;
    for my $p2 (0 .. 7) {
      my @level2;
      for my $i (0 .. 4) {
        my @level3;
        for my $j (0 .. 4) {
          my @level4;
          for my $k (0 .. 4) {
            my @level5;
            for my $l (0 .. 4) {
              if (defined($int22[$p1][$p2][$i][$j][$k][$l])) {
                push @level5, sprintf(" %6.0f", $int22[$p1][$p2][$i][$j][$k][$l]*100);
                warn "Unsymmetric entry [$p1][$p2][$i][$j][$k][$l]"
                  if $int22[$p1][$p2][$i][$j][$k][$l] !=
                     $int22[$p2][$p1][$k][$l][$i][$j];
              } else {
                push @level5, "    INF";
              }
            }
            push @level4, "{" . join(",", @level5) . "} /* $self->{_pnames}[$p1],$self->{_pnames}[$p2],$self->{_bnames}[$i],$self->{_bnames}[$j],$self->{_bnames}[$k] */";
          }
          push @level3, "{" . join("\n    ,", @level4) . "\n    }";
        }
        push @level2, "{" . join("\n   ,", @level3) . "\n   }";
      }
      push @level1, "{" . join("\n  ,", @level2) . "\n  }";
    }
    push @rows, "{" . join("\n ,", @level1) . "\n }";
  }
  $string .= "{" . join("\n,", @rows) . "\n};\n";
  return $string;
}

sub write_tloop {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my $tl   = shift;

  return undef if not $type eq "both";

  my %tl_e  = %{$tl->{'energies'}};
  my %tl_dh = %{$tl->{'enthalpies'}};
  my @common_keys = sort { $a cmp $b } grep { exists $tl_dh{$_} } keys %tl_e;

  my $string = "";
  $string .= "# " . $id;
  $string .= $self->{'suffix'};
  $string .= "\n";

  foreach my $i (@common_keys){
    $string .= sprintf("\t%.8s\t%6.0f\t%6.0f\n", $i, $tl_e{$i}*100, $tl_dh{$i}*100);
  }
  return $string;
}

sub write_tloop_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my $tl   = shift;

  return "" if not $type eq "both";

  my %tl_e  = %{$tl->{'energies'}};
  my %tl_dh = %{$tl->{'enthalpies'}};
  my $charlength;
  $charlength = 40*6 + 1 if $id eq "Triloops";
  $charlength = 40*7 + 1 if $id eq "Tetraloops";
  $charlength = 40*9 + 1 if $id eq "Hexaloops";

  my @common_keys = sort { $a cmp $b } grep { exists $tl_dh{$_} } keys %tl_e;

  my $string = "PUBLIC int ".$C_varnames{$id}."37".$self->{'suffix'}."[40] = ";
  my @vals = ();
  foreach my $i (@common_keys){
    push @vals, sprintf(" %6.0f", $tl_e{$i}*100);
  }
  $string .= "{" . join(",", @vals) . "};\n";

  $string .= "PUBLIC int ".$C_varnames{$id}."dH".$self->{'suffix'}."[40] = ";
  @vals = ();
  foreach my $i (@common_keys){
    push @vals, sprintf(" %6.0f", $tl_dh{$i}*100);
  }
  $string .= "{" . join(",", @vals) . "};\n\n";

  $string .= "PUBLIC char ".$id.$self->{'suffix'}."[".$charlength."] = \n";
  foreach my $i (@common_keys){
    $string .= sprintf("  \"%.8s \"\n", $i);
  }
  $string .= "  \"\";\n";

  return $string;
}

sub write_loop {
  my $self = shift;
  my $id   = shift;
  my $type = shift;

  my $string = "";
  $string .= "# " . $id . "_" . $type if $type eq "enthalpies";
  $string .= "# " . $id if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "\n";

  my @a = @{$_[0]};
  my $s=0;
  foreach (@a) {
    if ($_ eq '.') {
      $string .= '   INF';
    } else {
      $string .= sprintf(" %6.0f", $_*100);
    }
    $string .= "\n" if (++$s % 10) == 0;
  }
  $string .= "\n" if ($s % 10) !=0;

  return $string;
}

sub write_loop_C {
  my $self = shift;
  my $id   = shift;
  my $type = shift;
  my @a = @{$_[0]};

  my $string = "PUBLIC int ".$C_varnames{$id};
  $string .= "dH" if $type eq "enthalpies";
  $string .= "37" if $type eq "energies";
  $string .= $self->{'suffix'};
  $string .= "[31] = ";

  my @vars;
  foreach (@a) {
    if ($_ eq '.') {
      push @vars, "   INF";
    } else {
      push @vars, sprintf(" %6.0f", $_*100);
    }
  }
  $string .= "{" . join(",", @vars) . "};\n";

  return $string;
}

sub write_misc {
  my $self = shift;
  my $id    = shift;
  my $type  = shift;
  my $dat   = shift;

  return undef if not $type eq "both";

  my @ml_e      = @{$dat->{'energies'}{multi}};
  my $termAU_e  = $dat->{'energies'}{termAU};
  my @ml_dh     = @{$dat->{'enthalpies'}{multi}};
  my $termAU_dh = $dat->{'enthalpies'}{termAU};
  my @ninio_e   = @{$dat->{'energies'}{ninio}};
  my @ninio_dh  = @{$dat->{'enthalpies'}{ninio}};
  my $maxninio  = $dat->{'energies'}{max_ninio};
  my $duplex_e  = $dat->{'energies'}{init};
  my $duplex_dh = $dat->{'enthalpies'}{init};
  my $lxc_e     = $dat->{'energies'}{extrapolation};
  my $lxc_dh    = $dat->{'enthalpies'}{extrapolation};

  my $string = "";
  $string .= "# ML_params";
  $string .= $self->{'suffix'};
  $string .= "\n";
  $string .= "/* F = cu*n_unpaired + cc + ci*loop_degree (+TermAU) */\n";
  $string .= "/*\t  cu\t cu_dH\t    cc\t cc_dH\t    ci\t ci_dH  */\n";
  warn "wrong number of ML_params\n" unless  @ml_e == 3 or @ml_dh == 3;
  $string .= sprintf("\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%6.0f\n", 
    $ml_e[1]*100, $ml_dh[1]*100, $ml_e[0]*100, $ml_dh[0]*100, $ml_e[2]*100, $ml_dh[2]*100);

  $string .= "\n";
  $string .= "# NINIO";
  $string .= $self->{'suffix'};
  $string .= "\n";
  $string .= "/* Ninio = MIN(max, m*|n1-n2| */\n";
  $string .= "/*\t   m\t  m_dH    max  */\n";
  $string .= sprintf("\t%6.0f\t%6.0f   %4d\n", $ninio_e[1]*100, $ninio_dh[1]*100, $maxninio*100);

  $string .= "\n";
  $string .= "# Misc";
  $string .= $self->{'suffix'};
  $string .= "\n";
  $string .= "/* all parameters are pairs of 'energy enthalpy' */\n";
  $string .= "/*    DuplexInit     TerminalAU      LXC */\n";
  $string .= sprintf("\t%6.0f\t%6.0f\t%6.0f\t%6.0f\t%6.2f\t%6.2f\n",
    $duplex_e*100, $duplex_dh*100, $termAU_e*100, $termAU_dh*100, $lxc_e*100, $lxc_dh*100);

  return $string;
}

sub write_misc_C {
  my $self = shift;
  my $id    = shift;
  my $type  = shift;
  my $dat   = shift;

  return undef if not $type eq "both";

  my @ml_e      = @{$dat->{'energies'}{multi}};
  my $termAU_e  = $dat->{'energies'}{termAU};
  my @ml_dh     = @{$dat->{'enthalpies'}{multi}};
  my $termAU_dh = $dat->{'enthalpies'}{termAU};
  my @ninio_e   = @{$dat->{'energies'}{ninio}};
  my @ninio_dh  = @{$dat->{'enthalpies'}{ninio}};
  my $maxninio  = $dat->{'energies'}{max_ninio};
  my $duplex_e  = $dat->{'energies'}{init};
  my $duplex_dh = $dat->{'enthalpies'}{init};
  my $lxc_e     = $dat->{'energies'}{extrapolation};
  my $lxc_dh    = $dat->{'enthalpies'}{extrapolation};
  my $cslope_e  = $dat->{'energies'}{cslope};
  my $cslope_dh = $dat->{'enthalpies'}{cslope};
  my $cinter_e  = $dat->{'energies'}{cinter};
  my $cinter_dh = $dat->{'enthalpies'}{cinter};
  my $tripleC_e  = $dat->{'energies'}{ch3};
  my $tripleC_dh = $dat->{'enthalpies'}{ch3};

  my $string = "";
  $string .= "PUBLIC double lxc37".$self->{'suffix'}."     = ".sprintf("%6.2f",$lxc_e*100).";\n";
  $string .= "PUBLIC double lxcdH".$self->{'suffix'}."     = ".sprintf("%6.2f",$lxc_dh*100).";\n";
  $string .= "PUBLIC int ML_intern37".$self->{'suffix'}."  = ".sprintf("%6.0f",$ml_e[2]*100).";\n";
  $string .= "PUBLIC int ML_interndH".$self->{'suffix'}."  = ".sprintf("%6.0f",$ml_dh[2]*100).";\n";
  $string .= "PUBLIC int ML_closing37".$self->{'suffix'}." = ".sprintf("%6.0f",$ml_e[0]*100).";\n";
  $string .= "PUBLIC int ML_closingdH".$self->{'suffix'}." = ".sprintf("%6.0f",$ml_dh[0]*100).";\n";
  $string .= "PUBLIC int ML_BASE37".$self->{'suffix'}."    = ".sprintf("%6.0f",$ml_e[1]*100).";\n";
  $string .= "PUBLIC int ML_BASEdH".$self->{'suffix'}."    = ".sprintf("%6.0f",$ml_dh[1]*100).";\n";
  $string .= "PUBLIC int MAX_NINIO".$self->{'suffix'}."    = ".sprintf("%6.0f",$maxninio*100).";\n";
  $string .= "PUBLIC int ninio37".$self->{'suffix'}."      = ".sprintf("%6.0f",$ninio_e[1]*100).";\n";
  $string .= "PUBLIC int niniodH".$self->{'suffix'}."      = ".sprintf("%6.0f",$ninio_dh[1]*100).";\n";
  $string .= "PUBLIC int TerminalAU37".$self->{'suffix'}." = ".sprintf("%6.0f",$termAU_e*100).";\n";
  $string .= "PUBLIC int TerminalAUdH".$self->{'suffix'}." = ".sprintf("%6.0f",$termAU_dh*100).";\n";
  $string .= "PUBLIC int DuplexInit37".$self->{'suffix'}." = ".sprintf("%6.0f",$duplex_e*100).";\n";
  $string .= "PUBLIC int DuplexInitdH".$self->{'suffix'}." = ".sprintf("%6.0f",$duplex_dh*100).";\n\n";
  $string .= "PUBLIC int TripleC37".$self->{'suffix'}."    = ".sprintf("%6.0f",$tripleC_e*100).";\n";
  $string .= "PUBLIC int TripleCdH".$self->{'suffix'}."    = ".sprintf("%6.0f",$tripleC_dh*100).";\n";
  $string .= "PUBLIC int MultipleCA37".$self->{'suffix'}." = ".sprintf("%6.0f",$cslope_e*100).";\n";
  $string .= "PUBLIC int MultipleCAdH".$self->{'suffix'}." = ".sprintf("%6.0f",$cslope_dh*100).";\n";
  $string .= "PUBLIC int MultipleCB37".$self->{'suffix'}." = ".sprintf("%6.0f",$cinter_e*100).";\n";
  $string .= "PUBLIC int MultipleCBdH".$self->{'suffix'}." = ".sprintf("%6.0f",$cinter_dh*100).";\n\n";
  $string .= "PUBLIC int GQuadAlpha37".$self->{'suffix'}." = -1800;\n";
  $string .= "PUBLIC int GQuadAlphadH".$self->{'suffix'}." = -11934;\n";
  $string .= "PUBLIC int GQuadBeta37".$self->{'suffix'}."  = 1200;\n";
  $string .= "PUBLIC int GQuadBetadH".$self->{'suffix'}."  = 0;\n";

  return $string;
}

1;

__END__

=head1 AUTHOR

Ronny Lorenz (ronny@tbi.univie.ac.at)

=head1 NAME

RNA::Params - Convert mfold/RNAstructure RNA/DNA parameter files into
ViennaRNA Package format

=head1 MODULE

B<Use as Perl module:>

  use RNA::Params;
  my $params = RNA::Params->new();
  $params->addFiles("datafile_dir");
  $params->print_parameters("RNA");

=head1 SYNOPSIS

B<Use as an executable perl script:>

  Params.pm [options]

=head1 DESCRIPTION

This package provides various subroutines for parsing, extracting, and
converting different RNA/DNA energy paramter files. In particular, this
package provides functions to convert energy parameter files as shipped
with the RNAstructure program package into ViennaRNA Package formatted
parameter files, and C-source code.

=head1 METHODS

=over 4

=item B<new([%options])>

Create a new RNA::Params object.
The optional I<%options> hash may be used to set set some non-default
operating settings.
These are the default options:

  %default_options = (
    'verbose'           => 0,   # print some verbose output if != 0
    'suffix'            => ""   # Add this string as suffix for descriptor
                                # or variable names in the output
  );

=item B<addFiles(filelist)>

Add input files and read their content.
The argument I<filelist> can either be a scalar string that contains the
file to parse, or a string containing comma separated multiple file names.
The method also allows for multiple such scalars provided as an array.

=item B<print_parameters(type[,fh])>

Print a ViennaRNA Parameter file.
The I<type> parameter must be one of "RNA","DNA", or "HYBRID" to print a
specific set of corresponding parameters. For printing the output to a
file instead of I<stdout>, the optional parameter I<fh> may be used to
provide a filehandle.

=item B<print_parameters_C(type)>

Print ViennaRNA Parameter C - source code files
This method creates C source code files ready to be included in the
ViennaRNA Package source code tree. Again, the I<type> parameter must
be one of "RNA","DNA", or "HYBRID" to write a specific set of
corresponding parameter files.

=back

=head1 OPTIONS

=over 4

=item B<--help>
Show help message

=item B<-c>

Generate C source code

=item B<-p>

Generate Parameter file

=item B<-o> <filename>

Print output to file

=item B<-i> <string>[,<string]

Input files to be processed. The files can either be a comma separated
list of files, or even directory names. Directories will not be parsed
recursively.

=item B<--dna,-d>

Print DNA parameters

=item B<--rna,-r>

Print RNA parameters

=item B<--hybrid,-h>

Print RNA/DNA hybrid parameters

=item B<--suffix> <string>

Suffix for descriptor/variable names

=item B<--verbose,-v>

Be verbose

=back

=cut


