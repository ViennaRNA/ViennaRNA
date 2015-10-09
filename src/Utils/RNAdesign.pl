#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Pod::Usage;
use RNA::Design;
use Data::Dumper;

# our $VERSION = '0.1';

# ********************** #
# Initialize RNA::Design #
# ...................... #
my $ViennaDesign = RNA::Design->new();

# ******************* #
# Get Input Arguments #
# ................... #

my @structs;
my $constr;
my $length = 0;
my $cutpnt = -1;

# These are also the defaults in RNA::Design.pm
my $optfun = 'eos(1)+eos(2)-2*efe() + 0.3*(eos(1)-eos(2)+0.00)**2';
my ($avoid,$apen) = (['AAAAA','CCCCC','GGGGG','UUUUU'],5);
my $bprobs = ({A => 0.25, C => 0.25, G => 0.25, U => 0.25});

my ($n, $m) = (1, 1e10);
my $startseq =undef;
my $ParamFile=undef;
my $verbose = 0;

GetOptions(
  "o|optfun=s"  => sub{$optfun = check_input_costfunction($_[1])},
  "a|avoid=s"   => sub{($avoid,$apen) = check_input_avoid($_[1])},
  "bprobs=s"    => sub{$bprobs = check_input_base_probs($_[1])},
  "n|number=i"  => \$n, # number of independent runs
  "m|maxiter=i" => \$m, # maximal number of steps
  "s|start=s"   => \$startseq,

  "P|params=s"  => \$ParamFile,
  "v|verbose=i" => \$verbose,

  "h|help"      => sub{pod2usage(1)},
  "man"         => sub{pod2usage(-exitstatus => 1, -verbose => 2)},
) or pod2usage(1);

read_input(\@structs, \$constr, \$length, \$cutpnt);

$constr = 'N' x $length unless $constr;

# ******************* #
# Set Input Arguments #
# ................... #
$ViennaDesign->set_cut_point($cutpnt);
$ViennaDesign->set_structures(@structs);
$ViennaDesign->set_constraint($constr);
$ViennaDesign->set_optfunc($optfun);
$ViennaDesign->set_base_probs($bprobs);
$ViennaDesign->set_avoid_motifs(@$avoid);
$ViennaDesign->set_avoid_penalty($apen);
RNA::read_parameter_file($ParamFile) if ($ParamFile);
#alternative: $ViennaDesign->set_parameter_file($ParamFile) if $ParamFile;

if ($verbose) {
  print "Option-check:\n";
  print "Avoid: @{$ViennaDesign->get_avoid_motifs}".
  " with penalty: ".$ViennaDesign->get_avoid_penalty."\n";
  print "Base-prob $_ = $$bprobs{$_}\n" foreach keys %$bprobs;
  print "Objective Function: ".$ViennaDesign->get_optfunc."\n";
  $ViennaDesign->set_verbosity($verbose);
}

$ViennaDesign->find_dependency_paths;
$ViennaDesign->explore_sequence_space;
#$ViennaDesign->eval_sequence($ViennaDesign->find_a_sequence);

# **************** #
# Main Design Loop #
# ................ #

my %final;
my ($seq, $cost, $count);
for (1 .. $n) {
  $seq = ($startseq) ? $startseq : $ViennaDesign->find_a_sequence;
  $seq = $ViennaDesign->optimize_sequence($seq, $m);

  if (exists $final{$seq}) {
    $cost = $final{$seq}[0];
    $final{$seq}[1]++;
  } else {
    $cost = $ViennaDesign->eval_sequence($seq);
    $final{$seq}[0] = $cost;
    $final{$seq}[1]++;
  }
}

# ************* #
# Print Results #
# ............. #

my $id=1;
foreach $seq (sort {$final{$a}[0]<=>$final{$b}[0]} keys %final) {
  my $cseq = $seq;
  substr $cseq, $cutpnt-1,0,'&' if ($ViennaDesign->get_cut_point != -1);
  printf "%3d %s %6.2f %d\n", $id++, $cseq, $final{$seq}[0], $final{$seq}[1];
}

#############################

sub read_input {
  my ($str,$con,$len,$cut) = @_;

  while (<>) {
    chomp;
    next if $_ eq '';

    # Check if length is consistent
    # TODO treat sequeces with '$' different to allow for cotranscr
    $$len = length $_ unless $$len;
    if (length $_ != $$len) {
      die "All structures/constraints must have the same length!";
    }

    # Set cut-point
    my $tmp = index $_, '&';
    if ($tmp != -1 || $$cut != -1) { # it is a cofold-design!
      ++$tmp;
      if ($$cut == -1) {
        $$cut = $tmp;
      } elsif ($tmp != $$cut) {
        die "cut_points ('&') may not differ in input";
      }
      s/&//; 
      die "only one cut_point ('&') allowed" if s/&//;
    }

    if (m/[\(\)\.x\&]/g) { # is a structure constraint!
      if (m/[^\(\)\.\&x]/g) { # only valid characters!
        die "$_ contains forbidden character: $&";
      } 
      push @$str, $_;
    } elsif (m/[ACUGTURYSMWKVHDBN]/) { # is a sequence constraint!
      if (m/[^ACUGTURYSMWKVHDBN&]/g) { # only valid characters!
        die "Constraint contains a non-iupack caracter: $&";
      }
      die "Only one constraint allowed!" if ($$con);
      $$con = $_;
    }
  }
  return 1;
}

sub check_input_costfunction {
  # TODO: need to check this in a better way!
  # TODO: check for large integers (especially in **2)
  my $optfun = shift;

  # TODO: remove these lines once the web is up-to-date
  $optfun =~ s/gfe/efe/g;
  $optfun =~ s/pfc/efe/g;

  my @allowed = qw( + - / * . );
  foreach my $f1 (@allowed) {
    foreach my $f2 (@allowed) {
      my $forbidden = $f1.$f2;
      if ($forbidden eq '**') {
        die "'***' not allowed in objective function" if (index($optfun, '***') != -1);
      } else {
        die "\'$forbidden\' not allowed in objective function" if (index($optfun, $forbidden) != -1);
      }
    }
  }
  my $nakedfun = $optfun;
  foreach my $fb ('eos', 'efe', 'prob', '_circ', 'barr') {
    $nakedfun =~ s/$fb//g;
  }
  die "cannot interpret \'$&\' in objective function" if ($nakedfun =~ m/[^\-\+\/\(\)\.\*\d\,\s]+/g);
  die "unbalanced brackets in optimization function" if scalar($nakedfun =~ tr/\(//) != scalar($nakedfun =~ tr/\)//);

  return $optfun;
}

sub check_input_avoid {
  my $avoid = shift;
  my ($m,$penalty) = split ':', $avoid;
  my @motifs = split ',', $m;

  foreach my $a (@motifs) {
    $a = uc($a);
    die "String $a in 'avoid' list may only contain A,C,U/T or G" 
      if $a =~ m/[^ACUGT]/g;
  }

  return (\@motifs, $penalty);
}

sub check_input_base_probs {
  my $input = shift;
  my %bprobs;

  foreach my $a (split ',', $input) {
    my @tmp = (split ':', $a);
    $bprobs{$tmp[0]}=$tmp[1];
  }

  if (%bprobs) {
    foreach my $b (sort keys %bprobs) {
      die "cannot interpret base: $_ in baseprobabilties" if $b =~ m/[^ACUG]/g || length $b > 1;
      die "bad probability for $b: $bprobs{$b}!" unless $bprobs{$b} =~ m/^\d*\.?\d+$/g;
      die "probability for $b may not be 0!" unless $bprobs{$b};
    }
  }
  return (%bprobs) ? \%bprobs : undef;
}

