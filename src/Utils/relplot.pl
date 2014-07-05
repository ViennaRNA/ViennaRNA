#!/usr/bin/perl -w
#
# colorize a secondary structure plot with reliability annotation
# from positional entropy
#
use strict;
use Getopt::Std;
$main::VERSION = 1.3;
$Getopt::Std::STANDARD_HELP_VERSION=1;

our ($opt_p, $opt_a, $opt_s);
getopts('pas');

sub HELP_MESSAGE {
  print STDERR "\nusage: $0 [-p] [-a] FOO_ss.ps FOO_dp.ps > FOO_rss.ps\n";
  print STDERR "\nusage: $0 -s FOO_ss.ps FOO.shape > FOO_rss.ps\n";
  print STDERR "For more details run\n\tperldoc -F $0\n";
}

HELP_MESSAGE() unless $#ARGV >0;
my $macro_seen= 0;
my %mfe = ();         # hash of mfe pairs
my @ss_ps = ('','');  # head and tail of the ss.ps file
my $cut_point = -1;   # account for shift of base pair positions 
                      # between rna.ps and dot.ps upon cofolding

my $n = swallow_ss_ps();    # read ss plot
my @sp = $opt_s ? swallow_shape() : posent();   # read shape file or dot plot and compute entropies
my $Smax = ($opt_p || $opt_a || $opt_s) ? 1 : 0;
if (!$opt_p) {
  foreach (@sp) {
    $Smax = $_ if $_>$Smax;
  }
  $Smax = ($Smax>0.2) ? sprintf("%3.1f", $Smax) : 0.2 ;

}
print $ss_ps[0];     # print head
if (!$macro_seen) {
  print <<_E_O_F_
/range 0.8 def
/drawreliability {
  /Smax $Smax def
  0
  coor {
    aload pop
    S 3 index get
    Smax div range mul
    invert {range exch sub} if
    dup 0 ge
    {1 1 sethsbcolor} {pop 1 1 1 setrgbcolor} ifelse
    newpath
    fsize 2 div 0 360 arc
    fill
    1 add
  } forall
} bind def
/colorbar { % xloc yloc colorbar -> []
  /STR 8 string def
  gsave
    xmin xmax add size sub 2 div
    ymin ymax add size sub 2 div translate
    size dup scale
    translate
    0.015 dup scale
    /tics 64 def
    gsave
      10 tics div 1 scale
      0 1 tics
      {
          dup 0 moveto 0.5 add
          tics div range mul
          invert {range exch sub} if
          1 1 sethsbcolor
          1 0 rlineto 0 1 rlineto -1 0 rlineto closepath fill
      } for
    grestore
    0 setgray
    -0.1 1.01 moveto (0) gsave 0.1 dup scale show grestore
    10 1.01 moveto Smax STR cvs
    gsave 0.1 dup scale dup stringwidth pop -2 div 0 rmoveto show grestore
  grestore
} bind def
_E_O_F_
}
print "/S [\n";
foreach (@sp) {
  printf "  %7.5f\n", $_;
}
print "] def\n\n";
print "/invert ", $opt_p||$opt_a ? 'true' : 'false', " def\n";
print "drawreliability\n";
print "0.1 0.1 colorbar\n";
print $ss_ps[1];  # print tail

sub swallow_ss_ps {
  # read the secondary structure plot
  my $length=0;
  my $tail=0;
  while (<>) {
    $macro_seen=1 if /drawreliability /;
    if(/^\/sequence/ .. /^\) def/){
      $cut_point = $-[1] if /^[acgtunACGTUN]+(\s)[acgtunACGTUN]+/;
    }
    $length ++ if /^\/coor/ .. /^\] def/;
    if (/^\/pairs/ .. /^\] def/) {
      $mfe{$1,$2}=1 if /(\d+)\s+(\d+)/;
    }
    $tail=1 if /^drawoutline/;
    $ss_ps[$tail] .= $_;
    last if eof;
  }
  return $length-2;
}

sub posent {
  # compute positional entropy from pair probs in the dot plot file
  # or, with $opt_p, find pair probs corresponding to mfe pairs
  my @pp;
  my @sp;
  while (<>) {
    next unless /(\d+)\s+(\d+)\s+([0-9.Ee-]+)\s+ubox/;
    my ($i, $j, $p) = ($1, $2, $3);

    # account for position shift in rna.ps if input comes from cofold
    if($cut_point > 0){
      $i = $i + 1 if $i > $cut_point;
      $j = $j + 1 if $j > $cut_point;
    }

    $p *= $p;
    if ($opt_p) {
      $sp[$i] = $sp[$j] = $p if exists  $mfe{$i,$j};
    } else {
      if ($opt_a) {
        $sp[$i] = $sp[$j] = 1-$p if exists  $mfe{$i,$j};
      } else {
        my $ss = ($p>0)?$p*log($p):0;
        $sp[$i] += $ss;
        $sp[$j] += $ss;
      }
    }
    $pp[$i] += $p;
    $pp[$j] += $p;
  }
  my $log2 = log(2);
  for my $i (1..$n) {
    no warnings;  # $p[$i] may be undef
    if ($opt_p || $opt_a) {
      $sp[$i] = 1-$pp[$i] if !defined $sp[$i];
    } else {
      $sp[$i] += ($pp[$i]<1) ? (1-$pp[$i])*log(1-$pp[$i]) : 0;
      $sp[$i] /= -$log2;
    }
  }
  shift @sp; # get rid of [0] entry
  return @sp;
}

sub swallow_shape {
  my %shape;
  my @shapearray;
  my $max = 0;
  while(<>) {
    chomp;

    my @columns = split(/\s+/, $_);
    my $length = scalar @columns;
    next if not $length;

    my $index = int($columns[0]);
    my $value = $length > 1 ? $columns[-1] : -1;
    next if $value ne $value + 0;
    $shape{$index} = $value;
    $max = $index if $index > $max;
  }

  for(my $i = 1; $i <= $max; $i++)
  {
    push(@shapearray, exists $shape{$i} ? $shape{$i} : -1);
  }

  return @shapearray;
}

=head1 NAME

relplot - annotate a secondary structure plot with reliability information

=head1 SYNOPSIS

   relplot [-p] [-a] file_ss.ps file_dp.ps > file_rss.ps
   relplot -s file_ss.ps file.shape > file_rss.ps

=head1 DESCRIPTION

relplot reads an RNA secondary structure plot and a dot plot
containing pair probabilities, as produced by C<RNAfold -p> or C<RNAcofold -p>, and
writes a new secondary structure with reliability annotation to
stdout.  The anotation is used to colorize the plot and can use either
"positional entropy" (default), or pair probabilities (with -p).

Positional entropies are computed from the pair probabilities as
C<S(i) = - Sum_i p(ij) log(p(ij))>. Low entropy regions have little
structural flexibility and the reliability of the predicted structure
is high. High entropy implies many structural alternatives. While
these alternatives may be functionally important, they make structure
prediction more difficult and thus less reliable.

If the -p switch is given, the script colors base pairs by their pair
probability, unpaired bases use the probability of being unpaired.

If the -a switch is given, the script colors base pairs by their accessibility, i.e. the probability of being unpaired.

Entropy (repsectively probability) is encoded as color hue, ranging
from red for low entropy, well-defined regions, (high probability
pairs) via yellow and green to blue and violet for regions with very
high entropy (low probability).

If the -s switch is given, the second given file will be interpreted as SHAPE file
and reactivity values defined in the file will be used to annotate the secondary structure plot.

You may have to manually move the color legend to a convenient
position. Just edit the postscript file and change the two numbers in the
line reading C<0.1 0.1 colobar>. Or delete the line to remove the legend.

=head1 AUTHOR

Ivo L. Hofacker <ivo@tbi.univie.ac.at>

=cut

#  End of file
