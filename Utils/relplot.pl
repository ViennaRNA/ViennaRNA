#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2003-07-19 19:54:11 ivo>

# colorize a secondary structure plot with reliability annotation
# from positional entropy
use strict;

sub usage {
  printf STDERR "\nusage: $0 FOO_ss.ps FOO_dp.ps > FOO_rss.ps\n";
  exit(1);
}

usage() unless $#ARGV >0;
my $macro_seen= 0;
my @ss_ps = ('',''); # head and tail of the ss.ps file

my $n = swallow_ss_ps();    # read ss plot
my @sp = posent();   # read dot plot and compute entropies

print $ss_ps[0];     # print head
if (!$macro_seen) {
  print <<_E_O_F_
/drawreliability {
  /Smax {sequence length log 5.5 min} bind def
  0
  coor {
    aload pop
    S 3 index get Smax div
    0.9 min 1 1 sethsbcolor
    newpath
    fsize 2 div 0 360 arc
    fill
    1 add
  } forall
} bind def
_E_O_F_
}
print "/S [\n";
foreach (@sp) {
  printf "  %7.5f\n", $_;
}
print "] def\n\n";
print "drawreliability\n";
print $ss_ps[1];  # print tail

sub swallow_ss_ps {
  # read the secondary structure plot
  my $length=0;
  my $tail=0;
  while (<>) {
    $macro_seen=1 if /drawreliability /;
    $length ++ if /^\/coor/ .. /^\] def/;
    $tail=1 if /^drawoutline/;
    $ss_ps[$tail] .= $_;
    last if eof;
  }
  return $length-2;
}

sub posent {
  # compute positional entropy from pair probs in the dot plot file
  my @pp;
  my @sp;
  while (<>) {
    my ($i, $j, $p, $id) = split;
    next unless defined $id && $id eq 'ubox';
    $p *= $p;
    my $ss = $p*log($p);
    $sp[$i] += $ss;
    $sp[$j] += $ss;
    $pp[$i] += $p;
    $pp[$j] += $p;
  }
  my $log2 = log(2);
  for my $i (1..$n) {
    no warnings;  # $p[$i] may be undef
    $sp[$i] += (1-$pp[$i])*log(1-$pp[$i]);
    $sp[$i] /= -$log2;
  }
  shift @sp; # get rid of [0] entry
  return @sp;
}

=head1 NAME

relplot - annotate a secdonary structure plot with reliability information

=head1 SYNOPSIS

   relplot file_ss.ps file_dp.ps > file_rss.ps

=head1 DESCRIPTION

relplot reads an RNA secondary structure plot and a dot plot
containing pair probabilities, as produces by C<RNAfold -p>. From the
pair probabilities it computes the "positional entropy" C<S(i) = - Sum
p(ij) log(p(ij))> which is then used to colorize the secondary
structure plot. Low entropy regions have little structural flexibility
and the reliability of the predicted structure is high. High entropy
implies many structural alternatives. While these alternatives may be
functionally important, they make structure prediction more difficult
and thus less reliable. The new colorized postscript file is written
to stdout.

Entropy is encoded as color hue, ranging from red for low entropy,
well-defined regions, via yellow and green to blue and violet for
regions with very high entropy.

=head1 AUTHOR

Ivo L. Hofacker <ivo@tbi.univie.ac.at>

=cut

#  End of file

