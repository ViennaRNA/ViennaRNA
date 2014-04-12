#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2005-11-07 12:08:38 ivo>
# $Id: colorrna.pl,v 1.1 2005/11/07 12:42:27 ivo Exp $
# colorize a secondary structure plot with reliability annotation
# from positional entropy
use strict;
use Getopt::Std;
$main::VERSION = 1.0;
$Getopt::Std::STANDARD_HELP_VERSION=1;

our $opt_p;
getopts('p');

sub HELP_MESSAGE {
  print STDERR "\nusage: $0 FOO_ss.ps FOO_dp.ps > FOO_css.ps\n";
  print STDERR "For more details run\n\tperldoc -F $0\n";
}

HELP_MESSAGE() unless $#ARGV >0;
my $macro_seen= 0;
my %mfe = ();        # hash of mfe pairs
my @ss_ps = ('',''); # head and tail of the ss.ps file

my $n = swallow_ss_ps();    # read ss plot
my %hsb = read_dp();   # read dot plot and compute entropies

print $ss_ps[0];     # print head
if (!$macro_seen) {
  print <<_E_O_F_
/hsb {
    dup 0.3 mul 1 exch sub sethsbcolor
} bind def
/colorpair { % i j hue sat colorpair
   % draw basepair i,j in color
   % 1 index 0.00 ne {
   gsave
   newpath
   hsb
   fsize setlinewidth
   1 sub coor exch get aload pop moveto
   1 sub coor exch get aload pop lineto
   stroke
   grestore
   % } if
} bind def
_E_O_F_
}
foreach my $k (keys %hsb) {
  my ($i,$j) = split($;,$k);
  print "$i $j ", $hsb{$k}->[0], " ", $hsb{$k}->[1]," colorpair\n";
}
print $ss_ps[1];  # print main part
print $ss_ps[2];  # print showpage etc

sub swallow_ss_ps {
  # read the secondary structure plot
  my $length=0;
  my $tail=0;
  while (<>) {
    $macro_seen=1 if /colorpair/;
    $length ++ if /^\/coor/ .. /^\] def/;
    if (/^\/pairs/ .. /^\] def/) {
      $mfe{$1,$2}=1 if /(\d+)\s+(\d+)/;
    }
    next if /\d gmark$/;
    $tail++ if /^drawpairs/;
    s/^drawpairs/% drawpairs/;
    $tail++ if /^showpage/;
    $ss_ps[$tail] .= $_;
    last if eof;
  }
  return $length-2;
}

sub read_dp {
  # read dot plot and extract hsb values
  my %hsb;
  while (<>) {
    next unless /lbox$/;
    my @F = split;
    $hsb{$F[3],$F[4]} = [$F[0], $F[1]];
  }
  return %hsb;
}

=head1 NAME

colorrna - colorize an alirna.ps file

=head1 SYNOPSIS

   colorrna.pl alirna.ps alidot.ps > colorRNA.ps

=head1 DESCRIPTION

colorrna reads an RNA secondary structure plot and a dot plot
containing pair probabilities and covariance annotation as, produced
by C<RNAalifold -p>, and writes a new secondary structure plot with
color annotated seqence annotation to stdout.  The color annotation
is taken directly from the color dot plot file.

=head1 AUTHOR

Ivo L. Hofacker <ivo@tbi.univie.ac.at>
Stefan Washietl <wash@tbi.univie.ac.at>

=cut

#  End of file

