#!/usr/bin/perl -w
# -*-CPerl-*-
# Last changed Time-stamp: <2003-07-16 17:12:00 ivo>

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

# End of file
