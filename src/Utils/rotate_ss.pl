#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2004-08-09 16:50:14 ivo>
# $Id: rotate_ss.pl,v 1.4 2006/02/28 14:44:14 ivo Exp $

use Getopt::Long;
use strict;
use vars qw/$ss_ps $opt_a $opt_m/;
use constant PI => 3.14159265;

my $ss_ps = { Header  => [],
	      Coords  => [],
	      BPairs  => [],
	      Tailer  => [] };
$opt_a = 0;
usage() unless GetOptions("a=f" => \$opt_a, "m");

swallow_ss_ps();
my $ar = to_Array($ss_ps->{Coords});
my $mp = get_Midpt($ar);
$ar = flip_Array($mp, $ar) if $opt_m;
$ar = rot_Array($mp, $ar) if $opt_a;
print_ss_ps($ar);

#---
sub usage {
  printf STDERR "\nusage: $0 [-a ANGLE] [-m] FOO_ss.ps > BAR_ss.ps\n";
  exit(1);
}

#---
sub swallow_ss_ps {
  local $_;
  my $toggle = 0;

  while (<>) {
    $toggle = 1 if m/^\/coor/;
    push @{$ss_ps->{Header}}, $_ if $toggle == 0;
    push @{$ss_ps->{Coords}}, $_ if $toggle == 1;
    push @{$ss_ps->{BPairs}}, $_ if $toggle == 2;
    push @{$ss_ps->{Tailer}}, $_ if $toggle >  2;
    $toggle++ if m/\]\s+def/;
  }
}

#---
sub to_Array {
  my $aref = shift;
  return [ map {chomp; s/\[//; s/\]//; [split] } @$aref[$[+1 .. $#$aref-1] ];
}

#---
sub get_Midpt {
  local $_;
  my $aref = shift;
  my ($xl, $yl, $xu, $yu) = ( 100000, 100000, -100000, -100000 );
  for ( @$aref ) {
    $xl = ($_->[0] < $xl) ? $_->[0] : $xl;
    $yl = ($_->[1] < $yl) ? $_->[1] : $yl;
    $xu = ($_->[0] > $xu) ? $_->[0] : $xu;
    $yu = ($_->[1] > $yu) ? $_->[1] : $yu;
  }
  return [($xl+$xu)/2, ($yl+$yu)/2, $xl, $yl, $xu, $yu];
}

#---
sub rot_Array {
  my ($mp, $ar) = @_;
  my $a = $opt_a/360*2*PI;
  my ($ca, $sa) = (cos($a), sin($a));
  return
      [
       map { [ ($ca*($_->[0]-$mp->[0]) + (-$sa*($_->[1]-$mp->[1])))+$mp->[0],
	       ($sa*($_->[0]-$mp->[0]) +   $ca*($_->[1]-$mp->[1])) +$mp->[1]
	       ] } @$ar ];
}
#---
sub flip_Array {
  my ($mp, $ar) = @_;
  return
      [ map { [ ($mp->[0]-$_->[0]), $_->[1]] } @$ar ];
}

#---
sub print_ss_ps {
  local $_;
  my $ar = shift;
  for (qw/Header Coords BPairs Tailer/) {
    if ($_ eq 'Coords') {
      print "/coor [\n";
      for my $xy (@$ar) {
        print sprintf("[%7.3f %7.3f]\n", $xy->[0], $xy->[1]);
      }
      print "] def\n";
      next;
    }
    print join "", @{$ss_ps->{$_}};
  }
}

=head1 NAME

rotate_ss - rotate or mirror coordinates of secondary structure plot

=head1 SYNOPSIS

   rotate_ss.pl [-a angle] [-m] old_ss.ps > new_ss.ps

=head1 DESCRIPTION

B<rotate_ss> reads a PostScript RNA secondary structure plot, as
produced by B<RNAfold> or B<RNAplot>, rotates and/or mirrors the
coordinates, and writes the new plot to STDOUT.

=head1 OPTIONS

=over 4

=item B<-a> I<angle>

Rotate the plot counter-clockwise by I<angle> degrees.

=item B<-m>

Mirror the coordinates in the plot, i.e. convert from a
couter-clockwise to clockwise layout. Note that, if both B<-m> and
B<-a> are given, the plot is first mirrored, then rotated.

=back

=head1 AUTHORS

Ivo L. Hofacker <ivo@tbi.univie.ac.at>,
Christoph Flamm <xtof@tbi.univie.ac.at>

=cut

__END__
