#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2003-07-18 09:42:38 xtof>
# $Id: rotate_ss.pl,v 1.1 2003/07/18 07:53:20 xtof Exp $

use Getopt::Long;
use strict;
use vars qw/$ss_ps $opt_a/;
use constant PI => 3.14159265;

my $ss_ps = { Header  => [],
	      Coords  => [],
	      BPairs  => [],
	      Tailer  => [] };
$opt_a = 0;
usage() unless GetOptions("a=i" => \$opt_a);

swallow_ss_ps();
my $ar = to_Array($ss_ps->{Coords});
print_ss_ps( rot_Array( get_Midpt($ar), $ar ) );

#---
sub usage {
  printf STDERR "\nusage: $0 [-a ANGLE] FOO_ss.ps > BAR_ss.ps\n";
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

__END__
