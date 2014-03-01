#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2005-03-09 19:38:25 ivo> 
# produce Pauline's mountain representation from bracket notation
# use e.g. as  RNAfold < foo.seq | b2mt | xmgr -pipe
# definition: h=number of base pairs enclosing base
use strict;
while (<>) {
    print if (s/>/#/);

    next unless (/\.\.\./);
    next if (/\[/);   # don't process output from partition function
    chop;
    my @F=split(//,$_);
    my $p=0; my $h=0;
    foreach my $i (@F) {
	$h-- if ($i eq ')');
	$p++;
	printf("%4d %4d\n",$p,$h);
	$h++ if ($i eq '(');	# increase $h *after* printing
    }
    print "&\n";
}

=head1 NAME

b2mt - produce coordinates for a mountain plot from bracket notation

=head1 SYNOPSIS

  b2mt.pl < seq.fold > mountain.dat

=head1 DESCRIPTION

read a secondary structures in bracket notation as output by RNAfold,
and compute coordinates for a mountain plot as introduced by Pauline Hogeweg.
The output is suitable for graphing  with xmgrace, e.g.:
C< RNAfold < foo.seq | b2mt.pl | xmgrace -pipe>

=head1 AUTHOR

Ivo L. Hofacker <ivo@tbi.univie.ac.at>
