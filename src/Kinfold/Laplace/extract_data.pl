#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2006-10-04 12:11:12 xtof>
# $Id: extract_data.pl,v 1.1 2006/10/05 13:29:03 xtof Exp $

use strict;

while (<>) {
  next if m/^\#/;
  next if m/\)$/;
  next if m/\.$/;
  my @F = split;
  print "$F[-2] $F[-1]\n" if scalar @F > 3;
}

__END__
