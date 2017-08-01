#!/usr/bin/perl -w
# zoom into a PostScript dot plot produced by RNAfold, RNAalifold or alidot
# Last changed Time-stamp: <2003-03-24 12:47:26 ivo> 
sub usage {
  die "Usage: $0 [-f first] [-l last] [dp_file]\n";
}

$from=1;
$to = -1;
# parse command line
while ($_ = shift) {
  /^-[\?h]/ && usage();
  /^-f$/ && ($from = shift, next );
  /^-l$/ && ($to   = shift, next );
  unshift (@ARGV, $_);
  last;
}

usage() if eof();

while (<>) {
  if (/\/sequence \{ \((\S*)[\\\)]/) {
    $seq = $1;              # empty for new version
    while (!/\) \} def/) {  # read until end of definition
      $_ = <>;
      /(\S*)[\\\)]/;      # ends either with `)' or `\'
      $seq .= $1;
    }
    if ($to>$from) { $seq=substr($seq,$from-1,$to-$from+1) }
    else { $seq=substr($seq,$from-1) }
    print  "% Subsequence from $from to $to\n";
    print "/sequence { (\\\n";
    for ($p=0; $p<length($seq); $p+=250)  {
      print substr($seq,$p, 250), "\\\n";
    }
    print ") } def\n";
  }
  elsif (/(\d+) (\d+) (.*box$)/) {
    $i=$1; $j=$2;
    next if (($i<$from) || (($j>$to)&&($to>0)));
    $i -= ($from-1);
    $j -= ($from-1);
    printf "%d %d %s\n", $i, $j, $3;
  }
  else {print}
}

=head1 NAME

dpzoom - zoom into a PostScript dot plot produced by RNAfold, RNAalifold or alidot

=head1 SYNOPSIS

  dpzoom [OPTIONS] dot.ps > new_dot.ps

=head1 DESCRIPTION

Read a PostScript dot plot as produced by RNAfold and zoom in to a particular
region. Output is again a PostScript dot plot where only a region for the given
subsequence is shown. The user must provide the index of the first and last
position for the subsequence to be zoomed in.

=head1 OPTIONS

=over 4

=item B<-f> <int>

First nucleotide

=item B<-l> <int>

Last nucleotide

=back

=head1 EXAMPLE

Produce PostScript dot plot for subsequence from position 10 to 50:
C<dzoom.pl -f 10 -l 50 dot.ps>

=head1 AUTHOR

Ivo L. Hofacker <ivo@tbi.univie.ac.at>

=cut


