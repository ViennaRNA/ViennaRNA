#!/usr/local/bin/perl
# zoom into a PostScript dot plot produced by RNAfold 
# Last changed Time-stamp: <1998-06-27 22:00:20 ivo> 
sub usage {
    die "Usage: $0 [-f first] [-l last] [dp_file]\n";
}

$from=1;

# parse command line
while ($_ = shift) {
    /^-[h\?]/ && &usage;
    /^-f$/ && ($from = shift, next );
    /^-l$/ && ($to   = shift, next );
    unshift (@ARGV, $_);
    last;
}

while (<>) {
    if (/\/sequence \{ \((\S*)[\\\)]/) {  
	$seq = $1;              # empty for new version
	while (!/\) \} def/) {  # read until end of definition
	    $_ = <>;
	    /(\S*)[\\\)]/;      # ends either with `)' or `\' 
	    $seq .= $1;
	}
	if ($to>$from) { $seq=substr($seq,$from-1,$to-$from+1) }
	else { $seq=substr($seq,$from-1); }
	print  "% Subsequence from $from to $to\n";
	printf "/sequence { (%s) } def\n", $seq;
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
