#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2003-07-18 09:50:16 ivo>
# produce Pauline Hogeweg's mountain representation *_dp.ps files
# writes 3 sets of x y data separated by a "&"
# the first two sets are mountain representations from base pair probabilities
# and mfe structure, respectively.
# definition: mm[i],mp[i] = (mean) number of base pairs enclosing base i
# third set a measure of well-definedness: the entropy of the pair probs of
# base i, sp[i] = -Sum p_i * ln(p_i). Well-defined regions have low entropy.
#
# use e.g. as  mountain.pl dot.ps | xmgrace -pipe

use strict;
our (@mm, @mp, @pp, @sp, @p0, $i, $max, $length);  # perl5 only

my $sep = "&";   # xmgr uses & to separate data sets  

while (<>) {
    my ($seq,$i,$j,$p,$id);
    chop;
    if (/\/sequence \{ \((\S*)[\\\)]/) {  
        $seq = $1;              # empty for new version
        while (!/\) \} def/) {  # read until end of definition
            $_ = <>;
            /(\S*)[\\\)]/;      # ends either with `)' or `\' 
            $seq .= $1;
        }
	$length = length($seq);
	next;
    }

    next unless /box$/;    
    ($i, $j, $p, $id) = split;
    if ($id eq "ubox") {
	$p *= $p;           # square it to probability
	$mp[$i+1] += $p; 
	$mp[$j]   -= $p;
	my $ss = $p*log($p);
	$sp[$i] += $ss;
	$sp[$j] += $ss;
	$pp[$i] += $p;
	$pp[$j] += $p;
    }
    if ($id eq "lbox") {
	$mm[$i+1]++;
	$mm[$j]--;
    }
}
$mp[0] = $mm[0] = $max = 0;
for ($i=1; $i<=$length; $i++) {
    no warnings;
    $mp[$i]+=$mp[$i-1];
    $max = $mp[$i] if ($mp[$i]>$max);
    $mm[$i]+=$mm[$i-1];
    $max = $mp[$i] if ($mp[$i]>$max);
    $sp[$i] += (1-$pp[$i])*log(1-$pp[$i]);
}

# print the results for plotting
for ($i=1; $i<=$length; $i++) {
    printf("%4d  %7.5g\n", $i, $mp[$i]);
}			
print "$sep\n";

for ($i=1; $i<=$length; $i++) {
    printf("%4d  %4d\n", $i, $mm[$i]);
}	
print "&\n";
my $log2 = log(2);
for ($i=1; $i<=$length; $i++) {
    printf("%4d  %7.5g\n", $i, -$sp[$i]/$log2);
}	
