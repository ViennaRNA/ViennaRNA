#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <1999-09-29 13:24:23 ivo>
# produce dot bracket notation of an RNA secondary structure from
# Zuker's .ct file

while (<>) {
    @F = split;
    if ((/ = /) && ($#F!=5)) {
	if (defined($seq)) {
	    print "$seq\n" unless ($seq eq $oldseq);
	    $oldseq = $seq;
	    print "$s ($E)\n";
	}
	$E = $F[3];
	$s = $seq = "";
	next;
    }
    $seq .= $F[1];
    if ($F[4]==0) {$s .= "."; next}
    $s .= "(" if ($F[0]<$F[4]);
    $s .= ")" if ($F[0]>$F[4]);
}
print "$seq\n" unless ($seq eq $oldseq);
print "$s ($E)\n" if defined($s);

# End of file
