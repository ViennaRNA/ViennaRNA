#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <95/05/22 17:52:36 ivo> 
# produce Pauline's mountain representation from bracket notation
# use e.g. as  RNAfold < foo.seq | b2mt | xmgr -pipe
# definition: h=number of base pairs enclosing base
while (<>) {
    print if (s/>/#/);
	
    next unless (/\.\.\./);
    next if (/\[/);   # don't process output from partition function
    chop;
    @F=split(//,$_);
    $p=0; $h=0;
    foreach $i (@F) {
	$h-- if ($i eq ')');
	$p++;
	printf("%4d %4d\n",$p,$h);
	$h++ if ($i eq '(');	# increase $h *after* printing
    }
    print "&\n";
}
    
