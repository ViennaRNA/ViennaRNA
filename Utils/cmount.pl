#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2002-11-27 16:39:03 ivo>
# Produce coloured Hogeweg mountain representation in PostScript.
# Input is a colour _dp.ps file from alidot (aka read_ali) or dp_zoom
# definition: mm[i],mp[i]=number of base pairs enclosing base i

sub usage {
    die "Usage: $0 alidot.ps\n";
}

usage() if eof();

print "%!PS-Adobe-2.0 EPSF-1.2
%%Title: Mount Ali
%%Creator: in hiding
%%BoundingBox:66 209 518 686
%%%%Pages: 1
%%EndComments\n";


print <<EOF; # PS macros
/trapez {   % use as i j height prob hue sat
  dup 0.3 mul 1 exch sub sethsbcolor
  newpath
  3 index 0.5 sub 2 index moveto % i-0.5 h moveto
  dup 1 exch rlineto
  4 2 roll exch sub 2 sub 0 rlineto
  neg 1 exch rlineto
  pop closepath
  gsave fill grestore
  0 setgray stroke
} def

/centershow {
  gsave 1 xs div 1 ys div scale 60 rotate
% dup stringwidth pop 2 div neg 0 rmoveto
  show
  grestore
} def
EOF

$from=1;
while (<>) {
    chomp;
    if (/^% Subsequence from (\d+) to (\d+)/) { # get start and end 
	$from=$1; $to=$2;
    }
    if (/\/sequence \{ \((\S*)[\\\)]/) {
	print "$_\n";
        $seq = $1;              # empty for new version
        while (!/\) \} def/) {  # read until end of definition
            $_ = <>;
	    print $_;
            /(\S*)[\\\)]/;      # ends either with `)' or `\'
            $seq .= $1;
        }
	print "/len { sequence length } def\n";
        $length = length($seq);
        next;
    }

    ($i, $j, $p, $h, $s, $tok) = split;

    if ($tok eq "lbox") { # only read lbox entries
	$mp[$i+1]+=$p*$p; 
	$mp[$j]  -=$p*$p;
	$pair[$i] = $j;
	$hue[$i]  = $h;
	$sat[$i]  = $s;
	$pr[$i]   = $p*$p;
    }
}
for ($i=1; $i<=$length; $i++) { #find maximum for scaling
    $mp[$i]+=$mp[$i-1];
    $max = $mp[$i] if ($mp[$i]>$max);
    $mm[$i]+=$mm[$i-1];
    $max = $mp[$i] if ($mp[$i]>$max);    
}

# postscript scaleing etc
print "72 216 translate\n";
print "/xs {72 6 mul len div} def /ys {72 6 mul $max div} def xs ys scale\n";
print "0.03 setlinewidth
/Times-Roman findfont 1.8 scalefont setfont
0 1 len 1 sub {
    dup
    0.8 add -0.5 moveto
    sequence exch 1 getinterval
    gsave 1 $max len div scale show grestore
} for\n\n";

print "/Times-Roman findfont 10 scalefont setfont
0.01 setlinewidth
len log 0.7 sub cvi 10 exch exp  % grid spacing
gsave 0.5 0 translate
/temp 12 string def
0 exch len {
   dup dup
   0 moveto
   $max 1.03 mul lineto
   cvi $from 1 sub add temp cvs centershow
} for
stroke
grestore\n";
   

for ($i=1; $i<=$length; $i++) { # print pairs as coloured trapezes
    next unless ($pair[$i]);
    print "$i $pair[$i] ";
    printf "%6.4f %6.4f $hue[$i] $sat[$i] trapez\n", $mp[$i], $pr[$i] ;
}

print "showpage\n";

