#!/usr/bin/perl -w

use Getopt::Long;
use POSIX;
#use GD;
use strict;

my $ssFile='alirna.ps';
my $columnWidth=120;
my $formatT=0;
my $help=0;

GetOptions ('-s:s' => \$ssFile,
			'-c:i' => \$columnWidth,
			'-t' => \$formatT,
			"help"=>\$help,
			"h"=>\$help,
		   );

if ($help){
  print "\ncoloraln.pl [-s structure.ps] < file.aln \n\n";
  print " -s ... RNAalifold consensus structure file (default: alirna.ps) \n";
  print " -c ... maximum number of columns in a block (default: 120)\n";
  print " -t ... print Ts (default: print Us)\n\n";
  exit(1);
}

if (!-e $ssFile){
  print "No RNAalifold consensus secondary structure file found (use option -s with the correct filename).\n";
  exit(1);
}

my $aln=readClustal();

my @ss=();
for my $i (1..length($aln->[0]->{seq})){
  push @ss,'.';
}

my $consStruc='';

open(ALIRNA,"<$ssFile");
my $pairsFlag=0;
my $sequenceFlag=0;
while (<ALIRNA>){
  $pairsFlag=1 if (/\/pairs/);
  if ($pairsFlag and /\[(\d+) (\d+)\]/){
	$ss[$1-1]='(';
	$ss[$2-1]=')';
  }
  $pairsFlag=0 if ($pairsFlag and /def/);
}

$consStruc=join('',@ss);

for my $row (@$aln){
  $row->{seq}=uc($row->{seq});
  if ($formatT){
	$row->{seq}=~s/U/T/g;
  } else {
	$row->{seq}=~s/T/U/g;
  }
}


my $consSeq=consensusSeq($aln);

#print "$consSeq\n$consStruc\n";

plotAln($aln,$consSeq,$consStruc);



######################################################################
#
# consensusSeq(\@aln alnref)
#
# Returns consensus sequence of alignment
#
######################################################################

sub consensusSeq{

  my @aln=@{$_[0]};
  my $out='';

  for my $i (0..length($aln[0]->{seq})-1){

	my %countHash=('A'=>0,'C'=>0,'G'=>0,'T'=>0,'U'=>0,'-'=>0);
	#print %countHash,"\n";
	for my $j (0..$#aln){
	  my $c=substr($aln[$j]->{seq},$i,1);
	  $countHash{$c}++;
	}
	#print %countHash,"\n";
	my $maxCount=0;
	my $maxChar='';

	for my $c ('A','C','G','T','U','-'){
	  if ($countHash{$c}>=$maxCount){
		$maxChar=$c;
		$maxCount=$countHash{$c};
	  }
	}
	$out.=$maxChar;
  }
  return $out;
}


######################################################################
#
# readClustal(filehandle)
#
# Reads Clustal W formatted alignment file and returns it in list of
# hash references with keys "name" and "seq" for the name and the sequence,
# resp.
#
######################################################################

sub readClustal{
#  my $fh=shift;
  my @out=();
  my (%order, $order, %alignments);
  while(<>) {
	next if ( /^\s+$/ );
	my ($seqname, $aln_line) = ('', '');
	if( /^\s*(\S+)\s*\/\s*(\d+)-(\d+)\s+(\S+)\s*$/ ) {
	  # clustal 1.4 format
	  ($seqname,$aln_line) = ("$1/$2-$3",$4);
	} elsif( /^(\S+)\s+([A-Z\-]+)\s*$/ ) {
	  ($seqname,$aln_line) = ($1,$2);
	} else {
	  next;
  }
	if( !exists $order{$seqname} ) {
	  $order{$seqname} = $order++;
	}
	$alignments{$seqname} .= $aln_line;
  }

  foreach my $name ( sort { $order{$a} <=> $order{$b} } keys %alignments ) {
	if( $name =~ /(\S+):(\d+)-(\d+)/ ) {
	  (my $sname,my $start, my $end) = ($1,$2,$3);
	} else {
	  (my $sname, my $start) = ($name,1);
	  my $str  = $alignments{$name};
	  $str =~ s/[^A-Za-z]//g;
	  my $end = length($str);
	}
	my $seq=$alignments{$name};
	push @out, {name=>$name,seq=>$seq};
  }
  return [@out];
}

######################################################################
#
# getPairs(\@aln alnref, $ss string)
#
# Evalutates the pairing of an alignment according to a given
# consensus secondary structure
#
# Returns list of all base pairs which is a hash with the following
# keys:
#
#  open ... column in the alignment of the first base in the pair "("
#  close ... column in the alignment of the second base in the pair ")"
#  all ... list of all basepairs in the different sequences in the alignment
#  pairing ... list of all different pairing basepairs
#  nonpairing ... list of all incompatible basepairs
#
######################################################################

sub getPairs{

  my @inputAln=@{$_[0]};
  my $ss=$_[1];

  # return nothing if there are no pairs
  if (!($ss=~tr/(/(/)){
	return ();
  }

  my @aln=();
  foreach my $row (@inputAln){
	my $seq=$row->{seq};
	$seq=uc($seq);
	$seq=~s/T/U/g;
	my @tmp=split(//,$seq);
	push @aln,\@tmp;
  }
  my @ss=split(//,$ss);

  my @pairs=();
  my @stack=();

  foreach my $column (0..$#ss){

	my $currChar=$ss[$column];

	if ($currChar eq '('){
	  push @stack,$column;
	}

	if ($currChar eq ')'){
	  my $openedCol=pop @stack;
	  push @pairs,{open=>$openedCol,close=>$column};
	}
  }

  @pairs=sort {$a->{open} <=> $b->{open}} @pairs;

  foreach my $i (0..$#pairs){
	#print "$i: $pairs[$i]->{open} - $pairs[$i]->{close}\n";

	my @all=();
	my @pairing=();
	my @nonpairing=();

	for my $j (0..$#aln){
	  my $currPair=$aln[$j][$pairs[$i]->{open}].$aln[$j][$pairs[$i]->{close}];
	  push @all,$currPair;
	}

	for my $pair (@all){
	  if (($pair eq 'AU') or
		  ($pair eq 'UA') or
		  ($pair eq 'GC') or
		  ($pair eq 'CG') or
		  ($pair eq 'UG') or
		  ($pair eq 'GU')){

		push @pairing, $pair;
	  } elsif ($pair eq "--"){
		# do nothing
	  } else {
		push @nonpairing,$pair;
	  }
	}

	undef my %saw;
    my @uniquePairing = grep(!$saw{$_}++, @pairing);

	$pairs[$i]->{all}=[@all];
	$pairs[$i]->{pairing}=[@uniquePairing];
	$pairs[$i]->{nonpairing}=[@nonpairing];
  }

  return @pairs;

}


######################################################################
#
# plotAln(\@aln ref-to-alignment, $consensus string, $ss string)
#
# Creates a image of the alignment with various annotations
#
# \@aln ... alignment in list of hash format
# $consensus ... consensus sequence
# $ss ... consensus secondary structure in dot bracket notation
#
# Returns png as string.
#
######################################################################

sub plotAln{

  my @aln=@{$_[0]};
  my $consensus=$_[1];
  my $ss=$_[2];

  my $ps='';

  # Get important measures
  (my $fontWidth, my $fontHeight) = (6,6);

  my $length=length($aln[0]->{seq});
  my $maxName=0;
  foreach my $row (@aln){
	$maxName=length($row->{name}) if (length($row->{name})>$maxName);
  }

  # Custom Sizes
  my $lineStep=$fontHeight+2; # distance between lines
  my $blockStep=3.5*$fontHeight; # distance between blocks
  my $consStep=$fontHeight*0.5; # distance between alignment and conservation curve
  my $ssStep=2;
  my $nameStep=3*$fontWidth; # distance between names and sequences
  my $maxConsBar=2.5*$fontHeight; # Height of conservation curve
  my $startY=2; # "y origin"
  my $namesX=$fontWidth; # "x origin"
  my $seqsX=$namesX+$maxName*$fontWidth+$nameStep; # y-coord. where sequences start

  # Calculate width and height

  my $tmpColumns=$columnWidth;

  $tmpColumns=$length if ($length<$columnWidth);

  my $imageWidth=$namesX+($maxName+$tmpColumns)*$fontWidth+$nameStep+$fontWidth;
  my $imageHeight=$startY+(ceil($length/$columnWidth))*((@aln+1)*$lineStep+$blockStep+$consStep+$ssStep);



  my $white="1 1 1";
  my $black="0 0 0" ;
  my $grey="1 1 1";

  my $red="0.0 1";
  my $ocre="0.16 1";
  my $green="0.32 1";
  my $turq="0.48 1";
  my $blue="0.65 1";
  my $violet="0.81 1";


  my $red1="0.0 0.6";
  my $ocre1="0.16 0.6";
  my $green1="0.32 0.6";
  my $turq1="0.48 0.6";
  my $blue1="0.65 0.6";
  my $violet1="0.81 0.6";

  my $red2="0.0 0.2";
  my $ocre2="0.16 0.2";
  my $green2="0.32 0.2";
  my $turq2="0.48 0.2";
  my $blue2="0.65 0.2";
  my $violet2="0.81 0.2";


  my @colorMatrix=([$red,$red1,$red2],
				   [$ocre,$ocre1,$ocre2],
				   [$green,$green1,$green2],
				   [$turq,$turq1,$turq2],
				   [$blue,$blue1,$blue2],
				   [$violet,$violet1,$violet2]);

  my @pairs=getPairs(\@aln,$ss);

  foreach my $pair (@pairs){

	foreach my $column ($pair->{open},$pair->{close}){
	  my $block=ceil(($column+1)/$columnWidth);
	  my $effCol=$column-($block-1)*$columnWidth;
	  my $x=$seqsX+$effCol*$fontWidth;

	  foreach my $row (0..$#aln){

		my $pairing=@{$pair->{pairing}};
		my $nonpairing=@{$pair->{nonpairing}};

		my $color;

		if ($nonpairing <=2){
		  $color=$colorMatrix[$pairing-1][$nonpairing];

		  if (($pair->{all}->[$row] eq 'AU') or
			  ($pair->{all}->[$row] eq 'UA') or
			  ($pair->{all}->[$row] eq 'GC') or
			  ($pair->{all}->[$row] eq 'CG') or
			  ($pair->{all}->[$row] eq 'GU') or
			  ($pair->{all}->[$row] eq 'UG')){

			my $y=$startY+($block-1)*($lineStep*(@aln+1)+$blockStep+$consStep)+$ssStep*($block)+($row+1)*$lineStep;

			my $xtmp=$x+$fontWidth;
			my $ytmp1=$y-1;
			my $ytmp2=$y+$fontHeight+1;
			$ps.="$x $ytmp1 $xtmp $ytmp2 $color box\n";
		  }
		}
	  }
	}
 }
  # Calculate conservation scores as (M+1)/(N+1), where M is the
  # number of matches to the consensus and N is the number of
  # sequences in the alignment

  my @conservation=(); # each column one score

  for my $column (0..$length-1){

	my $consChar=substr($consensus,$column,1);

	# if consensus is gap, score=0
	if ($consChar eq "-" or $consChar eq "_"){
	  push @conservation, 0;
	  next;
	}

	my $match=0;
	for my $row (0..$#aln){
	  my $currChar=substr($aln[$row]->{seq},$column,1);
	  $match++ if ($currChar eq $consChar);
	}

	my $score=($match-1)/(@aln-1);
	push @conservation,$score;
  }

  # Draw the alignments in chunks
  my $currY=$startY;
  my $currPos=0;

  while ($currPos<$length){

	# secondary structure in first line
	#$out->string($font,$seqsX,$currY,substr($ss,$currPos,$columnWidth),$black);
	my $tmpSeq=substr($ss,$currPos,$columnWidth);
	$tmpSeq=~s/\(/\\\(/g;
	$tmpSeq=~s/\)/\\\)/g;
	$ps.="($tmpSeq) $seqsX $currY string\n";
	$currY+=$lineStep+$ssStep;

	# sequences labeled only with the organism-specifier
	foreach my $row (@aln){
	  #$out->string($font,$namesX,$currY,$row->{name},$black);
	  $ps.="($row->{name}) $namesX $currY string\n";
	  #$out->string($font,$seqsX,$currY,substr($row->{seq},$currPos,$columnWidth),$black);
	  my $tmpSeq=substr($row->{seq},$currPos,$columnWidth);
	  $ps.="($tmpSeq) $seqsX $currY string\n";
	  $currY+=$lineStep;
	}

	# conservation curve
	$currY+=$consStep;
	$ps .= "0.6 setgray\n";
	for my $col ($currPos..$currPos+$columnWidth-1){
	  my $score=shift @conservation;
	  last if (!defined $score);
	  my $barHeight=$maxConsBar*$score;
	  $barHeight=1 if ($barHeight==0);
	  my $x=$seqsX+($col-($columnWidth*int($currPos/$columnWidth)))*$fontWidth;
	 # $out->filledRectangle($x,$currY+$maxConsBar-$barHeight,$x+$fontWidth,$currY+$maxConsBar,$grey);
	  my $ytmp1=$currY+$maxConsBar-$barHeight;
	  my $xtmp2=$x+$fontWidth;
	  my $ytmp2=$currY+$maxConsBar;
	  $ps .= "$x $ytmp1 $xtmp2 $ytmp2 box2\n";
	}
	$currY+=$blockStep;
	$currPos+=$columnWidth;
  }

  my $BB="0 0 $imageWidth $imageHeight";

  while (<DATA>){
	s/!BB!/$BB/;
	s/!HEIGHT!/$imageHeight/;
	print;
  }
  print $ps;
  print "showpage\n";
}

=head1 NAME

coloraln.pl - colorize an alignment with consensus structure

=head1 SYNOPSIS

  coloraln.pl [-s structure.ps] file.aln

=head1 DESCRIPTION

colorrna.pl reads an alignment in CLUSTAL format, and a consensus
secondary structure (which it extracts from a postscript secondary
structure plot). It produces a postscript figure of the alignment, in
which compensatory mutation supporting the consensus structure are
marked by color. The color scheme is the same employed by RNAalifold
and alidot: Red marks pairs with no sequence variation; ochre, green,
turquoise, blue, and violet mark pairs with 2,3,4,5,6 different tpyes
of pairs, respectively.

=head1 OPTIONS

=over 4

=item B<-s> I<file>

read I<file> to extract the consensus structure (default: C<alirna.ps>)

=item B<-c> I<width>

break alignments into blocks of at most I<width> columns, (default: 120)

=item B<-t>

suppress conversion of C<T> to C<U>, i.e. do not convert DNA to
RNA, (default: convert to C<U>)

=back

=cut


__DATA__
%!PS-Adobe-3.0 EPSF-3.0
%%BoundingBox: !BB!
%%EndComments

% draws Vienna RNA like colored boxes
/box { % x1 y1 x2 y2 hue saturation
  gsave
  dup 0.3 mul 1 exch sub sethsbcolor
  exch 3 index sub exch 2 index sub rectfill
  grestore
} def

% draws a box in current color
/box2 { % x1 y1 x2 y2
  exch 3 index sub exch 2 index sub rectfill
} def

/string { % (Text) x y
  6 add
  moveto
  show
} def

0 !HEIGHT! translate
1 -1 scale
/Courier findfont
[10 0 0 -10 0 0] makefont setfont

%100 100 106 110 0.5 0.5 1 box
%(A) 100 100 string
%(This is a string) 200 200 string
