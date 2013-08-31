#!/usr/bin/perl -w
# -*-Perl-*-
# Last changed Time-stamp: <2006-08-15 22:04:27 ivo>
# after predicting a conensus structure using RNAalifold,
# refold the indivual sequences, using the consensus structure as constraint

# read a clustal fromat alignment, plus either the standart output of
# RNAalifold or a dot plot produced by "RNAalifold -p"

# usage refold.pl [-t thresh] clustal.aln alidot.ps
# or    refold.pl clustal.aln clustal.alifold

use Getopt::Long;
#use POSIX;
#use GD;
use strict;

my $thresh = 0.9;
my $help = 0;
my $TURN = 3;

GetOptions ('-t=f' => \$thresh,
            'turn' => \$TURN,
            "help"=>\$help,
            "h"=>\$help,
            );

if ($help){
  print "\  refold.pl [-t threshold] myseqs.aln alidot.ps | RNAfold -C\n";
  print "or refold.pl [--turn looplength] myseqs.aln myseqs.alifold | RNAfold -C\n";
  print " -t ... use only pairs with p>thresh (default: 0.9)\n";
  print " --turn ... use only pairs closing hairpins with length > looplength (default: 3)\n\n";
  exit(1);
}

my $aln = readClustal();
my $len = length($aln->[0]->{seq});

$_ = <>;
#print STDERR "$_\n";
my $constraint = (/^%\!PS/) ? read_dot() : read_alifold();

#print STDERR "$constraint\n";

foreach my $a (@$aln) {
    print "> $a->{name}\n";

    # remove gaps
    my $seq = $a->{seq};
    my $cons = $constraint;
    my @pt = make_pair_table($cons);
    my %pair = ("AU" => 5,
                "GC" => 1, 
                "CG" => 2, 
                "UA" => 6,
                "GU" => 3, 
                "GT" => 3, 
                "TG" => 4, 
                "UG" => 4,
                "AT" => 5, 
                "TA" => 6);


    for my $p (0..length($seq)-1) {
        # remove non-compatible pairs as well as pairs to a gap position
        my $c = substr($seq,$p,1);
        if ($c eq '-') {
            substr($cons,$p,1) = 'x'; # mark for removal
            substr($cons,$pt[$p],1) = '.' 
                if $pt[$p]>0  && substr($cons,$pt[$p],1) ne 'x';   # open pair
        }
        elsif ($pt[$p]>$p) {
            substr($cons,$p,1) = substr($cons,$pt[$p],1) = '.'
                unless exists $pair{$c . substr($seq,$pt[$p],1)};
        }
    }
#    print STDERR length($seq), length($cons), "\n";
    $cons =~ s/x//g;
    $seq  =~ s/-//g;

    # remove all hairpins with length < TURN
    @pt = make_pair_table($cons);
    for my $p (0..length($seq)-1) {
      next if $p > $pt[$p];
      substr($cons,$p,1) = substr($cons,$pt[$p],1) = '.'
        if $pt[$p] - $p - 1 < $TURN;
    }
#    print STDERR length($seq), length($cons), "\n";
    print "$seq\n$cons\n";
}


sub make_pair_table {
   #indices start at 0 in this version!
   my $structure = shift;
   my (@olist, @table);
   my ($hx,$i) = (0,0);

   foreach my $c (split(//,$structure)) {
       if ($c eq '.') {
          $table[$i]= -1;
     } elsif ($c eq '(') {
          $olist[$hx++]=$i;
     } elsif ($c eq ')') {
         my $j = $olist[--$hx];
         die ("unbalanced brackets in make_pair_table") if ($hx<0);
         $table[$i]=$j;
         $table[$j]=$i;
     }
      $i++;
   }
   die ("too few closed brackets in make_pair_table") if ($hx!=0);
   return @table;
}

sub read_alifold {
    while (<>) {
        return $1 if /^([(.)]+)/;
    }
}

sub read_dot {
    my $cons = '.' x $len;
    while(<>) {
        next if /^%/;
        next unless /lbox$/;
        my @F = split;
        next if $F[5] < $thresh;
        substr($cons,$F[3]-1,1) = '(';
        substr($cons,$F[4]-1,1) = ')';
    }
    return $cons;
}

######################################################################
#
# readClustal(filehandle)
#
# Reads Clustal W formatted alignment file and returns it in list of
# hash references with keys "name" and "seq" for the name and the sequence,
# resp.
# 
# Fixme: the code below assumes all uppercase sequences 
######################################################################

sub readClustal{
#  my $fh=shift;
  my @out=();
  my (%order, $order, %alignments);
  while(<>) {
        last if eof;
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


=head1 NAME

refold.pl - refold using consensus structure as constraint

=head1 SYNOPSIS

  refold.pl [-t thresh] file.aln alidot.ps
  refold.pl [--turn length] file.aln file.alifold

=head1 DESCRIPTION

refold.pl reads an alignment in CLUSTAL format, and a consensus
secondary structure, either in the form of the standard output of
C<RNAalifold>, or in the form of a dot plot as produced by
C<RNAalifold -p>. For each sequence in the alignment it writes the
name, sequence, and constraint structure to stdout in a format
suitable for piping into C<RNAfold -C>.

The constraint string is produced by removing from the consensus
structure all gaps, as well as all pairs not compatible with the
particular sequence. If the structure is read from a dot plot file,
only pairs with probability > then some threshold (default 0.9) are
used.

=head1 OPTIONS

=over 4

=item B<-t> I<thresh>

use only pairs with p>thresh as constraint. Only applicable when
reading consensus structure from a dot plot.

=item B<--turn> I<length>

use only pairs closing a hairpin with looplength >length. This option
allows for setting minimal loop length of hairpin loops (default: 3)

=back

=head1 EXAMPLE

  RNAalifold test.aln > test.alifold
  refold.pl test.aln test.alifold > test.cfold
  RNAfold -C < test.cfold

=cut
