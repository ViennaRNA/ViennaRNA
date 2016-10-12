package RNA::Utils;

use strict;
use Exporter;
use warnings;

use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION     = 1.00;
@ISA         = qw(Exporter);
@EXPORT      = ();
@EXPORT_OK   = qw(make_pair_table DB2PairList DB2Shape DB2Helix MergeHelices);
%EXPORT_TAGS = ( Structures => [qw(&make_pair_table &DB2PairList &DB2Shape &DB2Helix &MergeHelices)] );

sub make_loop_index_from_pt {
  # number loops and assign each position its loop-number
  # handy for checking which pairs can be added
  my @pt = @_;
  shift @pt; # get rid of length in first entry
  
  my (@loop, @olist);
  my ($l, $nl, $hx) = (0, 0, 0);

  foreach my $i (0 .. $#pt) {
    if ($i < $pt[$i]) { # '('
      $l = ++$nl;       # start a new loop
      $olist[$hx++]=$i;
    }
    $loop[$i] = $l;
    if ($pt[$i] && $pt[$i] < $i) { # ')'
      --$hx;            # i pairs with olist[--hx]
      if ($hx > 0) {
        my $idx = $olist[$hx-1];
        $l = $loop[$idx];
      } else {
        $l = 0;
      }
    }
  }

  unshift @loop, $nl;
  return @loop;
}

sub make_loop_index {
  use integer;
  # number loops and assign each position its loop-number
  # handy for checking which pairs can be added
  my $struc = shift;
  my($c, $j, @olist, @loop);
  my ($hx,$i,$l, $nl)=(0,0,0,0);
  foreach $c (split(//,$struc)) {
    if ($c eq '(') {
      $nl++; $l=$nl;       # start a new loop
      $olist[$hx++]=$i;
    }
    $loop[$i]=$l;
    if ($c eq ')') {
      --$hx;                         # i pairs with olist[--hx]
      if ($hx>0) {                   # we're still in a loop
        my $jj = $olist[$hx-1];    # jj started the previous loop 
        $l = $loop[$jj] if (defined($jj)); # previous loop index
      } else {
        $l = 0;                    # external loop has index 0
      }
    }
    $i++;
  }
  push @loop, $nl;
  return @loop;
}


sub make_pair_table {
   use integer;
   # let table[i]=j if (i,j) is pair, 0 if i unpaired
   # indices start at 1 in this version!
   my($str) = shift(@_);
   my($i,$j,$hx,$c,@olist,@table);
   $hx=0;
   $i=1;
   $table[0] = length($str);

   foreach $c (split(//,$str)) {
      if ($c eq '(') {
         $olist[$hx++]=$i;
      } elsif ($c eq ')') {
         $j = $olist[--$hx];
         die ("unbalanced brackets in make_pair_table") if ($hx<0);
         $table[$i]=$j;
         $table[$j]=$i;
      } else { # ($c eq '.') {
         $table[$i]=0;
      } 
      $i++;
   }
   die ("too few closed brackets in make_pair_table") if ($hx!=0);
   return \@table;
}

sub DB2PairList {
  my $struct  = shift;
  my $pt      = make_pair_table($struct);
  my %plist   = ();

  for(my $i = 1; $i <= $pt->[0]; $i++){
    next unless $i < $pt->[$i];
    $plist{$i,$pt->[$i]} = 1;
  }

  return \%plist;
}

sub DB2Shape{
  my $pt = shift;
  my $i  = shift;
  my $j  = shift;
  my $level = shift;
  my $components;
  my $ext;

  my $shape_string = "";

  # find out if there is more than one component
  for(my $p = $i, $components = 0; $p <= $j; $p++){
    if($p < $pt->[$p]){
      $components++;
      $p = $pt->[$p];
      last if $components > 1;
    }
  }
  # find out if we process the exterior loop
  for(my $p = $i - 1, $ext = 1; $p > 0; $p--){
    if($pt->[$p] > $j){ $ext = 0; last;}
  }

  for(; $i <= $j; $i++){
    if($i < $pt->[$i]){ # we have base pair (i,j)
      my $k = $i + 1;
      my $l = $pt->[$i] - 1;
      my $sub_string = "";
      my $sub5 = "";
      my $sub3 = "";
      while(1){
        my $u5 = 0;
        my $u3 = 0;
        for(; $pt->[$k] == 0; $k++){ $u5++;}; # skip unpaired 5' side
        for(; $pt->[$l] == 0; $l--){ $u3++;}; # skip unpaired 3' side
        if($k >= $l){ # hairpin loop
          if($level == 0){
            $sub_string .= "_" x $u5;
          }
          last;
        } else{
          if($pt->[$k] != $l){ # multi branch loop
            $sub_string = DB2Shape($pt, $k, $l, $level);
            if($level == 1){
              $sub5 .= "_" if $u5 > 0;
              $sub3 = "_".$sub3 if $u3 > 0;
            } elsif($level == 0){
              $sub5 .= "_" x $u5;
              $sub3 = ("_" x $u3).$sub3;
            }
            last;
          } else{ # interior loop with enclosed pair (k,l)
            if($level < 5){
              my $no_stack = ($u5 + $u3 > 0);
              my $bulge = ($u5 > 0) && ($u3 > 0);
              if($level == 4){
                if($bulge){
                  $sub5 .= "[";
                  $sub3 = "]".$sub3;
                }
              } elsif($level == 3){
                if($no_stack){
                  $sub5 .= "[";
                  $sub3 = "]".$sub3;
                }
              } elsif(($level == 2) || ($level == 1)){
                if($no_stack){
                  $sub5 .= "_" if $u5 > 0;
                  $sub3 = "_".$sub3 if $u3 > 0;
                  $sub5 .= "[";
                  $sub3 = "]".$sub3;
                }
              } elsif($level == 0){
                  $sub5 .= "_" x $u5;
                  $sub3 = ("_" x $u3).$sub3;
                  $sub5 .= "[";
                  $sub3 = "]".$sub3;
              }
            }
            $k++;$l--;
          }
        }
      }
      $shape_string .= "[".$sub5.$sub_string.$sub3."]";
      $i = $pt->[$i];
    } elsif($pt->[$i] == 0){
      if($level < 3){
        if($level == 0){
          $shape_string .= "_";
        } else {
          if(($level == 1) || (($components < 2) && ($ext == 0))){
            $shape_string .= "_";
            for(; ($i <= $j) && ($pt->[$i] == 0); $i++){;};
            $i--;
          }
        }
      }
    }
  }
  return $shape_string;
}

sub DB2Helix{
  my ($pt,$i,$end, $helices)        = @_;
  my ($h_start, $h_length, $h_end)  = (0, 0, 0);

  for(; $i < $end; $i++){
    next if $i > $pt->[$i];
    $h_start  = $i;
    $h_end    = $pt->[$i];
    $h_length = 1;
    $h_length++, $i++ while($pt->[$i+1] == $pt->[$i]-1);

    DB2Helix($pt, $i+1, $h_end, $helices) if $i < $h_end;

    if($h_length >= 1){
      my %entry;
      $entry{start}   = $h_start;
      $entry{end}     = $h_end;
      $entry{length}  = $h_length;
      push(@{$helices}, \%entry);
    }

    $i = $pt->[$h_start] - 1;
  }
}

#
# this function will fail if the input helix list is not sorted
# ascending by helix start position
#
sub MergeHelices{
  my ($helices, $max_dist) = @_;
  my $merged = 0;

  do{
    $merged = 0;
    for(my $i=1; $i < @{$helices}; $i++){
      # GOAL:
      # merge two consecutive helices i and i-1, if i-1
      # subsumes i, and not more than i

      # check whether the next helix may be a neighboring stem within the previous helix
      my $has_neighbors = 0;
      for(my $j = $i+1; $j < @{$helices}; $j++){
        last if !defined $helices->[$j];
        last if $helices->[$j]->{start} > $helices->[$i-1]->{end};
        next if $helices->[$j]->{start} < $helices->[$i]->{end};
        $has_neighbors = 1;
      }
      next if $has_neighbors;

      # check whether we can merge i with i-1
      if($helices->[$i]->{end} < $helices->[$i-1]->{end}){
      
        # compute and memorize unpaired nucleotides for 5' part
        if(!exists($helices->[$i-1]->{up5})){
          $helices->[$i-1]->{up5} = $helices->[$i]->{start} - $helices->[$i-1]->{start} - $helices->[$i-1]->{length};
        } else {
          $helices->[$i-1]->{up5} += $helices->[$i]->{start} - $helices->[$i-1]->{start} - $helices->[$i-1]->{length} - $helices->[$i-1]->{up5};
        }

        # compute and memorize unpaired nucleotides for 3' part
        if(!exists($helices->[$i-1]->{up3})){
          $helices->[$i-1]->{up3} = $helices->[$i-1]->{end} - $helices->[$i-1]->{length} - $helices->[$i]->{end};
        } else {
          $helices->[$i-1]->{up3} += $helices->[$i-1]->{end} - $helices->[$i-1]->{length} - $helices->[$i-1]->{up3} - $helices->[$i]->{end};
        }

        # add unpiared nucleotides gathered for the current helix
        $helices->[$i-1]->{up5} += $helices->[$i]->{up5} if exists($helices->[$i]->{up5});
        $helices->[$i-1]->{up3} += $helices->[$i]->{up3} if exists($helices->[$i]->{up3});

        # add number of bp from i to that of previous helix i-1
        $helices->[$i-1]->{length} += $helices->[$i]->{length};
        # remove helix i
        splice @{$helices}, $i, 1;
        $merged = 1;
        last;
      }
    }
  } while($merged);
}

1;
