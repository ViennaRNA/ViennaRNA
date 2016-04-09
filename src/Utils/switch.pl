#!/usr/bin/perl -w
# -*-Perl-*-
# tool for the design of bistable RNA molecules

# additional paths were perl looks for RNA.pm and RNA/barriers.pm
# use lib '/home/RNA/barriers';

use RNA;
use Getopt::Long;
use strict;
use vars qw/$opt_debug $opt_v $ParamFile $opt_circ/;
my $sE;
my %pair = ("AU", 5, "GC", 1, "CG", 2, "GU", 3, "UG", 4, "UA", 6);
my %mut = ('A' => 'G', 'G' => 'A', 'U' => 'C', 'C' => 'U');
my %wc = ('A' => 'U', 'G' => 'C', 'U' => 'A', 'C' => 'G');

my ($fist, $sest, $startseq, $cons);
my $noo = 1;
my $bar;
my $nom = 2000000;
my $border = 30000;
my $small = 0.3;
my $Temperature1 = 37;
my $Temperature2;
my ($optseq, $optcost, $e);
my @fibo = (0,1);
my %cos;
my $cpnt= -1; # Cofold
my $dG = 0;   # deltaG for normal Switches

Getopt::Long::config('no_ignore_case');
&usage() unless GetOptions("T=f" => \$Temperature1,
                           "T2=f" => \$Temperature2,
                           "4" => sub {$RNA::tetra_loop = 0},
                           "d|d0" => sub {$RNA::dangles=0},
                           "d2" => sub {$RNA::dangles=2},,
                           "debug",
                           "noGU" => \$RNA::noGU,
                           "noCloseGU" => \$RNA::no_closingGU,
                           "noLP" => \$RNA::noLonelyPairs,
                           "P=s" => \$ParamFile,
                           "n=i" => \$noo,
                           "trial=i" => \$nom,
                           "bar=f" => \$bar,
                           "g=f" => \$small,
                           "dG=f" => \$dG,
                           "circ", # circular RNA switches
                           "v");

$RNA::temperature = $Temperature1;
$Temperature2 = $Temperature1 unless defined $Temperature2;

$noo = ($noo >= 1) ? $noo : 1;
$nom = ($nom >= 1) ? $nom : 1000;
$border = ($border >=1) ? $border : 100;
my $interactiv = init_ia(0);
my $ia = 1 if -t STDIN && -t STDOUT && $#ARGV < 0;
RNA::read_parameter_file($ParamFile) if ($ParamFile);
RNA::init_rand();
srand();

for(;;) { # main loop
   $interactiv->() if $ia;
   last if !process_input();
   if ($bar && $opt_circ) {
     print STDERR "Cannot optimize barriers for circular species\n"; $bar=0;
   }

   my @fist = make_pair_table ($fist);
   my @sest = make_pair_table ($sest);

   fibo(length($fist)+2);

   my @clist = make_cycle_table (\@fist, \@sest);
   $border = estimate_border(\@clist);
   print_lol->(\@clist) if $opt_debug;
   my $nocs = calc_number_of_seq(\@clist);
   print "$nocs possible sequences\n" if ($opt_v);
   # make optimization
   for (1..$noo) {
      $startseq = make_compatible(@clist, ' ' x length($fist)); # unless $startseq;
      $optseq   = do_opti($nom, $startseq, \@clist);
      $optcost  = cost_function($optseq);
      substr($optseq, $cpnt-1, 0, '&') if ($cpnt != -1);
      $cos{$optseq}[1]++;
      $cos{$optseq}[0] = $optcost;
     printf "%s %8.4f\n", $optseq, $optcost;
   }

   print "the resulting sequences are:\n";
   foreach $e (sort {$cos{$b}[0]<=>$cos{$a}[0]} keys %cos) {
      printf "%s %8.4f",  $e, $cos{$e}[0];
      printf " %2d", $cos{$e}[1] if ($cos{$e}[1]>1);
      print "\n";
   }
   #print_statistics();
}

sub estimate_border {
   my @clist = @{$_[0]};
   my $max_cyc_len = 0;
   my $border = 0;
   for (@clist) {
      $max_cyc_len = $#{$_}+1 if  $#{$_}+1 > $max_cyc_len;
      $border += calc_number_of_seq([$_])-1;
   }
   printf "max cyclen is %4d border is %5d\n", $max_cyc_len, $border
       if $opt_debug;
   return $border;
}

sub print_lol {
   my @lol = @{$_[0]};
   for (@lol) {
      $, = ' ';
      print @{$_}, "\n";
      $, = '';
   }
}

sub fibo {
   my $length = shift;
   while ($#fibo < $length) {
      push @fibo, $fibo[-1] + $fibo[-2];
   }
   return $fibo[$length];
}

sub process_input {
    $_ = <>;
    return 0 if !defined $_;
    chomp;
    $fist = $_;
    return 0 if $fist eq '@';
    $cpnt = index($fist, "&")+1;
    chomp($_ = <>);
    $sest = $_;
    die "ERROR: Different Cut-Points set!\n" if ($cpnt != index($sest, "&")+1);
    die "ERROR: Structures have unequal length" if (length($fist) != length($sest));
    if ($cpnt) {
        die "ERROR: Cut-Point and circular?\n" if ($opt_circ);
        $fist =~ s/&//g;
        $sest =~ s/&//g;
        $RNA::cut_point = $cpnt;
    } else { $cpnt = -1 }
    # Sequence Constraints (no force bases supported)
    $_ = <>;
    if ($_) {
        chomp($_);
        $cons = uc $_;
        if (($cpnt != -1) && (length($cons) > $cpnt-1)) {
            die "ERROR: Different Cut-Points set!\n" if ($cpnt != index($cons, "&")+1);
            $cons =~ s/&//g;
        }
        $cons .= 'N' x (length($fist)-length); # pad with N
    } else {
        $cons .= 'N' x length($fist); # pad with N
    }
    #print "$cons\n";
    return 1;
}

sub init_ia {
   my $times = undef;
   return sub {
      if (!defined($times)) {
         print "Input structure1 structure2 & start string\n";
         print "    @ to quit, enter for random start string\n";
      }
      print "T1= $Temperature1,  T2= $Temperature2";
      $bar ? print "  Barrier=$bar\n": print "\n";
      for (1..8) { print '....,....', $_;}
      print "\n";
      $times=1;
   }
}


sub make_pair_table {
   use integer;
   # let table[i]=j if (i,j) is pair, -1 if i unpaired
   # indices start at 0 in this version!
   my($str) = shift(@_);
   my($i,$j,$hx,$c,@olist,@table);
   $hx=$i=0;
   foreach $c (split(//,$str)) {
      if ($c eq '.') {
         $table[$i]= -1;
      } elsif ($c eq '(') {
         $olist[$hx++]=$i;
      } elsif ($c eq ')') {
         $j = $olist[--$hx];
         die ("unbalanced brackets in make_pair_table") if ($hx<0);
         $table[$i]=$j;
         $table[$j]=$i;
      }
      $i++;
   }
   die ("too few closed brackets in make_pair_table") if ($hx!=0);
   return @table;
}

sub make_cycle_table {
   use integer;
   my @fist = @{$_[0]};
   my @sest = @{$_[1]};
   my (@clist,@seen);
   for my $i (0..$#fist) {
      my @cycle1 = ();
      next if ($seen[$i]);
      push @cycle1, $i;
      $seen[$i]=1;
      my $j = $fist[$i];
      while ($j>=0 && !$seen[$j]) {
         push @cycle1, $j;
         $seen[$j] = 1;
         $j = ($#cycle1 % 2) ? $sest[$j] : $fist[$j];
      }
      $j = $sest[$i];
      my @cycle2 = ();
      while ($j>=0 && !$seen[$j]) {
         unshift @cycle2, $j;
         $seen[$j] = 1;
         $j = ($#cycle2 % 2) ? $sest[$j] : $fist[$j];
      }
      # duplicate first element if closed cycle
      unshift @cycle2, $j if ($j!=-1);
      push @clist, [@cycle2, @cycle1];
   }
   return @clist;
}

sub make_compatible {
   # make random sequence or randomly mutate sequence
   # by replacing all @components by random paths
   my $startseq = pop;
   my @components = @_;
   foreach (@components) {
      my @comp = @{$_};
      my $l = @comp;
      if ($comp[$[] eq $comp[-1] && (@comp>1)) { # it's a cycle
         $l--; $l=-$l;
      }
      my @seq = make_pathseq ($l);
      for ($[..$#seq) {
         substr($startseq, $comp[$_], 1, $seq[$_]);
      }
   }
   return $startseq;
}

sub make_pathseq {
   # with one argument: return random cycle of length $l
   # two arguments: return $rand-th possible cycle sequence
   # if $l<0 assume closed cycle else path
   my ($l, $rand) = @_;
   my $ll = $l;
   if ($l<0) { # its a closed cycle
      $l = -$l;
      $ll = $l-1;
      die "cycles must have even length" if $l%2;
   }
   my $n = 2*($fibo[$l+1]+$fibo[$ll]);       # total number of seqs
   $rand = rand($n) unless (defined($rand));

   my @seq = ();
   # set first base in sequence
   if ($rand < $fibo[$ll]) {
      push @seq, 'A', 'U';
   } elsif ($rand < 2*$fibo[$ll]) {
      push @seq, 'C', 'G'; $rand -= $fibo[$ll];
   } else {
      $rand -= 2*$fibo[$ll]; $ll=$l;
      push @seq, ($rand>=$fibo[$l+1])?'U':'G';
      $rand -= $fibo[$l+1] if $rand >= $fibo[$l+1];
   }

   # grow sequence to length $l
   # if we have a cycle starting with A or C $ll=$l-1, else $ll=$l
   while (@seq < $l) {
      if ($rand < $fibo[$ll-@seq]) {
         push @seq, 'C','G' if ($seq[-1] eq 'G');
         push @seq, 'A','U' if ($seq[-1] eq 'U');
      } else {
         $rand -= $fibo[$ll-@seq];
         push @seq, ($seq[-1] eq 'G') ? 'U' : 'G';
      }
   }
   pop @seq if (@seq > $l); # in case we've added one base too many
   return @seq;
}

sub invert_cyclelist {
# foreach position store its cycle in @pos_list
   my @clist = @{$_[0]};
   my @pos_list;
   foreach my $lp (@clist) {
      my @l = @{$lp};
      for my $p (0..$#l) {
         $pos_list[$l[$p]] = [$lp, $p];
      }
   }
   return @pos_list;
}


sub do_opti {
   my $nom = shift;
   my $refseq = shift;
   my $clistp = shift;
   my $refcost = cost_function($refseq);
   my ($mutseq,$newcost);
   my $reject = 0;
   my %seen = ();
   my @poslist = invert_cyclelist($clistp);
   for my $d (1..$nom) {
      $mutseq = mutate_seq($refseq, $clistp, \@poslist);
      if (!exists $seen{$mutseq}) {
         $seen{$mutseq}=1;
         $newcost = cost_function($mutseq,$refcost);
         if ($opt_v) {
            printf '%4d %s %6.3f', $d, $mutseq, $newcost;
            (defined $sE) ? printf("%6.2f\n", $sE) : print "\n";
         }
         if ($newcost < $refcost) {
            $refseq = $mutseq;
            $refcost = $newcost;
            $reject = 0;
            print "New refseq!\n" if ($opt_v);
         }
      }
      else {
         $seen{$mutseq}++;
         $reject++;
         if ($reject >= $border) {
            last;
         }
      }
   }
   return $refseq;
}

sub mutate_seq {
   my $refseq = shift;
   my ($mutseq, $closed);
   my @clist = @{$_[0]};
   my @poslist = @{$_[1]};
   # my @cyc = @{$clist[int rand(@clist)]};
   my $pos = int rand(length $refseq);

   my ($cp, $cpos) = @{$poslist[$pos]};
   my @cyc = @{$cp};
   if ($#cyc>0 && $cyc[0] == $cyc[-1]) { #close cycle
      shift @cyc;
      $cpos--;
      $closed=1;
   }
   if (rand() < 1/@cyc) {
#      my $cp = $clist[int rand(@clist)];
      do { # prevent mutations to the same cyclesequence
         $mutseq = make_compatible($cp, $refseq);
      } while (substr($mutseq,$pos,1) eq substr($refseq,$pos,1));
   } else {
      $mutseq = $refseq;
      my $c = $mut{substr($mutseq, $pos, 1)};
      substr($mutseq, $pos, 1, $c);
      if ($cpos>0 || ($closed)) { # pos>0 or closed cyc
         my $cc = substr($mutseq, $cyc[$cpos-1], 1);
         substr($mutseq, $cyc[$cpos-1], 1, $wc{$c})
             if (!exists $pair{"$c$cc"});
      }
      if ($cpos+1<=$#cyc || $closed) {
         my $cp = ($cpos+1) % @cyc;
         my $cc = substr($mutseq, $cyc[$cp], 1);
         substr($mutseq, $cyc[$cp], 1, $wc{$c})
             if (!exists $pair{"$c$cc"});
      }
   }
   return $mutseq;
}

sub cost_function {
   my $seq = shift;
   my $refcost = shift;
   $RNA::temperature = $Temperature1;
   my ($f1, $e1, $e1s);
   if ($cpnt != -1) {
       $f1 = RNA::co_pf_fold($seq, undef);
       $e1 = RNA::energy_of_struct($seq, $fist);
       $e1s = RNA::energy_of_struct($seq, $sest);
   } elsif ($opt_circ) {
       $f1 = RNA::pf_circ_fold($seq, undef);
       $e1 = RNA::energy_of_circ_struct($seq, $fist);
       $e1s = RNA::energy_of_circ_struct($seq, $sest);
   } else {
       $f1 = RNA::pf_fold($seq, undef);
       $e1 = RNA::energy_of_struct($seq, $fist);
       $e1s = RNA::energy_of_struct($seq, $sest);
   } 
   my $cost;

   if ($Temperature1 != $Temperature2) {
      $RNA::temperature = $Temperature2;
      my ($f2, $e2, $e2s);
      if ($cpnt != -1) {
          $f2 = RNA::co_pf_fold($seq, undef);
          $e2 = RNA::energy_of_struct($seq, $fist);
          $e2s = RNA::energy_of_struct($seq, $sest);
      } elsif ($opt_circ) {
          $f2 = RNA::pf_circ_fold($seq, undef);
          $e2 = RNA::energy_of_circ_struct($seq, $fist);
          $e2s = RNA::energy_of_circ_struct($seq, $sest);
      } else {
          $f2 = RNA::pf_fold($seq, undef);
          $e2 = RNA::energy_of_struct($seq, $fist);
          $e2s = RNA::energy_of_struct($seq, $sest);
      } 
      
      $cost = ($e1-$f1)+($e2-$f2) +
          $small*(($e1-$e1s+$dG) + ($e2 - $e2s+$dG));
  } else {
      $cost = $e1+$e1s-2*$f1+$small*($e1-$e1s+$dG)*($e1-$e1s+$dG);
      if ($bar && ((!defined $refcost) || $cost<$refcost)) {
          $sE = RNA::find_saddle($seq, $fist, $sest, 20)/100.;
          $sE -= ($e1+$e1s)/2;
          printf "sE = %6.2f %6.2f\n", $sE,(0.1*$small*($sE - $bar)*($sE - $bar)) if ($opt_v);
          $cost += (0.1*$small*($sE - $bar)*($sE - $bar));
      }
  }

  my $d = class_dist($seq, $cons);
  $cost += $d*50*$small;

  return $cost;
}

sub calc_number_of_seq {
   my @allcyc = @{$_[0]};
   my $nos = 1;
   foreach (@allcyc) {
      my @cyc = @{$_};
      my $l = @cyc;
      my $nc;
      if ($cyc[0] eq $cyc[-1] && $#cyc>0) { # closed cycle
         $l--;
         $nos *= 2*(fibo($l+1)+fibo($l-1));
      }
      else { # open path
         $nos *= 2*(fibo($l+1)+fibo($l));
      }
   }
   return $nos;
}

sub class_dist {
   my %classes = ('A' => 'A',
                  'C' => 'C',
                  'G' => 'G',
                  'T' => 'U',
                  'U' => 'U',
                  'R' => 'AG',  # purine
                  'Y' => 'CT',  # pyrimidine
                  'S' => 'CG',
                  'M' => 'AC',
                  'W' => 'AU',
                  'K' => 'GU',
                  'V' => 'ACG', # not T
                  'H' => 'ACU', # not G
                  'D' => 'AGU', # not C
                  'B' => 'CGU', # not A
                  'N' => 'ACGU');

   my @seq = split(//, uc $_[0]);
   my @template = split(//, uc $_[1]);

   my $d=0;
   foreach (@seq) {
      my $c = $classes{shift @template} || 'N';
      $d++ unless /[$c]/;
   }
   return $d;
}

sub usage {
   print <<EOF;
usage: $0 [options] [infile]
program specific options:
 -n <int>     run $0 <int> times (default: $noo)
 -trial <int> max number of sequences tested per run (default: $nom)
 -g <float>   small positivie parameter for costfunction (default: $small)
 -bar <float> barrier hight, seperating the two states (default: $bar)
 -dG <float>  energy difference between the two states (default: $dG)
 -T2 <float>  for temperature sensitive switches, temperature at which
              2nd structure is ground-state (default: undef)
standard Vienna RNA options:
 -T <float>   sets temperature to <float> (default 37.0C)
 -4           no extrastable tetraloops
 -d           no dangling ends
 -d2          use double dangles
 -noGU        no GU or UG basepairs allowed
 -noClosingGU no closing GU or UG basepairs allowed
 -noLP        no lonely basepairs allowed
 -circ        design circular RNA
 -P <string>  read fold-parameters from file <string>
infile format:
1st-line structure one
2nd-line structure two
3rd-line sequence constraints or empty line

EOF
    exit;
}

=head1 NAME

switch.pl - design bistable RNA sequences

=head1 SYNOPSIS

  switch.pl [options] [infile]

=head1 DESCRIPTION

switch.pl designs RNA sequences that exhibit two secondary
structures of almost equal stability. For any two given structures
there always exist many sequences compatible with both
structures. If both structures are reasonable stable, we can find
sequences where both target structures have almost equal energy, and
all other structures have much higher energies.

For details of the algorithm see:  Flamm et al.,
"Design of Multi-Stable RNA Molecules", RNA 7:254-265 (2001)

Input consists of three lines, the first two containing the target
structures in dot bracket notations. Structures may contain a '&' to 
model bistability of two interacting RNAs (using RNAcofold algorithm). The third line may be used to
define sequence constraints: It contains a sequence string using
IUPAC codes for nucleotide classes (i.e. C<Y> for pyrimidine, C<R> for
purine, C<N> for anything...). If the line is empty or shorter than the
structures it is padded with C<N>s.

Sequence constraints are not strictly enforced, instead a constraint
violation term is added to the cost function. Thus it is legal to
specify sequence constraints that cannot be fulfilled.

switch.pl uses the Vienna RNA package for energy evaluation, the
Vienna RNA package and corresponding Perl module therefore have to be
installed. 

=head1 OPTIONS

=over 4

=item B<-n> <int>

number of independent optimization runs

=item B<-trial> <int>

maximum number of sequences to test per run, up to a million or so.

=item B<-g> <int>

Parameter of the cost function that weights the importance of
equal energies, and desired energy barriers.

The cost function primarily optimizes product of Boltzman
probabilities of the two structures C<P(S1)*P(S2)>, in addition it
contains a penalty proportional to C<[E(S1)-E(S2)]^2> that
enforces equal energies for both structures. The -g parameter
defines the weight of these additional cost function terms.

=item B<-T> <float>

Temperature in C for all energy calculations. Default 37C.

=item B<-T2> <float>

For temperature sensitive switches, use -T and -T2 to define the
temperatures at which structures S1 and S2 should be prefered.

=item B<-bar> <float>

Size of the desired energy barrier between the two structures in
kcal/mol. A fast heuristic that looks at shortest refolding paths is
used to estimate the barrier.

=item B<-dG> <float>

Specify an energy difference between the two optimized structures 
(kcal/mol). The energy is added to the second conformation.

=item ViennaRNA standard options

the -noLP, -P <paramfile>, -d, -d2, -4, -noGU, -noClosingGU, -circ,
should work as usual. See the RNAfold man page for details.

=back

=head1 AUTHORS

Ivo L. Hofacker, Christoph Flamm, Peter Stadler, Sebastian
Maurer-Stroh, Martin Zehl.
Send comments to <ivo@tbi.univie.ac.at>.

=cut

#  End of file
