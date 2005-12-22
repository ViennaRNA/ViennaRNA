#!/usr/bin/perl -T
# -*-Perl-*-
# Last changed Time-stamp: <2005-09-06 16:53:35 ivo>
use HTML::Entities;
use URI::Escape;
use Chart::Lines;
use RNA;
use Mail::Send;
use strict;
use warnings;
no warnings qw(uninitialized);

my $RNAdir = '/var/www/RNA/';   # where the RNAalifold executable lives
$ENV{PATH} = '/bin:/usr/bin:/usr/local/bin';
#print STDERR "RNAbatch options:", join(' : ', @ARGV), "  ", `pwd`, "\n";

my $email = pop;  # first argument is email to notify at completion
my $RNAz  = pop;  # second argument is "true" if RNAz should be run
my $options = join(' ', @ARGV);

$options =~ /([-\w .]+)/; $options= $1; # $options now untainted
sleep(10);
open(OUT, ">>index.html") or die "can't append to output file: $!";
my $old_stdout = select(OUT);
print "An equivalent RNAalifold command line would have been<BR>\n",
    "<kbd>RNAalifold $options alifold.aln</kbd><P>\n";
$options =~ s|dna.par|$RNAdir/dna.par|;

$| = 1;
print "<P><PRE>\n";
print `$RNAdir/RNAalifold $options alifold.aln 2>&1`;
print "</PRE><P>\n";

print "RNAalifold returned $?, maybe something went wrong.<br>\n" if ($?);
print "Other output files:<dl>\n",
    "<dd><a href=\"alirna.ps\">postscript drawing</a> ",
    "of the predicted structure<br>\n";
if ($options =~ /-p/) {
    print '<dd><a href="alidot.ps">postscript dot plot</a> ',
    "of pair probabilities<br>\n";
    
    print '<dd><a href="alifold.out">text output</a> with detailed ',
    "information about each base pair<br>\n";
    `$RNAdir/cmt.pl alifold.out > cmount.eps`;
    print "<dd>colored <a href=\"cmount.eps\">mountain plot</a> in postscript";
}
if ($RNAz =~ /true/) {
    $ENV{RNAZDIR} = "$RNAdir/RNAz_models";
    `$RNAdir/RNAz alifold.aln > RNAz.out`;
    print "<dd>The <a href=\"RNAz.out\">RNAz output</a> shows whether the alignment is likely to harbour a functional RNA structure.";
}
print "</dl>\n";
my ($user,$system,$cuser,$csystem) = times;
print "<BR>Time used for this call ",$user+$system, "+", $cuser+$csystem,
    " seconds<BR>\n";
print "Results will remain on this server for at least 1 day\n<HR>\n";
print <<END;
   We appreciate your feedback. Please send comments to
   <ADDRESS>Ivo Hofacker
   <A HREF="mailto:rna\@tbi.univie.ac.at">&lt;ivo\@tbi.univie.ac.at&gt;</A>
   </ADDRESS>
   <A HREF="http://www.tbi.univie.ac.at/~ivo/RNA/">Vienna RNA Home Page</A>
 </body></html>
END

select($old_stdout);
close(OUT);

my $dir = $1 if $ENV{PWD} =~ /(\w+)$/;
if ($dir) {			       
    my $msg = new Mail::Send Subject=>'ViennaRNA alifold job complete',To=>$email;
    $msg->set('Reply-To', 'rna@tbi.univie.ac.at');
    my $fh = $msg->open;
    print $fh "Your alifold job has completed.\n",
    "You should find the results at\n",
    "http://rna.tbi.univie.ac.at/RNAfold_dir/$dir/\n",
    "Your output will be deleted from our server after two days\n\n",
    "Contact rna\@tbi.univie.ac.at in case of problems\n";
    close $fh;
}
			       
sub mountain {
  my($name, $structure, $options) = @_;
  my $length = length $structure;
  my $width =  480;
  my $height = 300;

  # FIXME: legend_lables when doing mfe only
  my $skip = 10**(int (log($length)/log(10.) - 0.5));
  my $obj = Chart::Lines->new( $width, $height );
  $obj->set ('title' => $name,
	     'x_label' => 'Position',
	     'y_label' => 'Height',
	     'min_val' => 0,
	     'legend_labels' => ['mfe'],
	     'skip_x_ticks' => $skip);

  $obj->add_dataset ((0..$length));

  my @mp = (0); # valid indices 1..length
  my $h = 0;
  foreach (split(//, $structure)) {
    $h-- if $_ eq ')';
    push @mp, $h;
    $h++ if $_ eq '(';
  }
  warn "illegal structure $structure" if $mp[-1] != 0;
  $obj->add_dataset (@mp[0..$length]);  # @mp[1..$length] ??

  if ($options =~ /-p/) { # we have base pairing probs
    $obj->set('legend_labels' => ['mfe', 'pf']);
    open(PROB, '<', $name . '_dp.ps');
    @mp = (0);
    while (<PROB>) {
      my ($i, $j, $p, $id) = split;
      no warnings;
      if ($id eq "ubox") {
	$p *= $p;           # square it to probability
	$mp[$i+1] += $p;
	$mp[$j]   -= $p;
      }
    }
    close(PROB);
    for my $i (1..$length) {
      $mp[$i]+=$mp[$i-1];
    }
    $obj->add_dataset (@mp[0..$length]);
  }

  $name .= ".png";
  $obj->png("$name");
  my $png_escaped = uri_escape($name);
  print "<IMG ALT=\"Mountain Representation\" HEIGHT=\"$height\" ",
    "WIDTH=\"$width\" SRC=\"$png_escaped\"><P>\n";
}



# End of file
