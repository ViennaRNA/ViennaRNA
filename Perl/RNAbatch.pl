#!/usr/bin/perl -T
# -*-Perl-*-
# Last changed Time-stamp: <2004-08-19 10:15:27 ivo>
use HTML::Entities;
use URI::Escape;
use Chart::Lines;
use RNA;
use Mail::Send;
use strict;
use warnings;

my $paramdir = '/var/www/RNA';  # where parameters files etc are stored

$ENV{PATH} = '/bin:/usr/bin:/usr/local/bin';
#print STDERR "RNAbatch options:", join(' : ', @ARGV), "  ", `pwd`, "\n";
my $ssv_home = 'http://smi-web.stanford.edu/projects/helix/sstructview/home.html';
my $ssv_url  = '/~ivo/RNA/SStructView.zip';

my $options=shift;

$options =~ /([-\w,.]+)/; $options= $1; # $options now untainted
$options =~ s/,/ /g;
$options =~ s| dna.par| $paramdir/dna.par|;
open(FOLD, "RNAfold $options < job.in|") or die "can't start RNAfold: $!";
sleep(10);
open(OUT, ">>index.html") or die "can't append to output file: $!";
my $old_stdout = select(OUT);
print "An equivalent RNAfold command line would have been<BR>\n";
print "<kbd>RNAfold $options</kbd><P>\n";

$_ = <FOLD>; /^> ?(\S+)/; my $name = $1;
my $Sequence = <FOLD>; chomp($Sequence);
$Sequence =~ s/U/T/g  if ($options =~ /dna.par/);
$_ = <FOLD>; /([(.)]+) *\( *(-?\d+\.\d+)\)/;
my ($structure, $mfe) = ($1, $2);

print "The optimal secondary structure in bracket notation is given below\n";
my $eseq = $Sequence;
encode_entities($eseq);

print "<P><PRE>\n";
print "&gt; $name\n$eseq\n$structure";
printf " (%6.2f)\n", $mfe;
print "</PRE><P>\n";
print "The minimum free energy in kcal/mol is given in parenthesis.\n";

if (!($structure =~ /\(/)) {$RNA::rna_plot_type = 0};
sleep(1); # avoid race condition with RNAfold writing the _ss.ps file
RNA::PS_rna_plot($Sequence, $structure, $name .'_ss.ps') unless -f $name .'_ss.ps';
print "You may look at the PostScript drawing of the structure\n";
my $efname = $name. '_ss.ps'; encode_entities($efname);
my $ss_escaped = uri_escape($name . '_ss.ps');
print "in <A href=\"$ss_escaped\">$efname</a>.\n";
print "[For long sequences these plots get cluttered and ugly]<p>\n";

if ($options =~ /-p/) {
  $_ = <FOLD>; /(\S+) *\[ *(-?\d+\.\d+)\]/;
  my ($pf_struct, $energy) = ($1, $2);
  print "The free energy of the thermodynamic ensemble is ";
  printf "<b>%6.2f</b> kcal/mol\n<BR>", $energy;
  $_ = <FOLD>; print;
  my $fname_dp = $name . '_dp.ps';
  my $dp_escaped = uri_escape($fname_dp);
  print "<BR>The PostScript <em>dot plot</em> containing the base pair ",
    "probabilities is in\n<A href=\"$dp_escaped\">",
      encode_entities($fname_dp), "</a>.\n";
}
print "<P>";

mountain($name, $structure, $options) if grep /-mnt/, @ARGV;

if (grep /-ssview/, @ARGV) {
  my $fname_co = $name . '.coords';
  my $ret = RNA::ssv_rna_plot($Sequence, $structure, $fname_co);
  if ($ret == 0) {
    print "error: ssv_rna_plot plot failed, sorry<HR>\n";
    return;
  }
  my $escaped = uri_escape($fname_co);
  print "Please be patient while the applet is loading<br>\n";
  print "<APPLET ARCHIVE=\"$ssv_url\" CODE=\"SStructView.class\" ",
    "name=\"SSApplet\" WIDTH=500 HEIGHT=450>\n",
    "<PARAM name=\"structure-data-URL\" ",
    "value=\"$escaped\">\n",
    "<PARAM name=\"show-controls\" value=\"true\">\n",
#   "<PARAM name=\"result-display-frame\" value=\"_blank\">\n",
    "<PARAM name=\"default-background\" value=\"white\">\n",
    "</APPLET><HR>\n";
    print "For information on the java applet see the ",
    "<A href=\"$ssv_home\">SStructView</a> home page<p>\n";
}

if (grep /-svg/, @ARGV) {
  my $fname_svg = $name . '.svg';
  my $ret = RNA::svg_rna_plot($Sequence, $structure, $fname_svg);
  if ($ret == 0) {
    print "error: svg_rna_plot plot failed, sorry<HR>\n";
    return;
  }
  my $escaped = uri_escape($fname_svg);
  print "<iframe src=\"$escaped\" WIDTH=\"452\" HEIGHT=\"452\">\n",
    "your browser does not understand the &lt;iframe&gt; tag\n",
      "</iframe><br>\n";
  print 'Note: you need the <a href="http://www.adobe.com/svg">',
    "Adobe SVG plugin</a> or an SVG enabled browser to view the graphics<hr>\n";
}

print "Your output files will be deleted from this server after about one day.\n";

my ($user,$system,$cuser,$csystem) = times;
print "<BR>Time used for this call ",$user+$system+$cuser+$csystem,
  " seconds\n";
print "<HR>\n";
print <<END;
   We appreciate your feedback. Please send comments to
   <ADDRESS>Ivo Hofacker
   <A HREF="mailto:ivo\@tbi.univie.ac.at">&lt;ivo\@tbi.univie.ac.at&gt;</A>
   </ADDRESS><BR>
   <A HREF="http://www.tbi.univie.ac.at/~ivo/RNA/">Vienna RNA Home Page</A>
</body></html>
END

close(FOLD);
select($old_stdout);
close(OUT);

my $email = pop;
my $dir = $1 if $ENV{PWD} =~ /(\w+)$/;
if ($dir  && ($email ne 'rna@tbi.univie.ac.at')) {			       
    my $msg = new Mail::Send Subject=>'ViennaRNA fold job complete',To=>$email;
    $msg->set('Reply-To', 'rna@tbi.univie.ac.at');
    my $fh = $msg->open;
    print $fh "Your RNAfold job has completed.\n",
    "You should find the results at\n",
    "http://rna.tbi.univie.ac.at/RNAfold_dir/$dir/\n",
    "Your output will be deleted from our server after one day\n\n",
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
