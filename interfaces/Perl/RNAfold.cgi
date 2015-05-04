#!/usr/bin/perl
# -*-Perl-*-
# Last changed Time-stamp: <2002-11-19 15:40:22 ivo>
# CGI script for a Web-based RNA fold server
# you need to have the perl5 RNA module installed
# that comes as part of the Vienna RNA package

use strict;
use CGI qw/ :standard /;
use CGI::Carp qw(fatalsToBrowser);
use HTML::Entities;
use URI::Escape;
use RNA;
use Chart::Lines;
use vars qw/$hdir $maxlength $ssv_home $ssv_url $paramdir $serverroot/;

# please configure these variables for your site

$paramdir = '/var/www/RNA';  # directory where parameters files etc are stored
$serverroot = '/u/html';     # where your web pages lie
$hdir = "/RNAfold_dir";      # output directory (relative to $serverroot)
$maxlength = 500;            # only process sequences up to this length
$ssv_home = "http://smi-web.stanford.edu/projects/helix/sstructview/home.html";
$ssv_url  = "/~ivo/RNA/SStructView.zip";

my $query = new CGI;

if ($query->param) {
   &do_work($query);
} else {
   &print_header($query);
}

&print_tail;
print $query->end_html,"\n";

sub print_header {
   my($query, $cookie) = @_;
   if ($cookie) {
      print $query->header(-expires=>'+1m',
			   -cookie=>$cookie);
   } else {
      print $query->header(-expires=>'+1m');
   }
   print $query->start_html(-title=>"RNAfold input form",
			    -author=>'ivo@tbi.univie.ac.at',
			    -BGCOLOR=>'#f8f8ff');
   print "\n<H1 align=\"center\">Vienna RNA Secondary Structure Prediction</H1>\n";
   print "<H2 align=\"center\"><a name=\"top\">",
   "A web interface to the RNAfold programm</a></H2>\n";

   print "This server will predict secondary structures of single stranded\n",
   "RNA or DNA sequences. If the options look confusing <strong>read the\n",
   "<a href=\"http://www.tbi.univie.ac.at/~ivo/RNA/RNAcgi.html\">",
   "help page</a></strong><p>\n";

   print "<b><u>News:</u></b> based on ViennaRNA-1.5<br>\n",
     "Try the new SVG plot if your browser supports it!\n<p>\n";
   &print_prompt($query);
}

sub print_prompt {
   my($query) = @_;

   print $query->startform(-action=>'/cgi-bin/RNAfold.cgi#Results');
   print "<strong>Name of sequence</strong> ",
   "(choose a unique identifier)<BR>\n&gt; ";
   print $query->textfield('name');


   print "<P>\n<strong>Type in your sequence</strong>\n<kbd>T</kbd>s will be ",
   "automatically replaced by <kbd>U</kbd>s.\nAny symbols except ",
   "<kbd>AUCGTXKI</kbd> will be interpreted as nonbonding bases. ",
   "Maximum sequence length is currently $maxlength <BR>\n";

   print $query->textarea(-name=>'Sequence',
			  -rows=>5,
			  -columns=>70,
			  -default=>"");

   print "<P>\n<strong>Choose Fold Algorithm</strong><BR>\n";

   print $query->popup_menu(-name => 'pffold',
			    -values => ['mfe','pf'],
			    -default => 'pf',
			    -labels => {'mfe', 'minimum free energy only',
					'pf',
					'partition function and pair probabilities'}
			    );
   print $query->popup_menu(-name => 'Params',
			    -values => ['RNA','oldRNA','DNA'],
			    -default => 'RNA',
			    -labels => {'RNA', 'use RNA parameters',
					'oldRNA', 'old RNA parameters',
					'DNA', 'use DNA parameters'}
			    );

   print "<H4>Options to modify the fold algorithm</H4>\n";
   print "Rescale energy parameters to temperature\n",
     $query->textfield(-name=>'Temp',
		       -default=>37,
		       -size=>3);
   print "C<P>\n";
   print $query->checkbox_group(
				-name=>'toggles',
				-values=>['-4', '-d',  '-noGU', '-noCloseGU',
					  '-noLP'],
				-labels=>{'-4', 'no special tetraloops',
					  '-d', 'no dangling end energies',
					  '-noGU', 'no GU pairs',
					  '-noCloseGU',
					  'no GU pairs at the end of helices',
					  '-noLP', 'avoid isolated base pairs'
					   },
				-linebreak=>'yes',
				-defaults=>[]);

   print "\n<BR><strong>",
     "Should we produce a mountain plot of the structure?</strong> \n",
       $query->checkbox(-name=>'plot', -checked=>'ON');
   print "\n<BR>View a plot of the mfe structure inline using\n",
     "an SVG image (may requires plugin)\n", 
       $query->checkbox(-name=>'SVG', checked=>'ON');
   print "<br>\nor using the <code>SStructView</code> java applet?\n",
     $query->checkbox(-name=>'SSview');

   print "<P>\n",$query->reset, "\n";
   print $query->submit('Action','Fold it'), "\n";
   print $query->endform, "\n";

   if (!$query->param) {
      my @old_files = $query->cookie('old_files');
      if (@old_files) {
	 print "<P>The results from your last run should still be in\n";
	 foreach my $f (@old_files) {
	    print " <A href=\"$hdir/$f\">$f</a>";
	 }
	 print ".\n";
      }
   }
   print "<HR>\n";
}

sub do_work {
#   setpriority(PRIO_PROCESS, 0, 19);
   my($query) = @_;
   my($OK_CHARS) = '-a-zA-Z0-9_.@ ';
   my($RNAfold_id, $junk);
   my($WORK_DIR) = $serverroot . $hdir;
   my($sfact) = 1.02;
   chdir $WORK_DIR || die("can't change directory");

   # clean old files
   foreach my $f (< *.ps *.png *.coords *.svg>) {
      unlink($f) if (-M $f)>1.5;
   }

   my $Sequence = $query->param('Sequence');
   if (length($Sequence) < 1) {   # do nothing
      &print_header($query);
      return;
   }
   $Sequence =~ s/[\n\s]//gm;   # remove whitespace and linebreaks
   $Sequence = uc $Sequence;    # upper case
   if ($query->param('Params') eq 'DNA') {
      $Sequence =~ s/U/T/g;     # U -> T
   } else {
      $Sequence =~ s/T/U/g;     # T -> U
   }

   my $name = $query->param('name');

   my $options  = "";
   $options .= " -p" if ($query->param('pffold') eq 'pf');
   $options .= " -P dna.par" if $query->param('Params') eq 'DNA';
   $options .= " -P vienna13.par" if $query->param('Params') eq 'oldRNA';
   $options .= " " . join(' ', $query->param('toggles'));
   $RNA::tetra_loop = 0 if $options =~ /-4/;
   $RNA::dangles = 0 if $options =~ /-d/;
   $RNA::noGU = 1 if $options =~ /-noGU/;
   $RNA::no_closingGU = 1 if $options =~ /-noCloseGU/;
   $RNA::noLonelyPairs = 1 if $options =~ /-noLP/;

   if ($query->param('Temp') =~ /([-+]?\d+\.?\d*)/) {
     $RNA::temperature = $1;
     $options .= " -T $1" if $1 != 37;
   }

   # remove strange characters, so we can pass strings to shell
   $options =~ s/[^$OK_CHARS]//go;
   if ($name =~ /(\S+)/) {
      $name = $1;
   } else {
      $name = "s$$";
   }

   $name =~ s/\///g;    # no / allowed in file names
   my $fname_ss = $name . "_ss.ps";
   my $fname_dp = ($options =~ /-p /)?($name . "_dp.ps"):"";
   my $fname_co = $name . ".coords";
   my $fname_svg = $name . ".svg";
   my $ss_escaped = uri_escape($fname_ss);
   my $dp_escaped = uri_escape($fname_dp);
   my $png_escaped = uri_escape($name . ".png");

   my $the_cookie = $query->cookie(-name=>'old_files',
				-value=>[$ss_escaped,$dp_escaped,$png_escaped],
				-expires=>'+1d',
				-path=>'/cgi-bin');

   &print_header($query, $the_cookie);  # don't print anything before this line

   print "<H2><a name=\"Results\">Here are your results</a></H2>\n";

   if ($query->param('Params') eq 'DNA') {
      print "DNA parameters provided by courtesy of " .
	  "<b><a href=\"http://sun2.science.wayne.edu/~jslsun2/\">",
	  "John SantaLucia Jr.</a></b><br>\n",
	  "In any publication using these relsults, please cite:<br>\n",
	  "SantaLucia, J Jr (1998) ",
	  "\"A unified view of polymer, dumbbell, and oligonucleotide DNA ",
	  "nearest-neighbor thermodynamics.\"\n",
	  "<i>Proc. Natl. Acad. Sci. USA</i> <b>95</b>, 1460-1465.\n",
	  "<p><hr>\n";
   } else {
      print "RNA parameters are described in<br>\n" .
	  "D.H.  Mathews, J. Sabina, M. Zuker and H. Turner\n",
	  "\"Expanded Sequence Dependence of Thermodynamic Parameters ",
	  "Provides Robust Prediction of RNA Secondary Structure\",\n",
	  "JMB, 288, pp 911-940, 1999<hr>\n";
   }
#   print $query->query_string, "<br>\n";
   my $length = length($Sequence);
   if ($length > $maxlength) {
      print "Sorry, this interface is still being tested.\n",
      "Currently, sequences longer than $maxlength bases will currently ",
      "not be processed. Your sequence has length $length.<HR>\n";
      return;
   }

   print "An equivalent RNAfold command line would have been<BR>\n";
   print "<kbd>RNAfold $options</kbd><P>\n";
   if ($query->param('Params') eq 'DNA') {
    RNA::read_parameter_file("$paramdir/dna.par");
   } elsif ($query->param('Params') eq 'oldRNA') {
    RNA::read_parameter_file("$paramdir/vienna13.par");
   }

   my ($structure, $mfe) = RNA::fold($Sequence);
   print "The optimal secondary structure in bracket notation is given below\n";
   my $eseq = $Sequence;
   encode_entities($eseq);

   print "<P><PRE>\n";
   print "&gt; $name\n$eseq\n$structure";
   printf " (%6.2f)\n", $mfe;
   print "</PRE><P>\n";
   print "The minimum free energy in kcal/mol is given in parenthesis.\n";

   if (!($structure =~ /\(/)) {$RNA::rna_plot_type = 0};
   RNA::PS_rna_plot($Sequence, $structure, $fname_ss);
   print "You may look at the PostScript drawing of the structure\n";
   my $efname = $fname_ss; encode_entities($efname);
   print "in <A href=\"$hdir/$ss_escaped\">$efname</a>.<P>\n";

   if ($options =~ /-p/) {
     $RNA::dangles = 2 if $RNA::dangles;
     my $min_en = RNA::energy_of_struct($Sequence, $structure);
     my $kT = ($RNA::temperature+273.15)*1.98717/1000.; # in Kcal
     $RNA::pf_scale = exp(-($sfact*$min_en)/$kT/$length);
     my ($pf_struct, $energy) = RNA::pf_fold($Sequence);
     print "The free energy of the thermodynamic ensemble is ";
     printf "<b>%6.2f</b> kcal/mol\n", $energy;
     RNA::PS_dot_plot($Sequence, $fname_dp);
     print "<BR>The PostScript <em>dot plot</em> containing the base pair ",
       "probabilities is in\n<A href=\"$hdir/$dp_escaped\">",
	 encode_entities($fname_dp), "</a>.\n";
   }
   print "<P>";

   mountain($name, $structure, $options) if $query->param('plot') eq 'on';

   if ($query->param('SSview') eq 'on') {
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
      "value=\"$hdir/$escaped\">\n",
      "<PARAM name=\"show-controls\" value=\"true\">\n",
#	"<PARAM name=\"result-display-frame\" value=\"_blank\">\n",
      "<PARAM name=\"default-background\" value=\"white\">\n",
      "</APPLET><HR>\n";
      print "For information on the java applet see the ",
      "<A href=\"$ssv_home\">SStructView</a> home page<p>\n";
   }

   if ($query->param('SVG') eq 'on') {
      my $ret = RNA::svg_rna_plot($Sequence, $structure, $fname_svg);
      if ($ret == 0) {
	 print "error: svg_rna_plot plot failed, sorry<HR>\n";
	 return;
      }
      my $escaped = uri_escape($fname_svg);
      print "<iframe src=\"$hdir/$escaped\" WIDTH=\"452\" HEIGHT=\"452\">\n",
	"your browser does not understand the &lt;iframe&gt; tag\n",
	  "</iframe><br>\n";
      print 'Note: you need the <a href="http://www.adobe.com/svg">',
	"Adobe SVG plugin</a> or an SVG enabled browser to view the graphics<hr>\n";
   }

   print "Your output files will be deleted from this server after one day.\n";

#    $myself = $query->self_url;
#    print "<BR><A HREF=$myself#top>Back to the input form</A>";
   print "<BR>Scroll back to the top to submit another sequence.\n";

   my ($user,$system,$cuser,$csystem) = times;
   print "<BR>Time used for this call ",$user+$system, " seconds\n";
   print "<HR>\n";
}

sub print_tail {
   print <<END;
   We appreciate your feedback. Please send comments to
   <ADDRESS>Ivo Hofacker
   <A HREF="mailto:ivo\@tbi.univie.ac.at">&lt;ivo\@tbi.univie.ac.at&gt;</A>
   </ADDRESS><BR>
   <A HREF="http://www.tbi.univie.ac.at/~ivo/RNA/">Vienna RNA Home Page</A>
END
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

  #my $bp = $RNA::base_pair;
  #my $npair = RNA::bond_i_get($bp);
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
    my $bpp = RNA::Make_bp_profile($length);
    my @bpp = unpack("f*",RNA::cdata($bpp, ($length+1)*4*3));
    $mp[0] = $mp[1] = 0;  # use indices [1..length]
    for my $i (1..$length) {
      $mp[$i]  -= $bpp[$i*3+2]; # paired downstream
      $mp[$i+1] = $mp[$i] + $bpp[$i*3+1]; # upstream
    }

    $obj->add_dataset (@mp[0..$length]);
  }

  $name .= ".png";
  $obj->png("$name");
  my $png_escaped = uri_escape($name);
  print "<IMG ALT=\"Mountain Representation\" HEIGHT=\"$height\" ",
    "WIDTH=\"$width\" SRC=\"$hdir/$png_escaped\"><P>\n";
}
