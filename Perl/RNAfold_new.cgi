#!/usr/bin/perl -T
# -*-CPerl-*-
# Last changed Time-stamp: <2005-07-15 18:14:39 ivo>
# CGI script for a Web-based RNA fold server
# you need to have the perl5 RNA module installed
# that comes as part of the Vienna RNA package

use strict;
use CGI qw/ -no_xhtml :standard /;
use HTML::Entities;
use URI::Escape;
use RNA;
use Chart::Lines;
use Email::Valid;
use File::Temp qw/ tempdir/;
use warnings;
no warnings qw(uninitialized);
use vars qw/$hdir $maxlength1 $maxlength2 $maxlength_pf $batchscript
            $ssv_home $ssv_url $paramdir $serverroot $logfile/;

# please configure these variables for your site

$paramdir = '/var/www/RNA';  # directory where parameters files etc are stored
$serverroot = '/u/html';     # where your web pages live
$hdir = '/RNAfold_dir';      # output directory (relative to $serverroot)
$maxlength1 =  300;          # max length for immediate jobs
$maxlength2 = 5000;          # max length for queued jobs
$maxlength_pf = 4000;
$CGI::POST_MAX = 16*1024;    # prevent absurdly huge posts right away
$ssv_home = 'http://smi-web.stanford.edu/projects/helix/sstructview/home.html';
$ssv_url  = '/~ivo/RNA/SStructView.zip';
$batchscript = '/var/www/RNA/RNAbatch.pl'; # script to use for batch submissions
# for local test
#$serverroot='.';
#$maxlength1=100;

BEGIN {
  use CGI::Carp qw(carpout fatalsToBrowser);
  # redirect errors to a log file instead of the web server log
  my $log = '/var/www/RNA/RNAfold_cgi.log'; 
  no strict "subs";
  open(LOG, ">>$log") && carpout(LOG)
      or warn "Unable to open logfile $log: $!\n";
}

END {close(LOG);}

$ENV{PATH} = '/bin:/usr/bin:/usr/local/bin';
delete @ENV{'IFS', 'CDPATH', 'ENV', 'BASH_ENV'};


my $query = new CGI;

if ($query->cgi_error) {
  print $query->header(-status=>$query->cgi_error),
    $query->start_html('Problems'),
      $query->h2('Request not processed'),
        $query->strong($query->cgi_error);
  exit 0;
}

if ($query->param && $query->param('Action') =~ /Fold/) {
   do_work($query);
} else {
   print_form($query);
   print_tail();
}


sub print_form {
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
     "Try the new SVG plot if your browser supports it!<br>\n",
       "You can now submit sequences up to $maxlength2 as batch jobs.<p>\n";
   print $query->startform(-enctype=>'multipart/form');
   print "Name of sequence (optional, used to name output files)<BR>\n&gt; ";
   print $query->textfield('name');

   print "<P>\n<strong>Type in your sequence</strong>\n<kbd>T</kbd>s will be ",
   "automatically replaced by <kbd>U</kbd>s.\nAny symbols except ",
   "<kbd>AUCGTXKI</kbd> will be interpreted as nonbonding bases.\n",
   "Any non-alphabetic characters will be removed.<br>\n";

   print $query->textarea(-name=>'Sequence',
			  -rows=>5,
			  -columns=>70,
			  -default=>"");
   print "\n<br>Maximum sequence length for immediate jobs is $maxlength1.\n",
     "Sequences up to $maxlength2 (mfe only) or $maxlength_pf ",
       "(pair probabilities) will be queued as batch jobs\n";

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
				-values=>['-4', '-d',  '-noCloseGU',
					  '-noLP'],
				-labels=>{'-4', 'no special tetraloops',
					  '-d', 'no dangling end energies',
					  '-noCloseGU',
					  'no GU pairs at the end of helices',
					  '-noLP', 'avoid isolated base pairs'
					   },
				-linebreak=>'yes',
				-defaults=>[-noLP]);

   print "\n<BR><strong>",
     "Should we produce a mountain plot of the structure?</strong> \n",
       $query->checkbox(-name=>'plot', -checked=>'ON');
   print "\n<BR>View a plot of the mfe structure inline using\n",
     "an SVG image (may require plugin)\n", 
       $query->checkbox(-name=>'SVG', checked=>'ON');
   print "<br>\nor using the <code>SStructView</code> java applet?\n",
     $query->checkbox(-name=>'SSview');
   print "<br>\nEmail address. When the job has completed, we'll send a ",
         "mail containing a link to the results page, this is useful for ",
	 "long jobs that won't give results immediately. Please ",
	 "don't use fake addresses (just leave the filed as is, or empty).\n";
   print $query->textfield(-name => 'email',
			   -default => 'you@where.org',
			   -size => 30);
   print "<P>\n",$query->reset, "\n";
   print $query->submit('Action','Fold it'), "\n";
   print $query->endform, "\n";

   if (!$query->param) {
     my $last = $query->cookie('last');
     if ($last) {
       print "<P>The results from your last run should still be in\n";
       print " <A href=\"$last/\">$last</a>.\n";
     }
   }
   print "<HR>\n";
}

sub do_work {
   setpriority(0, 0, 10);
   my($query) = shift;
   my($OK_CHARS) = '-a-zA-Z0-9_.@ ';
   my($RNAfold_id, $junk);
   my($sfact) = 1.02;
   my $WORK_DIR = tempdir( "foldXXXXXX", DIR => "$serverroot$hdir");
   $hdir = substr($WORK_DIR, length($serverroot));
#   $hdir .= "/fold$$";
#   my($WORK_DIR) = $serverroot . $hdir;
#   mkdir $WORK_DIR or warn "can't create directory $WORK_DIR";

   chdir $WORK_DIR or die"can't change directory";

   open(STATE, ">cgistate.txt");
   $query->save(\*STATE);
   close(STATE);

   # clean old files
#   foreach my $f (< *.ps *.png *.coords *.svg>) {
#      unlink($f) if (-M $f)>1.5;
#   }
   foreach my $f (<../fold*>) {
     if ($f =~ /(\.\.\/fold\d+)/) {
       system('/bin/rm', '-r', $1) if (-d $1) && (-M $1)>3;
     }
   }

   my $Sequence = $query->param('Sequence');
   if (length($Sequence) < 1) {   # do nothing
      print_form($query);
      return;
   }
   $Sequence =~ s/[\n\s]//gm;   # remove whitespace and linebreaks
   $Sequence = uc $Sequence;    # upper case

   $Sequence =~ s/T/U/g;        # T -> U
   $Sequence =~ s/[^A-Z]//g;    # remove non-alpha characters
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

   $RNA::temperature=$1 if ($query->param('Temp') =~ /([-+]?\d+\.?\d*)/);

   # remove strange characters, so we can pass strings to shell
   $options =~ /([$OK_CHARS]+)/; $options= $1; # $options now untainted

   if ($name =~ /([-\@\w.]+)/) {
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

   my $the_cookie = $query->cookie(-name=>'last',
				-value=>"$hdir",
				-expires=>'+1d',
				-path=>'/cgi-bin');

   my $time=estimate(length($Sequence), $query->param('pffold'));
   print $query->header(-refresh=>"$time; URL=$hdir/",
		    -status=>'202 Accepted',
		    -cookie=>$the_cookie,
		    -type=>'text/html');
   print $query->start_html(-title=>"RNAfold in progress",
			-author=>'ivo@tbi.univie.ac.at',
			-BGCOLOR=>'#f8f8ff');

   print "<H2>Now folding... </H2>\n";
   my $length = length($Sequence);

   # do checks
   if ($length > $maxlength2) {
      print "Sorry, this interface is still being tested.\n",
      "Currently, sequences longer than $maxlength2 bases will currently ",
      "not be processed. Your sequence has length $length.<HR>\n";
      print_tail();
      return;
   }
   if ($RNA::temperature>1000 || $RNA::temperature<=-273.15 ||
       ($RNA::temperature<-200 && $options =~ /-p/)) {
     print "<b>Warning</b>: Temperature is out of range using 37C<p>\n";
     $RNA::temperature = 37;
   }
   $options .= " -T $RNA::temperature" if $RNA::temperature != 37;

   $maxlength1 += 100 unless $options =~ /-p/;
   my $email;
   if ($length >$maxlength1) {
       $email = undef if $email eq 'you@where.org';
       $email = Email::Valid->address($query->param('email'), -mxcheck => 1);
       if (!defined($email))  {
	   print "It seems you have no given a valid mail address ",
	   encode_entities($query->param('email')), ".<br>\n",
	   "Be sure you retain the link below so you can find your ",
	   "resluts when the calculation is finished<p>\n";
	   #print_tail();
	   #return;
       }
   }

   if ( ($Sequence =~ tr/AUGCT//) < 0.5*length($Sequence)) {
       print "ERROR: less than 50% of your sequence is normal nucleotides ",
       "(AUGCT). Surely, there must have been a mistake. Please resubmit ",
       "the correct sequence.<p>\n";
       print_tail();
       return;
   }
   if ( ($Sequence =~ tr/0-9//) > 3 + 0.01*length($Sequence)) {
     print "ERROR: There are lots of digits (0-9) in your sequence less than 50% of your sequence is normal nucleotides ",
       "(AUGCT). Surely, there must have been a mistake. Please resubmit ",
	 "the correct sequence.<p>\n";
     print_tail();
     return;
   }

   if (($length > $maxlength_pf) && ($options =~ /-p/)) {
     print "Won't calculate partition function for sequences longer than ",
       "$maxlength_pf, doein mfe part only...<p>\n";
     $options =~ s/-p +//;
   }

   # if get here we should have a valid request
   if ($length<=$maxlength1) {
     print "You should be automatically forwarded to your ",
       "results page after $time seconds.<br>\n",
	 "Otherwise, click <a href=\"$hdir/\">here</a> to view your ",
	   "results after a few seconds<p>\n";
   } else {
     open(FOLD_IN, ">job.in") or croak "can't write input file for batch job";
     print FOLD_IN ">$name\n$Sequence\n";
     close FOLD_IN;
     my $args = ',' . $options; # make sure $args is not empty;
     $args =~ s/ +/,/g;
     $args .= ' --mnt'    if $query->param('plot') eq 'on';
     $args .= ' --svg'    if $query->param('SVG') eq 'on';
     $args .= ' --ssview' if $query->param('SSview') eq 'on';
     $args .= ' ' .  $email;
#     system("$batchscript @args &");
     my $qu = 'RNA';
     $qu = 'RNAshort' if $length<2*$maxlength1;
     open(QU, '|-', "/usr/bin/qsub -q $qu -ln 15 -s /bin/sh") or
	 die "can't submit batch job: $!";
     print QU "cd $WORK_DIR; $batchscript $args\n";
     close QU;
     warn "submitted job $WORK_DIR of length $length for $email";
     print "<p>Your job has been submitted to our queuing system.<br>\n",
     "Assuming no other load on the system, expect results to be available ",
     "in about ";
     printf "%.2f %s", $time/60, " minutes.<br>\n";
     print "You should be notified by mail when the computation is complete. ",
       "<br>You can then follow the link to ",
	 "<a href=\"$hdir/\">your results page</a><p>\n";
   }

   print_tail();

   # start writing output file

   open(OUT, ">index.html");
   my $old_stdout = select(OUT);
   print $query->start_html(-title=>'RNAfold results',
                       -author=>'ivo@tbi.univie.ac.at',
                       -BGCOLOR=>'#f8f8ff');

   print "<H2><a name=\"Results\">Here are your RNAfold results</a></H2>\n";

   if ($query->param('Params') eq 'DNA') {
      print "DNA parameters provided by courtesy of " .
	  "<b><a href=\"http://ozone.chem.wayne.edu/\">",
	  "John SantaLucia Jr.</a></b><br>\n",
	  "In any publication using these relsults, please cite:<br>\n",
	  "SantaLucia, J Jr (1998) ",
	  "\"A unified view of polymer, dumbbell, and oligonucleotide DNA ",
	  "nearest-neighbor thermodynamics.\"\n",
	  "<i>Proc. Natl. Acad. Sci. USA</i> <b>95</b>, 1460-1465.\n",
	  "<p><hr>\n";
   } else {
      print "RNA parameters are described in<br>\n" .
	  "D.H.  Mathews, J. Sabina, M. Zucker and H. Turner\n",
	  "\"Expanded Sequence Dependence of Thermodynamic Parameters ",
	  "Provides Robust Prediction of RNA Secondary Structure\",\n",
	  "JMB, 288, pp 911-940, 1999<hr>\n";
   }

   if ($length>$maxlength1) {
     print "[If this file is incomplete (because the computation isn't finished) reload this page after a while]<p>\n";
     select($old_stdout);
     close(OUT);
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
   my $seq = $Sequence;
   $Sequence =~ s/U/T/g if ($query->param('Params') eq 'DNA');   
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
   print "in <A href=\"$ss_escaped\">$efname</a>.<P>\n";

   if ($options =~ /-p/) {
     $RNA::dangles = 2 if $RNA::dangles;
     my $min_en = RNA::energy_of_struct($seq, $structure);
     my $kT = ($RNA::temperature+273.15)*1.98717/1000.; # in Kcal
     $RNA::pf_scale = exp(-($sfact*$min_en)/$kT/$length);
     my ($pf_struct, $energy) = RNA::pf_fold($seq);
     print "The free energy of the thermodynamic ensemble is ";
     printf "<b>%6.2f</b> kcal/mol\n", $energy;
     RNA::PS_dot_plot($Sequence, $fname_dp);
     print "<BR>The PostScript <em>dot plot</em> containing the base pair ",
       "probabilities is in\n<A href=\"$dp_escaped\">",
	 encode_entities($fname_dp), "</a>.\n";
   }
   print "<P>";

   # compute Tm of mfe structure
   if (($RNA::temperature == 37.) && ($mfe < 0)) {
     $RNA::temperature=-273.15;
     my $enth = RNA::energy_of_struct($Sequence, $structure);
     my $Tm = 310.15 * $enth/($enth - $mfe) - 273.15;
     printf "The enthalpy of the mfe structure is %6.2f corresponding to a Tm of %5.1fC<p>\n", $enth, $Tm;
 }
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
      "value=\"http://rna.tbi.univie.ac.at$hdir/$escaped\">\n",
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
      print "<iframe src=\"$escaped\" WIDTH=\"452\" HEIGHT=\"452\">\n",
	"your browser does not understand the &lt;iframe&gt; tag\n",
	  "</iframe><small>\n",
	    "Click anywhere on the canvas to toggle display of ",
	      "the sequence</small><br>\n";
      print 'Note: you need the <a href="http://www.adobe.com/svg">',
	"Adobe SVG plugin</a> or an SVG enabled browser to view the graphics<hr>\n";
   }

   print "Your output files will be deleted from this server after about one day.\n";

   my $myself = $query->self_url;
   $myself =~ s/Action.*;//;
   print "<BR><A HREF=$myself>Back to the input form</A>"
     if length($myself)<500;

   my ($user,$system,$cuser,$csystem) = times;
   print "<BR>Time used for this call ",$user+$system, " seconds\n";
   print "<HR>\n";
   print_tail();
   select($old_stdout);
   close(OUT);
}

sub print_tail {
   print <<END;
   We appreciate your feedback. Please send comments to
   <ADDRESS>Ivo Hofacker
   <A HREF="mailto:ivo\@tbi.univie.ac.at">&lt;rna\@tbi.univie.ac.at&gt;</A>
   </ADDRESS><BR>
   <A HREF="http://www.tbi.univie.ac.at/~ivo/RNA/">Vienna RNA Home Page</A>
END
   print $query->end_html,"\n";
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
    $mp[$length] = 0;
    $obj->add_dataset (@mp[0..$length]);
  }

  $name .= ".png";
  $obj->png("$name");
  my $png_escaped = uri_escape($name);
  print "<IMG ALT=\"Mountain Representation\" HEIGHT=\"$height\" ",
    "WIDTH=\"$width\" SRC=\"$png_escaped\"><P>\n";
}

sub estimate {
  # estimate time for folding from seq length
  # adjust these parameters for your server machine!
  my ($l, $pf) = @_;
  my $t = 1 + 8e-6 * $l*$l + 5.5e-9 * $l*$l*$l;
  $t *= 4 if ($pf eq 'pf');
  $t += 100 if $l > $maxlength1; # extra time for batch submission
  return int($t+1);
}

