#!/usr/local/bin/perl 
# -*-Perl-*-
# Last changed Time-stamp: <2001-10-30 16:43:21 ivo>
# CGI script for a Web-based RNA fold server
# you need to have the perl5 RNA module installed
# that comes as part of the Vienna RNA package

use strict;
use CGI qw/ :standard /;
#use CGI::Carp;# qw(fatalsToBrowser);
use HTML::Entities;
use URI::Escape;
use RNA;
use Chart::Lines;

$CGI::POST_MAX = 20480;
# please configure these variables
use vars qw/$hdir $maxlength $ssv_home $ssv_url/;

$hdir = "/RNAfold_dir";      # were the output files are stored
$maxlength = 300;            # only process sequences up to this length

my $q = new CGI;

if ($q->cgi_error) {
   print $q->header(-status=>$q->cgi_error),
   $q->start_html('Problems'),
   $q->h2('Request not processed'),
   $q->strong($q->cgi_error);
   exit 0;
}

if ($q->param) {
   do_work($q);
} else {
   print_form($q);
}


sub print_form {
   my($q, $cookie) = @_;
   if ($cookie) {
      print $q->header(-expires=>'+1m',
		       -cookie=>$cookie);
   } else {
      print $q->header(-expires=>'+1m');
   }
   print $q->start_html(-title=>"RNAalifold input form",
			-author=>'ivo@tbi.univie.ac.at',
			-BGCOLOR=>'#f8f8ff');
   print "\n<H1 align=center>Vienna RNA Secondary Structure Prediction</H1>\n";
   print "<H2 align=center><a name=top>",
   "A Web Interface for the Prediction Consensus Structures of Aligned Sequences</a></H2>\n";
   
   print "This server will predict secondary structures of single stranded\n",
   "RNA or DNA sequences. If the options look confusing <strong>read the\n",
   "<a href=\"http://www.tbi.univie.ac.at/~ivo/RNA/RNAcgi.html\">",
   "help page</a></strong><p>\n";
   
   print "<b><u>News</u></b> This service is brand new and probably still has lots of bugs<p>\n",
      
   $q->start_multipart_form(-action=>'/cgi-bin/alifold.cgi');
   
   print "<P><b>Enter</b> the file containing your <b>aligned sequences</b> ",
   "for uploading.<br>\n",
   "The alignment file <b>must</b> be in <kbd>Clustalw</kbd> format.<br>\n",
   "Maximum file size is ", int($CGI::POST_MAX/1024), "KB.<BR>\n";
   
   print $q->filefield(-name=>'clustal_file',
		       -default=>'starting value',
		       -size=>50,
		       -maxlength=>80);
   
   print "<P>\n<strong>Choose Fold Algorithm</strong><BR>\n";
   
   print $q->popup_menu(-name => 'pffold',
			-values => ['mfe','pf'],
			-default => 'pf',
			-labels => {'mfe', 'minimum free energy only',
				    'pf',
				    'partition function and pair probabilities'}
			);
   print $q->popup_menu(-name => 'Params',
			-values => ['RNA','oldRNA','DNA'],
			-default => 'RNA',
			-labels => {'RNA', 'use RNA parameters',
				    'oldRNA', 'old RNA parameters',
				    'DNA', 'use DNA parameters'}
			);
   
   print "<H4>Options to modify the fold algorithm</H4>\n";
   print "Weight of covariance term (0..10)\n";
   print $q->textfield(-name=>'cv',
                           -default=>1.0,
                           -size=>3);

   print "&nbsp;Penalty for non-compatible sequences (0..10)\n";
   print $q->textfield(-name=>'nc',
                           -default=>1.0,
                           -size=>3);

   print "Rescale energy parameters to temperature\n";
   print $q->textfield(-name=>'Temp',
		       -default=>37,
		       -size=>3);
   print "C<P>\n";
   print $q->checkbox_group(
			    -name=>'toggles',
			    -values=>['-4',  '-noGU', '-noCloseGU',
				      '-noLP'],
			    -labels=>{'-4', 'no special tetraloops',
				      '-noGU', 'no GU pairs',
				      '-noCloseGU',
				      'no GU pairs at the end of helices',
				      '-noLP', 'avoid isolated base pairs'
				       },
			    -linebreak=>'yes',
			    -defaults=>[]);
   
   print "<P>\n",$q->reset;
   print $q->submit('Action','Fold it');
   print $q->endform;

   if (!$q->param) {
      my @old_files = $q->cookie('alifold_result');
      if (@old_files) {
	 print "<P>The results from your last run should still be in\n";
	 foreach my $f (@old_files) {
	    print " <A href=\"$hdir/$f\">$f</a>";
	 }
	 print ".\n";
      }
   }
   print "<HR>\n";
   print_tail();
}

sub do_work {
   setpriority(0, 0, 19);
   my $q  = shift;
   my $name = "alifold$$";
   my $OK_CHARS ='-a-zA-Z0-9_.@ ';  
   my($RNAfold_id, $junk);
   my $WORK_DIR ="/u/www" . $hdir . "/$name";
   # my $WORK_DIR = $name;
   my $sfact =1.02;
   mkdir $WORK_DIR;
   chdir $WORK_DIR || die("can't change to directory $WORK_DIR");
   
   # clean old files
   foreach my $f (<../alifold*>) {
      print STDERR "removing $f\n";
      system('/bin/rm', '-r', $f) if ((-M $f)>1.2);
   }
   
   my $clustal_file = $q->param('clustal_file');
   my $l=0;
   open(ALN, ">alifold.aln");
   while (<$clustal_file>) {
      print ALN $_;
      $l += length;
   }
   close(ALN); close($clustal_file);
   
   my $options  = '';
   $options .= ' -p' if ($q->param('pffold') eq 'pf');
   $options .= ' -P dna.par'
       if ($q->param('Params') eq 'DNA');
   $options .= ' -P vienna13.par'
       if ($q->param('Params') eq 'oldRNA');
   $options .= " " . join(' ', $q->param('toggles'));
   $RNA::tetra_loop = 0 if ($options =~ /-4/);
   $RNA::dangles = 0 if ($options =~ /-d/);
   $RNA::noGU = 1 if ($options =~ /-noGU/);
   $RNA::no_closingGU = 1 if ($options =~ /-noCloseGU/);
   $RNA::noLonelyPairs = 1 if ($options =~ /-noLP/);
   
   if ($q->param('Temp') =~ /([-+]?\d+\.?\d*)/) {
      $RNA::temperature = $1;
      $options .= " -T $1" if ($1 != 37);
   }

   if ($q->param('cv')) {
      my $cv = $q->param('cv')*1.0;
      $cv = 0 if $cv <0;
      $cv = 10 if $cv>10;
      $options .= " -cv $cv" if $cv != 1;
   }

   if ($q->param('nc')) {
      my $nc = $q->param('nc')*1.0;
      $nc = 0 if $nc <0;
      $nc = 10 if $nc>10;
      $options .= " -nc $nc" if $nc != 1;
   }

   
   # remove strange characters, so we can pass strings to shell
   $options =~ s/[^$OK_CHARS]//go;
   

   my $the_cookie = $q->cookie(-name=>'alifold_result',
			       -value=>$name,
			       -expires=>'+1d',
			       -path=>'/cgi-bin');
   my $time = int($l*$l/1000000);   # crude estimate
   $time += 3;
   print $q->header(-refresh=>"$time; URL=$hdir/$name",
		    -expires=>'+1m',
		    -cookie=>$the_cookie,
		    -type=>'text/html');
   print $q->start_html(-title=>"Alifold in progress",
			-author=>'ivo@tbi.univie.ac.at',
			-BGCOLOR=>'#f8f8ff');

   print "<H2>Now folding... </H2>\n";

   print "if all goes well you should be automatically forwarded to your",
   "results page after $time seconds.<br>\n",
   "Otherwise, click <a href=\"$hdir/$name/\">here</a> to view your ",
   "results after a few seconds<p>\n";

   print_tail();

   open(OUT, ">index.html");
   my $old_stdout = select(OUT);
   print $q->start_html(-title=>'Alifold results',
		       -author=>'ivo@tbi.univie.ac.at',
		       -BGCOLOR=>'#f8f8ff');

   print 
       "<H1>Your Alifold Results</H1>\n",
       "For more information on the alifold method see ???<br>\n",
       "If the output seems incomplete press relaod after a few seconds<p>\n";
   
   if ($q->param('Params') eq 'DNA') {
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
	  "D.H.  Mathews, J. Sabina, M. Zucker and H. Turner\n",
	  "\"Expanded Sequence Dependence of Thermodynamic Parameters ",
	  "Provides Robust Prediction of RNA Secondary Structure\",\n",
	  "JMB, 288, pp 911-940, 1999<hr>\n";
   }
   
#   print $q->query_string, "<br>\n";
#   my $length = length($Sequence);
#   if ($length>$maxlength) {
#      print "Sorry, this interface is still being tested.\n",
#      "Currently, sequences longer than $maxlength bases will currently ",
#      "not be processed. Your sequence has length $length.<HR>\n";
#      return;
#   }
   
   print "An equivalent RNAalifold command line would have been<BR>\n";
   print "<kbd>RNAalifold $options alifold.aln</kbd><P>\n";
   $| = 1;
   print "<P><PRE>\n";
#   print STDERR "/home/blini/ivo/RNA/ViennaRNA/AliFold/RNAalifold $options alifold.aln\n";
   print `/home/blini/ivo/RNA/ViennaRNA/AliFold/RNAalifold $options alifold.aln 2>&1`;
   print "</PRE><P>\n";

   print "Other output files:<dl>\n",
   "<dd><a href=\"alirna.ps\">postscript drawing</a> ",
   "of the predicted structure<br>\n";
   if ($options =~ /-p/) {
      print "<dd><a href=\"alidot.ps\">postscript dot plot</a> ",
      "of pair probabilities<br>\n";
      
      print '<dd><a href="alifold.out">text output</a> with detailed ',
      "information about each base pair<br>\n";
      `/home/blini/ivo/RNA/ViennaRNA/AliFold/cmt.pl alifold.out > cmount.eps`;
      print "<dd>colored <a href=\"cmount.eps\">mountain plot</a> in postscript";
   }
   print "</dl>\n";
   
   my ($user,$system,$cuser,$csystem) = times;
   print "<BR>Time used for this call ",$user+$system, " seconds<BR>\n";
   print "Results will remain on this server for approximatly 1 day\n<HR>\n";
   print_tail();
   select($old_stdout);
}

sub print_tail {
   print <<END;
   We appreciate your feedback. Please send comments to    
   <ADDRESS>Ivo Hofacker
   <A HREF="mailto:ivo\@tbi.univie.ac.at">&lt;ivo\@tbi.univie.ac.at&gt;</A>
   </ADDRESS><BR>
   <A HREF="http://www.tbi.univie.ac.at/~ivo/RNA/">Vienna RNA Home Page</A>
END
    print $q->end_html,"\n";
}
