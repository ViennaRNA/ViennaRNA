#!/usr/bin/perl  -T
# -*-CPerl-*-
# Last changed Time-stamp: <2005-09-06 17:13:49 ivo>
# CGI script for a Web-based RNA fold server
# you need to have the perl5 RNA module installed
# that comes as part of the Vienna RNA package

use strict;
use CGI qw/-no_xhtml :standard /;
use CGI::Carp;# qw(fatalsToBrowser);
use HTML::Entities;
use URI::Escape;
use Email::Valid;
use warnings;
no warnings qw(uninitialized);
#use RNA;
#use Chart::Lines;

use vars qw/$hdir $maxlength $maxlength1 $help_url $RNAalifold $ServerRoot $RNAdir $batchscript/;

# please configure these variables
$ServerRoot = '/u/www';
$hdir = '/RNAfold_dir';      # were the output files are stored
$RNAdir = '/var/www/RNA/';   # where the RNAalifold executable lives
$help_url = 'http://www.tbi.univie.ac.at/~ivo/RNA/alifoldcgi.html';
$CGI::POST_MAX = 12*1024;    # maximum filsize for the alignment
$maxlength  = 2000;          # only process sequences up to this length
$maxlength1 = 300;           # limit for immediate jobs
$batchscript = '/var/www/RNA/ALIbatch_new.pl'; # script for batch submissions

$ENV{PATH} = '/bin:/usr/bin:/usr/local/bin';
delete @ENV{'IFS', 'CDPATH', 'ENV', 'BASH_ENV'};

my $q = new CGI;

if ($q->cgi_error) {
  print $q->header(-status=>$q->cgi_error),
    $q->start_html('Problems'),
      $q->h2('Request not processed'),
	$q->strong($q->cgi_error);
  exit 0;
}

if ($q->param &&  $q->param('Action') =~ /Fold/) {
  do_work($q);
} else {
  print_form($q);
}


sub print_form {
  my($q, $cookie) = @_;
  if ($cookie) {
    print $q->header(-cookie=>$cookie);
  } else {
    print $q->header();
  }
  print $q->start_html(-title=>"RNAalifold input form",
		       -author=>'ivo@tbi.univie.ac.at',
		       -BGCOLOR=>'#f8f8ff');
  print "\n<H1 align=center>Vienna RNA Secondary Structure Prediction</H1>\n",
    "<H2 align=center><a name=top>",
    "A Web Interface for the Prediction Consensus Structures of Aligned Sequences</a></H2>\n";

  print "This server will predict the consensus secondary structures for a\n",
    "set of aligned single stranded RNA or DNA sequences.\n",
      "If the options look confusing <strong>read the\n",
	"<a href=\"http://www.tbi.univie.ac.at/~ivo/RNA/alifoldcgi.html\">",
	  "help page</a></strong><p>\n";

  print "<b><u>News</u></b>You can now upload larger alignments for " .
    "batch processing.\n Maximum upload size is now 10Kbyte, maximum " .
      "alignment length is $maxlength.<p>\n";

  print $q->start_multipart_form();

  if ($q->param('aln')) {
    my @aln = split(' ',$q->param('aln'));
    print "Alignment from cgi arguments is:<p><PRE>";
    print join "\n", @aln;
    print "</PRE>\n";
    print $q->hidden(-name=>'aln', -default=>$q->param('aln'));
  } else {
    print "<P><b>Enter</b> the file containing your <b>aligned sequences</b> ",
      "for uploading.<br>\n",
      "The alignment file <b>must</b> be in <kbd>Clustalw</kbd> format.<br>\n",
      "Maximum file size is ", int($CGI::POST_MAX/1024 -1.1), "KB.<BR>\n";

    print $q->filefield(-name=>'clustal_file',
			-default=>'starting value',
			-size=>50,
			-maxlength=>80);
  }
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
				     '-noLP', 'RNAz'],
			   -labels=>{'-4', 'no special tetraloops',
				     '-noGU', 'no GU pairs',
				     '-noCloseGU',
				     'no GU pairs at the end of helices',
				     '-noLP', 'avoid isolated base pairs',
				     'RNAz', 'run RNAz on alignment'
				    },
			   -linebreak=>'yes',
			   -defaults=>['-noLP', 'RNAz']);

  print "\n<br><small>",
    "The <a href='http://www.tbi.unive.ac.at/~wash/RNAz/>RNAz</a> ",
    "program decides whether an alignment is likely to contain a ",
      "functional RNA structure</small><p>\n";
  print "Email address. For batch jobs (over $maxlength1) this ",
    "is mandatory, so we can notify you when the job has completed.\n";
  print $q->textfield(-name => 'email',
		      -default => 'you@where.org',
		      -size => 30);
  print "<P>\n",$q->reset;
  print $q->submit('Action','Fold it');
  print $q->endform;
  
  if (!$q->param) {
    my $last = $q->cookie('alifold_result');
    if ($last) {
      print "<P>The results from your last run should still be in\n";
      print " <A href=\"$hdir/$last/\">$last</a>";
    }
    print ".\n";
  }
  print "<HR>\n";
  print_tail();
}

sub do_work {
  setpriority(0, 0, 19);
  my $q  = shift;
  my $name = "alifold$$";
  my $OK_CHARS ='-a-zA-Z0-9_.@ ';
  my $WORK_DIR = $ServerRoot . $hdir . "/$name";
  # my $WORK_DIR = $name;
  my $sfact =1.02;
  mkdir $WORK_DIR or warn "can't create directory $WORK_DIR";
  chdir $WORK_DIR or die "can't change to directory $WORK_DIR";

  open(STATE, ">cgistate.txt");
  $q->save(\*STATE);
  close(STATE);

  # clean old files
  foreach my $f (<../alifold*>) {
    if ($f =~ /(\.\.\/alifold\d+)/) {
      system('/bin/rm', '-r', $1) if (-d $1) && (-M $1)>3;
    }
  }

  open(ALN, ">alifold.aln") or croak("can't write alignment");
  my @lines;
  if ($q->param('aln')) {
    my @aln = split(' ', $q->param('aln'));
    print ALN "CLUSTAL W alignment from pmmatch\n\n";
    my $sn=1;
    foreach (@aln) {
      push @lines, sprintf "%-10s %s", $sn++, $_;
    }
    push @lines, '          ';
  } else {
    my $clustal_fh = $q->upload('clustal_file');
    error_page(0) unless defined($clustal_fh);
    undef($/); # slurp mode
    $_ = <$clustal_fh>; error_page('firstline') unless /^CLUSTAL/;
    close($clustal_fh);
    @lines =  split /\015\012?|\012/;
    print ALN shift(@lines), "\n";
  }
  my ($size, $num_seq, $length) = (0,0,0);
  my %l = ();
  foreach (@lines) {
    print ALN $_, "\n";
    if (/^\S+\s+\S+/) {
      my($name, $sseq) = split;
      $l{$name} += length $sseq;
    }
    $size += length;
  }
  close(ALN); 
  error_page($size) if $size<14;
  $num_seq = scalar(keys %l);
  error_page('numseq') if $num_seq<2;
  $length = $l{(keys %l)[0]};
  foreach (keys %l) {
    error_page('unequal') if $length != $l{$_};
  }
  error_page("length $length") if $length>$maxlength; 
  my $options  = '';
  $options .= ' -p' if ($q->param('pffold') eq 'pf');
  $options .= ' -P dna.par'
    if ($q->param('Params') eq 'DNA');
  $options .= ' -P vienna13.par'
    if ($q->param('Params') eq 'oldRNA');
  $options .= " " . join(' ', $q->param('toggles'));

  if ($q->param('Temp') =~ /([-+]?\d+\.?\d*)/) {
     my $T = ($1>-273.15 && $1<1000) ? $1 : 37; # restrict interval
     $options .= " -T $T" if $T != 37;
   }

  if ($q->param('cv') =~ /(\d+.?\d*)/) {
    my $cv = $1;
    $cv = 0 if $cv <0;
    $cv = 10 if $cv>10;
    $options .= " -cv $cv" if $cv != 1;
  }

  if ($q->param('nc') =~ /(\d+.?\d*)/) {
    my $nc = $1;
    $nc = 0 if $nc <0;
    $nc = 10 if $nc>10;
    $options .= " -nc $nc" if $nc != 1;
  }

  # remove strange characters, so we can pass strings to shell
  $options =~ /([$OK_CHARS]+)/; $options= $1; # $options now untainted

  my $the_cookie = $q->cookie(-name=>'alifold_result',
			      -value=>$name,
			      -expires=>'+1d',
			      -path=>'/cgi-bin');
  my $time = estimate($num_seq,$length,$options);
  print $q->header(-refresh=>"$time; URL=$hdir/$name/",
		   -status=>'202 Accepted',
		   -cookie=>$the_cookie,
		   -type=>'text/html');
  print $q->start_html(-title=>"Alifold in progress",
		       -author=>'ivo@tbi.univie.ac.at',
		       -BGCOLOR=>'#f8f8ff');

  print "<H2>Now folding... </H2>\n";

  my $RNAz = 'false';
  if ($options =~ s/RNAz//) {
    if (($length<=400) && ($num_seq<=6)) {
      $RNAz = 'true';
    } else {
      print "Warning: RNAz will not be run, since your alignment has too many sequences or is too long.<br>",
	"RNAz is trained only for alignments up to length 400 with 2 to 6 sequences.<p>\n";
    }
  }


  if ($length<$maxlength1) {
    print "if all goes well you should be automatically forwarded to your ",
      "results page after $time seconds.<br>\n",
	"Otherwise, click <a href=\"$hdir/$name/\">here</a> to view your ",
	  "results after a few seconds<p>\n";
  } else {
    my $email;
    $email = Email::Valid->address($q->param('email'), -mxcheck => 1);
    if (!defined($email) || $email eq 'you@where.org') {
      print "Sorry, the address ",
	encode_entities($q->param('email')),
	  " does not appear to be valid.<br>\n",
	    "For long jobs we require an email address so we can notify ",
	      "you when the calculation is finished<p>\n";
      print_tail();
      return;
    }

    my $args = $options;
    #$args .= ' alifold.aln';
    open(QU, '|-', '/usr/bin/qsub -q RNA -s /bin/sh') or
      die "can't submit batch job: $!";
    print QU "cd $WORK_DIR\n $batchscript $args $RNAz $email\n";
    close QU;
    print "<p>Your job has been submitted to our queuing system.<br>\n",
      "Assuming no other load on the system, expect results to be available ",
	"in about ";
    printf "%.2f %s", $time/60, " minutes.<br>\n";
    print "You should be notified by mail when the computation is complete. ",
      "<br>You can then follow the link to ",
	"<a href=\"$hdir/$name/\">your results page</a><p>\n";
  }
  print_tail();

  open(OUT, ">index.html");
  my $old_stdout = select(OUT);
  print $q->start_html(-title=>'Alifold results',
		       -author=>'ivo@tbi.univie.ac.at',
		       -BGCOLOR=>'#f8f8ff');

  print "<H1>Your Alifold Results</H1>\n",
    'For detailed information on the alifold method see <a href=',
      '"http://www.tbi.univie.ac.at/papers/Abstracts/01-11-067abs.html">',
	"this manuscript</a><br>\n",
	  "If the output seems incomplete press relaod after a few seconds<p>\n";

  if ($q->param('Params') eq 'DNA') {
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

  if ($length>=$maxlength1) {
    print "[If this file is incomplete (because the computation ",
      "isn't finished) reload this page after a while]<p>\n";
    select($old_stdout);
    close(OUT);
    return;
  }

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
  print "</dl>\n";

  my ($user,$system,$cuser,$csystem) = times;
  print "<BR>Time used for this call ",$user+$system, "+", $cuser+$csystem,
    " seconds<BR>\n";
  print "Results will remain on this server for approximatly 1 day\n<HR>\n";
  print_tail();
  select($old_stdout);
}

sub print_tail {
  print <<END;
   We appreciate your feedback. Please send comments to
   <ADDRESS>Ivo Hofacker
   <A HREF="mailto:rna\@tbi.univie.ac.at">&lt;rna\@tbi.univie.ac.at&gt;</A>
   </ADDRESS>
   <A HREF="http://www.tbi.univie.ac.at/~ivo/RNA/">Vienna RNA Home Page</A>
END
  print $q->end_html,"\n";
}

sub error_page {
  local $_ = shift;
  my $error = $q->cgi_error() ? $q->cgi_error() : '444 Unusable upload';

  print $q->header(-status=>$error);
  print $q->start_html(-title=>"Alifold error",
		       -author=>'ivo@tbi.univie.ac.at',
		       -BGCOLOR=>'#f8f8ff');

  print "<H2>Alifold error </H2>\n";

  if (/^\d+/ && $_<14) {
    print "The alignment file you uploaded was empty.<br>\n",
      "Alifold needs a multiple sequence alignment in CLUSTAL",
	"format as input.<p>\n";
  }
  if (/length (\d+)/) {
    print "Sorry, in order to provide a swift service we have to limit ",
      "both the maximum file size and alignment length.<br>\n",
	"At the moment we can only accept alignments up to ",
	  "$maxlength in length. Your alignment has $_.<br>\n";
  } else {
    print "First line of the alignment should start with CLUSTAL<br>\n"
      if /firstline/;
    print "There should be at least two sequences in your alignment<br>\n"
      if /numseq/;
    print "Sequences in the alignment should be equal length<br>\n" 
      if /unequal/;
    print <<END;
The input file you uploaded does not look like a sequence alignment
in CLUSTAL format. Perhaps you uploaded the wrong file, or your
alignment is in some other format. In the latter case you can probably
load it into clustalx and re-export it in CLUSTAL format.<p>
END
  }

  print "See also the <a href=\"$help_url\">alifold help page</a>.<p>\n",
    "Use your browsers 'Back' button to try again.<p>\n";

  print_tail();
  exit;
}

sub estimate {
  # rough estimate for folding time from seq length
  # adjust these parameters for your server machine!
  my ($N, $l, $pf) = @_;
  my $t = 1 + $N * 8e-6 * $l*$l + 9e-9 * $l*$l*$l;
#  $t *= $N;
  $t = 4*$t + $N * 9e-9 * $l*$l*$l if ($pf eq 'pf');
  $t += 100 if $l > $maxlength1; # extra time for batch submission
  return int($t+1);
}
