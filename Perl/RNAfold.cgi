#!/usr/local/bin/perl 
# -*-Perl-*-
# Last changed Time-stamp: <1998-05-19 13:04:46 ivo>
# CGI script for a Web-based RNA fold server
# you need to have the perl5 RNA module installed
# that comes as part of the Vienna RNA package

use CGI;
use CGI::Carp qw(fatalsToBrowser);
use HTML::Entities;
use URI::Escape;
use RNA;
use Chart::Lines;

# please configure these variables
$hdir = "/RNAfold_dir";      # were the output files are stored
$maxlength = 300;            # only process sequences up to this length
$ssv_home = "http://smi-web.stanford.edu/projects/helix/sstructview/home.html";
$ssv_url  = "/~ivo/RNA/SStructView.zip";


$query = new CGI;

if ($query->param) {
    &do_work($query);
} else {
    &print_header($query);
}

&print_tail;
print $query->end_html;

sub print_header {
    my($query) = @_;
    if ($the_cookie) {
	print $query->header(-expires=>'+1m',
			     -cookie=>$the_cookie);
#	print "setting cookie $the_cookie\n";
    } else {
	print $query->header(-expires=>'+1m');
#	print "not setting cookie\n";
    }
    print $query->start_html(-title=>"RNAfold input form",
			     -author=>'ivo@tbi.univie.ac.at',
			     -BGCOLOR=>'#f8f8ff');
    print "\n<H1 align=center>Vienna RNA Secondary Structure Prediction</H1>\n";
    print "<H2 align=center><a name=top>",
    "A web interface to the RNAfold programm</a></H2>\n";
#    print "<H3 align=center>Experimental version, not for heavy use</H3>\n";
    print "If this looks confusing <strong>read the ",
    "<a href=\"http://www.tbi.univie.ac.at/~ivo/RNA/RNAcgi.html\">",
    "help page</a></strong><p>\n";
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
			   -columns=>60,
			   -default=>"");

    print "<P>\n<strong>Choose Fold Algorithm</strong><BR>\n";

    print $query->popup_menu(-name => 'pffold',
			     -values => ['mfe','pf'],
			     -default => 'pf',
			     -labels => {'mfe', 'minimum free energy only',
					 'pf',
					 'partition function and pair probabilities'}
			     );
    

    print "<H4>Options to modify the fold algorithm</H4>\n";
    print "Rescale energy parameters to temperature\n";
    print $query->textfield(-name=>'Temp',
			    -default=>37,
			    -size=>3);
    print "C<P>\n";
    print $query->checkbox_group(
				 -name=>'toggles',
				 -values=>['-4', '-d',  '-noGU', '-noCloseGU'],
				 -labels=>{'-4', 'no special tetraloops',
					   '-d',
					   'no dangling end energies',
					   '-noGU', 'no GU pairs',
					   '-noCloseGU',
					   'no GU pairs at the end of helices'
					   },
				 -linebreak=>'yes',
				 -defaults=>[]);

    print "\n<BR><strong>",
    "Should we produce a mountain plot of the structure?</strong> \n";
    print $query->checkbox(-name=>'plot', -checked=>'ON');
    print "\n<BR>View the mfe structure in the SStructView java applet?\n";
    print $query->checkbox(-name=>'SSview');

    print "<P>\n",$query->reset;
    print $query->submit('Action','Fold it');
    print $query->endform;

    if (!$query->param) {
	@old_files = $query->cookie('old_files');
	if (@old_files) {
	    print "<P>The results from your last run should still be in\n";
	    foreach $f (@old_files) {
		print " <A href=\"$hdir/$f\">$f</a>";
	    }
	    print ".\n";
	}
    }
    print "<HR>\n";
}

sub do_work {
    my($query) = @_;
    my($OK_CHARS)='-a-zA-Z0-9_.@ ';  
    my($RNAfold_id, $junk);
    my($WORK_DIR)="/usr/local/etc/httpd/htdocs" . $hdir;
    my($sfact)=1.02;
    chdir $WORK_DIR || die("can't change directory");

    # clean old files
    foreach $f (<*.ps *.gif *.coords>) {
	unlink($f) if ((-M $f)>1.2);
    }
    
    $Sequence = $query->param('Sequence');
    if (length($Sequence)<1) {   # do nothing
	&print_header($query);
	return;
    }
    $Sequence =~ s/[\n\s]//gm;   # remove whitespace and linebreaks
    $Sequence =~ s/[tT]/U/g;     # T -> U
    $Sequence = uc $Sequence;    # upper case

    $name = $query->param('name');

    $options  = "";
    $options .= " -p" if ($query->param('pffold') eq 'pf');
    $options .= " " . join(' ', $query->param('toggles'));
    $RNA::tetra_loop = 0 if ($options =~ /-4/);
    $RNA::dangles = 0 if ($options =~ /-d/);
    $RNA::noGU = 1 if ($options =~ /-noGU/);
    $RNA::no_closingGU = 1 if ($options =~ /-noCloseGU/);

    if ($query->param('Temp') =~ /([-+]?\d+\.?\d*)/) {
	$RNA::temperature = $1;
	$options .= " -T $1";
    }
	
    # remove strange characters, so we can pass strings to shell
    $options =~ s/[^$OK_CHARS]//go; 
    if ($name =~ /(\S+)/) {
	$name = $1;
    } else {
	$name = "s$$";
    }

    $name =~ s/\///g;    # no / allowed in file names
    $fname_ss = $name . "_ss.ps";
    $fname_dp = ($options =~ /-p /)?($name . "_dp.ps"):"";
    $fname_co = $name . ".coords";
    $ss_escaped = uri_escape($fname_ss);
    $dp_escaped = uri_escape($fname_dp);
    $gif_escaped = uri_escape($name . ".gif")
	if ($query->param('plot') eq 'on');
    
    $the_cookie = $query->cookie(-name=>'old_files',
				 -value=>[$ss_escaped,$dp_escaped,$gif_escaped],
				 -expires=>'+1d',
				 -path=>'/cgi-bin/RNAfold.cgi');

    &print_header($query);   # don't print anything before this line
	
    print "<H2><a name=\"Results\">Here are your results</a></H2>\n";
#   print $query->query_string, "<br>\n";
    $length = length($Sequence);
    if ($length>$maxlength) {
	print "Sorry, this interface is still being tested.\n",
	"Currently, sequences longer than $maxlength bases will currently ",
	"not be processed. Your sequence has length $length.<HR>\n";
	return;
    }
    
    print "An equivalent RNAfold command line would have been<BR><kbd>\n";
    print "RNAfold $options</kbd><P>\n";
    
    $structure = $Sequence;  # wierd way of allocating space
    $mfe = RNA::fold($Sequence, $structure);
    print "The optimal secondary structure in bracket notation is given below\n";
    $eseq = $Sequence;
    encode_entities($eseq);
    
    print "<P><PRE>\n";
    print "&gt; $name\n$eseq\n$structure";
    printf " (%6.2f)\n", $mfe;
    print "</PRE><P>\n";
    print "The minimum free energy in kcal/mol is given in parenthesis.\n"; 
#    print `pwd`," $fname_ss";
    if (!($structure =~ /\(/)) {$RNA::rna_plot_type=0};
  RNA::PS_rna_plot($Sequence, $structure, $fname_ss);
    print "You may look at the PostScript drawing of the structure\n";
    $efname = $fname_ss; encode_entities($efname);
    print "in <A href=\"$hdir/$ss_escaped\">$efname</a>.<P>\n";
    
    if ($options =~ /-p/) {
	$RNA::dangles=2 if ($RNA::dangles);
	$min_en = RNA::energy_of_struct($Sequence, $structure);
	$kT = ($RNA::temperature+273.15)*1.98717/1000.; # in Kcal 
	$RNA::pf_scale = exp(-($sfact*$min_en)/$kT/$length);
	$pf_struct = $structure;
	$energy = RNA::pf_fold($Sequence, $pf_struct);
	print "The free energy of the thermodynamic ensemble is ";
	printf "<b>%6.2f</b> kcal/mol\n", $energy;
      RNA::PS_dot_plot($Sequence, $fname_dp);
	print "<BR>The PostScript <em>dot plot</em> containing the base pair ",
	"probabilities is in\n<A href=\"$hdir/$dp_escaped\">",
	encode_entities($fname_dp), "</a>.\n";
    }
    print "<P>";
    
    mountain($name, $length, $options) if ($query->param('plot') eq 'on');

    if ($query->param('SSview') eq 'on') {
	$ret = RNA::ssv_rna_plot($Sequence, $structure, $fname_co);
	if ($ret == 0) {
	    print "error: ssv_rna_plot plot failed, sorry<HR>\n";
	    return;
	}
	$co_escaped = uri_escape($fname_co);
	print "Please be patient while the applet is loading<br>\n";
	print "<APPLET ARCHIVE=\"$ssv_url\" CODE=\"SStructView.class\" ",
	"name=\"SSApplet\" WIDTH=500 HEIGHT=450>\n",
	"<PARAM name=\"structure-data-URL\" ",
	"value=\"http://www.tbi.univie.ac.at$hdir/$co_escaped\">\n",
	"<PARAM name=\"show-controls\" value=\"true\">\n",
	"<PARAM name=\"result-display-frame\" value=\"_blank\">\n",
	"</APPLET><HR>\n";
	print "For information on the java applet see the ",
	"<A href=\"$ssv_home\">SStructView</a> home page<p>\n";
    }

    print "Your output files will be deleted from this server after one day.\n";

#    $myself = $query->self_url;
#    print "<BR><A HREF=$myself#top>Back to the input form</A>";
    print "<BR>Scroll back to the top to submit another sequence.\n";

    ($user,$system,$cuser,$csystem) = times;
    print "<BR>Time used for this call ",$user+$system, " seconds\n"; 
    print "<HR>\n";
    
}

sub print_tail {
    print <<END;
We appreciate your feedback. Please send comments to    
<ADDRESS>Ivo Hofacker
<A HREF="mailto:ivo\@tbi.univie.ac.at">&lt;ivo\@tbi.univie.ac.at&gt;</A>
</ADDRESS><BR>
<A HREF="/~ivo/RNA/">Vienna RNA Home Page</A>
END
}

sub mountain {
    my($name, $length, $options) = @_;
    $width =  480;
    $height = 300;

    $skip = 10**(int (log($length)/log(10.) - 0.5));
    $obj = Chart::Lines->new( $width, $height );
    $obj->set ('title' => $name,
	       'x_label' => 'Position',
	       'y_label' => 'Height',
	       'legend' => 'false',
	       'skip_x_ticks' => $skip);
    
    $obj->add_dataset ((0..$length));
    $mp[0]=$mp[1] = 0;  # use indices [1..length]

    $bp = $RNA::base_pair;
    $npair = RNA::bond_i_get($bp);
    @mp = ();
    for (1..$npair) {
	$bp = RNA::ptradd($bp, 8); # sizeof(struct bond) == 8
	$mp[RNA::bond_i_get($bp)+1] =  1;
	$mp[RNA::bond_j_get($bp)  ] = -1;
    }
    for ($i=1; $i<=$length; $i++) {
	$mp[$i] += $mp[$i-1];
    }
    $obj->add_dataset (@mp[0..$length]);
    
    if ($options =~ /-p/) { # we have base pairing probs
	$bpp = RNA::Make_bp_profile($length);
	$mp[0]=$mp[1] = 0;  # use indices [1..length]
	for ($i=1; $i<=$length; $i++) {
	    $pp = RNA::deref_any($bpp, $i); 
	    $mp[$i] -= RNA::ptrvalue($pp, 2, 'float'); # paired downstream
	    $mp[$i+1] = $mp[$i] + RNA::ptrvalue($pp, 1, 'float'); # upstream
	}
	
	$obj->add_dataset (@mp[0..$length]);
    }

    $name .= ".gif";
    $obj->gif ("$name");
    print "<IMG ALT=\"Mountain Representation\" HEIGHT=\"$height\" ",
    "WIDTH=\"$width\" SRC=\"$hdir/$gif_escaped\"><P>\n";
}
