use ExtUtils::MakeMaker qw(prompt WriteMakefile);

# =====> PATHS: CHECK AND ADJUST <=====
my @INC     = qw(-I../src  -I../src/win32 -I../src/PS -I../src/GD); 
my @LIBPATH = qw(-L../Debug -L../../gd-1.8.4);
my @LIBS    = qw(-lg2);

# FEATURE FLAGS
warn "\nPlease choose the features that match how g2 was built:\n";

my $PS    	= lc prompt('Build PostScript support?','y') eq 'y';
my $GD    	= lc prompt('Build gd (Bitmap) support?','y') eq 'y';
my $WIN32   = lc prompt('Build Win32 support?','y') eq 'y';
my $X11     = lc prompt('Build X11 support?','y') eq 'y';

warn "\nIf you experience compile problems, please check the \@INC, \@LIBPATH and \@LIBS\n",
     "arrays defined in Makefile.PL and manually adjust, if necessary.\n\n";

#### no user-serviceable parts below #####

push @LIBS,'-lgd' if $GD;
push @LIBS, '-lm' unless $^O eq 'MSWin32';

# FreeBSD 3.3 with libgd built from ports croaks if -lXpm is specified 
if ($^O ne 'freebsd' && $^O ne 'MSWin32') { 
  push @LIBS,'-lX11','-lXpm' if $XPM; 
} 

my $CAPI = defined $ExtUtils::MakeMaker::CAPI_support ? 'TRUE' : 'FALSE';
my $DEFINES = '';
$DEFINES    .= ' -DDO_PS'     if $PS;
$DEFINES    .= ' -DDO_GD'     if $GD;
$DEFINES    .= ' -DDO_WIN32'  if $WIN32;
$DEFINES    .= ' -DDO_X11'    if $X11;

WriteMakefile(
    'NAME'	     => 'g2',
    'VERSION_FROM'	=> 'g2.pm',
    'dist'      => {'COMPRESS'=>'gzip -9f', 'SUFFIX' => 'gz',
	            'ZIP'=>'/usr/bin/zip','ZIPFLAGS'=>'-rl'},
    'OBJECT'    => 'g2_wrap.o',
    'OPTIMIZE'  => '-g',
    'LIBS'      => [join(' ',$ENV{'G2_LIBS'},@LIBPATH,@LIBS)],
    'INC'       => join(' ',$ENV{'G2_INC'},@INC), 
    'AUTHOR'    => 'Horst Wagner',
    'ABSTRACT'  => 'Perl interface to g2 Graphics Library',
    'CAPI'      => $CAPI,
    'DEFINE'	=> $DEFINES,
);

sub MY::postamble {
    my $postamble = <<'END';
html: g2.pm
	pod2html --outfile=g2.html g2.pm
END
    $postamble;
}

