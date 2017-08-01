#!/usr/bin/perl
#

use Getopt::Long;

my $input;
my $output;

GetOptions(
  "i=s" => \$input,
  "o=s" => \$output
);

sub readFile {
  my $file = shift;
  my $html = "";

  local $/=undef;

  open(HEAD, "<$file");
  $html .= <HEAD>;
  close(HEAD);

  return $html;
}


if(defined($input)){
  my $html = readFile($input);

  $html =~ s/sec(\d+)\.(\d+)/sec$1_$2/gsm;

  # store navigation
  my $navigation = $1 if $html =~ /<\/a>Contents<\/h3>\s*(<ul(.*?)<\/ul>)\s*<p>/sm;
  # remove navigation
  $html =~ s/<\/a>Contents<\/h3>\s*(<ul(.*?)<\/ul>)\s*<p>/<\/a>Contents<\/h3>\n<p>\n/sm;

  # remove TOC section
  $html =~ s/<h3 class="likesectionHead">(.*?)Contents<\/h3>//sm;

  # reinsert navigation at correct position
  $html =~ s/__SIDENAVLIST__/$navigation/sm;

  # remove stupid senseless empty tags
  while($html =~ s/<p>\s*<\/p>//sm){;}

  # remove stupid senseless empty tags
  while($html =~ s/\s+<\/pre>/<\/pre>/sm){;}

  print $html;
}
