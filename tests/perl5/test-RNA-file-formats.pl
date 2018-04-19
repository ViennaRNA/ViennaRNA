#!/usr/bin/perl -I.libs

######################### We start with some black magic to print on failure.
# (It may become useful if the test is moved to ./t subdirectory.)

use strict;
use warnings;
use Test::More tests => 28;

use RNA;
use RNAHelpers qw(:Messages :Paths);


my $datadir = getDataDirPath();
my $msa_format;
my $n_seq;
my $sequence_identifiers;
my $alignment;
my $alignment_id;
my $consensus_structure;
my $fh;
my $counter;


#
# test detecting MSA formats
#
MsgChecking("whether we are able to detect STOCKHOLM 1.0 alignment format");
$msa_format = RNA::file_msa_detect_format($datadir . "/rfam_seed_selected.stk");
is($msa_format, RNA::FILE_FORMAT_MSA_STOCKHOLM);

MsgChecking("whether we are able to detect CLUSTAL alignment format");
$msa_format = RNA::file_msa_detect_format($datadir . "/070313_ecoli_cdiff_16S_clustalw.aln");
is($msa_format, RNA::FILE_FORMAT_MSA_CLUSTAL);

MsgChecking("whether we are able to detect FASTA alignment format");
$msa_format = RNA::file_msa_detect_format($datadir . "/070313_ecoli_cdiff_16S_fasta.aln");
is($msa_format, RNA::FILE_FORMAT_MSA_FASTA);

MsgChecking("whether we are able to detect MAF alignment format");
$msa_format = RNA::file_msa_detect_format($datadir . "/test.maf");
is($msa_format, RNA::FILE_FORMAT_MSA_MAF);

MsgChecking("whether we are able to reject non-alignment file format");
$msa_format = RNA::file_msa_detect_format($datadir . "/rnafold.cmds");
is($msa_format, RNA::FILE_FORMAT_MSA_UNKNOWN);

#
# test reading MSA from different file types
#
MsgChecking("whether we are able to parse STOCKHOLM 1.0 alignment format");
($n_seq, $sequence_identifiers, $alignment, $alignment_id, $consensus_structure) = RNA::file_msa_read($datadir . "/rfam_seed_selected.stk",
                                                                                                      RNA::FILE_FORMAT_MSA_STOCKHOLM | RNA::FILE_FORMAT_MSA_SILENT);

# first alignment in test file has 712 sequences
ok($n_seq == 712);
ok(scalar(@{$sequence_identifiers}) == 712);
ok(scalar(@{$alignment}) == 712);
# test alignment contains an ID as well as a SS_cons line, so we must have retrieved both
ok(defined($alignment_id));
ok(defined($consensus_structure));
# length of consensus structure must be equal to length of the sequences in alignment (we check against 1st sequence)
ok(length($consensus_structure) == length($alignment->[0]));

# lets parse the entire STOCKHOLM file now, it contains multiple alignments
MsgChecking("whether we are able to parse multi STOCKHOLM 1.0 alignment format");
open $fh, "<" . $datadir . "/rfam_seed_selected.stk";

$counter = 0;
while ( (($n_seq, $sequence_identifiers, $alignment, $alignment_id, $consensus_structure) = RNA::file_msa_read_record($fh)) && ($n_seq != -1)) {

  # skip empty alignments, i.e. number of sequences == 0
  next if $n_seq == 0;

  $counter++;
}

close($fh);
# There are 4 alignments in the STOCKHOLM file
is($counter, 4);


MsgChecking("whether we are able to parse CLUSTAL alignment format");
($n_seq, $sequence_identifiers, $alignment, $alignment_id, $consensus_structure) = RNA::file_msa_read($datadir . "/070313_ecoli_cdiff_16S_clustalw.aln",
                                                                                                      RNA::FILE_FORMAT_MSA_CLUSTAL | RNA::FILE_FORMAT_MSA_SILENT);

# alignment in test file has 2 sequences
ok($n_seq == 2);
ok(scalar(@{$sequence_identifiers}) == 2);
ok(scalar(@{$alignment}) == 2);
# CLUSTAL format does not provide ID or consensus structure information
is($alignment_id, "");
is($consensus_structure, "");


MsgChecking("whether we are able to parse FASTA alignment format");
($n_seq, $sequence_identifiers, $alignment, $alignment_id, $consensus_structure) = RNA::file_msa_read($datadir . "/070313_ecoli_cdiff_16S_fasta.aln",
                                                                                                      RNA::FILE_FORMAT_MSA_FASTA | RNA::FILE_FORMAT_MSA_SILENT);

# alignment in test file has 2 sequences
ok($n_seq == 2);
ok(scalar(@{$sequence_identifiers}) == 2);
ok(scalar(@{$alignment}) == 2);
# FASTA format does not provide ID or consensus structure information
is($alignment_id, "");
is($consensus_structure, "");


MsgChecking("whether we are able to parse MAF alignment format");
($n_seq, $sequence_identifiers, $alignment, $alignment_id, $consensus_structure) = RNA::file_msa_read($datadir . "/test.maf",
                                                                                                      RNA::FILE_FORMAT_MSA_MAF | RNA::FILE_FORMAT_MSA_SILENT);

# first alignment in test file has 5 sequences
ok($n_seq == 5);
ok(scalar(@{$sequence_identifiers}) == 5);
ok(scalar(@{$alignment}) == 5);
# MAF format does not provide ID or consensus structure information
is($alignment_id, "");
is($consensus_structure, "");

# lets parse the entire MAF file now, it contains multiple alignments
MsgChecking("whether we are able to parse multi MAF alignment format");
open $fh, "<" . $datadir . "/test.maf";

$counter = 0;
while ( (($n_seq, $sequence_identifiers, $alignment, $alignment_id, $consensus_structure) = RNA::file_msa_read_record($fh, RNA::FILE_FORMAT_MSA_MAF)) && ($n_seq != -1)) {

  # skip empty alignments, i.e. number of sequences == 0
  next if $n_seq == 0;

  $counter++;
}

close($fh);
# There are 3 alignments in the MAF file
is($counter, 3);
