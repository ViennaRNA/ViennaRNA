import RNApath

RNApath.addSwigInterfacePath()


import RNA
import unittest

datadir = RNApath.getDataDirPath()

class file_utils_msa_Test(unittest.TestCase):

    def test_file_msa_detect_format_stk(self):
        print "test_file_msa_detect_STOCKHOLM_format"
        msa_format = RNA.file_msa_detect_format(datadir + "/rfam_seed_selected.stk")
        self.assertTrue(msa_format == RNA.FILE_FORMAT_MSA_STOCKHOLM)


    def test_file_msa_detect_format_clustal(self):
        print "test_file_msa_detect_CLUSTAL_format"
        msa_format = RNA.file_msa_detect_format(datadir + "/070313_ecoli_cdiff_16S_clustalw.aln")
        self.assertTrue(msa_format == RNA.FILE_FORMAT_MSA_CLUSTAL)


    def test_file_msa_detect_format_fasta(self):
        print "test_file_msa_detect_FASTA_format"
        msa_format = RNA.file_msa_detect_format(datadir + "/070313_ecoli_cdiff_16S_fasta.aln")
        self.assertTrue(msa_format == RNA.FILE_FORMAT_MSA_FASTA)


    def test_file_msa_detect_format_maf(self):
        print "test_file_msa_detect_MAF_format"
        msa_format = RNA.file_msa_detect_format(datadir + "/test.maf")
        self.assertTrue(msa_format == RNA.FILE_FORMAT_MSA_MAF)


    def test_file_msa_detect_format_unknown(self):
        print "test_file_msa_detect_non_MSA_format"
        msa_format = RNA.file_msa_detect_format(datadir + "/rnafold.cmds")
        self.assertTrue(msa_format == RNA.FILE_FORMAT_MSA_UNKNOWN)


    def test_file_msa_read_stk(self):
        print "test_file_msa_read_stk"
        n_seq, sequence_identifiers, alignment, alignment_id, consensus_structure = \
            RNA.file_msa_read(datadir + "/rfam_seed_selected.stk",
                              RNA.FILE_FORMAT_MSA_STOCKHOLM | \
                              RNA.FILE_FORMAT_MSA_SILENT)

        # first alignment in test file has 712 sequences
        self.assertTrue(n_seq == 712)
        self.assertTrue(len(sequence_identifiers) == 712)
        self.assertTrue(len(alignment) == 712)
        # test alignment contains an ID as well as a SS_cons line, so we must have retrieved both
        self.assertTrue(alignment_id != "")
        self.assertTrue(consensus_structure != "")
        # length of consensus structure must be equal to length of the sequences in alignment
        # (we check against 1st sequence)
        self.assertTrue(len(consensus_structure) == len(alignment[0]))


    def test_file_msa_read_clustal(self):
        print "test_file_msa_read_clustal"
        n_seq, sequence_identifiers, alignment, alignment_id, consensus_structure = \
            RNA.file_msa_read(datadir + "/070313_ecoli_cdiff_16S_clustalw.aln",
                              RNA.FILE_FORMAT_MSA_CLUSTAL | \
                              RNA.FILE_FORMAT_MSA_SILENT)

        # alignment in test file has 2 sequences
        self.assertTrue(n_seq == 2)
        self.assertTrue(len(sequence_identifiers) == 2)
        self.assertTrue(len(alignment) == 2)
        # CLUSTAL format does not provide ID or consensus structure information
        self.assertTrue(alignment_id == "")
        self.assertTrue(consensus_structure == "")


    def test_file_msa_read_fasta(self):
        print "test_file_msa_read_fasta"
        n_seq, sequence_identifiers, alignment, alignment_id, consensus_structure = \
            RNA.file_msa_read(datadir + "/070313_ecoli_cdiff_16S_fasta.aln",
                              RNA.FILE_FORMAT_MSA_FASTA | \
                              RNA.FILE_FORMAT_MSA_SILENT)

        # alignment in test file has 2 sequences
        self.assertTrue(n_seq == 2)
        self.assertTrue(len(sequence_identifiers) == 2)
        self.assertTrue(len(alignment) == 2)
        # FASTA format does not provide ID or consensus structure information
        self.assertTrue(alignment_id == "")
        self.assertTrue(consensus_structure == "")


    def test_file_msa_read_maf(self):
        print "test_file_msa_read_maf"
        n_seq, sequence_identifiers, alignment, alignment_id, consensus_structure = \
            RNA.file_msa_read(datadir + "/test.maf",
                              RNA.FILE_FORMAT_MSA_MAF | \
                              RNA.FILE_FORMAT_MSA_SILENT)

        # alignment in test file has 5 sequences
        self.assertTrue(n_seq == 5)
        self.assertTrue(len(sequence_identifiers) == 5)
        self.assertTrue(len(alignment) == 5)
        # MAF format does not provide ID or consensus structure information
        self.assertTrue(alignment_id == "")
        self.assertTrue(consensus_structure == "")


    def test_file_msa_read_multi_stk(self):
        print "test_file_msa_read_multi_stk"
        f = open(datadir + "/rfam_seed_selected.stk", 'r')
        counter = 0
        while True:
            n_seq, sequence_identifiers, alignment, alignment_id, consensus_structure = \
                RNA.file_msa_read_record(f)

            # stop parsing on error or EOF
            if n_seq == -1:
                break

            # skip empty alignments, i.e. number of sequences == 0
            if n_seq == 0:
                continue

            counter = counter + 1

        f.close()

        # There are 4 alignments in the STOCKHOLM file
        self.assertTrue(counter == 4)


    def test_file_msa_read_multi_maf(self):
        print "test_file_msa_read_multi_maf"
        f = open(datadir + "/test.maf", 'r')
        counter = 0
        while True:
            n_seq, sequence_identifiers, alignment, alignment_id, consensus_structure = \
                RNA.file_msa_read_record(f, RNA.FILE_FORMAT_MSA_MAF)

            # stop parsing on error or EOF
            if n_seq == -1:
                break

            # skip empty alignments, i.e. number of sequences == 0
            if n_seq == 0:
                continue

            counter = counter + 1

        f.close()

        # There are 4 alignments in the MAF file
        self.assertTrue(counter == 3)



if __name__ == '__main__':
    unittest.main()
